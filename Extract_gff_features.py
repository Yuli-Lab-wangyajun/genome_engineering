#!/usr/bin/env python3
import argparse
from collections import defaultdict

def parse_gff3(path):
    chr_len = {}
    genes = {}
    transcripts = {}
    pending = []  # 暂存exon/CDS，等转录本就位后处理

    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                if line.startswith('##sequence-region'):
                    p = line.split()
                    if len(p) >= 4:
                        chr_len[p[1]] = max(chr_len.get(p[1], 0), int(p[3]))
                continue
            p = line.strip().split('\t')
            if len(p) != 9:
                continue
            seqid, src, ftype, s, e, score, strand, phase, attr_s = p
            start, end = int(s), int(e)
            attrs = {}
            for item in attr_s.split(';'):
                if '=' in item:
                    k, v = item.split('=', 1)
                    attrs[k] = v
            rid = attrs.get('ID')
            parent = attrs.get('Parent')

            # 基因 / 假基因
            if ftype in ('gene', 'pseudogene') and rid:
                genes[rid] = {
                    'id': rid,
                    'name': attrs.get('Name', attrs.get('gene', rid)),
                    'seqid': seqid,
                    'strand': strand,
                    'start': start,
                    'end': end,
                    'type': ftype,
                    'transcripts': []
                }
            # 转录本 (任何Parent指向已知基因的行都视为转录本)
            elif parent and parent in genes:
                transcripts[rid] = {
                    'id': rid,
                    'parent': parent,
                    'seqid': seqid,
                    'strand': strand,
                    'exons': [],
                    'cds': []
                }
                genes[parent]['transcripts'].append(rid)
            # exon / CDS：若父转录本已存在，立即添加；否则暂存
            elif ftype in ('exon', 'CDS') and parent:
                if parent in transcripts:
                    if ftype == 'exon':
                        transcripts[parent]['exons'].append((start, end))
                    else:
                        transcripts[parent]['cds'].append((start, end))
                else:
                    pending.append((parent, ftype, start, end))

    # 处理暂存的 exon/CDS
    for parent, ftype, start, end in pending:
        if parent in transcripts:
            if ftype == 'exon':
                transcripts[parent]['exons'].append((start, end))
            else:
                transcripts[parent]['cds'].append((start, end))

    return chr_len, genes, transcripts

def clip(seqid, s, e, lengths):
    m = lengths.get(seqid)
    s = max(1, s)
    e = min(e, m) if m is not None else e
    return (s, e) if s <= e else None

def write_bed(path, data):
    with open(path, 'w') as f:
        for chrom, s, e, name, score, strand in sorted(data):
            f.write(f"{chrom}\t{s-1}\t{e}\t{name}\t{score}\t{strand}\n")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('gff3')
    ap.add_argument('-o', '--prefix', default='out')
    ap.add_argument('--chromosome')
    ap.add_argument('-up', '--upstream', type=int, default=2000)
    args = ap.parse_args()

    chr_len, genes, transcripts = parse_gff3(args.gff3)
    if args.chromosome:
        genes = {k: v for k, v in genes.items() if v['seqid'] == args.chromosome}

    gene_bed, exon_bed, cds_bed, intron_bed, up_bed, distal_bed = [], [], [], [], [], []

    gene_list = []
    for gid, g in genes.items():
        seqid, name, strand, s, e, gtype = g['seqid'], g['name'], g['strand'], g['start'], g['end'], g['type']
        gene_bed.append((seqid, s, e, name, gtype, strand))
        gene_list.append((s, e, strand, name, gid, seqid))

        if strand == '+':
            us, ue = s - args.upstream, s - 1
        else:
            us, ue = e + 1, e + args.upstream
        c = clip(seqid, us, ue, chr_len)
        if c:
            up_bed.append((seqid, c[0], c[1], name, 'upstream', strand))

    for gid, g in genes.items():
        seqid, name, strand = g['seqid'], g['name'], g['strand']
        for tid in g['transcripts']:
            t = transcripts[tid]
            exons = sorted(t['exons'])
            for s, e in exons:
                c = clip(seqid, s, e, chr_len)
                if c:
                    exon_bed.append((seqid, c[0], c[1], name, 'exon', strand))
            for s, e in t['cds']:
                c = clip(seqid, s, e, chr_len)
                if c:
                    cds_bed.append((seqid, c[0], c[1], name, 'CDS', strand))
            for i in range(len(exons)-1):
                ins = exons[i][1] + 1
                ine = exons[i+1][0] - 1
                c = clip(seqid, ins, ine, chr_len)
                if c:
                    intron_bed.append((seqid, c[0], c[1], name, 'intron', strand))

    # 远端基因间区
    chrom_genes = defaultdict(list)
    for g in gene_list:
        chrom_genes[g[5]].append(g)
    for seqid, group in chrom_genes.items():
        group.sort(key=lambda x: x[0])
        for i in range(len(group)-1):
            a, b = group[i], group[i+1]
            inter_s = a[1] + 1
            inter_e = b[0] - 1
            if inter_s > inter_e:
                continue
            if b[2] == '+':
                bu_s = b[0] - args.upstream
                bu_e = b[0] - 1
            else:
                bu_s = b[1] + 1
                bu_e = b[1] + args.upstream
            c_up = clip(seqid, bu_s, bu_e, chr_len)
            if c_up is None:
                c = clip(seqid, inter_s, inter_e, chr_len)
                if c:
                    distal_bed.append((seqid, c[0], c[1], f"{a[3]}/{b[3]}", 'intergenic', 'NA'))
            else:
                segs = [[inter_s, inter_e]]
                s1, e1 = c_up
                res = []
                for iv in segs:
                    if iv[1] < s1 or iv[0] > e1:
                        res.append(iv)
                    else:
                        if iv[0] < s1:
                            res.append([iv[0], s1-1])
                        if iv[1] > e1:
                            res.append([e1+1, iv[1]])
                for r in res:
                    c = clip(seqid, r[0], r[1], chr_len)
                    if c:
                        distal_bed.append((seqid, c[0], c[1], f"{a[3]}/{b[3]}", 'intergenic', 'NA'))

    write_bed(f"{args.prefix}_gene.bed", gene_bed)
    write_bed(f"{args.prefix}_exon.bed", exon_bed)
    write_bed(f"{args.prefix}_CDS.bed", cds_bed)
    write_bed(f"{args.prefix}_intron.bed", intron_bed)
    write_bed(f"{args.prefix}_upstream{args.upstream}.bed", up_bed)
    write_bed(f"{args.prefix}_distal_intergenic.bed", distal_bed)
    print("Done.")

if __name__ == '__main__':
    main()