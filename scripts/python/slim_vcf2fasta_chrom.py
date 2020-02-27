'''
slim_vcf2fasta_chrom.py - 
use a chlamy region as the 'ref' then a mutation
matrix to determine the alternate base at each variant site
'''

import argparse
from tqdm import tqdm
from cyvcf2 import VCF
from Bio import SeqIO
import numpy as np
import itertools
import ant
import random

def args():
    parser = argparse.ArgumentParser(
        description='slim vcf -> fasta w/ reference', 
        usage='python3.5 slim_vcf2fasta_chrom.py [options]')

    parser.add_argument('-v', '--vcf', required=True,
                        type=str, help='bgzipped + tabixed vcf')
    parser.add_argument('-t', '--table', required=True,
                        type=str, help='annotation table')
    parser.add_argument('-r', '--region', required=True,
                        type=str, help='samtools format region (1 index)')
    parser.add_argument('-m', '--mut_mat', required=True,
                        type=str, help='LDhelmet mut mat file')
    parser.add_argument('-d', '--downsample', required=False,
                        type=float, help='Percent of SNPs to downsample to')
    parser.add_argument('-o', '--outfile', required=True,
                        type=str, help='File to write to')

    args = parser.parse_args()

    return args.vcf, args.table, args.region, \
        args.mut_mat, args.downsample, args.outfile

def prep_samples(vcf, table, region):
    print('prepping samples...')
    vcf_in = VCF(vcf)
    samples = vcf_in.samples
    samples_phased = []
    p = ant.Reader(table)
    chrom, coords = region.split(':')
    start, end = [int(n) for n in coords.split('-')]
    ref_seq = ''.join(record.ref for record in tqdm(p.fetch(chrom, start-1, end)))
    assert len(ref_seq) == int(1e5)
    for sample in samples:
        samples_phased.extend([sample + 'A', sample + 'B'])
    # import reference sequence
    sequences = {}
    for sample in samples_phased:
        sequences[sample] = ref_seq
    return samples_phased, sequences

def convert_record(rec):
    bases = rec.gt_bases
    split_bases = [pair.split('|') for pair in bases]
    combined_bases = list(itertools.chain.from_iterable(split_bases))
    indices_to_convert = [i for i, base in enumerate(combined_bases) if base == 'T']
    return indices_to_convert

def parse_mut_mat(mut_mat):
    print('parsing mutation matrix...')
    mut_dict = {}
    bases = ['A', 'C', 'G', 'T']
    # prep dict
    for base in bases:
        mut_dict[base] = {}
        for alt_base in bases:
            if base != alt_base:
                mut_dict[base][alt_base] = 0.0
    # populate
    with open(mut_mat, 'r') as f:
        lines = f.readlines()
    for i, line in enumerate(lines):
        split = line.rstrip().split(' ')
        current_base = bases[i]
        for j, base in enumerate(bases):
            if i == j:
                continue
            else:
                mut_dict[current_base][base] = float(split[j])
    # rescale
    for base in mut_dict:
        denom = sum(mut_dict[base].values())
        for alt_base in mut_dict[base]:
            mut_dict[base][alt_base] = round(mut_dict[base][alt_base] / denom, 3)
    
    return mut_dict

def get_alt_allele(ref_base, mut_dict):
    possible_bases = list(mut_dict[ref_base].keys())
    weights = [mut_dict[ref_base][alt_base] for alt_base in possible_bases]
    alt_allele = np.random.choice(possible_bases, 1, p=weights)[0]
    return alt_allele

def write_fasta(vcf, table, region, mut_mat, downsample, outfile):
    with open(outfile, 'w') as f:
        vcf_in = VCF(vcf)
        mut_dict = parse_mut_mat(mut_mat)
        samples = vcf_in.samples
        samples_phased, sequences = prep_samples(vcf, table, region)
        total_count, skipped_count = 0, 0
        for record in tqdm(vcf_in):
            pos = record.POS - 1 # gives python index
            # if downsampling - don't assign some alternates
            if downsample:
                total_count += 1
                draw = random.random()
                if draw >= downsample:
                    skipped_count += 1
                    continue # skip alternate assignment
                else:
                    pass
            indices_to_convert = convert_record(record)
            for i in indices_to_convert:
                if i % 2 == 1:
                    sample = 'i' + str(int((i - 1) / 2)) + 'B'
                else:
                    sample = 'i' + str(int(i / 2)) + 'A'
                alt = get_alt_allele(record.REF, mut_dict)
                sequences[sample] = sequences[sample][:pos] + alt + sequences[sample][record.POS:int(1e5)]
        for sample in samples_phased:
            f.write('>' + sample + '\n')
            f.write(sequences[sample] + '\n')
        if downsample:
            print('{} skipped of {} SNPs.'.format(skipped_count, total_count))
            print('{}% of SNPs reverted'.format(round(skipped_count / total_count, 3)))

def main():
    write_fasta(*args())

if __name__ == '__main__':
    main()

        

