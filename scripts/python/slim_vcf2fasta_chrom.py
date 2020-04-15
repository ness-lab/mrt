'''
slim_vcf2fasta_chrom.py -
use a chlamy region as the 'ref' then a mutation
matrix to determine the alternate base at each variant site
'''

import argparse
import itertools

from tqdm import tqdm
from cyvcf2 import VCF
from Bio import SeqIO
import numpy as np


def args():
    parser = argparse.ArgumentParser(
        description='slim vcf -> fasta w/ reference',
        usage='python3.5 slim_vcf2fasta_chrom.py [options]')

    parser.add_argument('-v', '--vcf', required=True,
                        type=str, help='bgzipped + tabixed vcf')
    parser.add_argument('-f', '--fasta', required=True,
                        type=str, help='reference fasta')
    parser.add_argument('-r', '--region', required=True,
                        type=str, help='samtools format region (1 index)')
    parser.add_argument('-m', '--mut_mat', required=True,
                        type=str, help='LDhelmet mut mat file')
    parser.add_argument('-o', '--outfile', required=True,
                        type=str, help='File to write to')

    args = parser.parse_args()

    return args.vcf, args.fasta, args.region, args.mut_mat, args.outfile


def prep_samples(vcf: str, fasta: str, region: str): -> list, dict
    """Prepare samples and reference sequences to populate with alt alleles.

    Args:
        vcf: a bgzipped and tabixed SLiM VCF output file.
        fasta: a FASTA containing the reference genome (or at minimum, 
            the reference chromosome(s) of interest)
        region: samtools formatted region to extract from reference
            FASTA.

    Returns:
        - A list with phased sample names.
        - A dictionary with phased sample names as keys and the desired
              reference sequence as values.
    """
    print('prepping samples...')
    vcf_in = VCF(vcf)
    samples = vcf_in.samples
    samples_phased = []
    # read in reference sequence
    print('reading in reference sequence...')
    chrom, coords = region.split(':')
    start, end = [int(n) for n in coords.split('-')]
    seq_length = (end - start) + 1
    for record in SeqIO.parse(fasta, 'fasta'):
        if chrom in str(record.id):
            ref_seq = str(record.seq)[start-1:end]
            break
    assert len(ref_seq) == int(seq_length)
    for sample in samples:
        samples_phased.extend([sample + 'A', sample + 'B'])
    sequences = {}
    for sample in samples_phased:
        sequences[sample] = ref_seq
    return samples_phased, sequences


def convert_record(rec):
    """Obtain indices of alternate alleles at a given site.
    Used in tandem with get_alt_allele() to assign mutations.

    Args:
        rec: a VCF record, of type cyvcf2.cyvcf2.Variant

    Returns:
        A list of indices to convert to alternate alleles.
    """
    bases = rec.gt_bases
    split_bases = [pair.split('|') for pair in bases]
    combined_bases = list(itertools.chain.from_iterable(split_bases))
    indices_to_convert = [i for i, base in enumerate(combined_bases) if base == 'T']
    return indices_to_convert


def parse_mut_mat(mut_mat):
    """Parse space-separated mutation matrix.
    Mutation matrix should be ordered A, C, G, T in
    both dimensions.

    Args:
        mut_mat: Mutation matrix filename

    Returns:
        A nested dictionary, where top-level keys represent starting
        bases, second-level keys represent bases being mutated to, and
        values represent probabilities.
    """
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
    """Use mutation matrix to assign a alt base via a weighted draw.

    Args:
        ref_base: Reference base (str)
        mut_dict: Mutation matrix as dictionary, obtained from
            parse_mut_mat

    Returns:
        A single alternate base (str)
    """
    possible_bases = list(mut_dict[ref_base].keys())
    weights = [mut_dict[ref_base][alt_base] for alt_base in possible_bases]
    alt_allele = np.random.choice(possible_bases, 1, p=weights)[0]
    return alt_allele


def write_fasta(vcf, fasta, region, mut_mat, outfile):
    """Use above functions to convert a SLiM VCF to FASTA based off of a
    provided reference genome and region, assigning alternate bases
    where there are variants in the VCF.

    Args:
        vcf: SLiM output VCF filename
        fasta: FASTA containing reference sequence
        region: samtools format region to use as reference
        mut_mat: mutation matrix filename
        outfile: file to write to

    Returns:
        None
        Writes results to outfile. 
    """
    # Get start and end of sequence
    coords = region.split(':')[1]
    start, end = [int(n) for n in coords.split('-')]

    with open(outfile, 'w') as f:
        vcf_in = VCF(vcf)
        mut_dict = parse_mut_mat(mut_mat)
        samples_phased, sequences = prep_samples(vcf, table, region)
        for record in tqdm(vcf_in):
            pos = record.POS - 1  # gives python index
            indices_to_convert = convert_record(record)
            for i in indices_to_convert:
                if i % 2 == 1:
                    sample = 'i' + str(int((i - 1) / 2)) + 'B'
                else:
                    sample = 'i' + str(int(i / 2)) + 'A'
                alt = get_alt_allele(record.REF, mut_dict)
                sequences[sample] = sequences[sample][:pos] + alt + sequences[sample][record.POS:int(end + 1)]
        for sample in samples_phased:
            f.write('>' + sample + '\n')
            f.write(sequences[sample] + '\n')


def main():
    write_fasta(*args())


if __name__ == '__main__':
    main()
