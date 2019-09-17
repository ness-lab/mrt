#!/usr/bin/env python3.6

# Script takes a directory filled with VCFs exported from SLiM and
# for each one, uses the site-frequency spectrum to calculate Ne
# from both Waterson and Pi. Exports results as single CSV.

from cyvcf2 import VCF
import os
import tqdm
import SFS as SFS


def create_sfs_dict(inpath):

    # Initialize sfs dict
    sfs_dict = {}

    # Get files from inpath
    files = os.listdir(inpath)

    for filename in tqdm.tqdm(files):
        # print(filename)

        if filename.endswith(".vcf"):
            # print(filename)

            # Load in VCF file
            my_vcf = VCF(inpath + '/' + filename)

            # Split VCF filename and split to extract population size, bottleneck and generation
            split_filename = filename.split('_')
            N_sims = int(split_filename[1].split('N')[1])
            bot = float(split_filename[2].split('bot')[1])
            gen = int(split_filename[3].split('gen')[1].split('.')[0])
            # print(N_sims, bot, gen)
            N_samples = len(my_vcf.samples)

            # Initialize list for sfs
            sfs_list = [0] * ((2 * N_samples) + 1)

            # Add singletons, doubletons, ...n-tons
            for variant in my_vcf:
                AC = variant.INFO.get('AC')
                # print(AC)
                # idx = line.find(';AC=', 0, 150)
                # if(idx == -1):
                #     continue
                # AC = int(line[idx + 4:idx + 8].split(';')[0])
                sfs_list[AC] += 1

            # Add invariant sites to sfs list
            ch_length = 1e8
            sfs_list[0] = int(ch_length - sum(sfs_list))
            # print(sfs_list)

            # Use sfs list to instantiate SFS class (Rob's code)
            sfs = SFS.SFS(sfs_list)

            # Create string for dictionary key
            l1 = str(N_sims) + '-' + str(bot)

            # Add SFS class instances to appropriate dictionary keys.
            if l1 in sfs_dict:
                sfs_dict[l1][gen] = sfs
            else:
                sfs_dict[l1] = {gen: sfs}

    return(sfs_dict)


def write_thetaNe_values(sfs_dict, outpath):

    # Open csv to write resutls
    outfile = open(outpath, 'w+')

    # Write header
    outfile.write('N,bot,gen,theta_pi,theta_w,Ne_pi,Ne_w\n')

    mu = float(1e-8)

    # Interate through sfs dictionary and write to csv
    for key in sfs_dict.keys():

        # Get dict values
        d = sfs_dict[key]

        l = d.keys()
        l = sorted(l)

        # Iterate through nested dictionary values (generations)
        for g in l:
            pi = d[g].theta_pi()
            watersons = d[g].theta_w()
            Ne_pi = (pi / (4 * mu))
            Ne_w = (watersons / (4 * mu))
#             print(g, d[g].theta_pi(), d[g].theta_w(),Ne_pi,Ne_w)
            N_and_Gen = key.split('-')
            outfile.write('{0},{1},{2},{3},{4},{5},{6}\n'.format(N_and_Gen[0], N_and_Gen[1], g, pi,watersons, Ne_pi, Ne_w))
    outfile.close()


inpath = "../../data/test_SLiM_output/"
sfs_dict = create_sfs_dict(inpath = inpath)

outpath = '../../data/clean/theta_NeValues.csv'
write_thetaNe_values(sfs_dict = sfs_dict, outpath = outpath)

