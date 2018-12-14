import argparse
import json
import numpy as np

from collections import Counter

from utils import (calculate_aa_freq, calculate_asa_avg, calculate_distribution,
                   calculate_dssp, calculate_overal_rel_asa, calculate_ss_freq,
                   check_num_fixed_peptides, check_num_random_peptides,
                   filter_epitopes, get_b_factors, get_epitopes_locations,
                   get_pdb_sequence, get_uniprot_sequence, plot_aa_freq,
                   plot_ss_freq, check_average_exp_residues,
                   check_num_fixed_peptides_all_peaks,
                   check_average_exp_residues_all_peaks)

parser = argparse.ArgumentParser(description='FVIII Epitope Identification.')
parser.add_argument('--pdb_id', help="pdb_id to use", default='2r7e')
parser.add_argument(
    '--bfactors', help='calculate and print bfactors', action="store_true")
parser.add_argument(
    '--dssp',
    help='calculate relative asa and secondary structure',
    action="store_true")
parser.add_argument(
    '--dssp_random',
    help='calculate exposed residues with random sampling',
    action="store_true")
parser.add_argument(
    '--plots', help='plot aa_freq and ss_freq', action="store_true")

args = parser.parse_args()

pdb_id = args.pdb_id

try:
    pdb_seq = get_pdb_sequence(pdb_id)
    pdb_seq = pdb_seq[1:]
    uniprot_seq = get_uniprot_sequence(pdb_id)

    print('Uniprot Sequence >>> {}\n'.format(uniprot_seq))
    print('Pdb Sequence >>> {}\n'.format(pdb_seq))

except Exception as ex:
    print('Error retrieving files for PDB: {}\n Error: {}'.format(pdb_id, ex))

if args.bfactors:
    with open('local_data/epitopes.json', 'r') as epitopes_file:
        epitopes_ls = json.load(epitopes_file)
        locations = get_epitopes_locations(pdb_seq, epitopes_ls)
        for loc in locations:
            try:
                residues, avg = get_b_factors(pdb_id, loc['start'], loc['end'])
                print('Peptide location: {}-{} Average bfactor: {}'.format(
                    loc['start'], loc['end'], avg))
                for res in residues:
                    print('Residue Name: {} ID: {} Average bfactor: {}'.format(
                        res['name'], res['index'], res['avg_bfactor']))
            except Exception as ex:
                print('Error getting bfactor for epitope in location: {}-{}'.
                      format(loc['start'], loc['end']))

if args.dssp:
    dssp_ls = calculate_dssp(pdb_id)
    overal_rel_asa = calculate_overal_rel_asa(dssp_ls)
    print(overal_rel_asa)

if args.dssp_random:
    dssp_ls = calculate_dssp(pdb_id)

    sample_size = 10000

    print('Calculating exposed residues...')
    print('Sample size: {}'.format(sample_size))

    total_percent = 0
    total_percent_fixed = 0
    count = 0
    count_all = 0
    count_fixed = 0
    total_avg_exp_res_5 = 0
    total_avg_exp_res_10 = 0
    total_avg_exp_res_15 = 0
    total_avg_exp_res_20 = 0
    total_avg_exp_res_15_all = 0
    total_percent_fixed_all = 0

    percentages_random = {}
    for i in range(0, sample_size):
        percent, dist_ran = check_num_random_peptides(12, dssp_ls)
        total_percent += percent
        rounded_per = int(round(percent))

        if percent >= 91:
            count += 1
        if rounded_per not in percentages_random:
            percentages_random[rounded_per] = 1
        else:
            percentages_random[rounded_per] += 1

        if rounded_per >= 91:
            count += 1

    percentages_fixed = {}

    for i in range(0, sample_size):
        percent_fixed, dist = check_num_fixed_peptides(dssp_ls)
        avg_exp_res_5, avg_exp_res_10, avg_exp_res_15 = check_average_exp_residues(
            dssp_ls)

        total_avg_exp_res_5 += avg_exp_res_5
        total_avg_exp_res_10 += avg_exp_res_10
        total_avg_exp_res_15 += avg_exp_res_15

        # total_5 += d_5
        # total_10 += d_10
        # total_15 += d_15

        total_percent_fixed += percent_fixed
        rounded_per = int(round(percent_fixed))

        if rounded_per not in percentages_fixed:
            percentages_fixed[rounded_per] = 1
        else:
            percentages_fixed[rounded_per] += 1
        if rounded_per >= 91:
            count_fixed += 1

    # calculations for all the peaks

    for i in range(0, sample_size):
        percent_fixed_all = check_num_fixed_peptides_all_peaks(dssp_ls)
        avg_exp_res_15_all, avg_exp_res_20 = check_average_exp_residues_all_peaks(
            dssp_ls)

        total_percent_fixed_all += percent_fixed_all

        total_avg_exp_res_15_all += avg_exp_res_15_all
        total_avg_exp_res_20 += avg_exp_res_20
        if percent_fixed_all >= 97:
            count_all += 1

    # print('Frequencies of exposed residues in peptides with length 5: {}'.
    #       format(total_5))
    # print('Frequencies of exposed residues in peptides with length 10: {}'.
    #       format(total_10))
    # print('Frequencies of exposed residues in peptides with length 15: {}'.
    #       format(total_15))
    print(" ")

    print('Percentage: {0:.2f}%'.format(total_percent / sample_size))
    print("More or equal to 91% out of {}: {}".format(sample_size, count))
    for k, v in sorted(percentages_random.items()):
        print('{}%: {}'.format(k, v))

    print('Percentage fixed: {0:.2f}%'.format(
        total_percent_fixed / sample_size))
    print('Average of exposed residues in 5: {}'.format(
        (total_avg_exp_res_5 / sample_size)))
    print('Average of exposed residues in 10: {}'.format(
        (total_avg_exp_res_10 / sample_size)))
    print('Average of exposed residues in 15: {}'.format(
        (total_avg_exp_res_15 / sample_size)))
    print("More or equal to 91% out of {}: {}".format(sample_size, count_fixed))
    for k, v in sorted(percentages_fixed.items()):
        print('{}%: {}'.format(k, v))

#print for all the peaks
    print(" ")
    print("Calculations for all peaks")
    print('Percentage: {0:.2f}%'.format(total_percent_fixed_all / sample_size))
    print("More or equal to 97% out of {}: {}".format(sample_size, count_all))
    print('Average of exposed residues in 15: {}'.format(
        (total_avg_exp_res_15_all / sample_size)))
    print('Average of exposed residues in 20: {}'.format(
        (total_avg_exp_res_20 / sample_size)))

if args.plots:
    with open('local_data/dssp.json', 'r') as dssp_file, open(
            'local_data/epitopes.json', 'r') as epitopes_file:
        dssp_ls = json.load(dssp_file)
        epitopes_ls = json.load(epitopes_file)
        filtered_epitopes_ls = filter_epitopes(epitopes_ls, dssp_ls)

        # calculate ss frequencies
        overal_ss_freqs = calculate_ss_freq(dssp_ls)
        ep_ss_freqs = calculate_ss_freq(filtered_epitopes_ls)
        # calculate aa frequencies
        overal_aa_count, overal_aa_freqs = calculate_aa_freq(dssp_ls)
        ep_aa_count, ep_aa_freqs = calculate_aa_freq(filtered_epitopes_ls)

        print('ss')
        print(overal_ss_freqs)
        print(ep_ss_freqs)

        print('aa')
        print(overal_aa_freqs)
        print(ep_aa_freqs)

        print(calculate_distribution(dssp_ls))
        print(calculate_distribution(filtered_epitopes_ls))

        plot_ss_freq(overal_ss_freqs, ep_ss_freqs)
        plot_aa_freq(overal_aa_freqs, ep_aa_freqs)

        print(calculate_asa_avg(dssp_ls))
        print(calculate_asa_avg(filtered_epitopes_ls))
