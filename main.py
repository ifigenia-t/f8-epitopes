import argparse
import json

from utils import (calculate_aa_freq, calculate_asa_avg,
                   calculate_distribution, calculate_dssp,
                   calculate_overal_rel_asa, calculate_ss_freq,
                   check_num_random_peptides, filter_epitopes, get_b_factors,
                   get_epitopes_locations, get_pdb_sequence,
                   get_uniprot_sequence, plot_aa_freq, plot_ss_freq)

parser = argparse.ArgumentParser(description='FVIII Epitope Identification.')
parser.add_argument('--pdb_id', help="pdb_id to use", default='2r7e')
parser.add_argument('--bfactors', help='calculate and print bfactors',
                    action="store_true")
parser.add_argument('--dssp',
                    help='calculate relative asa and secondary structure',
                    action="store_true")
parser.add_argument('--dssp_random',
                    help='calculate exposed residues with random sampling',
                    action="store_true")
parser.add_argument('--plots', help='plot aa_freq and ss_freq',
                    action="store_true")

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
                        res['name'], res['index'], res['avg_bfactor']
                    ))
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
    for i in range(0, sample_size):
        percent = check_num_random_peptides(12, dssp_ls)
        total_percent += percent

    print('Percentage: {0:.2f}%'.format(total_percent/sample_size))

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
