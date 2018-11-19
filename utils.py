import operator
import os
import random
import matplotlib.pyplot as plt
import numpy as np
import PeptideBuilder
from Bio import SeqIO
from Bio.PDB import DSSP, PDBIO, PDBList, PDBParser
from bioservices import UniProt
from PeptideBuilder import Geometry

LOCAL_STORAGE = "./local_files"

pdbl = PDBList(verbose=False)
# quickfix for a bioservices bug
UniProt._url = "https://www.uniprot.org"
u = UniProt(verbose=False)


def download_pdb_file(pdb_id):
    # Download pdb file
    return pdbl.retrieve_pdb_file(
        pdb_code=pdb_id, pdir=LOCAL_STORAGE, file_format="pdb")


def download_cif_file(pdb_id):
    return pdbl.retrieve_pdb_file(
        pdb_code=pdb_id, pdir=LOCAL_STORAGE)


def get_pdb_sequence(pdb_id):
    pdb_file = download_pdb_file(pdb_id)
    pdb_seq = ''

    # Open file
    handle = open(pdb_file, "rU")

    # Parse to retrieve the sequences for each chain
    for record in SeqIO.parse(handle, "pdb-seqres"):
        pdb_seq += str(record.seq)

    handle.close()

    return pdb_seq


def get_uniprot_sequence(pdb_id):
    id_mapping = u.mapping(fr="PDB_ID", to="ACC", query=pdb_id.upper())
    uni_id = id_mapping[pdb_id.upper()][0]

    return u.get_fasta_sequence(uni_id)


# Comparing PDB and UniProt sequences
# in order to identify residues missing from the structure.
def get_missing_residues(uniprot_seq, pdb_seq):
    """
    Comparing the structural and UniProt sequence to identify the different
    residues.
    Arguments:
    parameter 1: The UniProt sequence.
    parameter 2: The PDB sequence.
    Returns:
    This is a description of what is returned.
    """
    # pdb_seq_index = 0

    # results = []

    # for uniprot_seq_index, l in enumerate(uniprot_seq):
    #     if pdb_seq[pdb_seq_index] is None:
    #         res = {"location": uniprot_seq_index,
    #                "residue": uniprot_seq[uniprot_seq_index]}
    #         continue
    #     if uniprot_seq[uniprot_seq_index] == pdb_seq[pdb_seq_index]:
    #         pdb_seq_index += 1
    #     else:
    #         res = {"location": uniprot_seq_index,
    #                "residue": uniprot_seq[uniprot_seq_index]}
    #         results.append(res)
    #         pdb_seq_index += 1

    # return results

    uniprot_seq_index = 0
    pdb_seq_index = 0

    results = []

    for _ in range(len(uniprot_seq)):
        # print('{}:{}'.format(uniprot_seq[uniprot_seq_index],
        # pdb_seq[pdb_seq_index]))
        # print('{}:{}'.format(uniprot_seq_index, pdb_seq_index))
        if uniprot_seq[uniprot_seq_index] == pdb_seq[pdb_seq_index]:
            # print('match')
            uniprot_seq_index += 1
            pdb_seq_index += 1
        else:
            d = {"location": uniprot_seq_index, "residue": uniprot_seq
                 [uniprot_seq_index]}
            results.append(d)
            uniprot_seq_index += 1

    return results


def generate_linear_peptides(pdb_seq):
    """
    Assembling all the linear Pepscan peptides into the protein sequence.
    parameter: The PDB sequence.
    Returns: A list containing all the looped peptides.
    """
    step = 5
    peptide_mer = 20
    linear_list = []

    for i in range(0, len(pdb_seq), step):
        peptide = pdb_seq[i:i+peptide_mer]
        if len(peptide) is 20:
            linear_list.append(peptide)

    return linear_list


def generate_looped_peptides(pdb_seq):
    """
    Assembling all the looped Pepscan peptides into the protein sequence,
    adding a C in the begging and in the end of every peptide respectfully.
    parameter: The PDB sequence.
    Returns: A list containing all the looped peptides.
    """
    step = 5
    peptide_mer = 15
    looped_list = []

    for i in range(0, len(pdb_seq), step):
        peptide = pdb_seq[i:i+peptide_mer]
        if len(peptide) is 15:
            looped_list.append('C{}C'.format(peptide))

    return looped_list


# Generate Pepscan Sequence:
def generate_pepscan_seq(peptides_list, step, skip_first=False):
    if skip_first is True:
        start = 1
    else:
        start = 0

    pepscan_seq_lp = ''

    for i in peptides_list:
        pepscan_seq_lp += i[start:step]

    return pepscan_seq_lp

# Generate .pdb files of the peptides with PeptideBuilder


def generate_peptide_pdb(peptide, filename):
    for i, res in enumerate(peptide):
        geo = Geometry.geometry(res)
        if i == 0:
            structure = PeptideBuilder.initialize_res(geo)
        else:
            structure = PeptideBuilder.add_residue(structure, geo)

    out = PDBIO()
    out.set_structure(structure)
    out.save(filename)


def peptide_list_to_pdb_files(peptide_list):
    OUTPUT_FOLDER = 'output_files'

    if not os.path.exists(OUTPUT_FOLDER):
        os.makedirs(OUTPUT_FOLDER)

    for i, peptide in enumerate(peptide_list):
        filename = '{}/peptide_{}.pdb'.format(OUTPUT_FOLDER, i)
        generate_peptide_pdb(peptide, filename)


# Calculate the aa frequencies
def count_aa(peptides, looped=False):
    count = {}
    freq = {}

    pep_len = 20

    for peptide in peptides:
        if looped is True:
            pep_len = 15
            # if list is looped_poly or looped_mono:
            peptide = peptide[1:16]

        for aa in peptide:
            if aa not in count:
                count[aa] = 1
            else:
                count[aa] += 1

    total_len = pep_len * len(peptides)

    sorted_count = sorted(
        count.items(), key=operator.itemgetter(1), reverse=True)

    freq = {}

    for k, v in sorted_count:
        freq[k] = v/total_len

    sorted_freq = sorted(
        freq.items(), key=operator.itemgetter(1), reverse=True)

    return sorted_count, sorted_freq


def calculate_dssp(pdb_code):
    pdbl = PDBList()
    parser = PDBParser()

    f = pdbl.retrieve_pdb_file(pdb_code=pdb_code, file_format="pdb")
    structure = parser.get_structure(pdb_code, f)

    dssp = DSSP(structure[0], f)

    dssp_ls = []
    for k in list(dssp.keys()):
        res = {
            "dssp_index": dssp[k][0],
            "amino_acid": dssp[k][1],
            "secondary_structure": dssp[k][2],
            "relative_asa": dssp[k][3]
        }
        dssp_ls.append(res)

    return dssp_ls


def calculate_overal_rel_asa(dssp_ls):
    count = 0

    for aa in dssp_ls:
        if aa["relative_asa"] >= 0.25:
            count += 1
    return count/len(dssp_ls)


def plot_dssp(dssp_ls):
    x = []
    y = []

    for aa in dssp_ls:
        x.append(aa['dssp_index'])
        y.append(aa['relative_asa'])

    plt.plot(x, y, 'ro')
    plt.show()


def calculate_aa_freq(dssp_ls):
    overal_freq = {}
    count = {}
    total_len = 0
    for aa in dssp_ls:
        if aa["relative_asa"] >= 0.25:
            total_len += 1

            if aa["amino_acid"] not in count:
                count[aa["amino_acid"]] = 1
            else:
                count[aa["amino_acid"]] += 1

    sorted_count = sorted(
        count.items())

    for k, v in sorted_count:
        overal_freq[k] = v/total_len

    sorted_overal_freq = sorted(
        overal_freq.items())

    return sorted_count, sorted_overal_freq


def plot_overal_freq(sorted_overal_freq):
    n_groups = len(sorted_overal_freq)
    freqs = [v for k, v in sorted_overal_freq]
    aa_names = [k for k, v in sorted_overal_freq]

    index = np.arange(n_groups)
    plt.bar(index, freqs, color='teal')
    plt.xticks(np.arange(len(aa_names)), aa_names)

    plt.title('Overal Frequency of Surface Amino Acids')
    ax = plt.subplot()
    ax.set_xlabel('Amino Acid')
    ax.set_ylabel('Frequency')

    plt.show()


def calculate_ss_freq(dssp_list):
    counts = {
        'loop': 0,
        'strand': 0,
        'helix': 0,
        'total': 0
    }

    for aa in dssp_list:
        if aa['relative_asa'] >= 0.25:
            counts['total'] += 1

            ss = aa['secondary_structure'].upper()

            if ss == '-' or ss == 'S' or ss == 'T':
                counts['loop'] += 1
            elif ss == 'E' or ss == 'B':
                counts['strand'] += 1
            elif ss == 'G' or ss == 'H' or ss == 'I':
                counts['helix'] += 1

    freqs = {
        'loop': counts['loop']/counts['total'],
        'strand': counts['strand']/counts['total'],
        'helix': counts['helix']/counts['total'],
    }

    return freqs


def get_seq(dssp_list):
    result = ''

    for aa in dssp_list:
        result += aa['amino_acid']

    return result


def plot_ss_freq(overal_freqs, ep_freqs):
    n_groups = 3
    o_ls = (overal_freqs['loop'], overal_freqs['strand'],
            overal_freqs['helix'])
    e_ls = (ep_freqs['loop'], ep_freqs['strand'], ep_freqs['helix'])

    fig, ax = plt.subplots()
    index = np.arange(n_groups)
    bar_width = 0.2

    plt.bar(index, o_ls, bar_width, color='teal', label='Overal')
    plt.bar(index + bar_width, e_ls, bar_width, color='indianred',
            label='Epitopes')

    plt.xlabel('Secondary Structure Elements')
    plt.ylabel('Frequency')

    plt.xticks(index + bar_width, ('Loop', 'Strand', 'Helix'))
    plt.legend()

    plt.tight_layout()
    plt.show()


def filter_epitopes(epitopes_ls, dssp_ls):
    ep_ls = []
    sequence_from_dssp = get_seq(dssp_ls)

    for ep in epitopes_ls:
        ep_index = sequence_from_dssp.find(ep)
        ep_end = ep_index + len(ep)

        for x in range(ep_index, ep_end):
            ep_ls.append(dssp_ls[x])

    return ep_ls


def get_epitopes_locations(pdb_seq, epitopes_ls):
    indexes = []

    for ep in epitopes_ls:
        ep_start = pdb_seq.find(ep)
        ep_end = ep_start + len(ep)

        indexes.append({'start': ep_start, 'end': ep_end})

    return indexes


def copy_freqs(overal, ep):
    res = []

    for o in overal:
        val = 0

        for e in ep:
            if o[0] == e[0]:
                val = e[1]
        res.append((o[0], val))

    return res


def plot_aa_freq(overal_freqs, ep_freqs):
    n_groups = 20

    ep_freqs_merged = copy_freqs(overal_freqs, ep_freqs)

    o_ls = tuple([v for k, v in overal_freqs])
    e_ls = tuple([v for k, v in ep_freqs_merged])
    ticks = tuple([k for k, v in overal_freqs])

    fig, ax = plt.subplots()
    index = np.arange(n_groups)
    bar_width = 0.3

    plt.bar(index, o_ls, bar_width, color='teal', label='Overal')
    plt.bar(index + bar_width, e_ls, bar_width, color='indianred',
            label='Epitopes')

    plt.xlabel('Amino Acids')
    plt.ylabel('Frequency')

    plt.xticks(index + bar_width, ticks)
    plt.legend()

    plt.tight_layout()
    plt.show()


def round_rate(x):
    return round(x / 0.05) * 0.05


def calculate_distribution(res_ls):
    dist = {}

    for res in res_ls:
        asa = res['relative_asa']

        if asa >= 0.25:
            k = "{0:.2f}".format(round_rate(asa))

            if k in dist:
                dist[k] += 1
            else:
                dist[k] = 1

    return dist


def calculate_asa_avg(res_ls):
    total = []

    for res in res_ls:
        if res['relative_asa'] >= 0.25:
            total.append(res['relative_asa'])

    return np.average(total)


def get_b_factors(pdb_id, index_start, index_end):
    # Retrieve the pdb file
    pdbl = PDBList()
    f = pdbl.retrieve_pdb_file(pdb_code=pdb_id, file_format="pdb")
    parser = PDBParser()
    structure = parser.get_structure(pdb_id, f)

    residues = []

    for res in structure.get_residues():
        if res.get_id()[1] > index_start and res.get_id()[1] <= index_end:
            r = {}
            r['name'] = res.get_resname()
            r['index'] = res.get_id()[1]
            r['atoms'] = []

            atom_counter = 0
            for atom in res.get_atoms():
                if atom_counter < 3:
                    r['atoms'].append({
                        'name': atom.get_name(),
                        'bfactor': atom.get_bfactor(),
                    })
                    atom_counter += 1

            r['avg_bfactor'] = np.average([a['bfactor'] for a in r['atoms']])
            residues.append(r)

    bfactors_ls = [res['avg_bfactor'] for res in residues]
    if len(bfactors_ls) > 0:
        avg = np.average(bfactors_ls)
    else:
        avg = 0
    return residues, avg


def check_random_peptide(start, end, peptide_len, dssp_ls):
    index =  random.randint(start, end-peptide_len)
    index_stop = index + peptide_len
    exposed_res = 0
    for res in dssp_ls:
        # find the index in dssp_ls
        if res["dssp_index"] >= index and res["dssp_index"] <= index_stop:
            # calculate the rel_asa of these res
            if res['relative_asa'] >= 0.25:
                # see how many of these res have rel_asa >= 0.25
                exposed_res += 1        
            
        # see if this peptide have >= 3 res with >= 0.25
    if exposed_res >= 3:
        return True
    
    return False


def check_num_random_peptides(num, dssp_ls):
    total = 0
    for i in range(1, num):
        if check_random_peptide(1, 1340, random.choice([5, 10, 15]), dssp_ls):
            total += 1

    return (total/num)*100

        