import json

from utils import (generate_pepscan_seq, get_missing_residues,
                   get_pdb_sequence, get_uniprot_sequence,
                   peptide_list_to_pdb_files)

# pdb_id = input("Enter PDB ID: ")

# if pdb_id is '':
#     pdb_id = "2r7e"

# try:
#     pdb_seq = get_pdb_sequence(pdb_id)
#     pdb_seq=pdb_seq[1:]
#     uniprot_seq = get_uniprot_sequence(pdb_id)

#     print('Uniprot Sequence >>> {}\n'.format(uniprot_seq))
#     print('Pdb Sequence >>> {}\n'.format(pdb_seq))

#     residues = get_missing_residues(uniprot_seq, pdb_seq)

#     print('# Residues >>> {}'.format(len(residues)))
#     # for res in residue


# except Exception as ex:
#     print(ex)


try:
    with open('local_data/lists.json') as f:
        data = json.load(f)

    # linear = generate_pepscan_seq(data['pepscan_linear'], 5, False)
    # looped = generate_pepscan_seq(data['pepscan_looped'], 5, True)
    print('Generating PDB Structure Files...')
    peptide_list_to_pdb_files(data['pepscan_linear'])

except Exception as e:
    print(e)
