from utils import get_pdb_sequence, get_uniprot_sequence, get_missing_residues

pdb_id = input("Enter PDB ID: ")

if pdb_id is '':
    pdb_id = "2r7e"

try:
    pdb_seq = get_pdb_sequence(pdb_id)
    uniprot_seq = get_uniprot_sequence(pdb_id)

    print('Uniprot Sequence >>> {}\n'.format(uniprot_seq))
    print('Pdb Sequence >>> {}\n'.format(pdb_seq))

    residues = get_missing_residues(uniprot_seq, pdb_seq)

    print('# Residues >>> {}'.format(len(residues)))
    # for res in residues:
    #     print(res)

except Exception as ex:
    print(ex)
