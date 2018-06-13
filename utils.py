from Bio.PDB import PDBList
from Bio import SeqIO
from bioservices import UniProt

LOCAL_STORAGE = "./local_files"

pdbl = PDBList(verbose=False)
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
