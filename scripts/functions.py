import pandas as pd
import subprocess
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SearchIO
import numpy
from Bio import PDB
import os
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning

def read_foldseek_output(json_file, chain_index=0):
    df = pd.read_json(json_file)
    df = pd.json_normalize(df['results'].iloc[0], record_path='alignments')
    return df


def write_dict_to_fasta(dict_with_sequences_as_values, file_name,textulating):
    with open(file_name, 'w+') as f:
        for name, seq in dict_with_sequences_as_values.items():
            f.write('>' + name + '\n')
            f.write(str(seq) + '\n')
    if textulating:
        print('%s SeqRecord objects were saved to "%s"' %(len(dict_with_sequences_as_values), file_name))
    return None


def get_Twist_price(len_nt, number_of_fragments=1, clonal=True, 
    linear_nt_fee=0.07, clonal_nt_fee=0.09, cloning_fee=25, VAT=1.2):
    nt_price = clonal_nt_fee if clonal else linear_nt_fee
    if clonal:
        return (len_nt * nt_price + cloning_fee * number_of_fragments) * VAT
    else:
        # does not check for lenth limits & number of fragments
        return (len_nt * nt_price) * VAT



def build_blast_database(path, dbtype='nucl', blast_command='makeblastdb'):
    command = '%s -in %s -dbtype %s' % (blast_command, path, dbtype)
    return subprocess.check_output(command, shell=True)


def run_blast(blast_algorithm, query_file, blastdb, output_file, textulating=False):
    assert ' ' not in query_file and ' ' not in blastdb and ' ' not in output_file
    command = '%s -db %s -query %s -outfmt 5 -out %s' % (blast_algorithm, blastdb, query_file, output_file)
    # print 'Command for BLAST:\n%s' % command
    subprocess.check_output(command, shell=True)
    if textulating:
        print('BLAST results have been saved to\n%s' %output_file)


def run_blast_for_all_files_in_folder(blast_algorithm, folder_with_query_files, blast_db, folder_for_output_files):
    # function not tested
    counter = 0
    for f in os.listdir(folder_with_query_files):
        counter += 1
        assert ' ' not in f
        output_file = f.rstrip('.fa') + '__AGAINST__' + blast_db.split('/')[-1].rstrip('.fa') + '.txt'
        run_blast(blast_algorithm,
                  os.path.join(folder_with_query_files, f),
                  blast_db,
                  os.path.join(folder_for_output_files, output_file))
    print('%s files blasted using %s against %s' % (counter, blast_algorithm, blast_db.split('/')[-1].rstrip('.fa')))
    return counter


def blast_protein_at_NCBI(protein_sequence, remote_db='nr'):
    '''
    Runs a BLAST online using blastp algorhytm. Returns BlastRecord object
    :param protein_sequence: protein sequence as string
    :param remote_db: remote database at NCBI, i.e. 'nr' ot 'pdb'
    :return: BlastRecord object
    '''
    result_handle = NCBIWWW.qblast('blastp', remote_db, protein_sequence)
    blast_record = NCBIXML.read(result_handle)
    return blast_record


def get_blast_top_hit(blast_record):
    return blast_record.alignments[0].hsps[0]


def get_difference(blast_record_hsp):
    return blast_record_hsp.align_length - blast_record_hsp.identities


def parse_blast_output(file_with_blast_output, blast_format='blast-xml'):
    """
    Parses blast output file and returns an iterator over Biopython's BlastRecord objects
    :param file_with_blast_output: path to file
    :param blast_format: format of blast output (i.e. 'blast-text', 'blast-tab', 'blast-xml')
    :return: Iterator over BlastRecord objects
    """
    return SearchIO.parse(file_with_blast_output, blast_format.lower())


def get_hit_evalues(blast_record, ln=False):
    """
    Returns list of consisting of E-values of the first alignment for every hit in blast_record
    :param blast_record: BlastRecord or QueryRecord Biopython object
    :return: list of E-values
    """
    if ln:
        return [numpy.log(hit.hsps[0].evalue) for hit in blast_record]
    else:
        return [hit.hsps[0].evalue for hit in blast_record]

def GetBestBlastAlignmentBitscoreEff(sequences,known_toxins_fn,blast_output_folder,blast_output_fn):
    
    with open("tmp.fasta", "w") as fasta_file:
        for index, value in enumerate(sequences):
            # Write each element of the column to the FASTA file
            seq_record = f">Seq{index+1}\n{value}\n"
            fasta_file.write(seq_record)
    run_blast('blastp', 'tmp.fasta', 
            known_toxins_fn, blast_output_fn, textulating=False)
    os.system("rm tmp.fasta")
    
def calculate_protein_length(pdb_file_path):
    warnings.simplefilter('ignore', PDBConstructionWarning)

    structure = PDB.PDBParser().get_structure('protein', pdb_file_path)

    protein_length = 0

    for model in structure:
        for chain in model:
            for residue in chain:
                if PDB.is_aa(residue):
                    protein_length += 1

    return protein_length
    
def filter_hits(blast_record, evalue_cutoff=0.001):
    return [hit for hit in blast_record if hit.hsps[0].evalue < evalue_cutoff]


def list_files_in_directory(directory):
    files_list = []
    for filename in os.listdir(directory):
        full_path = os.path.join(directory, filename)
        if os.path.isfile(full_path):
            files_list.append(full_path)
    return files_list


def list_folders_in_directory(directory):
    folders_list = []
    for filename in os.listdir(directory):
        full_path = os.path.join(directory, filename)
        if os.path.isdir(full_path):
            folders_list.append(full_path)
    return folders_list
def check_existence(directory):
    if not os.path.exists(directory):
        # Create the directory (including intermediate directories)
        os.makedirs(directory)
        print(f"Directory '{directory}' was created.")
    else:
        print(f"Directory '{directory}' already exists.")

def check_and_assert_existence(directory):
    assert os.path.exists(directory), f"Directory does not exist: {directory}"
    print(f"Directory '{directory}' already exists.")

def parse_blast(blast_output_fn):
    blast_format='blast-xml'
    blast_results = SearchIO.parse(blast_output_fn, blast_format.lower())
    #bitscore, bitscore_raw = 'hi','there'  # Default values in case there are no hits
    bitscores = []
    raw_bitscores =[]
    e_values=[]
    aln_lens=[]
    names = []
    positives = []
    identities = []
    gaps = []
    for query_result in blast_results:
        bitscore = 0
        bitscore_raw = 0
        e_value = 10
        aln_len = 0
        hits = query_result.hits
        name = 'no protein'
        positive = 0
        identity = 0
        gap = 0
        if len(hits)>0:
            bitscore = hits[0].hsps[0].bitscore
            bitscore_raw = hits[0].hsps[0].bitscore_raw
            e_value = hits[0].hsps[0].evalue
            aln_len = hits[0].hsps[0].aln.get_alignment_length()
            name = hits[0].description 
            positive = hits[0].hsps[0].pos_num
            identity = hits[0].hsps[0].ident_num
            gap = hits[0].hsps[0].gap_num
        bitscores.append(bitscore)
        raw_bitscores.append(bitscore_raw)
        e_values.append(e_value)
        aln_lens.append(aln_len)
        names.append(name)
        positives.append(positive)
        identities.append(identity)
        gaps.append(gap)
        
    return pd.Series(bitscores), pd.Series(raw_bitscores), pd.Series(e_values), pd.Series(aln_lens), pd.Series(names), pd.Series(positives), pd.Series(identities), pd.Series(gaps)

def add_blast(df,BLASTXML,BLAST):
    bitscores, raw_bitscores, evals, aln_lens,names,positives,identities,gaps  = parse_blast(BLASTXML)
    df[f'{BLAST} Blast bitscores'] = bitscores
    df[f'{BLAST} Blast raw_bistscores'] = raw_bitscores
    df[f'{BLAST} Blast eval'] = evals
    df[f'{BLAST} Blast aln len'] = aln_lens
    df[f'{BLAST} Blast name'] = names
    df[f'{BLAST} Blast positives'] = positives
    df[f'{BLAST} Blast identities'] = identities
    df[f'{BLAST} Blast gaps'] = gaps
    return df 