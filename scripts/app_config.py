from pipeline_functions import *
from functions import *
import time

class AppConfig:
    def __init__(self, UserChoice, CurrentSearchDirectory):
        
        # Initialise attributes that do not depend on user input, these are effectively constants, contingent only on the current organisation of the directory
        ## Attributes relating to blast databases
        self.possible_blasts = list_folders_in_directory('BLASTDatabases')
        self.valid_blast_inputs = [entry.split('/')[1].split('_')[0] for entry in self.possible_blasts]
        
        ## Attributes relating to FoldSeek databases
        self.database_path = 'databases/' # The databases that FoldSeek will look through are in this directory - there are 5
        Databases = list_folders_in_directory(self.database_path)
        self.databases = [f"{element.rsplit('/', 1)[0]}/{element.rsplit('/', 1)[1]}/{element.rsplit('/', 1)[1].lower()}" if '/' in element else element for element in Databases]
        assert (len(self.databases)==5), f'Five databases should be present in the database directory!'
        
        ## Attributes relating to FoldSeek output and command line commands
        self.format_output = ("--format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,tseq,evalue,bits,cigar,"
                "alntmscore,qtmscore,ttmscore,rmsd,u,t,lddt,lddtfull,prob,taxid,taxname,taxlineage") 

        self.columns = ['query','target','fident','alnlen','mismatch','gapopen','qstart','qend','tstart','tend','tseq','evalue','bits','cigar',
               'alntmscore','qtmscore','ttmscore','rmsd','u','t','lddt','lddtfull','prob','taxid','taxname','taxlineage','database']
    
        # The ESMA database has no taxonomic data, and so requires its own handling
        self.format_output_esma = ("--format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,tseq,evalue,bits,cigar,"
                    "alntmscore,qtmscore,ttmscore,rmsd,u,t,lddt,lddtfull,prob") 
        
        
        # Initialise two core attributes common to all functions which are user defined
        self.user_choice = UserChoice # attributes of a class are in snake case to differentiate them from any variable counterparts
        self.current_search_directory = CurrentSearchDirectory
        
        # Initialise other core common attributes that are contingent upon the prior user defined variables
        self.query_path = f'{CurrentSearchDirectory}/ProteinQueries'
        assert os.path.exists(self.query_path), f"Directory does not exist: {QueryPath}" # Assert that the folder with protein queries already exists - this folder should always exist otherwise how could FoldSeek have run

        self.finished_foldseek_directory = f"{self.current_search_directory}/FinishedFoldSeekProteins" # CSV files that are generated go here
        self.csv_file_directory = f"{self.current_search_directory}/FoldSeekExtractionOutputCSVs" # CSV files that are generated go here
        self.blast_file_directory = f'{self.current_search_directory}/BLAST_XMLs' # BLASTs will be saved and parsed as .xml files

        check_existence(self.finished_foldseek_directory)
        check_existence(self.csv_file_directory)
        check_existence(self.blast_file_directory)

        
        self.setup_config()

    def setup_config(self):
        if self.user_choice in [1, 4]:
            self.setup_for_running_foldseek()
        if self.user_choice in [2, 4]:
            self.setup_for_running_blasts()
        if self.user_choice in [3, 4]:
            self.setup_for_applying_filters()


    def setup_for_running_foldseek(self):
        # Configuration setup for function 1
        self.query_proteins = list_files_in_directory(self.query_path)
        assert (len(self.query_proteins)>0), f'No proteins in directory: {self.query_path}' # FoldSeek gots to get those protein files!

        self.aln_tmp_file_directory = f'{self.current_search_directory}/FoldSeekAlignments' # BLASTs will be saved and parsed as .xml files
        check_existence(self.aln_tmp_file_directory)

        CreateHTMLs = input("Do you want to generate HTMLs to visualise alignments (this will double runtime)? Type 'Yes' if so.\n")
        if CreateHTMLs.strip().lower() == 'yes':
            print('HTML files will be created.')
            self.create_htmls = True
            self.html_file_directory = f"{self.current_search_directory}/FoldSeekExtractionOutputHTMLs"
            check_existence(self.html_file_directory)
        
    def setup_for_running_blasts(self):

        # Initialise empty lists to be filled with user inputs
        self.blast_list = []
        self.blast_bitscore_thresholds = []
        
        # Ask for additional user input specific to running BLASTs 
        print("Choose which BLAST databases you'd like to run downstream (choose at least 1) from the list below, and type out EXACTLY as written here")
        for BLAST in self.possible_blasts:
            print(BLAST.split('/')[1].split('_')[0])
            
        while True:
            # Ask the user for a single file directory
            user_input = input("Enter a BLAST database (or 'done' to finish): \n")
            
            # Check if the user has finished entering directories
            if user_input.lower() == 'done':
                break
            else:
                # Add the entered directory to the list
                assert user_input.strip() in self.valid_blast_inputs
                self.blast_list.append(user_input)  # Remove any leading/trailing whitespace
                user_input = input("Enter the bitscore threshold you would like to apply to this BLAST (you can put 0 if you don't want any filtering to be  done)\n")
                self.blast_bitscore_thresholds.append(float(user_input.strip()))
    
    def setup_for_applying_filters(self):
        # Get user input for thresholds and filters
        TMScoreThreshold = input("How high would you like the TM-score threshold to be?\n")
        self.tm_score_threshold = float(TMScoreThreshold)
        Methionine = input("Would you like to filter out proteins that do not start with a Methionine?\n")
        if Methionine.strip().lower() == 'yes':
            self.must_start_with_methionine = True
            print('A CSV containing only hits beginning with Methionines will be generated after structural and BLAST filtering')
        else:
            self.must_start_with_methionine = False
        
        Duplicates = input('Would you like to remove duplicate hits at the end of filtering?\n')
        if Duplicates.strip().lower() == 'yes':
            self.remove_duplicates = True
            print('Duplicates will be removed after structural and BLAST filtering')
        else:
            self.remove_duplicates = False
            

