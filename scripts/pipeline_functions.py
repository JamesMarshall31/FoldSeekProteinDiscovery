import os
import numpy as np
import pandas as pd
import time
import glob
from functions import *

def run_foldseek(config):

    print("FoldSeek will commence in 10 seconds") # Just to give some suspense to the whole thing
    time.sleep(10)

    print("Running FoldSeek")
    FirstFoldSeekdf = pd.DataFrame(columns=config.columns) # This is the first major dataframe, it will be saved at the end of the function
    # Loop through all the Query proteins
    for QueryProtein in config.query_proteins:
        PDBidentifier = QueryProtein.split('/')[-1].split('.')[-2] # This should hold true for any path to any .pdb file, this variable is used for naming of files
    
        for Database in config.databases:     # Loop through all the databases
            
            #First make the HTMLS (This doesn't depend on if the database is ESMA or not)
            if config.create_htmls:
                    os.system(f"foldseek easy-search {QueryProtein} {Database} {config.html_file_directory}/{PDBidentifier}_{Database.split('/')[-1]}.html {config.aln_tmp_file_directory}/tmp --format-mode 3")
                    os.system(f"rm -r {config.aln_tmp_file_directory}/tmp") # remove the temporary file (for storage and sanity)
                
            # Generate the alignment file - with a different output format if the database is ESMA    
            if Database.split('/')[-1] == 'esma':
                os.system(f"foldseek easy-search {QueryProtein} {Database} {config.aln_tmp_file_directory}/aln_{PDBidentifier}_{Database.split('/')[-1]} {config.aln_tmp_file_directory}/tmp {config.format_output_esma}")
                df = pd.read_csv(f"{config.aln_tmp_file_directory}/aln_{PDBidentifier}_{Database.split('/')[-1]}",sep='\t', header=None) # Read the current alignment file (it is a tab separated file)
                df.columns = config.columns[:-4] # we have to remove the last 4 (taxonomic and databases)
            else:
                os.system(f"foldseek easy-search {QueryProtein} {Database} {config.aln_tmp_file_directory}/aln_{PDBidentifier}_{Database.split('/')[-1]} {config.aln_tmp_file_directory}/tmp {config.format_output}")
                df = pd.read_csv(f"{config.aln_tmp_file_directory}/aln_{PDBidentifier}_{Database.split('/')[-1]}",sep='\t', header=None) # Read the current alignment file (it is a tab separated file)
                df.columns = config.columns[:-1] # Databases is yet to be added
            
            df['database'] = Database # Add in what database the current alignments are from
            df = df.reindex(columns=config.columns) # To rearrange everything in its final order - this gives ESMA data three NaN columns
            
            os.system(f"rm -r {config.aln_tmp_file_directory}/tmp") # remove the temporary directory (for storage and sanity)
            FirstFoldSeekdf = pd.concat([df,FirstFoldSeekdf]) # Add the current alignment to the main dataframe

        os.system(f"mv {QueryProtein} {config.finished_foldseek_directory}") # Move the query protein to the 'done folder once all databases are run with it
    
    FirstFoldSeekdf.to_csv(f'{config.csv_file_directory}/FirstFoldSeekdf.csv',index=False) # Save this dataframe to a csv for future manipulation
    
    FirstFoldSeekdf = pd.read_csv(f'{config.csv_file_directory}/FirstFoldSeekdf.csv') # Just to make sure we're working with it and nothing went wrong in the writing of the csv - also has the nice effect of resetting the index
    
    ## Obtain other variables like query protein length and target protein length that will probably be helpful!
    FirstFoldSeekdf['QueryProteinLength'] = 0 # Initialise the column to be filled in later
    FirstFoldSeekdf['TargetProteinLength'] = FirstFoldSeekdf['tseq'].apply(lambda x: len(x))
    
    for query in list(FirstFoldSeekdf['query'].unique()): # we are looping through each unique query to calculate their length
        queryUnique = query.split('_')[0] # pdb files can contain multiple chains
        NumChains = sum(query.split('_')[0] in pdbfile for pdbfile in list(FirstFoldSeekdf['query'].unique())) # count chains - need to divide by this after as it will count amino acids from all chains!
        ProteinLength = calculate_protein_length(f'{config.finished_foldseek_directory}/{queryUnique}')/NumChains # at this point all PDBs will be in the finished folder
        FirstFoldSeekdf.loc[FirstFoldSeekdf['query'] == query, 'QueryProteinLength'] = ProteinLength

    FirstFoldSeekdf.to_csv(f'{config.csv_file_directory}/FirstFoldSeekdf.csv',index=False) # Final writing including all variables from this stage

    print(f'FoldSeek has finished running and the resulting csv file has been saved to {config.csv_file_directory}/FirstFoldSeekdf.csv')












def run_blasts(config):
    print("Running BLASTs")

    FirstFoldSeekdf = pd.read_csv(f'{config.csv_file_directory}/FirstFoldSeekdf.csv')
    ## We make this df so we run every •unique• sequence through BLAST and join back - rather than running the same sequence multiple times!
    UniqueBLASTDF = pd.DataFrame()
    UniqueBLASTDF['tseq'] = FirstFoldSeekdf['tseq'].unique()

    ## Actually running BLAST
    for BLAST in config.blast_list: # Loop through each BLAST we need to do and BLAST the whole 'tseq' column
        BLASTInput = glob.glob(f'BLASTDatabases/{BLAST}_BLAST_DB/*.pdb')[0].split('.pdb')[0] # this will be the .fa file for those rebuilt, or for those not rebuilt from fasta e.g. patent - just the string 'pataa' which is what BLAST wants - there should be a better way of achieving this
        BLASTOutput = f'{config.blast_file_directory}/{BLAST}_Output.xml' # Specific output .xml file for this BLAST
        GetBestBlastAlignmentBitscoreEff(UniqueBLASTDF['tseq'], # This function BLASTs the whole column at once
                                      BLASTInput,
                                      config.blast_file_directory,
                                      BLASTOutput)
        print(f'{BLAST} BLAST is finished and has been saved to {config.blast_file_directory}/{BLAST}_Output.xml')
    ## Add this BLAST info (via parsing) to the unique df we made - this could be done in the loop above but this splitting is nice conceptually
    for BLAST in config.blast_list:
        UniqueBLASTDF = add_blast(UniqueBLASTDF,f'{config.blast_file_directory}/{BLAST}_Output.xml',BLAST)
    BLASTeddf = pd.merge(FirstFoldSeekdf, UniqueBLASTDF, on='tseq', how='inner')
    BLASTeddf.to_csv(f'{config.csv_file_directory}/BLASTedFoldSeekdf.csv',index=False)
    print(f'All BLASTs have been run and their data added to the FoldSeek data in {config.csv_file_directory}/BLASTedFoldSeekdf.csv')

    df = pd.DataFrame({
    'BLAST List': config.blast_list,
    'BLAST Thresholds': config.blast_bitscore_thresholds})

    df.to_csv(f'{config.csv_file_directory}/BLASTFilterSettings.csv', index=False)

    print(f'Desired BLAST threshold settings have been saved to {config.csv_file_directory}/BLASTFilterSettings.csv')








def apply_filters(config):
    print("Applying filters")
    BLASTeddf = pd.read_csv(f'{config.csv_file_directory}/BLASTedFoldSeekdf.csv')
    BLASTThresholddf = pd.read_csv(f'{config.csv_file_directory}/BLASTFilterSettings.csv')
     ## Structural filtering

    StructFiltdf = BLASTeddf[(BLASTeddf['qtmscore']>=config.tm_score_threshold) & (BLASTeddf['ttmscore']>=config.tm_score_threshold)]
    StructFiltdf.to_csv(f'{config.csv_file_directory}/StructFiltdf.csv',index=False)

    ## BLAST filtering
    BLASTFiltdf = StructFiltdf
    for BLAST, thresh in zip(list(BLASTThresholddf['BLAST List']), list(BLASTThresholddf['BLAST Thresholds'])): # Make sure the lists are the same length!
        BLASTFiltdf = BLASTFiltdf[BLASTFiltdf[f'{BLAST} Blast bitscores'] <= thresh]
    BLASTFiltdf.to_csv(f'{config.csv_file_directory}/BLASTFiltdf.csv',index=False)

    ## Methionine and/or RemoveDuplicates filtering as specified

    if config.must_start_with_methionine:
        MethionineFiltdf = BLASTFiltdf[BLASTFiltdf['tseq'].str.startswith('M')]
        MethionineFiltdf.to_csv(f'{config.csv_file_directory}/MethionineFiltdf.csv',index=False)

    if config.remove_duplicates and config.must_start_with_methionine:
        DuplicateRemoveddf = MethionineFiltdf.drop_duplicates(subset='tseq')
        DuplicateRemoveddf.to_csv(f'{config.csv_file_directory}/DuplicateRemoveddf.csv',index=False)
    elif config.remove_duplicates:
        DuplicateRemoveddf = BLASTFiltdf.drop_duplicates(subset='tseq')
        DuplicateRemoveddf.to_csv(f'{config.csv_file_directory}/DuplicateRemoveddf.csv',index=False)