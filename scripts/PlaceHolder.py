def RunFoldSeek(CurrentSearchDirectory,BLASTList,TMScoreThreshold,BLASTBitscoreThresholds,MustStartWithMethionine,RemoveDuplicates):

    BLASTDirectories = [f'BLASTDatabases/{BLAST}_BLAST_DB' for BLAST in BLASTList]
    for directory in BLASTDirectories: # make sure they exist!
        check_and_assert_existence(directory)

    ### Actually Running FoldSeek 
    
    FirstFoldSeekdf = pd.DataFrame(columns=Columns) # This is the first major dataframe, it will be the largest, and saved at the end of the cell
    # Loop through all the Query proteins
    for QueryProtein in QueryProteins:
        PDBidentifier = QueryProtein.split('/')[-1].split('.')[-2] # This should hold true for any path to any .pdb file, this variable is used for naming of files
    
        for Database in Databases:     # Loop through all the databases
            
            #First make the HTMLS (This doesn't depend on if the database is ESMA or not)
            if CreateHTMLs:
                    os.system(f"foldseek easy-search {QueryProtein} {Database} {HTMLFileDirectory}/{PDBidentifier}_{Database.split('/')[-1]}.html {AlignmentAndTempFileDirectory}/tmp --format-mode 3")
                    os.system(f"rm -r {AlignmentAndTempFileDirectory}/tmp") # remove the temporary file (for storage and sanity)
                
            # Generate the alignment file - with a different output format if the database is ESMA    
            if Database.split('/')[-1] == 'esma':
                os.system(f"foldseek easy-search {QueryProtein} {Database} {AlignmentAndTempFileDirectory}/aln_{PDBidentifier}_{Database.split('/')[-1]} {AlignmentAndTempFileDirectory}/tmp {FormatOutputESMA}")
                df = pd.read_csv(f"{AlignmentAndTempFileDirectory}/aln_{PDBidentifier}_{Database.split('/')[-1]}",sep='\t', header=None) # Read the current alignment file (it is a tab separated file)
                df.columns = Columns[:-4]
            else:
                os.system(f"foldseek easy-search {QueryProtein} {Database} {AlignmentAndTempFileDirectory}/aln_{PDBidentifier}_{Database.split('/')[-1]} {AlignmentAndTempFileDirectory}/tmp {FormatOutput}")
                df = pd.read_csv(f"{AlignmentAndTempFileDirectory}/aln_{PDBidentifier}_{Database.split('/')[-1]}",sep='\t', header=None) # Read the current alignment file (it is a tab separated file)
                df.columns = Columns[:-1]
            
            df['database'] = Database # Add in what database the current alignments are from
            df = df.reindex(columns=Columns) # Just to make sure everything is as it should be - the ESMA df should have three NaN columns
            
            os.system(f"rm -r {AlignmentAndTempFileDirectory}/tmp") # remove the temporary directory (for storage and sanity)
            FirstFoldSeekdf = pd.concat([df,FirstFoldSeekdf]) # Add the current alignment to the main dataframe
    
    FirstFoldSeekdf.to_csv(f'{CSVFileDirectory}/FirstFoldSeekdf.csv',index=False) # Save this dataframe to a csv for future manipulation
    
    FirstFoldSeekdf = pd.read_csv(f'{CSVFileDirectory}/FirstFoldSeekdf.csv') # Just to make sure we're working with it and nothing went wrong in the writing of the csv - also has the nice effect of resetting the index
    
    ## Obtain other variables like query protein length and target protein length that will probably be helpful!
    FirstFoldSeekdf['QueryProteinLength'] = 0 # Lets initialise these columns and fill them in, this data may be useful downstream
    FirstFoldSeekdf['TargetProteinLength'] = FirstFoldSeekdf['tseq'].apply(lambda x: len(x))
    
    for query in list(FirstFoldSeekdf['query'].unique()): # we are looping through each unique query to calculate their length
        queryUnique = query.split('_')[0] # pdb files can contain multiple chains
        NumChains = sum(query.split('_')[0] in pdbfile for pdbfile in list(FirstFoldSeekdf['query'].unique())) # count chains - need to divide by this after as it will count amino acids from all chains!
        ProteinLength = calculate_protein_length(f'{QueryPath}/{queryUnique}')/NumChains
        FirstFoldSeekdf.loc[FirstFoldSeekdf['query'] == query, 'QueryProteinLength'] = ProteinLength

    FirstFoldSeekdf.to_csv(f'{CSVFileDirectory}/FirstFoldSeekdf.csv',index=False) # Final writing including all variables from this stage
    
    ### Running the BLASTS

    

    ### Now do all the filtering

    ## Structural filtering

    StructFiltdf = BLASTeddf[(BLASTeddf['qtmscore']>=TMScoreThreshold) & (BLASTeddf['ttmscore']>=TMScoreThreshold)]
    StructFiltdf.to_csv(f'{CSVFileDirectory}/StructFiltdf.csv',index=False)

    ## BLAST filtering
    BLASTFiltdf = StructFiltdf
    for BLAST, thresh in zip(BLASTList, BLASTBitscoreThresholds): # Make sure the lists are the same length!
        BLASTFiltdf = BLASTFiltdf[BLASTFiltdf[f'{BLAST} Blast bitscores'] <= thresh]
    BLASTFiltdf.to_csv(f'{CSVFileDirectory}/BLASTFiltdf.csv',index=False)

    ## Methionine and/or RemoveDuplicates filtering as specified

    if MustStartWithMethionine:
        MethionineFiltdf = BLASTFiltdf[BLASTFiltdf['tseq'].str.startswith('M')]
        MethionineFiltdf.to_csv(f'{CSVFileDirectory}/MethionineFiltdf.csv',index=False)

    if RemoveDuplicates:
        DuplicateRemoveddf = MethionineFiltdf.drop_duplicates(subset='tseq')
        DuplicateRemoveddf.to_csv(f'{CSVFileDirectory}/DuplicateRemoveddf.csv',index=False)