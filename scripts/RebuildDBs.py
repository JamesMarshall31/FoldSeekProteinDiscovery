def RebuildDBs(): # defunct holding function until permanent solution to database building is found

    ## At some point BLAST databases may need to be rebuilt - this will require some further update and probably be moved to a separate function
    
    RebuildRequired = input("Will BLAST databases need to be rebuilt - respond Yes if so\n")
    if RebuildRequired == 'Yes':
        RebuildDBs = True # If this is true, an if statement will be entered where the BLASTable database is built up from a single fasta file
    else:
        RebuildDBs = False
    DBsToRebuild = ['BPPRC','AllergenOnline','DRAVPPatents','DRAVPProteins'] # What DBs need to be rebuilt in such a way - maybe hardcoding not great 
    if RebuildDBs:
        for DB in DBsToRebuild:
            PathToDB = f'BLASTDatabases/{DB}_BLAST_DB'
            DBFASTA = glob.glob(f'{PathToDB}/*.fa')[0] # finds the fasta in the directory - obviously make sure its there!!!
            build_blast_database(DBFASTA, dbtype='prot') # Builds out the necessary files for BLASTing