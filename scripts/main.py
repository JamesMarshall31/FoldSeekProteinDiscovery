# Example script: main_script.py
from app_config import AppConfig
from pipeline_functions import run_foldseek, run_blasts, apply_filters



def main():
    
    UserChoice = input("""Which of the below options would you like to run? Please give your answer as a single number.
    1.) Run FoldSeek on starting proteins.    Input: .pdb files -> Output: Alignments, FirstFoldSeekdf.csv, (HTMLs)
    2.) Run BLASTs on the hits returned by FoldSeek.    Input: FirstFoldSeekdf.csv -> Output: BLAST .xml files, BLASTedFoldSeekdf.csv
    3.) Filter based on structural or BLAST alignments.    Input: BLASTedFoldSeekdf.csv -> Output: .csv files
    4.) All of the above\n""")

    UserChoice = int(UserChoice)
    assert 1 <= UserChoice <= 4, "Number must be between 1 and 4 inclusive."


    CurrentSearchDirectory = input("Enter The path to the directory of the current search e.g FoldSeekSearches/Lead_To_Gold_Proteins_JM_20240304 \nRemember, this directory should already contain a folder called ProteinQueries, which is already populated with at least 1 pdb file! \n")
    
    config = AppConfig(UserChoice, CurrentSearchDirectory) # AppConfig is a class not a function, script or directory, hence the camel case
    
    if config.user_choice in [1, 4]:
        run_foldseek(config)
    if config.user_choice in [2, 4]:
        run_blasts(config)
    if config.user_choice in [3, 4]:
        apply_filters(config)
    
    
 
   

if __name__ == "__main__":
    from functions import *
    main()