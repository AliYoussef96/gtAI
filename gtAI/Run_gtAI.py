#Note it will not take genes with N,R,Y,M,k,S,W in the analysis
#################################################################

from gtAI import new_flow
from gtAI import bygaft
from Bio import SeqIO
from gtAI import gtAI
import pandas as pd 
from pandas import DataFrame
#from scipy import stats

def gtai_analysis(main_fasta, GtRNA, ref_fasta = "", genetic_code_number=11, size_pop=60, generation_number=100, bacteria=True):
    """

    Args:

        main_fasta (str): A main fasta file contains genes want to be analyzed (CDS)
        GtRNA (dict): the tRNA genes count
        ref_fasta (str): A reference genes with the highest gene expression in a genome (CDS)
        genetic_code_number (int): default = 1, The Genetic Codes number described by NCBI (https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)
        size_pop (int): A parameter for the genetic algorithm to identify the population size containing the possible solutions to optimize Sij-values
        generation_number (int): A parameter for the genetic algorithm to identify the generation number
        bacteria (bool): True If the tested organism is prokaryotic or archaea, else equal to False ( default = False )

    Returns:

         df_tai (dataframe): Contains each gene id and its gtAI value 
         final_dict_wi (dict): Contains each codon and its absolute adaptiveness value
         rel_values (dict): Contains each codon and its relative adaptiveness values

    Note: 
        All genes expected to be CDS.
        This function will generate Two files one for ENc values ( a CSV file ) and the other one named best_fit (Contains the optimization information for Sij-values)

    """



    if ref_fasta == "":
        # 1- Calculate the ENc values, and extract the reference set of genes.
        Enc_low_values = new_flow.ENc_calc(main_fasta)

        # 2- Calculate the RSCU values for the reference set.
        RSCU_values = new_flow.RSCU_calc(main_fasta,Enc_low_values,genetic_code_number)
    
    elif ref_fasta != "":
         # 1- Calculate the ENc values, and extract the reference set of genes.
        Enc_low_values = new_flow.ENc_calc_ref(ref_fasta)

        # 2- Calculate the RSCU values for the reference set.
        RSCU_values = new_flow.RSCU_calc(ref_fasta,Enc_low_values,genetic_code_number)
    
    #3- Optimize Sij-values.
    bygaft.gene_algo_corr(dict_tGCN=GtRNA , genetic_code_number = genetic_code_number, RSCU_df=RSCU_values, size_pop=size_pop, generation_number=generation_number, bacteria=bacteria )


    #thie function will return a file named best_fit.py ( must be in the same path )
    while True:
        try:
            import best_fit; best_fit_ = best_fit.best_fit
            steps, variants, fits = list(zip(*best_fit_))
            Sij = variants[-1]
            Sug = Sij[0]
            Sci = Sij[1]
            Sai = Sij[2]
            Sgu = Sij[3]
            Sal = Sij[4]
            if best_fit:
                break
        except:
            continue
    

    # 4- get the tRNA gene copy number ( from GtRNAdb database for example ).
    dict1 = gtAI.dict_codon_anticodon(GtRNA)

    dict2 = gtAI.dict_codon_anticodon_count(dic_codon_anticodon = dict1, dict_tGCN_main = GtRNA, bacteria= bacteria)

    # 5- Calculate the absolute adaptiveness values.
    final_dict_wi = gtAI.abs_Wi(Sug=Sug,Sci=Sci,Sai=Sai,Sgu=Sgu,Sal=Sal,dict_anticodon_number=dict2,bacteria=bacteria)

    # 6- Calculate the relative adaptiveness values.
    rel_values = gtAI.rel_Wi(final_dict_wi,genetic_code_number=genetic_code_number)

    ##########################################################################################################

    
    #7- Calculate the gtAI for all genes.
    tai_data = {}
    for dna in SeqIO.parse(main_fasta, "fasta"):
        dna.seq = dna.seq.upper()
        dna_seq = str(dna.seq)
        if "W" not in dna.seq and "S" not in dna.seq and "N" not in dna.seq and "R" not in dna.seq and "Y" not in dna.seq and "M" not in dna.seq and "K" not in dna.seq:
            TAi = gtAI.calc_Tai(dna_seq, rel_values, genetic_code_number=genetic_code_number)
            tai_data[dna.description] = TAi
        else:
            pass

    df_tai = pd.DataFrame(list(tai_data.items()))
    df_tai.columns = ["id","tai"]

    return df_tai , final_dict_wi, rel_values

