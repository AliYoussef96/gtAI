import ENc
from Bio import SeqIO
import pandas as pd 
from pandas import DataFrame
from pathlib import Path 
import math
import CA_RSCU  
import gtAI

def ENc_calc(fasta_file):
    #######################################################################
    #calculate the Enc values for all genes
    ###############################

    ENc.ENcfilename(fasta_file)

    Enc_file = pd.read_csv(fasta_file + ".enc", delimiter="\t")


    #get the 5% of lwo enc values ( high CUB )
    Enc_low_values = round(math.floor(len(Enc_file) * (5/100)),0)

    if Enc_low_values <= 0:
        Enc_low_values = 1
    else:
        pass
    Enc_low_values = Enc_file.nsmallest(Enc_low_values, 'eq5sun')

    return Enc_low_values


def ENc_calc_ref(fasta_file):
    #######################################################################
    #calculate the Enc values for all genes
    ###############################

    ENc.ENcfilename(fasta_file)

    Enc_file = pd.read_csv(fasta_file + ".enc", delimiter="\t")


    #get the 5% of lwo enc values ( high CUB )
    Enc_low_values = round(math.floor(len(Enc_file) * (100/100)),0)

    if Enc_low_values <= 0:
        Enc_low_values = 1
    else:
        pass
    Enc_low_values = Enc_file.nsmallest(Enc_low_values, 'eq5sun')

    return Enc_low_values
###########################################################################
###########################################################################

def RSCU_calc(fasta_file,Enc_low_values,genetic_code_number):  
    ###########################################################################
    # RSCU values for low and high ENC
    #############################

    #read fasta file
    dna_low_enc = ""
    dna_high_enc = ""
    for dna in SeqIO.parse(fasta_file, "fasta"):

        if str(dna.id) in list(Enc_low_values["id"]):
            dna_low_enc = dna_low_enc + str(dna.seq)



    RSCU_dataframe_low_enc = CA_RSCU.RSCU(dna_low_enc,"low_ENc_values",genetic_code_number)

    #join both dataframe

    return RSCU_dataframe_low_enc
###########################################################################
###########################################################################

def wi_tai_calc(dict_tGCN,Sug, Sci, Sai, Sgu, Sal, genetic_code_number=1,bacteria=False):
    ############################################# ##############################
    # calcularte the Wi wight for the codons and return the corr. between the high/low Enc with the Wi of each codon
    #################################

    #1. get the tRNA number

    #dict_tGCN = tai.get_trna_info(tRNA_url, url_start_from)

    dict1 = gtAI.dict_codon_anticodon(dict_tGCN)

    dict2 = gtAI.dict_codon_anticodon_count(dic_codon_anticodon = dict1, dict_tGCN_main = dict_tGCN, bacteria= bacteria)

    final_dict_wi = gtAI.abs_Wi(Sug=Sug,Sci=Sci,Sai=Sai,Sgu=Sgu,Sal=Sal,dict_anticodon_number=dict2,bacteria=bacteria)

    rel = gtAI.rel_Wi(final_dict_wi,genetic_code_number=genetic_code_number)
    new_rel = {}
    for i_codon in rel: # cuz RSCU value is T not U
        i_codon_new = i_codon.replace("U","T")
        new_rel[i_codon_new] = rel[i_codon]

    #convert the rel to dataframe
    df_rel = pd.DataFrame(list(new_rel.items()))
    df_rel.set_index(df_rel.columns[0],drop=True,inplace=True)
    df_rel.columns = ["Wi"]

    return df_rel
    ###########################################################################
    ###########################################################################


def corr_result(df_rel,RSCU_dataframe_low_enc):
    #merg the Wi.rel and the RSCU values
    RSCU_dataframe_low_enc = RSCU_dataframe_low_enc.merge(df_rel, how='outer', left_index=True, right_index=True)

    RSCU_wai_corr = RSCU_dataframe_low_enc.dropna(axis=0)

    return (RSCU_wai_corr.corr(method="spearman"))

     ###########################################################################
     ###########################################################################

  

