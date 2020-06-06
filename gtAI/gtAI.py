import itertools 
import pandas as pd
import urllib
import pandas as pd
from pandas import DataFrame
import numpy as np


def tRNADB_CE(url):
    """
    Get the tRNA genes count from tRNADB-CE database

    Args:

        url (string): a url to anticodon table for organism from tRNADB_CE database (http://trna.ie.niigata-u.ac.jp/cgi-bin/trnadb/index.cgi)

    Returns:

        A dictionary of each anticodon and its gene count
    
    Raises:

        ValueError if the URL, not a valid for tRNADB-CE database
    
    Example:

        > tRNADB_CE("trna.ie.niigata-u.ac.jp/cgi-bin/trnadb/whole_anticodon.cgi?GID=|CP001631&DTYPE=CMP&VTYPE=1")
        # Return an anticodon table of Acidimicrobium ferrooxidans DSM 10331

    """

    if "http://trna.ie.niigata-u.ac.jp" not in url:
        raise ValueError ("url is not includes in tRNADB-CE database")
    
    all_trna_df = DataFrame()

    read = pd.read_html(url)
    read1 = read[1].iloc[3:,3:5]
    read1.columns = ["anti_codon", "count"]

    read2 = read[1].iloc[3:,6:8]
    read2.columns = ["anti_codon", "count"]

    read3 = read[1].iloc[3:,9:11]
    read3.columns = ["anti_codon", "count"]

    read4 = read[1].iloc[3:,12:14]
    read4.columns = ["anti_codon", "count"]


    all_trna_df = all_trna_df.append(read1,ignore_index=True, sort=False)
    all_trna_df = all_trna_df.append(read2,ignore_index=True, sort=False)
    all_trna_df = all_trna_df.append(read3,ignore_index=True, sort=False)
    all_trna_df = all_trna_df.append(read4,ignore_index=True, sort=False)


    A_to_I_anticodon = ["AGA","AAG","AGG","ACG","AAU","AGU","AAC","AGC","ACC"]  #codons which will convert its A to I in Wobble pos.
    anti_codon_dict = {}
    for codon, count in zip( list(all_trna_df["anti_codon"]), list(all_trna_df["count"]) ):
        codon = codon.replace("T","U")
        if codon in A_to_I_anticodon:
            codon = list(codon)
            codon[0] = "I"
            codon = "".join(codon)
        else:
            pass
        anti_codon_dict[codon] = int(count)

    #remove anti_codons with A in the first postion # caz it not recognize ani codon
    for i_anticodon in list(anti_codon_dict.keys()):
        if i_anticodon[0] == "A":
            del anti_codon_dict[i_anticodon]

    if len(anti_codon_dict) == 0:
        raise ValueError ("No tables were found in this URL")
    else:
        return anti_codon_dict

########################################################
              #########################
########################################################

def GtRNAdb(url):
    """
        Get the tRNA genes count from GtRNAdb database

    Args:

        url (string): a url to anticodon table for organism from GtRNAdb database (http://gtrnadb.ucsc.edu/)

    Returns:

        A dictionary of each anticodon and its gene count
    
    Raises:

        ValueError if the URL, not a valid for GtRNAdb database

    Example:

        #example 1

        > GtRNAdb("http://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi19/")
      
    """

    amino_acid_list = ['Val', 'Ile', 'Leu', 'Glu', 'Gln', 'Asp', 'Asn', 'His', 'Trp', 'Phe', 'Tyr', 'Arg', 'Lys', 'Ser', 'Thr', 'Met', 'Ala', 'Gly', 'Pro', 'Cys']

    for start_from_table in range(1,4):

        if "gtrnadb.ucsc.edu" not in url:
            raise ValueError ("url is not includes in GtRNAdb database")

        tables = pd.read_html(url)
        anti_codon_dict = {}

        def SPLIT(x):
            return x.split(" ")

        A_to_I_anticodon = ["AGA","AAG","AGG","ACG","AAU","AGU","AAC","AGC","ACC"]  #codons which will convert its A to I in Wobble pos.
        cat_anti_codon = ()
        for i in range(start_from_table,len(tables)):
            tables_new = tables[i]
        
            for j in range(0, len(tables_new.index)):
                overview = tables_new.iloc[j,1:] 
                codon_anticodon = list ( map(str, overview) )    
                codon_anticodon = list ( map(SPLIT, codon_anticodon) ) 
                
                #anti_codon with A that will convert to I
                for i_codon in codon_anticodon:
                    if not i_codon[0].isdigit() and i_codon[0] != "&nbsp":
                        i_codon_rna = i_codon[0].replace("T","U")
                        if i_codon_rna in A_to_I_anticodon:
                            i_codon_rna = list(i_codon_rna)
                            i_codon_rna[0] = "I"
                            i_codon_rna = "".join(i_codon_rna)

                        try:
                            if i_codon_rna == "CAU" and "/" in i_codon[1]:
                                cat_anti_codon = cat_anti_codon + ( i_codon[1].split("/")[0] , )
                                cat_anti_codon = cat_anti_codon + ( i_codon[1].split("/")[1] , )
                            if i_codon_rna == "CAU" and "/" not in i_codon[1]:
                                cat_anti_codon = cat_anti_codon + ( i_codon[1] ,)
                            
                            
                            
                            anti_codon_dict[i_codon_rna] = int(i_codon[1])

                        except:
                            anti_codon_dict[i_codon_rna] = 0

        def seeifnumber(x):
            # take list and return a tuble with only number element
            result = ()
            for i in x:
                try:
                    int(i)
                    result = result + (int(i),)
                except:
                    pass
            return (result)


        try:
            anti_codon_dict["CAU"] = sum(seeifnumber(cat_anti_codon))
        except:
            pass

        dict_anti_codons = {}
        
        #remove anti_codons with A in the first postion # caz it not recognize ani codon
        for i_anticodon in list(anti_codon_dict.keys()):
            if i_anticodon[0] == "A":
                del anti_codon_dict[i_anticodon]

        if len(anti_codon_dict) == 0:
            raise ValueError ("No tables were found in this URL")
        else:
            see = sum([1 for i in amino_acid_list if i in anti_codon_dict])
            if see == 0:           
                 #if len(anti_codon_dict) == 58 or len(anti_codon_dict) == 57:
                return anti_codon_dict
                break
            else:
                pass

########################################################
              #########################
########################################################

def dict_codon_anticodon(anti_codon_dict):
    """
    Identify all potential anticodons for each codon 

    Args:

        anti_codon_dict (dict): a dictionary of all anticodons for an organism ( returned from ( tRNADB_CE() ) or ( GtRNAdb() ) )

    Returns:

        A dictionary of all potential anticodons for each codon 
    
    Raises:

        ValueError if the length of anti_codon_dict equal to zero
    
    Example:

        > anticodon_dict = tRNADB_CE("trna.ie.niigata-u.ac.jp/cgi-bin/trnadb/whole_anticodon.cgi?GID=|CP001631&DTYPE=CMP&VTYPE=1")
        # Return an anticodon table of Acidimicrobium ferrooxidans DSM 10331
        > dict_codon_anticodon ( anticanticodon_dictodon )

    """
    if len(anti_codon_dict) == 0:
        raise ValueError ("The dictionary can not be with length equal to zero")

    antiWC = { "A":"U", "U":"A", "C":"G", "G":"C" }
    antiWP = { "U":"G", "C":"I", "A":"I", "G":"U" }

    codon = ( "".join(i) for i in itertools.product('AUGC', repeat=3) )

    dic_codon_anticodon ={}

    #convert the codon to anticodon if anticodon in anti_codon_dict
    for i in codon:
        anti_codon_tuble = ( )

        if i[2] == "U":
            wc = antiWC[i[0]] + antiWC[i[1]] + "I"
            wc = wc[::-1]
        else:
            wc = antiWC[i[0]] + antiWC[i[1]] + antiWC[i[2]]
            wc = wc[::-1]
            
        if wc in anti_codon_dict:
            anti_codon_tuble = anti_codon_tuble + (wc, )
        else:
            pass
        #######

        wp = antiWC[i[0]] + antiWC[i[1]] + antiWP[i[2]] 
        wp = wp[::-1]

        if wp in anti_codon_dict:
            anti_codon_tuble = anti_codon_tuble + (wp, )
        else:
            pass

        if len(anti_codon_tuble) != 0: 
            dic_codon_anticodon[i] = ( anti_codon_tuble )
        else:
            pass

    return dic_codon_anticodon

########################################################
              #########################
########################################################

def dict_codon_anticodon_count(dic_codon_anticodon,dict_tGCN_main, bacteria = False):

    """ 
    Merge anticodon-codon dictionary with each anticodon tRNA gene copy number.

    Args:

        dic_codon_anticodon (dict): A dictionary of all potential anticodons for each codon returned from dict_codon_anticodon() function.

        dict_tGCN_main (dict): a dictionary of all anticodons for an organism ( returned from ( tRNADB_CE() function ) or ( GtRNAdb() function ) )

        bacteria (bool): True If the tested organism is prokaryotic or archaea, else equal to False ( default = False )

    Returns:
        
        a merged dictionary of anticodon-codon with each anticodon tRNA gene copy number.
    
    Raises:

        ValueError if the length of anti_codon_dict equal to zero
        ValueError if the length of dict_tGCN_main equal to zero

    Example:

        > anticodon_dict = tRNADB_CE("trna.ie.niigata-u.ac.jp/cgi-bin/trnadb/whole_anticodon.cgi?GID=|CP001631&DTYPE=CMP&VTYPE=1")
        # Return an anticodon table of Acidimicrobium ferrooxidans DSM 10331
        > anticodon_codon = dict_codon_anticodon ( anticanticodon_dictodon )
        > dict_codon_anticodon_count(anticodon_codon,anticodon_dict,bacteria = True)
    """

    if len(dic_codon_anticodon) == 0 or len(dict_tGCN_main) == 0:
        raise ValueError ("Both dic_codon_anticodon or dict_tGCN_main can not be with length equal to zero")

    dict_anticodon_number = {}
    for i in dic_codon_anticodon:
        temp = {}
        for j in range(len(dic_codon_anticodon[i])):
            
            
            temp[ dic_codon_anticodon[i][j] ]  = dict_tGCN_main[dic_codon_anticodon[i][j]] 

        dict_anticodon_number[i] = temp

    if bacteria:
        try:
            dict_anticodon_number["AUA"].update( { "CAU":dict_tGCN_main["CAU"] } )
        except:
            pass

    return dict_anticodon_number

########################################################
              #########################
########################################################
   
def abs_Wi(dict_anticodon_number, Sug,Sci,Sai,Sgu,Sal, bacteria=False):
    """
    Calculate the absolute adaptiveness values for each codon.

    Args:

        dict_anticodon_number (dict): a merged dictionary of anticodon-codon with each anticodon tRNA gene copy number returned from dict_codon_anticodon_count() function.

        Sug (int): the S-value for codon with (U) in the third position and (G) in first anticodon position

        Sci (int): the S-value for codon with (C) in the third position and (I) in first anticodon position

        Sai (int): the S-value for codon with (A) in the third position and (I) in first anticodon position

        Sgu (int): the S-value for codon with (G) in the third position and (U) in first anticodon position

        Sal (int): the S-value for codon with (A) in the third position and (L) in first anticodon position ( if bacteria = True)

        bacteria (bool): True If the tested organism is prokaryotic or archaea, else equal to False ( default = False )

    Returns:
        
        A dictionary of the absolute adaptiveness values for each codon

    Note:

        All Sij values (as Sug) should be a number from 0 to 1

    Raises:

        ValueError if the length of dict_anticodon_number equal to zero
        ValueError if any of Sij values are less than 0 or greater than 1

    Example:

        > anticodon_dict = tRNADB_CE("trna.ie.niigata-u.ac.jp/cgi-bin/trnadb/whole_anticodon.cgi?GID=|CP001631&DTYPE=CMP&VTYPE=1")
        # Return an anticodon table of Acidimicrobium ferrooxidans DSM 10331
        > anticodon_codon = dict_codon_anticodon ( anticanticodon_dictodon )
        > dict_codon_anticodon_count = dict_codon_anticodon_count(anticodon_codon,anticodon_dict,bacteria = True)
        > abs_Wi(dict_codon_anticodon_count, Sug=1, Sci=1, Sai=1, Sgu=1, Sal=1, bacteria=True) 
    """

    if len(dict_anticodon_number) == 0:
        raise ValueError ("dict_anticodon_number can not be with length equal to zero")

    if  0 > Sug > 1 or 0 > Sci > 1 or 0 > Sai > 1 or 0 > Sgu > 1 :
        raise ValueError ("All Sij values (as Sug) should be a number from 0 to 1") 

    if bacteria:
        if 0 > Sal > 1:
            raise ValueError ("All Sij values (as Sug) should be a number from 0 to 1") 


    # Sij = Scodon:anticodon
    # Sij of WC interaction is equal to 0
    Sij = { "Sui":0, "Scg":0, "Sau":0, "Sgc":0, 
            "Sug":Sug, "Sci":Sci, "Sai":Sai, 
            "Sgu":Sgu }

    if bacteria: #Lysidine (L) is a bacterial RNA modification of the DNA nucleotide cytidine interact with (A)
        try:
            if "CAU" in dict_anticodon_number["AUA"]:
                Sij["Sal"] = Sal
            else:
                dict_anticodon_number["AUA"].update( { "CAU":0 } )
        except:
            pass

    Wi_dict = {}
    for i_codon in dict_anticodon_number:
        if i_codon[2] == "U":

            #get the anticodon with I in first pos. for the codon with U in first pos.
            key_I = [j for j in dict_anticodon_number[i_codon] if j[0] == "I"]
            key_I = "".join(key_I)

            #get the anticodon with G in first pos. for the codon with U in first pos. (wobble pos)
            key_G = [j for j in dict_anticodon_number[i_codon] if j[0] == "G"]
            key_G = "".join(key_G)

            if key_I != "" and key_G != "":
                Wu = ( ( 1 - Sij["Sui"] ) * dict_anticodon_number[i_codon][key_I] ) +  ( ( 1 - Sij["Sug"] ) *  dict_anticodon_number[i_codon][key_G] ) 
            elif key_I != "" and key_G == "":
                Wu = ( ( 1 - Sij["Sui"] ) * dict_anticodon_number[i_codon][key_I] ) 
            elif key_I == "" and key_G != "":
                Wu = ( ( 1 - Sij["Sug"] ) *  dict_anticodon_number[i_codon][key_G] ) 
            else:
                Wu = 0

            Wi_dict[i_codon] = Wu

        elif i_codon[2] == "C":

            #get the anticodon with G in first pos. for the codon with C in first pos.
            key_G = [j for j in dict_anticodon_number[i_codon] if j[0] == "G"]
            key_G = "".join(key_G)
            #get the anticodon with I in first pos. for the codon with C in first pos. (wobble pos)
            key_I = [j for j in dict_anticodon_number[i_codon] if j[0] == "I"]
            key_I = "".join(key_I)

            if key_I != "": # key_G must be found
                Wc = ( ( 1 - Sij["Scg"] ) * dict_anticodon_number[i_codon][key_G] ) +  ( ( 1 - Sij["Sci"] ) *  dict_anticodon_number[i_codon][key_I] ) 
            else: #if key)I == ""
                Wc = ( ( 1 - Sij["Scg"] ) * dict_anticodon_number[i_codon][key_G] ) 

            Wi_dict[i_codon] = Wc
            
        elif i_codon[2] == "A":

            if i_codon == "AUA":

                if bacteria: 

                    #get the anticodon with U in first pos. for the codon with A in first pos.
                    key_U = [j for j in dict_anticodon_number[i_codon] if j[0] == "U"]
                    key_U1 = "".join(key_U[0])
                    try:
                        key_U2 = "".join(key_U[1])
                    except:
                        key_U2 = False

                    #get the anticodon with I in first pos. for the codon with A in first pos. (wobble pos)
                    key_I = [j for j in dict_anticodon_number[i_codon] if j[0] == "I"]   
                    key_I = "".join(key_I)


                    if key_U2 == True and key_I != "" :
                        Wa = ( ( 1 - Sij["Sau"] ) * dict_anticodon_number[i_codon][key_U1] ) +  ( ( 1 - Sij["Sai"] ) *  dict_anticodon_number[i_codon][key_I] ) + ( ( 1 - Sij["Sal"] ) *  dict_anticodon_number[i_codon][key_U2] )
                    elif key_U2 == False and key_I != "" :
                        Wa = ( ( 1 - Sij["Sau"] ) * dict_anticodon_number[i_codon][key_U1] ) +  ( ( 1 - Sij["Sai"] ) *  dict_anticodon_number[i_codon][key_I] ) 
                    elif key_U2 == True and key_I == "":
                        Wa = ( ( 1 - Sij["Sau"] ) * dict_anticodon_number[i_codon][key_U1] ) + ( ( 1 - Sij["Sal"] ) *  dict_anticodon_number[i_codon][key_U2] )
                    elif key_U2 == False and key_I == "":
                        Wa = ( ( 1 - Sij["Sau"] ) * dict_anticodon_number[i_codon][key_U1] ) 
                
                    Wi_dict[i_codon] = Wa
               

                if not bacteria:
                    #get the anticodon with U in first pos. for the codon with A in first pos.
                    key_U = [j for j in dict_anticodon_number[i_codon] if j[0] == "U"]
                    key_U = "".join(key_U)
                    #get the anticodon with I in first pos. for the codon with A in first pos. (wobble pos)
                    key_I = [j for j in dict_anticodon_number[i_codon] if j[0] == "I"]
                    key_I = "".join(key_I)
                    
                    if key_I != "":
                        Wa = ( ( 1 - Sij["Sau"] ) * dict_anticodon_number[i_codon][key_U] ) +  ( ( 1 - Sij["Sai"] ) *  dict_anticodon_number[i_codon][key_I] ) 
                    else:
                        Wa = ( ( 1 - Sij["Sau"] ) * dict_anticodon_number[i_codon][key_U] ) 

                    Wi_dict[i_codon] = Wa

            else:
                #get the anticodon with U in first pos. for the codon with A in first pos.
                key_U = [j for j in dict_anticodon_number[i_codon] if j[0] == "U"]
                key_U = "".join(key_U)
                #get the anticodon with I in first pos. for the codon with A in first pos. (wobble pos)
                key_I = [j for j in dict_anticodon_number[i_codon] if j[0] == "I"]
                key_I = "".join(key_I)
                
                if key_I != "":
                    Wa = ( ( 1 - Sij["Sau"] ) * dict_anticodon_number[i_codon][key_U] ) +  ( ( 1 - Sij["Sai"] ) *  dict_anticodon_number[i_codon][key_I] ) 
                else:
                    Wa = ( ( 1 - Sij["Sau"] ) * dict_anticodon_number[i_codon][key_U] ) 

                Wi_dict[i_codon] = Wa

            
        elif i_codon[2] == "G":

            #get the anticodon with C in first pos. for the codon with G in first pos.
            key_C = [j for j in dict_anticodon_number[i_codon] if j[0] == "C"]
            key_C = "".join(key_C)
            #get the anticodon with U in first pos. for the codon with G in first pos. (wobble pos)
            key_U = [j for j in dict_anticodon_number[i_codon] if j[0] == "U"]
            key_U = "".join(key_U)

            if key_U != "":
                Wg = ( ( 1 - Sij["Sgc"] ) * dict_anticodon_number[i_codon][key_C] ) +  ( ( 1 - Sij["Sgu"] ) *  dict_anticodon_number[i_codon][key_U] ) 
            else:
                Wg = ( ( 1 - Sij["Sgc"] ) * dict_anticodon_number[i_codon][key_C] ) 

            Wi_dict[i_codon] = Wg

    return Wi_dict

########################################################
              #########################
########################################################

def rel_Wi(dict_abs_Wi,genetic_code_number=1):
    """
    Calculate the relative adaptiveness values for each codon.

    Args:

        dict_abs_Wi (dict): dictionary of the absolute adaptiveness values for each codon returned from abs_Wi() function

        genetic_code_number (int): default = 1, The Genetic Codes number described by NCBI (https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)

    Returns:
        
        a dictionary of the relative adaptiveness values for each codon

    Raises:

        ValueError if the length of dict_abs_Wi equal to zero

    Example:

        > anticodon_dict = tRNADB_CE("trna.ie.niigata-u.ac.jp/cgi-bin/trnadb/whole_anticodon.cgi?GID=|CP001631&DTYPE=CMP&VTYPE=1")
        # Return an anticodon table of Acidimicrobium ferrooxidans DSM 10331
        > anticodon_codon = dict_codon_anticodon ( anticanticodon_dictodon )
        > dict_codon_anticodon_count = dict_codon_anticodon_count(anticodon_codon,anticodon_dict,bacteria = True)
        > abs_Wi = abs_Wi(dict_codon_anticodon_count, Sug=1, Sci=1, Sai=1, Sgu=1, Sal=1, bacteria=True) 
        > rel_Wi(abs_Wi, 11) 
    """

    if len(dict_abs_Wi) == 0:
        raise ValueError ("dict_abs_Wi can not be with length equal to zero")



    codons = ["".join(i) for i in itertools.product('AUGC', repeat=3)]

    dict_rel_Wi = {}

    # rel_wi if Wi != 0
    max_Wi = max (list(dict_abs_Wi.values()))
    for i_codon in codons:
        if i_codon in dict_abs_Wi and dict_abs_Wi[i_codon] != 0:
            dict_rel_Wi[i_codon] = dict_abs_Wi[i_codon] / max_Wi


    #if stop codons or  remove it
    def get_stop_codons(genetic_code_number):
        from Bio.Data import CodonTable
        genetic_table = CodonTable.unambiguous_dna_by_id[genetic_code_number]
        stop_codons = genetic_table.stop_codons
        stop_M_codons = [i.replace("T","U") for i in stop_codons]
        return stop_M_codons
    

    #remove stop codon from the analysis brefore calc. the geo_mean
    for i_codons_remove in get_stop_codons(genetic_code_number):
        if i_codons_remove in dict_rel_Wi:
            del dict_rel_Wi[i_codons_remove]
        else:
            pass

    #rel_wi if Wi == 0
    def geo_mean(iterable):
    
        #function to calculate the geometric mean
        #iterable = list of number
    
        import numpy as np
        a = np.log(iterable)
        return np.exp(a.sum()/len(a))

    #add the geo_mean if Wi == 0
    #calc geo_mean for rel_wi with Wi != 0
    wi_not_zero = [i for i in list(dict_rel_Wi.values()) if i != 0]
    geo_rel_wi = geo_mean(wi_not_zero)
    for i_codon in codons:
        if i_codon not in dict_rel_Wi or dict_rel_Wi[i_codon] == 0 :
            dict_rel_Wi[i_codon] = geo_rel_wi


    #remove stop codon  from the analysis after calculting the geo_mean
    for i_codons_remove in get_stop_codons (genetic_code_number):
        if i_codons_remove in dict_rel_Wi:
            del dict_rel_Wi[i_codons_remove] 
        else:
            pass

            
    return dict_rel_Wi

########################################################
              #########################
########################################################

def calc_Tai(DNA,rel_dict_wi,genetic_code_number=1):
    """
    Calculate the tRNA adaptation index of a gene.

    Args:

        DNA (str): a coding sequence of DNA ( should only contain A, C, T, and G )

        rel_dict_wi (dict): dictionary of the relative adaptiveness values for each codon returned from rel_Wi() function

        genetic_code_number (int): default = 1, The Genetic Codes number described by NCBI (https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)

    Returns:
        
        the tRNA adaptation index of a gene

    Raises:

        ValueError if the length of dict_anticodon_number equal to zero

    Example:

        > anticodon_dict = tRNADB_CE("trna.ie.niigata-u.ac.jp/cgi-bin/trnadb/whole_anticodon.cgi?GID=|CP001631&DTYPE=CMP&VTYPE=1")
        # Return an anticodon table of Acidimicrobium ferrooxidans DSM 10331
        > anticodon_codon = dict_codon_anticodon ( anticanticodon_dictodon )
        > dict_codon_anticodon_count = dict_codon_anticodon_count(anticodon_codon,anticodon_dict,bacteria = True)
        > abs_Wi = abs_Wi(dict_codon_anticodon_count, Sug=1, Sci=1, Sai=1, Sgu=1, Sal=1, bacteria=True) 
        > rel_Wi(abs_Wi, 11) 
    """

    if len(rel_dict_wi) == 0:
        raise ValueError ("rel_dict_wi can not be with length equal to zero")

    if DNA.count("A") + DNA.count("T") + DNA.count("C") + DNA.count("G") != len(DNA):
        raise ValueError("The DNA can not include any nucleotides expect (A, T, C, and G)")

    codon_DNA = [ DNA[i:i+3].replace("T","U") for i in range(0,len(DNA),3) if len( DNA[i:i+3] ) == 3  ]
    
    if codon_DNA[0] == "AUG": #start codon remove it
        codon_DNA = codon_DNA[1:]
    else:
        pass

    #if stop codons  remove it
    def get_stop_codons(genetic_code_number):
        from Bio.Data import CodonTable
        genetic_table = CodonTable.unambiguous_dna_by_id[genetic_code_number]
        stop_codons = genetic_table.stop_codons
        stop_M_codons = [i.replace("T","U") for i in stop_codons]
        return stop_M_codons


    stop_codons  = get_stop_codons(genetic_code_number)
    codon_wi = [ rel_dict_wi[i_codon] for i_codon in codon_DNA if i_codon not in stop_codons   ]

    def geo_mean(iterable):
        
        #function to calculate the geometric mean
        #iterable = list of number
        
        import numpy as np
        a = np.log(iterable)
        return np.exp(a.sum()/len(a))
    if sum(codon_wi) != 0:
        Tai_result = geo_mean(codon_wi)
    else:
        Tai_result = 0

    return Tai_result
