def RSCU(allseq,allseq_name,The_Genetic_Codes_number):
    """calculate RSCU values

    Args:
        
        allseq (str): DNA sequence
        allseq_name (str) : gene name
        The_Genetic_Codes_number (int) : default = 1, The Genetic Codes number described by NCBI (https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)

        
    Returns:
        DataFrame: DataFrame contains codons and RSCU values
    """
    from Bio.Data import CodonTable
    from Bio.Seq import Seq
    import re
    from Bio.Alphabet import generic_dna
    import pandas as pd
    from pandas import DataFrame
    from itertools import tee

    xcodontable = CodonTable.unambiguous_dna_by_id[1]
    ycodontable = xcodontable.forward_table
    zcodontable = [ycodontable[i] for i in ycodontable]
    qcodontable = [i for i in ycodontable ]


    for i in zcodontable:
        if zcodontable.count(i) == 1:
            zcodontable.remove(i)
    RSCU = {}

    sequ = str(allseq)

    allseqstr, allseqstr_1  = tee(sequ[i: i+3] for i in range(0, len(sequ), 3) if len(sequ[i: i+3]) == 3)

    qcodontable = ( i for i in qcodontable)
    dic2 = {}
    allseqstr2 = Seq('', generic_dna)
    for i in allseqstr:
        allseqstr2 += i
    aminoacid2 = allseqstr2.translate(table = The_Genetic_Codes_number , stop_symbol ='')
    aminoacid2 = str(aminoacid2)
    RSCUall = {}

    for ii in allseqstr_1:
        dic2[ii] = dic2.get(ii,0) + 1
    for k in  qcodontable:
        RSCUall[k] = 0
        if k in dic2:
            try:
                rscu2 = dic2[k] / ((1/ zcodontable.count(ycodontable[k]))*(aminoacid2.count(ycodontable[k])))
                RSCUall[k] = round(rscu2,6)

            except ZeroDivisionError:
                pass

    df = pd.DataFrame(index=pd.Series([i for i in RSCUall]))
    df[allseq_name] = [RSCUall[r] for r in RSCUall ]
    df.drop(['ATG'],0,inplace= True)
    df.sort_index(inplace=True)
    return df
