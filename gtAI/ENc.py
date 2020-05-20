def ENcfilename(file):
    import sys
    import argparse
    from Bio import SeqIO
    from pathlib import Path

   
    file_out_put = file + ".enc"



    parser = argparse.ArgumentParser(usage='%(prog)s [options]',
                                     description='Extracts codonusage from CDS input FASTA file. Output will be raw codon counts (.codoncnt), global ACTG counts (.actgcnt), first (.firstcnt), second (.secondcnt), third (.third) codon position counts and Relative Synonymous Codon Usage (.rscucnt). Optional different methods can be applied to calculate Effective Number of Codons (.enc).')
    parser.add_argument("-v", "--verbose", help="increase output verbosity", action="store_true")
    parser.add_argument('-i', help='specify CDS input file in FASTA format' , default = None)
    parser.add_argument('-o', help='specify output prefix',default = None)
    parser.add_argument('-r', action='store_true',
                        help='specify if CDS sequences with length modulo 3 unequal to 0 should be removed and reported to std.out')
    parser.add_argument('-enc', choices=['eq4Wright', 'eq2Sun', 'eq5Sun', 'all'], default = "eq5Sun",
                        help='specify equation to calculate ENC. Either equation (4) [eq4Wright] of (Wright. (1990) Gene 87:23-29) or equation (2) [eq2Sun] or equation (5) [eq5Sun] of (Sun et al. (2012) Mol. Biol. Evol. 30:191-196) or [all].')
    parser.add_argument('-six2fourtwo', default='False', choices=['True', 'False'],
                        help='specify if sixfold codons should be grouped into one fourfold and one twofold group [default: False]. This will only affect calculation of ENC values.')
    args = parser.parse_args()

############
    infile = file
############

#############
    outfile_enc = file_out_put 
############

    six2fourtwo = 'False'
    if args.six2fourtwo == 'True':
        six2fourtwo = True
    if args.six2fourtwo == 'False':
        six2fourtwo = False
    codontable = codontable1
    if six2fourtwo:
        codontable = codontable2

    print("read fasta")
    original_fasta = SeqIO.parse(infile, "fasta")
    print("extract codon counts")
    global_codons = codontable()
    global_actg = actgtable()
    global_first = actgtable()
    global_second = actgtable()
    global_third = actgtable()
    global_rscu = codontable()
    ids_mo3 = []

    if args.enc is not None:
        enchandle = open(outfile_enc, "w")
        if args.enc == 'eq4Wright':
            enchandle.write("id\tlen\tmo3\teq4Wright\n")
        if args.enc == 'eq2Sun':
            enchandle.write("id\tlen\tmo3\teq2Sun\n")
        if args.enc == 'eq5Sun':
            enchandle.write("id\tlen\tmo3\teq5sun\n")
        if args.enc == 'all':
            enchandle.write("id\tlen\tmo3\teq4Wright\teq2Sun\teq5Sun\n")

        c = 0
        cmo3 = 0
        for record in original_fasta:
            c += 1
            tmp_counts = codontable()
            tmp_actgcounts = actgtable()
            tmp_first = actgtable()
            tmp_second = actgtable()
            tmp_third = actgtable()
            tmp_id = record.id
            tmp_len = len(record)
            tmp_mo3 = 0
            if len(record) % 3 != 0:
                tmp_mo3 = 1
                cmo3 += 1
                ids_mo3.append(tmp_id)
                if args.r:
                    continue
            for i in range(0, len(record), 3):
                codon = str(record[i:i + 3].seq)
                for nuc in codon:
                    if nuc in global_actg:
                        global_actg[nuc][1] += 1
                        tmp_actgcounts[nuc][1] += 1
                if len(codon) != 3:
                    continue
                if codon in global_codons:
                    global_codons[codon][2] += 1
                    tmp_counts[codon][2] += 1
                    global_first[codon[0]][1] += 1
                    tmp_first[codon[0]][1] += 1
                    global_second[codon[1]][1] += 1
                    tmp_second[codon[1]][1] += 1
                    global_third[codon[2]][1] += 1
                    tmp_third[codon[2]][1] += 1
                if codon not in global_codons:
                    global_codons['XXX'][2] += 1
                    tmp_counts['XXX'][2] += 1
                    global_first['X'][1] += 1
                    tmp_first['X'][1] += 1
                    global_second['X'][1] += 1
                    tmp_second['X'][1] += 1
                    global_third['X'][1] += 1
                    tmp_third['X'][1] += 1
            if args.enc is not None:
                tmp_enc = {}
                tmp_gcbypos = gcbypos(tmp_counts, six2fourtwo)
                if args.enc == 'eq4Wright':
                    enchandle.write(tmp_id + "\t" + str(tmp_len) + "\t" + str(tmp_mo3) + "\t" + str(
                        calc_eq4wright(tmp_gcbypos[2])) + "\n")
                if args.enc == 'eq2Sun':
                    enchandle.write(tmp_id + "\t" + str(tmp_len) + "\t" + str(tmp_mo3) + "\t" + str(
                        calc_eq2sun(tmp_counts, six2fourtwo)) + "\n")
                if args.enc == 'eq5Sun':
                    enchandle.write(tmp_id + "\t" + str(tmp_len) + "\t" + str(tmp_mo3) + "\t" + str(
                        calc_eq5sun(tmp_counts, six2fourtwo)) + "\n")
                if args.enc == 'all':
                    enchandle.write(tmp_id + "\t" + str(tmp_len) + "\t" + str(tmp_mo3) + "\t" + str(
                        calc_eq4wright(tmp_gcbypos[2])) + "\t" + str(
                        calc_eq2sun(tmp_counts, six2fourtwo)) + "\t" + str(
                        calc_eq5sun(tmp_counts, six2fourtwo)) + "\n")

    print("finished writing")
    if args.r:
        print ('\n'.join(ids_mo3))


def inversetable():
    inversetablecount = {
        'A': {'GCT': 0, 'GCC': 0, 'GCA': 0, 'GCG': 0},
        'R': {'CGT': 0, 'CGC': 0, 'CGA': 0, 'CGG': 0, 'AGA': 0, 'AGG': 0},
        'N': {'AAT': 0, 'AAC': 0},
        'D': {'GAT': 0, 'GAC': 0},
        'C': {'TGT': 0, 'TGC': 0},
        'Q': {'CAA': 0, 'CAG': 0},
        'E': {'GAA': 0, 'GAG': 0},
        'G': {'GGT': 0, 'GGC': 0, 'GGA': 0, 'GGG': 0},
        'H': {'CAT': 0, 'CAC': 0},
        'I': {'ATT': 0, 'ATC': 0, 'ATA': 0},
        'M': {'ATG': 0},
        'L': {'TTA': 0, 'TTG': 0, 'CTT': 0, 'CTC': 0, 'CTA': 0, 'CTG': 0},
        'K': {'AAA': 0, 'AAG': 0},
        'F': {'TTT': 0, 'TTC': 0},
        'P': {'CCT': 0, 'CCC': 0, 'CCA': 0, 'CCG': 0},
        'S': {'TCT': 0, 'TCC': 0, 'TCA': 0, 'TCG': 0, 'AGT': 0, 'AGC': 0},
        'T': {'ACT': 0, 'ACC': 0, 'ACA': 0, 'ACG': 0},
        'W': {'TGG': 0},
        'Y': {'TAT': 0, 'TAC': 0},
        'V': {'GTT': 0, 'GTC': 0, 'GTA': 0, 'GTG': 0},
        '*': {'TAA': 0, 'TGA': 0, 'TAG': 0}
    }
    return inversetablecount


def codontable1():
    codontablecount = {
        'GCT': ['A', 'four', 0],
        'GCC': ['A', 'four', 0],
        'GCA': ['A', 'four', 0],
        'GCG': ['A', 'four', 0],
        'CGT': ['R', 'six', 0],
        'CGC': ['R', 'six', 0],
        'CGG': ['R', 'six', 0],
        'CGA': ['R', 'six', 0],
        'AGA': ['R', 'six', 0],
        'AGG': ['R', 'six', 0],
        'AAT': ['N', 'two', 0],
        'AAC': ['N', 'two', 0],
        'GAT': ['D', 'two', 0],
        'GAC': ['D', 'two', 0],
        'TGT': ['C', 'two', 0],
        'TGC': ['C', 'two', 0],
        'CAA': ['Q', 'two', 0],
        'CAG': ['Q', 'two', 0],
        'GAA': ['E', 'two', 0],
        'GAG': ['E', 'two', 0],
        'GGT': ['G', 'four', 0],
        'GGC': ['G', 'four', 0],
        'GGA': ['G', 'four', 0],
        'GGG': ['G', 'four', 0],
        'CAT': ['H', 'two', 0],
        'CAC': ['H', 'two', 0],
        'ATT': ['I', 'three', 0],
        'ATC': ['I', 'three', 0],
        'ATA': ['I', 'three', 0],
        'ATG': ['M', 'one', 0],
        'TTA': ['L', 'six', 0],
        'TTG': ['L', 'six', 0],
        'CTT': ['L', 'six', 0],
        'CTC': ['L', 'six', 0],
        'CTA': ['L', 'six', 0],
        'CTG': ['L', 'six', 0],
        'AAA': ['K', 'two', 0],
        'AAG': ['K', 'two', 0],
        'TTT': ['F', 'two', 0],
        'TTC': ['F', 'two', 0],
        'CCT': ['P', 'four', 0],
        'CCC': ['P', 'four', 0],
        'CCA': ['P', 'four', 0],
        'CCG': ['P', 'four', 0],
        'TCT': ['S', 'six', 0],
        'TCC': ['S', 'six', 0],
        'TCA': ['S', 'six', 0],
        'TCG': ['S', 'six', 0],
        'AGT': ['S', 'six', 0],
        'AGC': ['S', 'six', 0],
        'ACT': ['T', 'four', 0],
        'ACC': ['T', 'four', 0],
        'ACA': ['T', 'four', 0],
        'ACG': ['T', 'four', 0],
        'TGG': ['W', 'one', 0],
        'TAT': ['Y', 'two', 0],
        'TAC': ['Y', 'two', 0],
        'GTT': ['V', 'four', 0],
        'GTC': ['V', 'four', 0],
        'GTA': ['V', 'four', 0],
        'GTG': ['V', 'four', 0],
        'TAA': ['*', 'three', 0],
        'TGA': ['*', 'three', 0],
        'TAG': ['*', 'three', 0],
        'XXX': ['_missing', 'none', 0]
    }
    return codontablecount


def codontable2():
    codontablecount = {
        'GCT': ['A', 'four', 0],
        'GCC': ['A', 'four', 0],
        'GCA': ['A', 'four', 0],
        'GCG': ['A', 'four', 0],
        'CGT': ['R', 'four', 0],
        'CGC': ['R', 'four', 0],
        'CGG': ['R', 'four', 0],
        'CGA': ['R', 'dour', 0],
        'AGA': ['R_', 'two', 0],
        'AGG': ['R_', 'two', 0],
        'AAT': ['N', 'two', 0],
        'AAC': ['N', 'two', 0],
        'GAT': ['D', 'two', 0],
        'GAC': ['D', 'two', 0],
        'TGT': ['C', 'two', 0],
        'TGC': ['C', 'two', 0],
        'CAA': ['Q', 'two', 0],
        'CAG': ['Q', 'two', 0],
        'GAA': ['E', 'two', 0],
        'GAG': ['E', 'two', 0],
        'GGT': ['G', 'four', 0],
        'GGC': ['G', 'four', 0],
        'GGA': ['G', 'four', 0],
        'GGG': ['G', 'four', 0],
        'CAT': ['H', 'two', 0],
        'CAC': ['H', 'two', 0],
        'ATT': ['I', 'three', 0],
        'ATC': ['I', 'three', 0],
        'ATA': ['I', 'three', 0],
        'ATG': ['M', 'one', 0],
        'TTA': ['L_', 'two', 0],
        'TTG': ['L_', 'two', 0],
        'CTT': ['L', 'four', 0],
        'CTC': ['L', 'four', 0],
        'CTA': ['L', 'four', 0],
        'CTG': ['L', 'four', 0],
        'AAA': ['K', 'two', 0],
        'AAG': ['K', 'two', 0],
        'TTT': ['F', 'two', 0],
        'TTC': ['F', 'two', 0],
        'CCT': ['P', 'four', 0],
        'CCC': ['P', 'four', 0],
        'CCA': ['P', 'four', 0],
        'CCG': ['P', 'four', 0],
        'TCT': ['S', 'four', 0],
        'TCC': ['S', 'four', 0],
        'TCA': ['S', 'four', 0],
        'TCG': ['S', 'four', 0],
        'AGT': ['S_', 'two', 0],
        'AGC': ['S_', 'two', 0],
        'ACT': ['T', 'four', 0],
        'ACC': ['T', 'four', 0],
        'ACA': ['T', 'four', 0],
        'ACG': ['T', 'four', 0],
        'TGG': ['W', 'one', 0],
        'TAT': ['Y', 'two', 0],
        'TAC': ['Y', 'two', 0],
        'GTT': ['V', 'four', 0],
        'GTC': ['V', 'four', 0],
        'GTA': ['V', 'four', 0],
        'GTG': ['V', 'four', 0],
        'TAA': ['*', 'three', 0],
        'TGA': ['*', 'three', 0],
        'TAG': ['*', 'three', 0],
        'XXX': ['_missing', 'none', 0]
    }
    return codontablecount


def actgtable():
    actgtablecount = {
        'A': ['A', 0],
        'C': ['C', 0],
        'T': ['T', 0],
        'G': ['G', 0],
        'X': ['_missing', 0]
    }
    return actgtablecount


# Wright 1990 eqautions
def calc_fa(tmp):
    counts = [x[2] for x in tmp]
    na = sum(counts)
    return [((na * sum([(x / float(na)) ** 2 for x in counts])) - 1) / (na - 1), sum(counts)]


def calc_frc(fa, nrc):
    return 1 / float(nrc) * float(sum(fa))


# Latex function
# \widehat{N}_{c} = 2 + GC_{(3)} + (\frac{29}{GC^{2}_{(3)} + (1 - GC^{2}_{(3)})})
def calc_eq4wright(fgc):
    return 2 + float(fgc) + (29 / ((fgc ** 2) + ((1 - fgc) ** 2)))


# Sun 2012 equations
# FCF according to equation (2) without pseudocounts for eq2Sun
# F_{CF} = \sum_{i=1}^{m} (\frac {n_{i}}{n})^{2}
def calc_fcf(tmp):
    counts = [x[2] for x in tmp]
    na = sum(counts)
    if na == 0:
        return 'NA', 0
    return sum([(x / float(na)) ** 2 for x in counts]), sum(counts)


# FCF according to equation (3) using pseudocounts for eq5Sun
# F_{CF} = \sum_{i=1}^{m} (\frac {n_{i}+1}{n+m})^{2}
def calc_fcf_(tmp):
    counts = [x[2] for x in tmp]
    pseudocounts = [x + 1 for x in counts]
    na = sum(pseudocounts)
    return sum([(x / float(na)) ** 2 for x in pseudocounts]), sum(pseudocounts)


def calc_eq2sun(codoncounts, six2fourtwo):
    six_fcf = []
    six_nrc = []
    four_fcf = []
    four_nrc = []
    three_fcf = []
    three_nrc = []
    two_fcf = []
    two_nrc = []
    one_fcf = []
    one_nrc = []
    four_aa = []
    three_aa = []
    two_aa = []
    one_aa = []
    ncsix = None
    nc = None
    if not six2fourtwo:
        six_aa = ['S', 'R', 'L']
        four_aa = ['A', 'P', 'T', 'G', 'V']
        three_aa = ['I']
        two_aa = ['C', 'E', 'D', 'F', 'H', 'K', 'N', 'Q', 'Y']
        one_aa = ['M', 'W']
        for aa in six_aa:
            tmp_six = [codoncounts[x] for x in codoncounts.keys() if codoncounts[x][0] == aa]
            six = calc_fcf(tmp_six)
            six_fcf.append(six[0])
            six_nrc.append(six[1])
        # setting codons with no counts to average value 1/6
        six_fcf = [float(1) / float(6) if x == 'NA' else x for x in six_fcf]
        ncsix = len(six_fcf) / (sum(six_fcf) / len(six_fcf))
    if six2fourtwo:
        four_aa = ['A', 'P', 'T', 'G', 'V', 'S', 'R', 'L']
        three_aa = ['I']
        two_aa = ['C', 'E', 'D', 'F', 'H', 'K', 'N', 'Q', 'Y', 'S_', 'R_', 'L_']
        one_aa = ['M', 'W']
    for aa in four_aa:
        tmp_four = [codoncounts[x] for x in codoncounts.keys() if codoncounts[x][0] == aa]
        four = calc_fcf(tmp_four)
        four_fcf.append(four[0])
        four_nrc.append(four[1])
    # setting codons with no counts to average value 1/4
    four_fcf = [float(1) / float(4) if x == 'NA' else x for x in four_fcf]
    ncfour = len(four_fcf) / (sum(four_fcf) / len(four_fcf))
    for aa in three_aa:
        tmp_three = [codoncounts[x] for x in codoncounts.keys() if codoncounts[x][0] == aa]
        three = calc_fcf(tmp_three)
        three_fcf.append(three[0])
        three_nrc.append(three[1])
    # setting codons with no counts to average value 1/3
    three_fcf = [float(1) / float(3) if x == 'NA' else x for x in three_fcf]
    ncthree = len(three_fcf) / (sum(three_fcf) / len(three_fcf))
    for aa in two_aa:
        tmp_two = [codoncounts[x] for x in codoncounts.keys() if codoncounts[x][0] == aa]
        two = calc_fcf(tmp_two)
        two_fcf.append(two[0])
        two_nrc.append(two[1])
    # setting codons with no counts to average value 1/2
    two_fcf = [float(1) / float(2) if x == 'NA' else x for x in two_fcf]
    nctwo = len(two_fcf) / (sum(two_fcf) / len(two_fcf))
    for aa in one_aa:
        tmp_one = [codoncounts[x] for x in codoncounts.keys() if codoncounts[x][0] == aa]
        one = calc_fcf(tmp_one)
        one_fcf.append(one[0])
        one_nrc.append(one[1])
    # setting codons with no counts to average value 1/1
    one_fcf = [float(1) / float(1) if x == 'NA' else x for x in one_fcf]
    ncone = len(one_fcf) / (sum(one_fcf) / len(one_fcf))
    if not six2fourtwo:
        nc = ncone + nctwo + ncthree + ncfour + ncsix
    if six2fourtwo:
        nc = ncone + nctwo + ncthree + ncfour
    return nc


# N_{c} = \frac {K_{1} \times \sum_{i}^{K_{1}} n_{i}}{\sum_{i=1}^{K_{i}}(n_{i} \times F_{CF_{i}})}
#  + \frac {K_{2} \times \sum_{i}^{K_{2}} n_{i}}{\sum_{i=1}^{K_{2}}(n_{i} \times F_{CF_{i}})}
#  + \frac {K_{3} \times \sum_{i}^{K_{3}} n_{i}}{\sum_{i=1}^{K_{3}}(n_{i} \times F_{CF_{i}})}
#  + \frac {K_{4} \times \sum_{i}^{K_{4}} n_{i}}{\sum_{i=1}^{K_{4}}(n_{i} \times F_{CF_{i}})}
#  + \frac {K_{6} \times \sum_{i}^{K_{6}} n_{i}}{\sum_{i=1}^{K_{6}}(n_{i} \times F_{CF_{i}})}
def calc_eq5sun(codoncounts, six2fourtwo):
    six_fcf_nrc = []
    four_fcf_nrc = []
    three_fcf_nrc = []
    two_fcf_nrc = []
    one_fcf_nrc = []
    four_aa = []
    three_aa = []
    two_aa = []
    one_aa = []
    ncsix = None
    nc = None
    if not six2fourtwo:
        six_aa = ['S', 'R', 'L']
        four_aa = ['A', 'P', 'T', 'G', 'V']
        three_aa = ['I']
        two_aa = ['C', 'E', 'D', 'F', 'H', 'K', 'N', 'Q', 'Y']
        one_aa = ['M', 'W']
        for aa in six_aa:
            tmp_six = [codoncounts[x] for x in codoncounts.keys() if codoncounts[x][0] == aa]
            six = calc_fcf_(tmp_six)
            six_fcf_nrc.append(six)
        ncsix = len(six_fcf_nrc) / (sum([x[0] * x[1] for x in six_fcf_nrc]) / (sum([x[1] for x in six_fcf_nrc])))
    if six2fourtwo:
        four_aa = ['A', 'P', 'T', 'G', 'V', 'S', 'R', 'L']
        three_aa = ['I']
        two_aa = ['C', 'E', 'D', 'F', 'H', 'K', 'N', 'Q', 'Y', 'S_', 'R_', 'L_']
        one_aa = ['M', 'W']
    for aa in four_aa:
        tmp_four = [codoncounts[x] for x in codoncounts.keys() if codoncounts[x][0] == aa]
        four = calc_fcf_(tmp_four)
        four_fcf_nrc.append(four)
    ncfour = len(four_fcf_nrc) / (sum([x[0] * x[1] for x in four_fcf_nrc]) / (sum([x[1] for x in four_fcf_nrc])))
    for aa in three_aa:
        tmp_three = [codoncounts[x] for x in codoncounts.keys() if codoncounts[x][0] == aa]
        three = calc_fcf_(tmp_three)
        three_fcf_nrc.append(three)
    ncthree = len(three_fcf_nrc) / (sum([x[0] * x[1] for x in three_fcf_nrc]) / (sum([x[1] for x in three_fcf_nrc])))
    for aa in two_aa:
        tmp_two = [codoncounts[x] for x in codoncounts.keys() if codoncounts[x][0] == aa]
        two = calc_fcf_(tmp_two)
        two_fcf_nrc.append(two)
    nctwo = len(two_fcf_nrc) / (sum([x[0] * x[1] for x in two_fcf_nrc]) / (sum([x[1] for x in two_fcf_nrc])))
    for aa in one_aa:
        tmp_one = [codoncounts[x] for x in codoncounts.keys() if codoncounts[x][0] == aa]
        one = calc_fcf_(tmp_one)
        one_fcf_nrc.append(one)
    ncone = len(one_fcf_nrc) / (sum([x[0] * x[1] for x in one_fcf_nrc]) / (sum([x[1] for x in one_fcf_nrc])))
    if not six2fourtwo:
        nc = ncone + nctwo + ncthree + ncfour + ncsix
    if six2fourtwo:
        nc = ncone + nctwo + ncthree + ncfour
    return nc


# Relative Synonymous Codon Usage
# RSCU_{i,j} = \frac{X_{i,j}}{\frac{1}{n_{i}} \times \sum_{j=1}^{n_{i}} X_{i,j}}
def calc_rscu(codoncounts, codontable):
    rscutable = codontable()
    for codon in rscutable.keys():
        aa = rscutable[codon][0]
        tmp_aa = [[x, codoncounts[x][2]] for x in codoncounts.keys() if codoncounts[x][0] == aa]
        counts = [x[1] for x in tmp_aa]
        na = sum(counts)
        if na == 0:
            rscutable[codon][2] = 'NA'
        else:
            tmp_codon_na = [x[1] for x in tmp_aa if x[0] == codon][0]
            tmp_codon_freq = tmp_codon_na / float(na)
            tmp_codon_rscu = tmp_codon_freq * len(counts)
            rscutable[codon][2] = tmp_codon_rscu
    return rscutable


def gcbypos(codoncounts, six2fourtwo):
    aminoacids = ['A', 'C', 'E', 'D', 'G', 'F', 'I', 'H', 'K', 'M', 'L', 'N', 'Q', 'P', 'S', 'R', 'T', 'W', 'V', 'Y']
    if six2fourtwo:
        aminoacids = ['A', 'C', 'E', 'D', 'G', 'F', 'I', 'H', 'K', 'M', 'L', 'L_', 'N', 'Q', 'P', 'S', 'S_', 'R', 'R_',
                      'T', 'W', 'V', 'Y']
    gcone_gc = []
    gcone_sum = []
    gctwo_gc = []
    gctwo_sum = []
    gcthree_gc = []
    gcthree_sum = []
    gcone = None
    gctwo = None
    gcthree = None
    for aa in aminoacids:
        tmp_aa_one = [[x[0], codoncounts[x][2]] for x in codoncounts.keys() if codoncounts[x][0] == aa]
        tmp_aa_two = [[x[1], codoncounts[x][2]] for x in codoncounts.keys() if codoncounts[x][0] == aa]
        tmp_aa_three = [[x[2], codoncounts[x][2]] for x in codoncounts.keys() if codoncounts[x][0] == aa]
        tmp_aa_one_sum = sum([x[1] for x in tmp_aa_one])
        tmp_aa_two_sum = sum([x[1] for x in tmp_aa_two])
        tmp_aa_three_sum = sum([x[1] for x in tmp_aa_three])
        gc_one = sum([x[1] for x in tmp_aa_one if x[0] == 'G' or x[0] == 'C'])
        gcone_sum.append(tmp_aa_one_sum)
        gcone_gc.append(gc_one)
        gc_two = sum([x[1] for x in tmp_aa_two if x[0] == 'G' or x[0] == 'C'])
        gctwo_sum.append(tmp_aa_two_sum)
        gctwo_gc.append(gc_two)
        gc_three = sum([x[1] for x in tmp_aa_two if x[0] == 'G' or x[0] == 'C'])
        gcthree_sum.append(tmp_aa_three_sum)
        gcthree_gc.append(gc_three)
    if float(sum(gcone_gc)) == 0 and float(sum(gcone_sum)) == 0:
        print ('GCone_gc and GCone_sum zero\n')
        gcone = 0
    if float(sum(gcone_gc)) != 0 or float(sum(gcone_sum)) != 0:
        gcone = float(sum(gcone_gc)) / float(sum(gcone_sum))
    if float(sum(gctwo_gc)) == 0 and float(sum(gctwo_sum)) == 0:
        print ('GCtwo_gc and GCtwo_sum zero\n')
        gctwo = 0
    if float(sum(gctwo_gc)) != 0 or float(sum(gctwo_sum)) != 0:
        gctwo = float(sum(gctwo_gc)) / float(sum(gctwo_sum))
    if float(sum(gcthree_gc)) == 0 and float(sum(gcthree_sum)) == 0:
        print ('GCthree_gc and GCthree_sum zero\n')
        gcthree = 0
    if float(sum(gcthree_gc)) != 0 or float(sum(gcthree_sum)) != 0:
        gcthree = float(sum(gcthree_gc)) / float(sum(gcthree_sum))
    return [gcone, gctwo, gcthree]
