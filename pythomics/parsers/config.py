BASE_PAIR_COMPLEMENTS = {'a': 't', 't': 'a', 'c': 'g', 'g': 'c', 'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'n': 'n', 'N': 'N'}

CODON_TABLE = {'AAA': 'K', 'AAC': 'N', 'AAG': 'K', 'AAT': 'N', 'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
               'AGA': 'R', 'AGC': 'S', 'AGG': 'R', 'AGT': 'S', 'ATA': 'I', 'ATC': 'I', 'ATG': 'M', 'ATT': 'I',
               'CAA': 'Q', 'CAC': 'H', 'CAG': 'Q', 'CAT': 'H', 'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
               'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R', 'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
               'GAA': 'E', 'GAC': 'D', 'GAG': 'E', 'GAT': 'D', 'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
               'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G', 'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
               'TAA': '*', 'TAC': 'Y', 'TAG': '*', 'TAT': 'Y', 'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
               'TGA': '*', 'TGC': 'C', 'TGG': 'W', 'TGT': 'C', 'TTA': 'L', 'TTC': 'F', 'TTG': 'L', 'TTT': 'F',
               'NNN': 'X'}

for i in 'ACTG':
    for j in 'ACTG':
        CODON_TABLE['%s%sN' % (i,j)] = 'X'
        CODON_TABLE['%sN%s' % (i,j)] = 'X'
        CODON_TABLE['N%s%s' % (i,j)] = 'X'
    CODON_TABLE['%sNN' % i] = 'X'
    CODON_TABLE['N%sN' % i] = 'X'
    CODON_TABLE['NN%s' % i] = 'X'
            