# -*- coding: utf-8 -*-
"""
Determines amino acid sequence from section of DNA.

@author: Becca Getto

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('T')
    'A'
    >>> get_complement('C')
    'G'
    >>> get_complement('G')
    'C'
    """

    if nucleotide == 'A':
        return 'T'
    elif nucleotide == 'T':
        return 'A'
    elif nucleotide == 'G':
        return 'C'
    elif nucleotide == 'C':
        return 'G'
    else:
        return None 


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    # assigning the complement nucleotide string to new variable called complement
    complement = ''

    for i in range(0, len(dna)):
        complement += get_complement(dna[i])

    # reversing the order to the complement string
    return complement[::-1]


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """

    # look for stop codons by traversing string three nucleotides at a time
    for i in range(0, len(dna), 3):
        codon = dna[i:i+3]                      # assigning codon as set of three nucleotides
        if codon in ['TAG', 'TAA', 'TGA']:      # if stop codon is found, return dna up to that point
            return dna[:i]

    return dna


def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    >>> find_all_ORFs_oneframe("GTCATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """

    nnORFs = []
    i = 0

    while (i < len(dna)):
        codon = dna[i:i+3]
        if codon == 'ATG':                      # look for start codons
            orf = rest_of_ORF(dna[i:])          # collect nucleotides until stop codon is reached 
            nnORFs.append(orf)
            i += len(orf)
        else:                                   # keep looking if no start codon found
            i += 3

    return nnORFs 


def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    
    allORFs = []

    # checks for possible reading frames starting at the 1st, 2nd, then 3rd element
    for i in range(3):                          
        ORFs = find_all_ORFs_oneframe(dna[i:])  
        allORFs += ORFs

    return allORFs


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    allORFs_both_strands = []

    allORFs_forward = find_all_ORFs(dna)

    reverse_complement = get_reverse_complement(dna)
    allORFs_reverse_complement = find_all_ORFs(reverse_complement)

    allORFs_both_strands = allORFs_forward + allORFs_reverse_complement

    return allORFs_both_strands


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    allORFs = find_all_ORFs_both_strands(dna)
    longest_ORF = ''

    # check each ORF to determine if it is longer than the previous longest ORF
    for ORF in allORFs:
        if len(ORF) > len(longest_ORF):
            longest_ORF = ORF

    return longest_ORF



def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    
    all_the_longest_ORFs = []
    max_longest_ORF = ''
    i = 0

    while i <= num_trials:
        shuffled_dna = ''.join(random.sample(dna, len(dna)))    # shuffles the dna
        shuffled_longest_ORF = longest_ORF(shuffled_dna)        # find longest ORF in the shuffled dna
        
        if len(shuffled_longest_ORF) > len(max_longest_ORF):    # determines longest ORF from all of the shuffles
            max_longest_ORF = shuffled_longest_ORF
        i += 1

    return len(max_longest_ORF)


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment
                 if an invalid sequence (dna sequence that is not a multiple fo 3) returns False

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCGACGACGA")
        'MRRR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        Last element of dna is not a full codon
        False
    """

    # check to ensure codons are 3 nucleotides long
    if len(dna) % 3 != 0:
        print 'Last element of dna is not a full codon'
        return False

    amino_acid = []
    i = 0
    
    # converts each codon to the corresponding amino acid over the length of dna
    while (i < len(dna)):
        codon = dna[i:i+3]
        amino_acid.append(aa_table[codon])
        i += 3                                  # iterate by 3 as codons are 3 nucleotides in length

    amino_acid_string = ''.join(amino_acid)
    
    return amino_acid_string


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    
    # determines reasonable guess for how long a coding region should be
    threshold = longest_ORF_noncoding(dna, 1500)

    all_amino_acid_sequences = []
    allORFs_both_strands = find_all_ORFs_both_strands(dna)
    
    # converts dna sequences that pass the threshold test into amino acid sequences
    for ORF in allORFs_both_strands:
        if len(ORF) >= threshold:
            all_amino_acid_sequences.append(coding_strand_to_AA(ORF))
        # elif len(ORF) < threshold:
        #     pass

    return all_amino_acid_sequences

if __name__ == "__main__":
    import doctest
    doctest.testmod()

    from load import load_seq
    dna = load_seq("./data/X73525.fa")

    print gene_finder(dna)