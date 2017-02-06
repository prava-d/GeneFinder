
"""
Finds amino acid sequences that are coded by some DNA input

@author: Prava Dhulipalla

"""

import random
from load import load_seq
from amino_acids import aa, codons, aa_table   # you may find these useful
test = load_seq("./data/X73525.fa")


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide

        No additional unit tests are needed because all possible cases were
        implemented the same way, so if one works, all work.
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    """
    if nucleotide == 'A':
        return 'T'
    elif nucleotide == 'C':
        return 'G'
    elif nucleotide == 'G':
        return 'C'
    elif nucleotide == 'T':
        return 'A'
    return ''


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string

        No additional unit tests are needed because everything goes through the
        same code and since there are no conditionals, everything executes the
        same way. Two doctests are sufficient.
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    sequence = ''
    i = len(dna) - 1
    while(i >= 0):
        sequence += get_complement(dna[i])
        i -= 1
    return sequence


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string

        I added a doctest for the condition that if there is no stop codon in
        frame, it returns the original string. I am including this because it
        is an important feature of the purpose of the function. I also added a
        unit test for if the last codon is a stop codon, since that seems like
        a reasonable possibility for a genome sequence.
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    >>> rest_of_ORF("ATGCAGA")
    'ATGCAGA'
    >>> rest_of_ORF("ATGCAGTAG")
    'ATGCAG'
    """
    i = 3
    while (i < len(dna) - 2):
        if dna[i:i+3] == 'TAG' or dna[i:i+3] == 'TAA' or dna[i:i+3] == 'TGA':
            stop = i
            return dna[:stop]
        i += 3
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

        I added a doctest in the case that an ATG is within the reading frame
        within another sequence. I included this because am important feature
        of this function is that it only returns non-nested ORFs.
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    >>> find_all_ORFs_oneframe("ATGCGAATGTAGCAGATGGTATAG")
    ['ATGCGAATG', 'ATGGTA']
    >>> find_all_ORFs_oneframe("ATGATGATGTAG")
    ['ATGATGATG']
    """
    oneframe_ORFs = []
    i = 0
    stopped = True
    #Boolean tag for end codon
    while(i < len(dna) - 2):
        if dna[i:i+3] == 'ATG' and stopped:
            oneframe_ORFs.append(rest_of_ORF(dna[i:]))
            stopped = False
        elif dna[i:i+3] == 'TAG' or dna[i:i+3] == 'TAA' or dna[i:i+3] == 'TGA':
            stopped = True
        i += 3
    return oneframe_ORFs


def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

        I added a doctest with ATGs in all the reading frames. This
        is because it is important that the function goes through all the
        reading frames correctly. I also added a doctest for if there is two
        ORFs in one reading frame, because that also seemed possible with an
        actual genome and I wanted to make sure it would do it right.
    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    >>> find_all_ORFs("ATGAATGTATGGATGATGATGTAG")
    ['ATGAATGTATGGATGATGATG', 'ATGTATGGA', 'ATGGATGATGATGTAG']
    >>> find_all_ORFs("ATGCATGCCTAGCTGACCATGCCCTAG")
    ['ATGCATGCC', 'ATGCCC', 'ATGCCTAGC']
    """
    all_ORFs = []
    for i in range(3):
        all_ORFs += find_all_ORFs_oneframe(dna[i:])
    return all_ORFs


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

        I did not create another doctest because all this function is doing is
        applying find_all_ORFs to two strands, and I felt that the doctest
        given sufficed.
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    both_ORFs = []
    both_ORFs += find_all_ORFs(dna)
    both_ORFs += find_all_ORFs(get_reverse_complement(dna))
    return both_ORFs


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    longest = find_all_ORFs_both_strands(dna)
    longest.sort(key=len)
    return longest[len(longest) - 1]
    #sort and find the last one


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    long_ORF = ''
    for i in range(num_trials):
        shuffle_string(dna)
        if len(longest_ORF(dna)) > len(long_ORF):
            long_ORF = longest_ORF(dna)
    return len(long_ORF)


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    protein = ''
    i = 0
    while i < len(dna)-2:
        amino = aa_table[dna[i:i+3]]
        protein += amino
        i += 3
    return protein


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    #threshold = longest_ORF_noncoding(dna, 1500)
    #Ninja said that my code was fine but the threshold calculated was always
    #too long, so he set it to 200 and it worked (Protein BLAST confirmed that
    #it was Salmonella)
    threshold = 200
    sequence = []
    ORFs = find_all_ORFs_both_strands(dna)
    for orf in ORFs:
        if len(orf) > threshold:
            sequence.append(coding_strand_to_AA(orf))
    return sequence


if __name__ == "__main__":
    print(gene_finder(test))
    import doctest
    doctest.run_docstring_examples(get_complement, globals(), verbose=True)
    doctest.run_docstring_examples(get_reverse_complement, globals(),
                                   verbose=True)
    doctest.run_docstring_examples(rest_of_ORF, globals(),
                                   verbose=True)
    doctest.run_docstring_examples(find_all_ORFs_oneframe, globals(),
                                   verbose=True)
    doctest.run_docstring_examples(find_all_ORFs, globals(), verbose=True)
    doctest.run_docstring_examples(find_all_ORFs_both_strands, globals(),
                                   verbose=True)
    doctest.run_docstring_examples(longest_ORF, globals(), verbose=True)
    doctest.run_docstring_examples(coding_strand_to_AA, globals(),
                                   verbose=True)
