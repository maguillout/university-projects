# -*- coding: utf-8 -*-

import getopt
import sys
import pickle
from collections import OrderedDict
from index import Fmi

def kmer_list(read, k):
    """
    Return all the kmers of the read.

    Args_
        read (str): The read
        k (int): Length of kmers

    Returns_
        k_mer_list (list of strings):
            List of all the kmers of the reads.
            List indexes are equivalent to kmer position in the read.
    """
    k_mer_list = []
    for i in range(len(read)-k+1):
        k_mer = read[i:i+k]
        k_mer_list.append(k_mer)
    return k_mer_list

def reverse_comp(sequence):
    """
    Return the reversed brand of the complementary brand of a sequence.

    Args_
        sequence (str): A nucleotidic sequence

    Returns_
        str: The complementary brand of the sequence reversed

    """
    reverse_comp = ""
    for char in sequence:
        if char == 'A':
            reverse_comp += 'T'
        elif char == 'C':
            reverse_comp += 'G'
        elif char == 'G':
            reverse_comp += 'C'
        elif char == 'T':
            reverse_comp += 'A'
    return reverse_comp[::-1]

def mismatches(s1, s2, pos):
    """
    All existing substitutions if s2 is aligned on s1 at position pos.

    Args_
        s1 (str): sequence
        s2 (str): read
        pos (int): position of sequence s1 in the genome

    Returns_
        snp (Dictionnary): Keys are positions. Values are substitutions.

    """
    if len(s1) != len(s2):
        ValueError(f"{s1} and {s2} do not have the same size.")

    else:
        snp = {}
        for i in range(len(s1)):
            ref = s1[i]
            alt = s2[i]

            if ref != alt:
                snp[i + pos] = alt

        return snp

def mapping_read(read, k, genome, max_hamming):
    """
    Align the read on the reference genome. Returns the best position.

    Args_
        read (string): the read which will be aligned
        k (int): length of kmers
        genome (string): reference genome
        max_hamming (int): maximal number of substitutions authorized

    Returns_
        position_mapping (int): read position on the genome,
        in such a way that the alignment is optimal
        nb_min: the number of substitutions if the read is aligned at
        this position
        snp (dictionnary): it stores all the substitutions.
        Keys are positions and values are alternative alleles.

    """
    list_kmers = kmer_list(read, k)
    positions_testees = []
    # if there are substitutions, nb_min_substitions cannot be higher than max_hamming
    nb_min_substitions = max_hamming + 1
    # position_mapping cannot be higher than genome length
    position_mapping = len(genome)
    snp = ""

    for i in range(len(list_kmers)):
        k_mer = list_kmers[i]
        occurences = fm_index.get_occurences(k_mer)

        for occ in occurences:
            pos = occ - i  # kmer position on the genome, must shift by i

            if pos not in positions_testees and pos <= (len(genome)-len(read)):
                positions_testees.append(pos)
                partie_g = genome[pos:pos+len(read)]
                snp = mismatches(partie_g, read, pos)

                if snp is None:  # if there are no substitutions
                    nb_substitutions = 0
                    snp = ""

                else:
                    nb_substitutions = len(snp)

                if nb_substitutions < nb_min_substitions:
                    nb_min_substitions = nb_substitutions
                    position_mapping = pos

                elif (nb_substitutions == nb_min_substitions) and (pos < position_mapping):
                    nb_min_substitions = nb_substitutions
                    position_mapping = pos
                    
    return(position_mapping, nb_min_substitions, snp)



def multiple_mapping(list, k, genome, max_hamming):
    """
    Align several reads on the genome, and select the best brand.

    Args_
        list (list of strings): list of reads
        k (int): length of kmers
        genome (string): reference genome
        max_hamming (int): maximal number of substitutions authorized

    Returns_
        read_positions (dictionnary): Mapping dictionnary.
        Each entry corresponds to one read.
        Keys are reads, values are mapping position, brand ('+' or '-')
        and substitutions (if necessary)
        vcf (dictionnary): Variants dictionnary. Keys are positions in the
        genome. Values are alternative alleles and abundances.
    """
    read_positions = {}
    vcf = {}

    for read in list:

        r_rev_comp = reverse_comp(read)
        pos_brand_plus, subs_brand_plus, snp_plus = mapping_read(read, k, genome, max_hamming)
        # Mapping of the reverseved complementary read:
        pos_brand_minus, subs_brand_minus, snp_minus = mapping_read(r_rev_comp, k, genome, max_hamming)

        # Check if the read can map in both directions
        if (subs_brand_plus == max_hamming + 1) and  (subs_brand_minus == max_hamming + 1):
            # The read r cannot map with this sequence, for the two brands
            continue

        # Select the best alignment, among brand plus and brand minus

        if (subs_brand_plus < subs_brand_minus) or (subs_brand_plus == subs_brand_minus and pos_brand_plus < pos_brand_minus):
            read_positions[read] = [pos_brand_plus, "+", snp_plus]
            snp = snp_plus
        else:
            read_positions[read] = [pos_brand_minus, "-", snp_minus]
            snp = snp_minus

        for position, alt in snp.items():

            if position not in vcf.keys():
                vcf[position] = [alt, 1]

            else:
                # if this position has already be reported in the vcf dictionnary
                variant_found = False
                for i in range(0, len(vcf[position]), 2):
                    # searching for the current alternative allele among registered variants    
                    if vcf[position][i] == alt:
                        vcf[position][i+1] += 1
                        variant_found = True
                        break

                if not variant_found:
                    vcf[position] += [alt, 1]

    return (read_positions, vcf)


def variant_calling(list, k, genome, max_hamming, m_value):
    """
    Map the reads and call all the variants. Only variants with an abundance higher than m_value are reported.

    Args_
        list (list of strings): list of reads
        k (int): length of kmers
        genome (string): reference genome
        max_hamming (int): maximal number of substitutions authorized
        m_value (int): minimal abundance for a variant to be reported

    Returns_
        info_variant (string): Each line of this string corresponds to
        one variant
    """
    mapping, vcf = multiple_mapping(list, k, genome, max_hamming)
    # remove snp with low abundance
    for position, variant_list in vcf.items():
        alleles_final = []
        for i in range(1, len(variant_list), 2):
            if variant_list[i] >= m_value:
                alleles_final += [variant_list[i-1], variant_list[i]]
        vcf[position] = alleles_final

    # sorting variants by position in genome
    vcf_sorted = OrderedDict(sorted(vcf.items(), key=lambda t: t[0]))
    info_variant = ""
    for pos in vcf_sorted.keys():
        variant_list = vcf_sorted[pos]
        ref = genome[pos]
        for i in range(0, len(variant_list), 2):
            info_variant += f"{pos}\t{ref}\t{variant_list[i]}\t{variant_list[i+1]}\n"
    return info_variant

def usage():
    """Print usage of the program."""
    print(f"\nUSAGE\npython {sys.argv[0]} --ref [genome_file.fa] --index [dumped_index.dp] --reads [reads.fa] -k [k_value] --max_hamming [h_value] --min_abundance [m_value] --out snps.vcf")
    print("Map the reads on the reference genome and do a snp caller")
    print("\nOPTIONS")
    print("\t -g --ref FILE.fa : file with reference genome")
    print("\t -i --index FILE.dp : file with fm index")
    print("\t -r --reads FILE.fa : file with all the reads")
    print("\t -k --kmer int value : length of the k-mer")
    print("\t -d --max_hamming int value : max value of substitution in the read mapping")
    print("\t -a --min_abundance int value : minimum abundance of a variant")
    print("\t -o --out FILE.vcf : vcf file that you can assign")

if __name__ == "__main__":
    ref = None
    index = None
    reads = None
    k = 45
    max_hamming = 5
    min_abundance = 5
    file_vcf = "snps.vcf"

    try:
        opts, _ = getopt.getopt(sys.argv[1:],
                                "g:i:r:k:d:a:o:h",
                                ["ref=", "index=", "reads=", "kmer=", "max_hamming=", "min_abundance=", "out=", "help"])
    except getopt.GetoptError as err:
        print(err)
        usage()
        sys.exit(2)

    for option, arg in opts:
        if option in ("-h", "--help"):
            usage()
            sys.exit(0)
        elif option in ("-g", "--ref"):
            ref = arg
        elif option in ("-i", "--index"):
            index = arg
        elif option in ("-r", "--reads"):
            reads = arg
        elif option in ("-k", "--kmer"):
            k = int(arg)
        elif option in ("-d", "--max_hamming"):
            max_hamming = int(arg)
        elif option in ("-a", "--min_abundance"):
            min_abundance = int(arg)
        elif option in ("-o", "--out"):
            file_vcf = arg
            if file_vcf[-4:] != ".vcf":
                file_vcf = file_vcf + ".vcf"

    if not ref:
        usage()
        sys.exit(0)

    if not index:
        usage()
        sys.exit(0)

    if not reads:
        usage()
        sys.exit(0)

    with open(ref, "r") as genome_file:
        for line in genome_file:
            if line[0] == ">":
                continue
            else:
                genome = line.replace("\n", "")
                break

    fm_index = pickle.load(open(index, "rb"))

    with open(reads, "r") as reads_file:
        list_reads = []
        for line in reads_file:
            if line[0] == ">":
                continue
            list_reads.append(line.replace("\n", ""))

    with open(file_vcf, "w") as vcf_file:
        vcf_file.write(f"#REF: {ref.split('/')[-1]}\n")
        vcf_file.write(f"#READS: {reads.split('/')[-1]}\n")
        vcf_file.write(f"#K: {k}\n")
        vcf_file.write(f"#MAX_SUBST: {max_hamming}\n")
        vcf_file.write(f"#MIN_ABUNDANCE: {min_abundance}\n")
        vcf_file.write(variant_calling(list_reads, k, genome, max_hamming, min_abundance))
