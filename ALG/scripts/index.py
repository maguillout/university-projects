# -*- coding: utf-8 -*-

import getopt, sys
import pickle 
import tools_karkkainen_sanders as tks


class Fmi:
    """Create an object with all information to use the indexation algorithm fm-index."""
    
    def __init__(self, genome):
        self.genome = genome + "$"
        self.bwt = ""
        self.sa = []
        self.n = {}
        self.r = []
        
        self.set_sa()        
        self.set_bwt()
        self.set_n()
        self.set_r()
        
    def set_sa(self):
        """Set the suffixe array of the genome sequence."""
        self.sa = tks.simple_kark_sort(self.genome)
        
    def set_bwt(self):
        """Set the bwt of the genome sequence."""
        bwt = ""
        for i in range(len(self.sa)):
            if self.sa[i] == 0:
                bwt += "$"
            else:
                bwt += self.genome[self.sa[i]-1]
        self.bwt = bwt
    
    def set_n(self) -> {}:
        """Set the number of each iterations of A, C, G, T and $."""
        n = {"$": 0, "A": 0, "C":0, "G": 0, "T":0}
        for letter in self.bwt: 
            n[letter] += 1
        self.n = n
    
    def set_r(self) -> []:
        """Set the list of r[i] = rank in the bwt of the bwt[i] character."""
        n = {}
        r = []
        for letter in self.bwt: 
            if letter not in n: 
                n[letter] = 0
            n[letter] += 1
            r.append(n[letter])
        self.r = r
    
    def left_first(self, alpha: chr, k: int) -> int:
        """
        Return the k-ième alpha character in the bwt.
        
        Args_
            alpha (chr): The character which we want in the bwt.
            k (int): The rank of the character in the bwt.
        
        Raises_
            ValueError: if character alpha not in the bwt
        
        Returns_
            int value of the k-ième alpha character
        """
        assert k <= self.n[alpha], f"Cannot ask for the {k}^th {alpha}, it does not exist."
        if alpha == "$": return 0
        if alpha == "A": return self.n["$"] + k - 1
        if alpha == "C": return self.n["$"] + self.n["A"] + k - 1
        if alpha == "G": return self.n["$"] + self.n["A"] + self.n["C"] + k - 1
        if alpha == "T": return self.n["$"] + self.n["A"] + self.n["C"] + self.n["G"] + k - 1
        raise ValueError(f"Character {alpha} not in the bwt")

    def get_first(self, alpha: chr, i: int, j: int) -> int:
        """
        Detect the first occurence of alpha in bwt for position in [i,j].
        
        Args_
            alpha (chr): The character we need to find in the bwt
            i (int): Start of the interval
            j (int): Stop of the interval
        Returns_
            The position of the first occurence in bwt
        """
        pos_bwt = i
        for pos_bwt in range(i,j+1): #on parcourt bwt pour trouver la 1re ligne x (entre i et j) telle que BWT[x] = current_char
            #+1 parce que range s'arrête une itération avant
            if self.bwt[pos_bwt] == alpha:
                return pos_bwt
        return -1

    def get_last(self, alpha: chr, i: int, j: int) -> int:
        """
        Detect the last occurence of alpha in bwt for position in [i,j].
        
        Args_
            alpha (chr): The character we need to find in the bwt
            i (int): Start of the interval
            j (int): Stop of the interval
        Returns_
            The position of the last occurence in bwt
        """
        pos_bwt = j
        for pos_bwt in range(j,i-1,-1): #on parcours bwt à l'envers pour trouver la dernière ligne x (entre i et j) telle que BWT[x] = current_char
            #-1 parce que range s'arrête une itération avant
            if self.bwt[pos_bwt] == alpha:
                return pos_bwt
        
    def get_occurences(self, pattern: chr) -> []:
        """
        Return occurences list of pattern in the genome.
        
        Args_
            pattern (chr): The pattern who need to find the genome
            
        Returns_
            A list of position in the genoma
        """
        # Initialisation :
        i = 0
        j = len(self.bwt)-1
        # Parcourir la query : Chercher le plus grand intervalle [new_i, new_j] avec i ≤ new_i ≤ new_j ≤ j, 
        # tel que BWT[new_i] et BWT[new_j] aient pour valeur le caractère précédent dans Q
        for posP in range(len(pattern)-1,-1,-1):
            current_char = pattern[posP]
            
            new_i = self.get_first(current_char, i, j)
            if new_i == -1: #on sort si on a parcouru tout bwt sans trouver current_char
                    return []
            new_j = self.get_last(current_char, i, j)
            i = self.left_first(self.bwt[new_i], self.r[new_i])
            j = self.left_first(self.bwt[new_j], self.r[new_j])
        occurences = []
        for l in range(i,j+1):
            occurences.append(self.sa[l])
        return occurences


def usage():
    """Print usage of the program."""
    print(f"\nUSAGE\npython {sys.argv[0]} --ref [genome_file.fa] --out [dumped_index.dp]")
    print("Create the FMI of the reference genome for read mapping")
    print("\nOPTIONS")
    print("\t -g --ref FILE.fa : fasta file name of the reference genome")
    print("\t -o --out FILE.dp : dumped file that you can assign")

if __name__ == "__main__":
    
    file_name = None
    file_dp = "dumped_index.dp"
    
    try:
        opts, _ = getopt.getopt(sys.argv[1:],"g:o:h",
                                ["ref=", "help", "out="])
    except getopt.GetoptError as err:
        print(err)
        usage()
        sys.exit(2)
        
    for option, arg in opts:  
        if option in ("-h", "--help"):
            usage()
            sys.exit(0)
            
        elif option in ("-g", "--ref"):
            file_name = arg
            
        elif option in ("-o", "--out"):
            file_dp = arg
            if file_dp[-3:] != ".dp":
                file_dp = file_dp + ".dp"
    
    if not file_name:
        usage()
        sys.exit(0)
        
    with open(file_name, "r") as genome_file:
        for line in genome_file:
            if line[0] == ">":
                continue
            else:
                fm_index = Fmi(line.replace("\n", ""))
                break
    pickle.dump(fm_index, open(file_dp, "wb"))


