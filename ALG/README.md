# To do list

## But du projet

Implémenter un mappeur-SNPcaller de reads sur génome de référence. Il n'y aura pas d'insertion ou de délétion, que des substitutions.

## Différentes étapes

- [x] Indexation

    Programme `index.py` qui créé un FM-index     (FMI) du génome en entrée.

    - input : genome.fa
    - Tableau des suffixe-array (programme 'tools_karkkainen_sanders.py')
    - Fonction donnant la BWT `get_bwt`.
    - Fonction donnant le nombre d'occurences de chaque base `get_n`.
    - Fonction indiquant le rang des occurences dans la BWT `get_r`.
    - Fonction LF `left_first`.
    - Stocker les informations dans un fichier output : dumped_index.dp

    -> `index.py` est spécifique aux méthodes du FM-index. (si jamais on veut utiliser une autre méthode d'indexation, on utilise un autre fichier)


- [x] Mapping

    Programme `map.py` mapper l'ensemble des reads sur un génome de référence.

    - input : genome.fa, reads.fa, dumped_index.dp
    - prend en argument : k : longueur des k-mers, max_hamming : nbre max de substitutions possibles et min d'abundance : permet de s'intéresser qu'au SNP  avec une abundance supérieur à la valeur.
    - fonction `get_occurence` de read avec l'utilisation du FMI, peut partir de celle vu en TP en ajoutant qu'elle peut prendre en compte les substitutions et en choississant la position idéale pour le mapping du read. stocker les positions pour chaque read. (A MODIFIER)

    - Etapes du mapping pour chaque read:
        - prend un mot de taille k au début du read.
        - cherche les occurences de ce k-mer sur le génome de référence.
        - Aligne le read sur les positions retenues en comptant le nombre de substitutions.
        - Retenir le mapping avec le moins de substitutions, si égalité prend celui le plus à gauche.

    - fonction `print_mapping` qui va afficher le mapping d'un ou des reads par rapport au génome de référence.
    - fonction `SNP_caller` qui va énumérer le nombre de SNP avec un minimum d'abondance en comparant les reads mappés et mettre en output un fichier vcf.

- [x] Fichier VCF

    Au début nous devrons indiquer plusieurs informations : le nom du fichier du génome de référence, nom du fichier contenant les reads, le k, le max hamming, le minimum d'abondance.
    Devra contenir comme information pour chacun des SNP:
    - la position dans le génome
    - l'allèle de référence
    - l'allèle muté 
    - le nombre de reads possédant l'allèle muté (abondance)


