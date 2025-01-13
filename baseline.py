import src.algorithm.folding_functions as fold
from src.algorithm.folding_functions import random_iterated

def main():
    sequences = ["HHPHHHPHPHHHPH",
                 "HPHPPHHPHPPHPHHPPHPH",
                 "PPPHHPPHHPPPPPHHHHHHHPPHHPPPPHHPPHPP",
                 "HHPHPHPHPHHHHPHPPPHPPPHPPPPHPPPHPPPHPHHHHPHPHPHPHH",
                 "PPCHHPPCHPPPPCHHHHCHHPPHHPPPPHHPPHPP",
                 "CPPCHPPCHPPCPPHHHHHHCCPCHPPCPCHPPHPC",
                 "HCPHPCPHPCHCHPHPPPHPPPHPPPPHPCPHPPPHPHHHCCHCHCHCHH",
                 "HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH"]

    for protein_sequence in sequences:
        random_iterated(protein_sequence, fold.random_fold)

if __name__ == "__main__":
    main()