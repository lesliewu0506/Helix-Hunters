#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

// --------------------------------------------------
// Hulpfunctie om de volgende absolute richting te bepalen
// op basis van de huidige absolute richting én de 'keuze' (0, 1, 2).
// Dit komt overeen met jouw direction_map in Python.
//
// current_dir kan 1, -1, 2 of -2 zijn.
// choice kan 0, 1, 2.
//
// We geven de nieuwe absolute richting terug.
// --------------------------------------------------
int get_next_direction(int current_dir, int choice) {
    // direction_map in Python:
    // {
    //   1:  [ 2,  1, -2],
    //  -1: [-2, -1,  2],
    //   2:  [-1,  2,  1],
    //  -2: [ 1, -2, -1]
    // }
    if (current_dir == 1) {
        if      (choice == 0) return  2;
        else if (choice == 1) return  1;
        else                  return -2;
    } else if (current_dir == -1) {
        if      (choice == 0) return -2;
        else if (choice == 1) return -1;
        else                  return  2;
    } else if (current_dir == 2) {
        if      (choice == 0) return -1;
        else if (choice == 1) return  2;
        else                  return  1;
    } else { // current_dir == -2
        if      (choice == 0) return  1;
        else if (choice == 1) return -2;
        else                  return -1;
    }
}

// --------------------------------------------------
// We vertalen de (sequence_length - 2) keuzes (ieder 0,1,2)
// naar een array van absolute richtingen door steeds
// get_next_direction() te gebruiken.
//
// De eerste richting zetten we hard op 1 (zoals in jouw Python),
// de laatste op 0 (maar we voegen die 0 pas toe als 'laatste stap').
// --------------------------------------------------
int* direction_translator(const int* choices, int seq_len) {
    // seq_len is de lengte van je protein_sequence
    // choices heeft (seq_len - 2) elementen
    // De resulterende folding_sequence heeft dezelfde lengte als het eiwit (seq_len),
    // want in Python return je: [1] + ... + [0] -> lengte seq_len.
    int* folding_sequence = (int*) malloc(seq_len * sizeof(int));
    if (!folding_sequence) {
        fprintf(stderr, "Error: could not allocate folding_sequence.\n");
        return NULL;
    }

    // Eerste absolute richting is 1 (zoals in je Python-code)
    folding_sequence[0] = 1;

    // Bouw de rest via get_next_direction
    for (int i = 1; i < seq_len - 1; i++) {
        // i-1 in folding_sequence is de vorige absolute richting
        // choices[i-1] is de 'relatieve' keuze (0,1,2)
        folding_sequence[i] = get_next_direction(folding_sequence[i-1], choices[i-1]);
    }

    // Laatste is 0
    folding_sequence[seq_len - 1] = 0;

    return folding_sequence;
}

// --------------------------------------------------
// Check of een gegeven folding geldig is door te kijken
// of je ooit twee keer hetzelfde coordinaat bezoekt.
//
// We modelleren de eiwit-structuur 2D door (x, y) te updaten
// op basis van direction_map:
//    1  -> x+1
//   -1  -> x-1
//    2  -> y+1
//   -2  -> y-1
//    0  -> geen beweging (eindpunt)
// --------------------------------------------------
bool check_valid_folding(const int* folding, int seq_len) {
    // We houden alle bezochte coördinaten bij in een simpel array van
    // (x, y) zodat we overlap kunnen detecteren. 
    // In de praktijk zou je een hashset of iets efficiënters gebruiken
    // als seq_len heel groot is.

    // x_coords[i], y_coords[i] = coördinaat van amino i
    int* x_coords = (int*) calloc(seq_len, sizeof(int));
    int* y_coords = (int*) calloc(seq_len, sizeof(int));
    if (!x_coords || !y_coords) {
        fprintf(stderr, "Error: could not allocate coordinate arrays.\n");
        free(x_coords);
        free(y_coords);
        return false;
    }

    // Startcoördinaat is (0,0)
    x_coords[0] = 0;
    y_coords[0] = 0;

    for (int i = 1; i < seq_len; i++) {
        // Bepaal nieuwe positie obv folding[i-1] 
        // (Let op: i-th amino ligt op index i, maar om van i-1 naar i te gaan, gebruik je folding[i-1].)
        int dir = folding[i-1];
        x_coords[i] = x_coords[i-1];
        y_coords[i] = y_coords[i-1];

        if (dir == 1) {
            x_coords[i]++;
        } else if (dir == -1) {
            x_coords[i]--;
        } else if (dir == 2) {
            y_coords[i]++;
        } else if (dir == -2) {
            y_coords[i]--;
        } 
        // dir == 0 => geen beweging (laatste amino)

        // Check overlap met eerdere coördinaten
        for (int j = 0; j < i; j++) {
            if (x_coords[i] == x_coords[j] && y_coords[i] == y_coords[j]) {
                // Overlap gevonden => niet geldig
                free(x_coords);
                free(y_coords);
                return false;
            }
        }
    }

    free(x_coords);
    free(y_coords);
    return true;
}

// --------------------------------------------------
// Genereer alle combinaties van {0,1,2}^(seq_len - 2)
// en check of ze geldig zijn. Zo ja, schrijf ze naar CSV-bestand.
//
// protein_sequence is hier bijvoorbeeld een string "HHHHHH..." van lengte seq_len.
// We gebruiken alleen de lengte, de letters zelf doen we in deze demo verder niks mee.
// --------------------------------------------------
void generate_all_foldings(const char* protein_sequence) {
    int seq_len = strlen(protein_sequence);
    if (seq_len < 2) {
        fprintf(stderr, "Sequence te kort!\n");
        return;
    }

    // Open het CSV-bestand
    char filename[256];
    snprintf(filename, 256, "%s.csv", protein_sequence);

    FILE* fp = fopen(filename, "w");
    if (!fp) {
        fprintf(stderr, "Kon %s niet openen voor schrijven.\n", filename);
        return;
    }

    // Aantal combinaties = 3^(seq_len - 2)
    long long num_choices = 1;
    for (int i = 0; i < seq_len - 2; i++) {
        num_choices *= 3; 
    }

    // We maken een array van (seq_len - 2) voor alle keuzes
    int* choices = (int*) malloc((seq_len - 2) * sizeof(int));
    if (!choices) {
        fprintf(stderr, "Kon choices niet alloceren.\n");
        fclose(fp);
        return;
    }
    for (int i = 0; i < seq_len - 2; i++) {
        choices[i] = 0; // initialize
    }

    // Itereer over alle combinaties in "3-tallig stelsel"
    // (Een simpele "tellus op" in basis 3)
    for (long long idx = 0; idx < num_choices; idx++) {
        // Zet idx om naar basis-3 in 'choices'
        long long temp = idx;
        for (int c = seq_len - 3; c >= 0; c--) {
            choices[c] = temp % 3;  // 0,1,2
            temp /= 3;
        }

        // Maak de absolute folding via direction_translator
        int* folding = direction_translator(choices, seq_len);
        if (!folding) {
            continue; // eventueel wat netter afhandelen
        }

        // Check of deze folding geldig is
        if (check_valid_folding(folding, seq_len)) {
            // Schrijf naar CSV (bijv. komma-gescheiden)
            // Let op: de folding heeft lengte seq_len
            for (int i = 0; i < seq_len; i++) {
                fprintf(fp, "%d", folding[i]);
                if (i < seq_len - 1) {
                    fprintf(fp, ",");
                }
            }
            fprintf(fp, "\n");
        }
        free(folding);
    }

    free(choices);
    fclose(fp);
}

int main(int argc, char** argv) {
    if (argc < 2) {
        printf("Gebruik: %s <protein_sequence>\n", argv[0]);
        return 1;
    }
    generate_all_foldings(argv[1]);
    return 0;
}
