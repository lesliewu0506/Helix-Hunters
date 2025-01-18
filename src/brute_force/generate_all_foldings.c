#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

// Function prototypes
void generate_all_foldings(const char* protein_sequence);

// Main
int main(int argc, char** argv)
{
    // For invalid commands
    if(argc != 2)
    {
        printf("Usage: %s <protein_sequence>\n", argv[0]);
        return 1;
    }
    
    generate_all_foldings(argv[1]);
    return 0;
}

void generate_all_foldings(const char* protein_sequence)
{
    // Get length of sequence
    int length = strlen(protein_sequence);
    
    // Create CSV file
    char filename[256];
    sprintf(filename, "%s.csv", protein_sequence);
    FILE * file = fopen(filename, "w");

    // Check for error in opening file
    if (!file)
    {
        printf("Error while opening file.\n");
        return;
    }


}