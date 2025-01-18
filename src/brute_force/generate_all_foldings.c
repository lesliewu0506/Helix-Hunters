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
    FILE* file = fopen(filename, "w");

    // Check for error in opening file
    if (!file)
    {
        printf("Error while opening file.\n");
        return;
    }

    // Calculate total amount of combinations
    unsigned long long int total = 1;
    for (int i = 0; i < length - 2; i++)
    {
        total *= 3;
    }
    printf("%lli\n", total);

    // Create buffer to store relative directions
    int* relative_directions = (int*)malloc(sizeof(int) * (length - 2));

    // Close file and free memory
    fclose(file);
    free(relative_directions);
}