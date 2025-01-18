#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

// Direction Mapping
static const int direction_map[4][3] = 
{
    { 2,  1, -2},
    {-2, -1,  2},
    {-1,  2,  1},
    { 1, -2, -1}
};

// Function prototypes
void generate_all_foldings(const char* protein_sequence);
int* check_folding(int* rel_dir, int protein_length);
void direction_translator(int* rel_dir, int* absolute, int relative_length);
int direction_to_row(int direction);

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
    int protein_length = strlen(protein_sequence);
    
    // Create CSV file
    char filename[60];
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
    for (int i = 0; i < protein_length - 2; i++)
    {
        total *= 3;
    }

    // Create buffer to store relative directions
    int* relative_directions = (int*)malloc(sizeof(int) * (protein_length - 2));

    // Turn every value into base 3 to get all relative direction combinations
    for (unsigned long long int value = 0; value < total; value++)
    {
        unsigned long long int temp = value;
        for (int i = protein_length - 3; i >= 0; i--)
        {
            relative_directions[i] = temp % 3;
            temp /= 3;
        }

        // Get absolute directions
        int* result = check_folding(relative_directions, protein_length);

        // Write to CSV
        if (result != NULL)
        {
            for (int i = 0; i < protein_length; i++)
            {
                fprintf(file, "%i", result[i]);
                if (i < protein_length - 1)
                {
                    fprintf(file, ",");
                }
            }
            fprintf(file, "\n");
        }
    }

    // Close file and free memory
    fclose(file);
    free(relative_directions);
}

// Translates relative directions into absolute directions and then checks if it is valid
int* check_folding(int* relative_directions, int protein_length)
{
    // Create array for only in this function
    static int absolute_directions[60];

    int relative_length = protein_length - 2;
    direction_translator(relative_directions, absolute_directions, relative_length);

    return absolute_directions;
}

// Translates array from relative directions to absolute directions
void direction_translator(int* relative_directions, int* absolute_directions, int relative_length)
{
    // Start with first direction
    int current_direction = 1;
    absolute_directions[0] = current_direction;

    for (int i = 0; i < relative_length; i++)
    {
        int row = direction_to_row(current_direction);
        int relative_direction = relative_directions[i];
        int absolute_direction = direction_map[row][relative_direction];

        absolute_directions[i + 1] = absolute_direction;
        current_direction = absolute_direction;
    }

    // Add final direction
    absolute_directions[relative_length + 1] = 0;
}

// Helper function for translating absolute direction to index row for mapping
int direction_to_row(int direction) 
{
    // Use switch instead of if else statements
    switch(direction) 
    {
        case 1:  return 0;
        case -1: return 1;
        case 2:  return 2;
        case -2: return 3;
    }
}