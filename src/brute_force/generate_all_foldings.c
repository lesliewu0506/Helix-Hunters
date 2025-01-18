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

static const int DX[5] = {0, 1, -1, 0, 0};
static const int DY[5] = {0, 0, 0, 1, -1};

// Function prototypes
void generate_all_foldings(const char* protein_sequence);
int* check_folding(int* rel_dir, int protein_length);
void direction_translator(const int* rel_dir, int* absolute, int relative_length);
int direction_to_row(int direction);
bool check_valid_folding(const int* absolute_directions, int protein_length);
int direction_to_index(int direction);

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

    if (check_valid_folding(absolute_directions, protein_length))
    {
        return absolute_directions;
    }

    return NULL;
}

// Translates array from relative directions to absolute directions
void direction_translator(const int* relative_directions, int* absolute_directions, int relative_length)
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

// Helper function for translating absolute direction to row index for mapping
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

// Check if there are collisions in protein fold
bool check_valid_folding(const int* absolute_directions, int protein_length)
{
    // Use two arrays to save coordinates
    int x_coordinates[60] = {};
    int y_coordinates[60] = {};

    int x_current = 0;
    int y_current = 0;
    int visited_count = 0;

    for (int i = 0; i < protein_length; i++)
    {
        // First check if already visited
        for (int j = 0; j < visited_count; j++)
        {
            if (x_coordinates[j] == x_current && y_coordinates[j] == y_current)
            {
                return false;
            }
        }
        // Add coordinates to visited
        x_coordinates[i] = x_current;
        y_coordinates[i] = y_current;
        visited_count++;

        // Update coordinates
        int direction_index = direction_to_index(absolute_directions[i]);
        x_current += DX[direction_index];
        y_current += DY[direction_index];
    }

    // Check if last positions is also valid
    // printf("Last\n");
    // printf("%i, %i\n", x_current, y_current);
    for (int j = 0; j < visited_count - 1; j++)
    {
        // printf("%i, %i\n", x_coordinates[j],y_coordinates[j]);
        if (x_coordinates[j] == x_current && y_coordinates[j] == y_current)
        {
            // printf("False\n");
            return false; 
        }
    }
    return true;
}

// Helper function for translating absolute directions to index for mapping
int direction_to_index(int direction)
{    
    // Use switch instead of if else statements
    switch(direction) 
    {
        case 0:  return 0;
        case 1:  return 1;
        case -1: return 2;
        case 2:  return 3;
        case -2: return 4;
    }
}