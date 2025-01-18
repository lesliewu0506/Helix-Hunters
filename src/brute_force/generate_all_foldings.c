#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

// Main
int main(int argc, char** argv)
{
    // For invalid commands
    if(argc < 2)
    {
        printf("Usage: %s <protein_sequence>\n", argv[0]);
        return 1;
    }
    
    generate_all_foldings(argv[1]);
    return 0;
}

void generate_all_foldings(const char* protein_sequence)
{
    
}