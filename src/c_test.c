#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

#define DIRECTIONS_COUNT 4
#define MAX_SEQUENCE_LENGTH 100
#define BATCH_SIZE 1000  // Aantal resultaten per batch

// Mogelijke richtingen
int directions[DIRECTIONS_COUNT] = {-2, -1, 1, 2};

// Struct voor het doorgeven van argumenten aan threads
typedef struct {
    int thread_id;
    int sequence_length;
    FILE *file;
    pthread_mutex_t *file_mutex;
} ThreadArgs;

// Functie om een geldige folding te controleren
int check_valid_sequence(int *folding, int length) {
    if (folding[0] == -1) {
        return 0;
    }
    for (int i = 1; i < length; i++) {
        if (folding[i] == -folding[i - 1]) {
            return 0;
        }
    }
    return 1;
}

// Iteratieve functie om foldings te genereren en te controleren
void *generate_foldings(void *args) {
    ThreadArgs *thread_args = (ThreadArgs *)args;
    int sequence_length = thread_args->sequence_length;
    FILE *file = thread_args->file;
    pthread_mutex_t *file_mutex = thread_args->file_mutex;

    int folding[MAX_SEQUENCE_LENGTH - 2];
    int batch[BATCH_SIZE][MAX_SEQUENCE_LENGTH];
    int batch_index = 0;

    // Bereken start- en eindindex voor deze thread
    long start_index = thread_args->thread_id * (1L << (2 * (sequence_length - 2))) / 4;
    long end_index = (thread_args->thread_id + 1) * (1L << (2 * (sequence_length - 2))) / 4;

    for (long index = start_index; index < end_index; index++) {
        long temp = index;
        for (int i = 0; i < sequence_length - 2; i++) {
            folding[i] = directions[temp % DIRECTIONS_COUNT];
            temp /= DIRECTIONS_COUNT;
        }

        if (check_valid_sequence(folding, sequence_length - 2)) {
            batch[batch_index][0] = 1;
            for (int i = 0; i < sequence_length - 2; i++) {
                batch[batch_index][i + 1] = folding[i];
            }
            batch[batch_index][sequence_length - 1] = 0;
            batch_index++;

            if (batch_index == BATCH_SIZE) {
                pthread_mutex_lock(file_mutex);
                for (int i = 0; i < batch_index; i++) {
                    fprintf(file, "%d,", batch[i][0]);
                    for (int j = 1; j < sequence_length; j++) {
                        fprintf(file, "%d,", batch[i][j]);
                    }
                    fprintf(file, "0\n");
                }
                pthread_mutex_unlock(file_mutex);
                batch_index = 0;
            }
        }
    }

    // Schrijf resterende batch naar bestand
    if (batch_index > 0) {
        pthread_mutex_lock(file_mutex);
        for (int i = 0; i < batch_index; i++) {
            fprintf(file, "%d,", batch[i][0]);
            for (int j = 1; j < sequence_length; j++) {
                fprintf(file, "%d,", batch[i][j]);
            }
            fprintf(file, "\n");
        }
        pthread_mutex_unlock(file_mutex);
    }

    return NULL;
}

int main() {
    char protein_sequence[MAX_SEQUENCE_LENGTH];
    printf("Enter protein sequence: ");
    scanf("%s", protein_sequence);

    int sequence_length = 0;
    while (protein_sequence[sequence_length] != '\0') {
        sequence_length++;
    }

    FILE *file = fopen("foldings.csv", "w");
    if (file == NULL) {
        printf("Error opening file!\n");
        return 1;
    }

    pthread_mutex_t file_mutex;
    pthread_mutex_init(&file_mutex, NULL);

    int num_threads = 16;  // Aantal threads (afhankelijk van je CPU)
    pthread_t threads[num_threads];
    ThreadArgs thread_args[num_threads];

    // Start threads
    for (int t = 0; t < num_threads; t++) {
        thread_args[t].thread_id = t;
        thread_args[t].sequence_length = sequence_length;
        thread_args[t].file = file;
        thread_args[t].file_mutex = &file_mutex;

        pthread_create(&threads[t], NULL, generate_foldings, &thread_args[t]);
    }

    // Wacht tot alle threads klaar zijn
    for (int t = 0; t < num_threads; t++) {
        pthread_join(threads[t], NULL);
    }

    pthread_mutex_destroy(&file_mutex);
    fclose(file);

    printf("All valid foldings have been saved to 'foldings.csv'.\n");
    return 0;
}
