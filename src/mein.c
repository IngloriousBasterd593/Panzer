#include <stdio.h>
#include <string.h>

#define MAX_LENGTH 1000

int main() {
    char virkne[MAX_LENGTH];
    char samazinata[MAX_LENGTH];
    int count = 1;
    int k = 0;

    // Lietotāja ievade
    printf("Ievadiet simbolu virkni: ");
    if (scanf("%s", virkne) != 1) {
        printf("Nepareiza ievade!\n");
        return 1;
    }

    int len = strlen(virkne);

    for(int i = 1; i <= len; i++) {
        if (virkne[i] == virkne[i-1]) {
            count++;
        } else {
            samazinata[k++] = virkne[i-1];
            samazinata[k++] = count + '0';
            count = 1;
        }
    }

    samazinata[k] = '\0';

    // Salīdzināšana
    if (strlen(samazinata) < len) {
        printf("Samazinātā virkne: %s\n", samazinata);
    } else {
        printf("Samazināšana nav iespējama. Oriģinālā virkne: %s\n", virkne);
    }

    return 0;
}
