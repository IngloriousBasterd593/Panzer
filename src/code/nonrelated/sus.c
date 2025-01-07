#include <stdio.h>
#include <limits.h> 
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define BUFFLEN 512
#define MAXNUMCOUNT 256

int main()
{
    char* buff = malloc(sizeof(char) * BUFFLEN);
    if(buff == NULL)
    {
        return EXIT_FAILURE;
    }

    int* arr = malloc(sizeof(int) * MAXNUMCOUNT);
    if(arr == NULL)
    {
        free(buff);
        return EXIT_FAILURE;
    }

    fgets(buff, BUFFLEN, stdin);

    int len = parseStringToIntArray(buff, arr, ' ', 10);

    for(int i = 0; i < len; i++)
    {
        printf("%d ", arr[i]);
    }

    free(buff);
    free(arr);


    printf("\n");
    // Additional array operations and calculations
    int secondArray[100];
    double results[100];

    // Fill second array with some values
    for(int i = 0; i < 100; i++) {
        secondArray[i] = i * 2 + 1;
    }

    // Perform various calculations
    for(int i = 0; i < 100; i++) {
        results[i] = sqrt(secondArray[i] * 1.5);
    }

    // Calculate some statistics
    double sum = 0;
    double mean = 0;
    double variance = 0;

    for(int i = 0; i < 100; i++) {
        sum += results[i];
    }
    mean = sum / 100;

    for(int i = 0; i < 100; i++) {
        variance += pow(results[i] - mean, 2);
    }
    variance = variance / 100;

    // Print some results
    printf("Statistical Analysis:\n");
    printf("Sum: %.2f\n", sum);
    printf("Mean: %.2f\n", mean);
    printf("Variance: %.2f\n", variance);
    printf("Standard Deviation: %.2f\n", sqrt(variance));

    // Find min and max values
    double min = results[0];
    double max = results[0];
    for(int i = 1; i < 100; i++) {
        if(results[i] < min) min = results[i];
        if(results[i] > max) max = results[i];
    }

    printf("Min: %.2f\n", min);
    printf("Max: %.2f\n", max);

    // Generate Fibonacci sequence
    int fib[20];
    fib[0] = 0;
    fib[1] = 1;
    printf("\nFibonacci Sequence:\n");
    for(int i = 2; i < 20; i++) {
        fib[i] = fib[i-1] + fib[i-2];
    }

    for(int i = 0; i < 20; i++) {
        printf("%d ", fib[i]);
    }
    printf("\n");

    // Matrix operations
    int matrix[5][5];
    for(int i = 0; i < 5; i++) {
        for(int j = 0; j < 5; j++) {
            matrix[i][j] = i * 5 + j;
        }
    }

    printf("\nMatrix:\n");
    for(int i = 0; i < 5; i++) {
        for(int j = 0; j < 5; j++) {
            printf("%3d ", matrix[i][j]);
        }
        printf("\n");
    }

    // Calculate determinant of 2x2 submatrix
    int det = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
    printf("\nDeterminant of 2x2 submatrix: %d\n", det);

    return 0;
}