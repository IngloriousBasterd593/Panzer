#include <stdio.h>
#include <limits.h> 
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

#define BUFFLEN 512
#define MAXNUMCOUNT 256

int maxOfArray(int* arr, int len) {
    int max = arr[0];

    for(int i = 0; i < len; i++) {
        if(i + 1 >= len) {
            return max;
        }

        if(arr[i + 1] > arr[i]) {
            max = arr[i + 1];
        } else {
            max = arr[i];
        }
    }

    return max;
}

int sumOfArray(int* arr, int len) {
    int sum = 0;

    for(int i = 0; i < len; i++) {
        sum += arr[i];
    }

    return sum;
}

int parse(char* buff, int* arr, char delimeter, int maxNumSize) {
    if(buff == NULL) {
        return EXIT_FAILURE;
    }

    int lenght = 1;
    int len = strlen(buff);

    for(int i = 0; i < len; i++) {
        if(buff[i] == delimeter) {
            lenght++;
        }
    }
    
    char num[maxNumSize];
    int numOfInts = 1;
    int number = 0;
    int numCount = 0;
    int numSize = 0;
    int power = 0;

    memset(num, '\0', maxNumSize);

    for(int i = 0; i < len; i++) {
        if(i + 1 >= len) {
            for(int j = numSize - 1; j >= 0; j--) {
                number += (num[j] - 48) * ((int)pow(10, power));
                power++;
            }

            arr[numCount] = number;
    
            return lenght;
        }

        if(buff[i] != delimeter) {
            num[numSize] = buff[i];
            numSize++;
        }

        if(buff[i] == delimeter) {
            for(int j = numSize - 1; j >= 0; j--) {
                number += (num[j] - 48) * ((int)pow(10, power));
                power++;
            }

            memset(num, '\0', maxNumSize);
            arr[numCount] = number;
            power = 0;
            number = 0;
            numCount++;
            numSize = 0;
        }
    }

    return EXIT_SUCCESS;
}

int maxSubarraySum(int arr[], int n) {
    int maxSum = INT_MIN; 
    int currSum = 0;   

    for(int i = 0; i < n; i++) {
        currSum += arr[i];
        if(currSum > maxSum) {
            maxSum = currSum; 
        }
        if(currSum < 0) {
            currSum = 0; 
        }
    }

    return maxSum;
}

int main() {
    // 1 uzdevums
    int size;
    printf("Ievadiet masīva izmēru: ");
    scanf("%d", &size);

    int arr[size];
    printf("Ievadiet %d veselos skaitļus:\n", size);
    for(int i = 0; i < size; i++) {
        scanf("%d", &arr[i]);
    }

    int result = maxSubarraySum(arr, size);
    printf("Maksimālā summa: %d\n", result);

    // 2 uzdevums
    int n;

    printf("Ievadiet balonu skaitu: ");
    scanf("%d", &n);

    int balloons[n];
    printf("Ievadiet balonu augstumus:\n");
    for(int i = 0; i < n; i++) {
        scanf("%d", &balloons[i]);
    }

    int arrows[1000] = {0}; 
    int arrow_count = 0;    

    for(int i = 0; i < n; i++) {
        int height = balloons[i];
        int shot = 0;

        for(int j = 0; j < arrow_count; j++) {
            if(arrows[j] == height) {
                arrows[j]--; 
                shot = 1;    
                break;
            }
        }
        
        if(!shot) {
            arrows[arrow_count++] = height - 1;
        }
    }

    printf("Nepieciešamais bultu skaits: %d\n", arrow_count);

    return 0;
}
