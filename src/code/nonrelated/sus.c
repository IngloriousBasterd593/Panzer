#include <stdio.h>
#include <stdlib.h>

typedef struct matL
{
    int sizeX;
    int sizeY;
    long** m;
} matL;

matL nullMat = {0, 0, {{0}}};

matL createMatrix(int sizeX, int sizeY)
{
    matL res;

    res.m = (long**) malloc(sizeX * sizeof(long*));

    for(int i = 0; i < sizeY; i++)
    {
        res.m[i] = (long*) malloc(sizeY * sizeof(long));
    }
    
    for(int i = 0; i < sizeX; i++)
    {
        for(int j = 0; j < sizeY; j++)
        {
            res.m[i][j] = 0;
        }
    }

    res.sizeX = sizeX;
    res.sizeY = sizeY;

    return res;
}

matL multiplyMatrices(matL A, matL B)
{
    if(A.sizeY != B.sizeX)
    {
        printf("ahh");
        return nullMat;
    }

    matL res = createMatrix(A.sizeX, B.sizeY);

    for(int i = 0; i < A.sizeX; i++)
    {
        for(int j = 0; j < B.sizeY; j++)
        {
            res.m[i][j] = 0;
            for(int k = 0; k < A.sizeY; k++)
            {
                res.m[i][j] += A.m[i][k] * B.m[k][j];
            }
        }
    }

    return res;
}

void populateMatrix(matL* mat, long* data)
{
    for(int i = 0; i < mat->sizeX; i++)
    {
        for(int j = 0; j < mat->sizeY; j++)
        {
            mat->m[i][j] = data[i * mat->sizeY + j];
        }
    }
}

long* getMatrixInput()
{
    int n;
    scanf("%d", &n);
    long* res = (long*) malloc(n * sizeof(long));

    for(int i = 0; i < n; i++)
    {
        scanf("%li", &res[i]);
    }

    return res;

}

void printMatrix(matL mat)
{
    for(int i = 0; i < mat.sizeX; i++)
    {
        for(int j = 0; j < mat.sizeY; j++)
        {
            printf("%li ", mat.m[i][j]);
        }
        printf("\n");
    }
}

int main()
{
    int x, y;

    printf("Enter the size of the matrices: ");

    scanf("%d %d", &x, &y);

    matL A = createMatrix(x, y);
    matL B = createMatrix(x, y);

    // printMatrix(A);

    long* data = getMatrixInput();

    populateMatrix(&A, data);

    free(data);

    data = getMatrixInput();

    populateMatrix(&B, data);

    free(data);

    printMatrix(A);

    matL res = multiplyMatrices(A, B);

    printMatrix(res);






    return 0;
}