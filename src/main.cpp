#include <iostream>
#include <vector>
#include <climits>  

void printMatrix(int m, int n, const std::vector<std::vector<int>>& arr) 
{
    for (int i = 0; i < m; i++) 
    {
        for (int j = 0; j < n; j++) 
        {
            std::cout << arr[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

int sumMatrix(int m, int n, const std::vector<std::vector<int>>& arr) 
{
    int sum = 0;
    for (int i = 0; i < m; i++) 
    {
        for (int j = 0; j < n; j++) 
        {
            if(arr[i][j] == INT_MAX)
            {
                continue;
            }

            sum += arr[i][j];
        }
    }
    return sum;
}

int algo(int n, std::vector<std::vector<int>> arr) 
{
    for (int k = 0; k < n; k++) 
    {
        for (int i = 0; i < n; i++) 
        {
            for (int j = 0; j < n; j++) 
            {
                if(arr[i][k] + arr[k][j] < arr[i][j]) 
                {
                    arr[i][j] = arr[i][k] + arr[k][j];
                }
            }
        }
    }

    return sumMatrix(n, n, arr);
}

void deleteRowCol(int rc, int n, std::vector<std::vector<int>> &arr) 
{
    for (int i = 0; i < n; i++) 
    {
        arr[rc][i] = INT_MAX; 
        arr[i][rc] = INT_MAX; 
    }
}

int main() 
{
    int n;
    std::cin >> n;

    std::vector<std::vector<int>> arr(n, std::vector<int>(n, 0));
    



    std::vector<int> v(n);

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            int x;
            std::cin >> x;
            arr[i][j] = x;  
        }
    }

   
    for (int i = 0; i < n; i++)
    {
        std::cin >> v[i];
        v[i]--;
    }

   
    for (int i = 0; i < n; i++)
    {
        deleteRowCol(v[i], n, arr);
        std::cout << algo(n, arr) << std::endl;

    }
    
    return 0;
}