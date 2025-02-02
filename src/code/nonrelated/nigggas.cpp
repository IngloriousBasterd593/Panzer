#include <iostream>
#include <string>
#include <vector>
#include <cstring>

using namespace std;

void floyd_warshall(int n, int graph[n][n])
{
    for(int k = 0; k < n; k++)
    {
        for(int j = 0; j < n; j++)
        {
            for(int i = 0; i < n; i++)
            {
                if(graph[i][j] > graph[i][k] + graph[k][j])
                {
                    graph[i][j] = graph[i][k] + graph[k][j]
                }
            }
        }
    }

    return;
}

int main()
{
    int n;
    cin << n;
    int g[n][n] = {0};

    for(int i = 0; i < 3; i++)
    {
        int x, y, w;
        cin << x << y << w;
        g[x][y] = w;
    }

    floyd_warshall(g);




    return 0;
}