#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../headers/RCM.h"

int** AdjacencyMatrix(double** matrix, int size) {
    int** adjMatrix = malloc(size * sizeof(int*));

    for (int i = 0; i < size; i++) {
        adjMatrix[i] = calloc(size, sizeof(int));

        for (int j = 0; j < size; j++) {
            if (matrix[i][j] != 0.0) {
                adjMatrix[i][j] = 1;
            }
        }
    }

    return adjMatrix;
}

void freeAdjacencyMatrix(int** adjMatrix, int size) {
    for (int i = 0; i < size; i++) {
        free(adjMatrix[i]);
    }

    free(adjMatrix);
}

typedef struct Node {
    int index;
    int degree;
} Node;

int compare(const void* a, const void* b) {
    return ((Node*)a)->degree - ((Node*)b)->degree;
}

void reverseCuthillMcKee(double** matrix, int size, int* perm) {
    Node* nodes = (Node*)malloc(size * sizeof(Node));
    int* queue = (int*)malloc(size * sizeof(int));
    int* visited = (int*)calloc(size, sizeof(int));

    // Calculate degrees
    //Based adjacency on wether or not the value between 2 coordinates is not null
    for (int i = 0; i < size; i++) {
        nodes[i].index = i;
        nodes[i].degree = 0;
        for (int j = 0; j < size; j++) {
            if (matrix[i][j] != 0.0) {
                nodes[i].degree++;
            }
        }
    }

    // Sort nodes by degree
    qsort(nodes, size, sizeof(Node), compare);

    // BFS
    int start = 0;
    int end = 1;
    queue[0] = nodes[0].index;
    visited[nodes[0].index] = 1;

    while (start < end) {
        Node* neighbors = (Node*)malloc(size * sizeof(Node));
        int count = 0;

        // Find all unvisited neighbors
        for (int i = 0; i < size; i++) {
            if (matrix[queue[start]][i] != 0.0 && !visited[i]) {
                neighbors[count].index = i;
                neighbors[count].degree = nodes[i].degree;
                count++;
            }
        }

        // Sort unvisited neighbors by degree
        qsort(neighbors, count, sizeof(Node), compare);

        // Add neighbors to the queue
        for (int i = 0; i < count; i++) {
            queue[end++] = neighbors[i].index;
            visited[neighbors[i].index] = 1;
        }

        free(neighbors);
        start++;
    }

    // Generate the perm array
    for (int i = 0; i < size; i++) {
        perm[i] = queue[size - i - 1];
    }

    free(nodes);
    free(queue);
    free(visited);
}