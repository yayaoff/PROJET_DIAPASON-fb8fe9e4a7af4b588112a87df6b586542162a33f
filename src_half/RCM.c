#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../headers/RCM.h"

double distance(Point p1, Point p2) {
    double dx = p1.x - p2.x;
    double dy = p1.y - p2.y;
    return sqrt(dx * dx + dy * dy);
}

int compare(const void* a, const void* b) {
    Level* nodeA = (Level*)a;
    Level* nodeB = (Level*)b;
    return nodeA->level - nodeB->level;
}

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

void reverseCuthillMcKee(int** adjMatrix, int size, int* perm) {
    Level* nodes = malloc(size * sizeof(Level));
    int* visited = calloc(size, sizeof(int));
    int queue[size];
    int front = 0;
    int rear = 0;
    int level = 0;

    // Initialize the starting node
    int startLevel = 0;
    nodes[startLevel].node = startLevel;
    nodes[startLevel].level = level;
    visited[startLevel] = 1;
    queue[rear++] = startLevel;

    while (front < rear) {
        int currentLevel = queue[front++];
        level++;

        // Find neighbors of the current node
        for (int i = 0; i < size; i++) {
            if (adjMatrix[currentLevel][i] && !visited[i]) {
                // Check if the neighbor is already in the queue
                int alreadyInQueue = 0;
                for (int j = front; j < rear; j++) {
                    if (queue[j] == i) {
                        alreadyInQueue = 1;
                        break;
                    }
                }

                // Enqueue the neighbor only if it is not already in the queue
                if (!alreadyInQueue) {
                    nodes[i].node = i;
                    nodes[i].level = level;
                    visited[i] = 1;
                    queue[rear++] = i;
                }
            }
        }
    }

    // Sort the nodes based on the levels in increasing order
    qsort(nodes, size, sizeof(Level), compare);

    // Generate the permutation array
    for (int i = 0; i < size; i++) {
        perm[i] = nodes[i].node;
    }

    free(nodes);
    free(visited);
}