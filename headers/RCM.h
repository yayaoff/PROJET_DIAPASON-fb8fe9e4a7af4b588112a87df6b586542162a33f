typedef struct {
    int node;
    int level;
} Level;

typedef struct {
    double x;
    double y;
} Point;

double distance(Point p1, Point p2) ;
Point *coord_to_point(double * coord, int n_nodes);
int** AdjacencyMatrix(double** matrix, int size);
void freeAdjacencyMatrix(int** adjMatrix, int size);
int compare(const void* a, const void* b) ;


void reverseCuthillMcKee(int** adjMatrix, int size, int* perm);