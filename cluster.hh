
# ifndef _AST523_CLUSTER_
# define _AST523_CLUSTER_

float dist(int D, float * pos1, float * pos2);
int indexofSmallestElement(float * array, int size);
int KMeans(int D, int N, int k, float * points, int * grouping);
int TestKMeans();

# endif
