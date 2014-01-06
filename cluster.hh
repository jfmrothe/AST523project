
# ifndef _AST523_CLUSTER_
# define _AST523_CLUSTER_

double dist(int D, gsl_vector * pos1, gsl_vector * pos2);
int KMeans(gsl_vector ** points, int D, int N, int k, int * grouping, int * groupsize);
int TestKMeans();
double ClusteringQuality(std::vector<Ellipsoid*>& clustering, double Xtot);
int TestClusteringQuality();
void EllipsoidalPartitioning(gsl_vector ** coors, int D, int N, double Xtot, std::vector<Ellipsoid*>& clustering);
int TestEllipsoidalPartitioning();
# endif
