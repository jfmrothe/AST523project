
# ifndef _AST523_CLUSTER_
# define _AST523_CLUSTER_

double dist(int D, double * pos1, double * pos2);
double gsldist(int D, gsl_vector * pos1, gsl_vector * pos2);
int indexofSmallestElement(double * array, int size);
int KMeans(int D, int N, int k, double * points, int * grouping);
int TestKMeans();
int gslKMeans(gsl_vector ** points, int D, int N, int k, int * grouping, int * groupsize);
int TestgslKMeans();
double ClusteringQuality(int D, int N, int i, Ellipsoid * clustering, int Nell);
int TestClusteringQuality();
Ellipsoid ** EllipsoidalPartitioning(gsl_vector ** coors, int D, int N, int * Nell);
void EllipsoidalPartitioningNoEll(gsl_vector ** coors, int D, int N, int * Nell, gsl_vector ** centers, gsl_matrix ** covMats, double * enlFacs);
void EllipsoidalPartitioningVec(gsl_vector ** coors, int D, int N, double Xtot, std::vector<Ellipsoid*>& clustering);
int TestEllipsoidalPartitioning();
int TestEllipsoidalPartitioningNoEll();
int TestEllipsoidalPartitioningVec();
int JustEllipsoidalPartitioning();
# endif
