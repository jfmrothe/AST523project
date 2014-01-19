#ifndef SAMPLERS_H
#define SAMPLERS_H

#include "Data.h"
#include "Ellipsoid.h"

// **** Samplers Utilities 
double dist(int, Point *, gsl_vector *);
int KMeans(vector<Point *>&, int, int, int *);
Ellipsoid FindEnclosingEllipsoid(vector<Point *>&, int);
void SelectFromGrouping(vector<Point *>&, int, int *, int, vector<Point *>&);
bool u_in_hypercube(gsl_vector *, int);
// ****

class Samplers
{
    private:
        int D;
	int N;
        double Vtot;
        double e;
	double logZ; 
	double H;

        vector<Ellipsoid *> clustering;
	list<Point *> discard_pts;
	gsl_vector * newcoor_;
	FILE * file3;
    public:
        Samplers(int Dim, int Npts, Data data_obj, vector<string> prior_types, vector<double> min_vals, vector<double> max_vals, double eff);
	~Samplers();
        int DrawSample();
	double ResetWorstPoint(int nest, Data data_obj, int * nLeval);
	void printout();
	void EllipsoidalPartitioning(vector<Point *>& pts, double Xi);
	void EllipsoidalRescaling(double Xi);
        int get_NumEll() {return clustering.size();}
        void mcmc(Point*, Data, double);
        void CalcVtot(); 
        double ClusteringQuality(double Xtot) {return e*Vtot/Xtot;} 
        void ClearCluster();
	void EraseFirst() {clustering.erase(clustering.begin());}
	void printClustering();

	int getN() {return N;}
	double get_logZ(){return logZ;}
	gsl_vector * get_newcoor() {return newcoor_;}
};
#endif
