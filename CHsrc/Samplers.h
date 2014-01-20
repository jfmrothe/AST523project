#ifndef SAMPLERS_H
#define SAMPLERS_H
#include<string>
//#include "Data.h"
#include "Ellipsoid.h"
using namespace std;

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
        int D_;
	      int N;
        double Vtot;
        double e_;
	      double logZ; 
	      double H;
        double logLmax_,logLmin_;
        gsl_vector * newcoor_;
        int ellworst_,ptworst_,ellnew_;
        vector<Ellipsoid *> clustering;
	      list<Point *> discard_pts;
    public:
        //Samplers(int Dim, int Npts, vector<string> prior_types, vector<double> min_vals, vector<double> max_vals, double eff);
        Samplers(double* min_vals, int nmin,double *max_vals, int nmax, double eff, int Npts, string const& prior_types);
	      ~Samplers();
        void DrawSample();
        void SetAllPoint(double * logL,int nL);
	      void DisgardWorstPoint(int nest);
        void ResetWorstPoint(double *theta, int nt);
        void ResetWorstPointLogL(double logL);
	      //void printout();
	      void EllipsoidalPartitioning(vector<Point *>& pts, double Xi);
        int get_NumEll() {return clustering.size();}
        void CalcVtot(); 

        double ClusteringQuality(double Xtot) {return Vtot/Xtot;} 
        void ClearCluster();
	      void EraseFirst() {clustering.erase(clustering.begin());}
        void getAlltheta(double *Alltheta, int nx, int ny);
	      double getlogLmax(){return logLmax_;}
	      double getlogLmin(){return logLmin_;}
        void Recluster(double X_i);
        int getN() {return N;}
        int countTotal();
        void getlogZ(double *logzinfo, int nz);
        void getPosterior(double * posterior, int nx, int ny, double * prob, int np);
        void EllipsoidalRescaling(double Xi);
};
#endif
