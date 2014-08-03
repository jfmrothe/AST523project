#ifndef SAMPLERS_H
#define SAMPLERS_H
#include<string>
//#include "Data.h"
#include "Ellipsoid.h"
#include "Point.h"
#include <time.h>
using namespace std;



class Samplers
{
    private:
        struct Ran{
  	   unsigned long long int u,v,w;
  	   //Ran(unsigned long long int j) : v(4101842887655102017LL),w(1){
           //u = j^v; int64();
           //w = v; int64();
           //}
           inline unsigned long long int int64() {
           u = u*2862933555777941757LL + 7046029254386353087LL;
           v ^=v >> 17; v^=v<<31; v^=v >>8;
           w = 4294957665U*(w & 0xffffffff) + (w>>32);
           unsigned long long int x = u^(u<<21); x^=x>>35; x^=x<<4;
              return (x +v) ^w;
           }
           inline double doub(){return 5.42101086242752217E-20 * int64();}
          inline unsigned int int32() {return (unsigned int) int64();}
        };

	Ran myrand;
	int D_;
	int N_;
        double Vtot;
        double e_;
	double logZ; 
	double H;
        double logLmax_,logLmin_;
        gsl_vector * newcoor_;
        int ellworst_,ptworst_,ellnew_;
        vector<Ellipsoid *> clustering;
	list<Point *> discard_pts;
	// **** Samplers Utilities 
	double dist(int, Point *, gsl_vector *);
	int KMeans(vector<Point *>&, int, int, int *);
	Ellipsoid FindEnclosingEllipsoid(vector<Point *>&, int);
	void SelectFromGrouping(vector<Point *>&, int, int *, int, vector<Point *>&);
	bool u_in_hypercube(gsl_vector *, int);

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

        double ClusteringQuality(double Xtot) {return e_*Vtot/Xtot;} 
        void ClearCluster();
	void EraseFirst(); 
        void getAlltheta(double *Alltheta, int nx, int ny);
	double getlogLmax(){return logLmax_;}
	double getlogLmin(){return logLmin_;}
	void FullRecluster(double X_i);
        int Recluster(double X_i, double qualthresh);
        int getN() {return N_;}
        int countTotal();
        void OutputClusters();
	void getlogZ(double *logzinfo, int nz);
        void getPosterior(double * posterior, int nx, int ny, double * prob, int np);
        void EllipsoidalRescaling(double Xi);

	int TestUniformity() {return 0;}
};
#endif
