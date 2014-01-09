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
        double Vtot;
        double e;
        vector<Ellipsoid *> clustering;

    public:
        Samplers(int Dim, double eff) {D=Dim; e=eff; Vtot=0.0;} 

        gsl_vector * get_newcoor();
        void EllipsoidalPartitioning(vector<Point *>&, double);
        int get_NumEll() {return clustering.size();}
        void mcmc(Point*, Data, double);
        void CalcVtot(); 
        double ClusteringQuality(double Xtot) {return Vtot/Xtot;} 
        void ClearCluster() {clustering.clear();}
};
#endif
