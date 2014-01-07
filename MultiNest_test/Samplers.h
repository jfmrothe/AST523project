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
        vector<Ellipsoid *> clustering;

    public:
        Samplers(int Dim) {D=Dim;} 

        gsl_vector * get_newcoor();
        void EllipsoidalPartitioning(vector<Point *>&, double);
        int get_NumEll() {return clustering.size();}
        void mcmc(Point*, Data, double);
        void CalcVtot(); 
        void ClearCluster() {clustering.clear();}
};
#endif
