#ifndef POINT_H
#define POINT_H
#include <iostream>
#include <list>
#include <fstream>
#include <string>
#include <float.h>
#include <math.h>
#include <vector>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
using namespace std;

#define THRESH 1.0e-5
#define UNIFORM ((rand() + 0.5)/(RAND_MAX+1.0))
#define PLUS(x,y) (x>y ? x+log(1+exp(y-x)) : y+log(1+exp(x-y)))

class Point
{
    private:
        int D; 
        double logL;
        double logWt;
        struct Parameter
        {
            std::string prior;
            double u;
            double theta;
            double min;
            double max;
        };
        vector<Parameter> myparams;

    public:

        Point(int);
        Point(Point &pt);
        void get_params(string & prior_types, double *min_vals, double *max_vals);
        void set_params(string const& prior_types,double* min_vals,double* max_vals); //use the same prior_type right now
        void set_u_single(int i, double x){myparams[i].u = x;}
        //void set_u( * coor, int nt){for(int i = 0; i < D; i++) {myparams[i].u = coor[i];}}
        void set_u(gsl_vector * coor){for(int i = 0; i < D; i++) {myparams[i].u = gsl_vector_get(coor,i);}}
        void set_theta(int i, double x){myparams[i].theta = x;}
        void get_theta(double *theta, int nt);
        double get_u(int i){return myparams[i].u;}
        void set_logL(double L){logL = L;}
        double get_logL(){return logL;}
        void set_logWt(double W){logWt = W;}
        double get_logWt(){return logWt;}
        int get_D();
        void hypercube_prior();
        void transform_prior();
};
#endif
