#ifndef SAMPLERS_H
#define SAMPLERS_H

#include "Data.h"

float boxmuller();
float quadr();
void unisphere(float * coor, int D);

class Samplers
{
    private:
        int D;
        double f;
        gsl_vector * center;
        gsl_vector * coor;
        gsl_matrix * C;

    public:
        Samplers(int); 
        ~Samplers();

        void set_vectors_zero();
        double get_f() {return f;}
        void set_f_factor(double x) {f = f*x;}
        gsl_vector * get_coor() {return coor;}
        void SampleEllipsoid();
        void FindEnclosingEllipsoid(int N, Point *pts[]);
        void mcmc(Point*, Data, double);
};
#endif
