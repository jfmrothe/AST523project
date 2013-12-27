#ifndef ELLIPSOIDAL_SAMP_H
#define ELLIPSOIDAL_SAMP_H
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include "Point.h"
#include <vector>
#include <cmath>

float uniform(float min, float max);
float boxmuller();
float quadr();
void unisphere(float * coor, int D);
void SampleEllipsoid(int D, gsl_vector * center, gsl_matrix * C, double f, gsl_vector * coor);
void FindEnclosingEllipsoid(int D, int N, Point *pts[],  gsl_vector * center, gsl_matrix * C, double * f);
#endif
