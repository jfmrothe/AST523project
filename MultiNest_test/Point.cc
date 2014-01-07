#include "Point.h"

Point::Point(int Dim) : myparams(Dim)
{
    D=Dim;
    u = gsl_vector_alloc(D);
    theta = gsl_vector_alloc(D);
}

Point::~Point()
{
    gsl_vector_free(u);
    gsl_vector_free(theta);
}

void Point::set_params(vector<string> param_names, vector<string> prior_types, double min_vals[], double max_vals[])
{
    for(int i = 0; i<D; i++)
    {
        myparams[i].name = param_names[i];
        myparams[i].prior = prior_types[i];
        myparams[i].min = min_vals[i];
        myparams[i].max = max_vals[i];
    }
}

void Point::hypercube_prior()
{
    for(int i = 0; i<D; i++) {gsl_vector_set(u, i, UNIFORM);}
}

void Point::transform_prior()
{
    double r;
    double min, max;

    for(int i = 0; i<D; i++)
    {
        r = gsl_vector_get(u, i);
        min = myparams[i].min;
        max = myparams[i].max;
        
        if(myparams[i].prior == "uniform") {gsl_vector_set(theta, i, min + r*(max-min));}
    }
}
