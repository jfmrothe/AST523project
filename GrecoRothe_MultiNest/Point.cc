#include "Point.h"

void Point::set_params(vector<string> param_names, vector<string> prior_types, double min_vals[], double max_vals[], double guess_vals[])
{
    for(int i = 0; i<D; i++)
    {
        myparams[i].name = param_names[i];
        myparams[i].prior = prior_types[i];
        myparams[i].min = min_vals[i];
        myparams[i].max = max_vals[i];
        myparams[i].init_guess = guess_vals[i];
    }
}

void Point::hypercube_prior()
{
    for(int i = 0; i<D; i++) {myparams[i].u = UNIFORM;}
}

void Point::transform_prior()
{
    double r;
    double min, max;

    for(int i = 0; i<D; i++)
    {
        r = myparams[i].u;
        min = myparams[i].min;
        max = myparams[i].max;
        
        if(myparams[i].prior == "uniform") {myparams[i].theta = min + r*(max-min);}
    }
}
