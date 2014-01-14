/***************************************************

File: Point.cc

Description:
Source code for the Point class. 

To do:
Add more prior distribution transformations. 

Programmers: Johnny Greco & Johannes Rothe
Contacts: greco@princeton.edu, jrothe@princeton.edu

****************************************************/
#include "Point.h"

Point::Point(int Dim) : myparams(Dim) {D=Dim;}

void Point::set_params(vector<string> prior_types, vector<double> min_vals, vector<double> max_vals)
// set parameter values needed for prior distributions
{

    for(int i = 0; i<D; i++)
    {
        myparams[i].prior = prior_types[i];
        myparams[i].min = min_vals[i];
        myparams[i].max = max_vals[i];
    }
}

void Point::hypercube_prior()
// set hypercube coordinates -> uniform(0,1)
{
    for(int i = 0; i<D; i++) {myparams[i].u = UNIFORM;}
}

void Point::transform_prior()
// transform hypercube prior to physical prior
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
