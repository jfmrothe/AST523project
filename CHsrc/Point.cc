#include "Point.h"

Point::Point(int Dim) : myparams(Dim) {D=Dim;}

void Point::set_params(string const& prior_types, vector<double> min_vals, vector<double> max_vals)
{
    for(int i = 0; i<D; i++)
    {
        myparams[i].prior = prior_types;
        myparams[i].min = min_vals[i];
        myparams[i].max = max_vals[i];
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

void Point::get_theta(double *theta,int nt){
  if(nt!=D){
    cout<<"Point::get_theta:the dimention of parameters are not correct\n"<<endl;
    exit(1);
  }
  for (int i=; i<D; i++){
    theta[i] = myparams[i].theta;
  }
  return
}

int Point::get_D(){return D;}
