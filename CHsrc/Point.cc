#include "Point.h"

Point::Point(int Dim):myparams(Dim){
  D=Dim;
  logL=0;
  logWt=0;
}
Point::Point(Point &pt):myparams(pt.get_D()){
  int tempD=pt.get_D();
  D = tempD;
  set_logL(pt.get_logL());
  set_logWt(pt.get_logWt());
  double *theta = new double [D];
  for (int i =0; i< D; i++){
    double tempu=0;
    tempu=pt.get_u(i);
    myparams[i].u=tempu;
  }
  pt.get_theta(theta,D);
  for (int i=0; i<D; i++){
    set_theta(i,theta[i]);
  }
  string prior_types;
  double *min_vals = new double[D];
  double *max_vals = new double [D];
  pt.get_params(prior_types,min_vals,max_vals);
  for(int i = 0; i<D; i++)
  {
        myparams[i].prior = prior_types;
        myparams[i].min = min_vals[i];
        myparams[i].max = max_vals[i];
  }
  delete [] theta;
  delete [] min_vals;
  delete [] max_vals;
}
void Point::get_params(string & prior_types, double *min_vals, double *max_vals){
    prior_types = myparams[0].prior;
    for(int i = 0; i<D; i++)
    {
        min_vals[i] = myparams[i].min; 
        max_vals[i] = myparams[i].max; 
    }
}
void Point::set_params(string const& prior_types, double *min_vals, double *max_vals)
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
    printf("Point::get_theta:the dimention of parameters are not correct\n");
  }
  for (int i=0; i<D; i++){
    //printf("i=%d,myparams[i]=%f\n",i,myparams[i].theta);
    theta[i] = myparams[i].theta;
  }
  return;
}

int Point::get_D(){return D;}
