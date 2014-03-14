/***************************************************

File: Data.cc 

Description:
Source code for the methods of the Data class.

To do:
At the moment, the likelihood evaluations are quite
suboptimal; the likelihood function is determined by 
the file name given by the user. This will be changed 
to a pointer to a user defined function that is passed
to the Data constructor at run time. 

Programmers: Johnny Greco & Johannes Rothe
Contacts: greco@princeton.edu, jrothe@princeton.edu

****************************************************/
#include "Data.h"

Data::Data(int Dim, int n_col, string filename) : data(n_col)
{
    D = Dim;
    num_col = n_col;
    data_filename = filename;
}

void Data::get_data()
// Get the data for likelihood evaluation. Currently, 
// only the lighthouse problem uses this method.
{
    double dat;
    string tmpname = "DataFiles/"+data_filename;
    ifstream datafile(tmpname.c_str());
    while(!datafile.eof())
    {
        for(int i=0; i<num_col; i++)
        {
            datafile >> dat;
            data[i].push_back(dat);
        }
    }
    datafile.close();
}

void Data::Get_L(double *theta, double *logLarr, int nl)
// Evaluate likelihood function. In a future version, this will become 
// a very general function, which handles a user defined likelihood. 
{
    for (int i =0; i<nl; i++){
      if(data_filename == "lighthouse.dat") // lighthouse likelihood
      {
          double x, y;
          double logL = 0;
          x = theta[D*i];
          y = theta[D*i+1];
          for(datum=data[0].begin(); datum!=data[0].end(); datum++)
          {
              logL += log((y/M_PI) / ((*datum-x)*(*datum-x) + y*y));
          }
          logLarr[i]=logL;
      }
      else if (data_filename == "eggbox") // egg-box likelihood
      {
          double x = 1.0;
          for(int j = 0; j < D; j++) {x *= cos(theta[D*i+j]/2.0);}
          logLarr[i]=pow(2.0 + x, 5.0);

      }
      else if (data_filename == "gauss_shells") // Gaussian shells likelihood
      {
          vector<double> c1(D), c2(D);
          double r = 2.0;
          double ww = 0.1*0.1;
          double logL, arg1, arg2;
          double dist1_sq = 0.0, dist2_sq = 0.0;
          c1[0] = -3.0;
          c2[0] = 3.0;
          double Norm = 1.0/sqrt(2*M_PI*ww);
          for(int j = 1; j<D; j++){c1[j] = 0.0;}

          for(int j = 0; j<D; j++)
          {
              dist1_sq += pow(theta[D*i+j] - c1[i], 2);
              dist2_sq += pow(theta[D*i+j] - c2[i], 2);
          }

          arg1 = -pow(sqrt(dist1_sq) - r, 2)/(2*ww);
          arg2 = -pow(sqrt(dist2_sq) - r, 2)/(2*ww);

          logL = log(Norm) + PLUS(arg1, arg2);

          logLarr[i]=logL;
      }
      else{cout << "cannot find likelihood function!" << endl;}
    }
}
