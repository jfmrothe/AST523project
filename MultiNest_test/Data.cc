#include "Data.h"

void Data::get_data(string filename)
{
    double dat;
    ifstream datafile(filename.c_str());

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

void Data::get_results(list<Point *>& samples, double logZ)
{
    vector<double> x(D, 0.0);
    vector<double> xx(D, 0.0);

    double w;
    list<Point *>::iterator s;
    for(s=samples.begin(); s!=samples.end(); s++)
    {
      w = exp((*s)->get_logWt() - logZ);
        for(int i = 0; i<D; i++)
        {
            x[i] += w*((*s)->get_theta(i));
            xx[i] += w*((*s)->get_theta(i))*((*s)->get_theta(i));
        }
    }

    for(int i = 0; i<D; i++)
    {
        cout << "<x" << i+1 << "> = " << x[i] << " +/- " << sqrt(xx[i] - x[i]*x[i]) << endl;
    }
}

void Data::lighthouse_logL(Point* pt)
{
    double x, y;
    double logL = 0;
    x = pt->get_theta(0);
    y = pt->get_theta(1);
    for(datum=data[0].begin(); datum!=data[0].end(); datum++)
    {
        logL += log((y/M_PI) / ((*datum-x)*(*datum-x) + y*y));
    }
    pt->set_logL(logL);
}

void Data::eggbox_logL(Point* pt)
{
    double x = 1.0;
    for(int i = 0; i < D; i++) {x *= cos(pt->get_theta(i)/2.0);}
    pt->set_logL(pow(2.0 + x, 5.0));
}

