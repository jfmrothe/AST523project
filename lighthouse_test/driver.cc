#include <iostream>
#include <list>
#include <fstream>
#include <float.h>
#include <time.h> 
#include "EllSamp.h"

#define MAXITER 10000
#define UNIFORM ((rand() + 0.5)/(RAND_MAX+1.0))
#define PLUS(x,y) (x>y ? x+log(1+exp(y-x)) : y+log(1+exp(x-y)))

// global variables
ifstream infile("lighthouse.dat");
list<double> data;
list<double>::iterator datum;
// ****

// function to calculate log(L)
void logLhood(Point* pt)
{
    double x, y;
    double logL = 0; 
    x = pt->get_theta(0);
    y = pt->get_theta(1);
    for(datum=data.begin(); datum!=data.end(); datum++)
        logL += log((y/M_PI) / ((*datum-x)*(*datum-x) + y*y));
    pt->set_logL(logL);
}

// uniform prior
void prior(Point* pt)
{
    pt->set_u(0, UNIFORM);
    pt->set_u(1, UNIFORM);
    pt->set_theta(0, 4.0*pt->get_u(0) - 2.0);
    pt->set_theta(1, 2.0*pt->get_u(1)); 
}

// mcmc algorithm to find knew point 
void mcmc(Point* pt, double logLmin, int D)
{
    double step = 0.1, new_u, new_v;
    double u = pt->get_u(0);
    double v = pt->get_u(1);
    int accept = 0;
    int reject = 0;
    Point *trial; 
    trial = new Point(D);

    for(int j = 20; j > 0; j--)
    {
        new_u = u + step*(2.0*UNIFORM - 1.0);
        new_v = v + step*(2.0*UNIFORM - 1.0);
        new_u -= floor(new_u);
        new_v -= floor(new_v);
        trial->set_u(0, new_u);
        trial->set_u(1, new_v);
        trial->set_theta(0, 4.0*new_u - 2.0);
        trial->set_theta(1, 2.0*new_v);
        logLhood(trial);

        if(trial->get_logL() > logLmin){*pt = *trial; accept++;}
        else reject++;

        if(accept > reject) step *= exp(1.0 / accept);
        else if(accept < reject) step /= exp(1.0 / reject);
    }

    delete trial;
}

// function to display the results
void results(list<Point> samples, double logZ)
{
    double x = 0.0, xx = 0.0; 
    double y = 0.0, yy = 0.0; 
    double w;
    list<Point>::iterator s;
    for(s=samples.begin(); s!=samples.end(); s++)
    {
        w = exp(s->get_logWt() - logZ);
        x += w*(s->get_theta(0));
        xx += w*(s->get_theta(0))*(s->get_theta(0));
        y += w*(s->get_theta(1));
        yy += w*(s->get_theta(1))*(s->get_theta(1));
    }
 
    cout << "<x> = " << x << " +/- " << sqrt(xx - x*x) << endl; 
    cout << "<y> = " << y << " +/- " << sqrt(yy - y*y) << endl; 
}

// main code starts here
int main()
{
    const int D = 2;
    int N = 1000;
    srand(time(NULL)); // seed random number generator
    double logZ = -DBL_MAX;
    double logwidth, logZnew, logZ_err; 
    double logLmin, H=0.0;
    int j, nest, worst, copy;
    double dat;
    while (infile >> dat) {data.push_back(dat);}
    Point *pts[N]; // create N active points
    list<Point> samples;

    // JRothe vars
    gsl_vector * mycenter = gsl_vector_alloc(D);
    gsl_vector * newcoor = gsl_vector_alloc(D);
    gsl_matrix * C = gsl_matrix_calloc(D,D);
    double f;


    for(j=0; j<N; j++)
    {
        pts[j] = new Point(D); 
        prior(pts[j]);
        logLhood(pts[j]);
    }

    logwidth = log(1.0 - exp(-1.0/N));

    for(nest=0; nest<MAXITER; nest++)
    {
        worst = 0; 
        for(j=1; j<N; j++)
            if(pts[j]->get_logL() < pts[worst]->get_logL()) worst = j; 

        pts[worst]->set_logWt(logwidth + pts[worst]->get_logL());

        logZnew = PLUS(logZ, pts[worst]->get_logWt()); 
        H = exp(pts[worst]->get_logWt() - logZnew)*(pts[worst]->get_logL()) + exp(logZ - logZnew)*(H + logZ) - logZnew;
        logZ = logZnew; 

        do copy = (int)(N*UNIFORM) % N; 
        while (copy == worst);
       
        samples.push_back(*pts[worst]);
        logLmin = pts[worst]->get_logL();
        *pts[worst] = *pts[copy];
        
        gsl_vector_set_zero(mycenter);
        gsl_vector_set_zero(newcoor);
        
        FindEnclosingEllipsoid(D, N, pts, mycenter, C, &f);

        do 
        {
            f = 1.5;
            SampleEllipsoid(D, mycenter, C, f, newcoor);
            pts[worst]->set_theta(newcoor, D);
            logLhood(pts[worst]);
        }
        while(logLmin > pts[worst]->get_logL());  

//        mcmc(pts[worst], logLmin, D);
    
        logwidth -= 1.0/N;
    }

    gsl_vector_free(mycenter);
    gsl_vector_free(newcoor);
    gsl_matrix_free(C);

    for(j=0; j<N; j++){delete pts[j];}

    logZ_err = sqrt(H/N);
    results(samples, logZ);
    cout << "information: H =  " << H/log(2.0) << " bits " << endl;
    cout << "evidence: logZ = " << logZ << " +/- " << logZ_err << endl;

    ofstream outfile;
    outfile.open("posterior_pdfs.dat");
    list<Point>::iterator s;
    for(s=samples.begin(); s!=samples.end(); s++)
        outfile << s->get_theta(0) << " " << s->get_theta(1) << " " << exp(s->get_logL() - logZ) << endl;

    return 0;
}
