#include <iostream>
#include <list>
#include <fstream>
#include <vector>
#include <cmath>
#include <float.h>
#include <time.h> 
using namespace std;

#define N 100
#define MAXITER 1000
#define UNIFORM ((rand() + 0.5)/(RAND_MAX+1.0))
#define PLUS(x,y) (x>y ? x+log(1+exp(y-x)) : y+log(1+exp(x-y)))

// This is a class for the active points
class Point 
{   
    private:
        vector<double> theta; // the dimension is set in main
        vector<double> u;
        double logL;
        double logWt;

    public:
        Point(int D) : theta(D), u(D){}
        void set_u(int i, double x){u[i] = x;}
        void set_theta(int i, double x){theta[i] = x;}
        double get_theta(int i){return theta[i];}
        vector<double> get_theta(){return theta;}
        double get_u(int i){return u[i];}
        void set_logL(double L){logL = L;} 
        double get_logL(){return logL;} 
        void set_logWt(double W){logWt = W;} 
        double get_logWt(){return logWt;} 
};   

// global variables
const int D = 2;
const double PI = 3.1415927;
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
        logL += log((y/PI) / ((*datum-x)*(*datum-x) + y*y));
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
void mcmc(Point* pt, double logLmin)
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
    srand(time(NULL)); // seed random number generator
    double logZ = -DBL_MAX;
    double logwidth, logZnew, logZ_err; 
    double logLmin, H=0.0;
    int j, nest, worst, copy;
    double dat;
    while (infile >> dat) {data.push_back(dat);}
    Point *pt[N]; // create N active points
    list<Point> samples;

    for(j=0; j<N; j++)
    {
        pt[j] = new Point(D); 
        prior(pt[j]);
        logLhood(pt[j]);
    }

    logwidth = log(1.0 - exp(-1.0/N));

    for(nest=0; nest<MAXITER; nest++)
    {
        worst = 0; 
        for(j=1; j<N; j++)
            if(pt[j]->get_logL() < pt[worst]->get_logL()) worst = j; 

        pt[worst]->set_logWt(logwidth + pt[worst]->get_logL());

        logZnew = PLUS(logZ, pt[worst]->get_logWt()); 
        H = exp(pt[worst]->get_logWt() - logZnew)*(pt[worst]->get_logL()) + exp(logZ - logZnew)*(H + logZ) - logZnew;
        logZ = logZnew; 

        do copy = (int)(N*UNIFORM) % N; 
        while (copy == worst);
       
        samples.push_back(*pt[worst]);
        logLmin = pt[worst]->get_logL();
        *pt[worst] = *pt[copy];
        
        mcmc(pt[worst], logLmin);
    
        logwidth -= 1.0/N;
    }

    for(j=0; j<N; j++){delete pt[j];}

    logZ_err = sqrt(H/N);
    results(samples, logZ);
    cout << "information: H =  " << H/log(2.0) << " bits " << endl;
    cout << "evidence: logZ = " << logZ << " +/- " << logZ_err << endl;

    ofstream outfile;
    outfile.open("posterior_pdfs.dat");
    list<Point>::iterator s;
    for(s=samples.begin(); s!=samples.end(); s++)
        outfile << s->get_theta(0) << " " << s->get_theta(1) << " " << exp(s->get_logL()-logZ) << endl;

    return 0;
}
