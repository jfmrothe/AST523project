#ifndef POINT_H
#define POINT_H
#include <vector>
using namespace std;

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
        void set_theta(gsl_vector * coor, int D){for(int i = 0; i < D; i++){theta[i] = gsl_vector_get(coor,i);}}
        double get_theta(int i){return theta[i];}
        vector<double> get_theta(){return theta;}
        double get_u(int i){return u[i];}
        void set_logL(double L){logL = L;}
        double get_logL(){return logL;}
        void set_logWt(double W){logWt = W;}
        double get_logWt(){return logWt;}
};
#endif
