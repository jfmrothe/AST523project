#ifndef DATA_H
#define DATA_H

#include "Point.h"

class Data
{
    private:
        int D;
        int num_col; 
        vector< list<double> > data;
        list<double>::iterator datum;

    public:
        Data();
        Data(int Dim, int n_col) : data(n_col) {D = Dim; num_col = n_col;}
        void get_data(string);
        void get_results(list<Point *>&, double);
        void lighthouse_logL(Point*);
        void eggbox_logL(Point*);
};
#endif
