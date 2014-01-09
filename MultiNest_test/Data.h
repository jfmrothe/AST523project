#ifndef DATA_H
#define DATA_H

#include "Point.h"

class Data
{
    private:
        int D;
        int num_col; 
        string data_filename;
        vector< list<double> > data;
        list<double>::iterator datum;

    public:
        Data(int Dim, int n_col, string filename);
        void get_data();
        void get_results(list<Point *>&, double);
        void lighthouse_logL(Point*);
        void logL(Point*);
};
#endif
