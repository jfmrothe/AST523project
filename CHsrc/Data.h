/***************************************************

File: Data.h 

Description:
Header file for the Data class. Data objects handle
all interaction with data and likelihood evaluations.

Programmers: Johnny Greco & Johannes Rothe
Contacts: greco@princeton.edu, jrothe@princeton.edu

****************************************************/
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
        //void get_results(list<Point *>&, double);
       void Get_L(double *theta, double *logLarr, int nl);
};
#endif
