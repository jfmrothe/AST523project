#ifndef Parameter_H
#define Parameter_H

class Parameter
{
    private:
        std::string name;
        std::string prior;
        double min; 
        double max; 
        double guess;

    public:
        Parameter(std::string, std::string, double, double, double);  
        double get_min(){return min;}
        double get_max(){return max;}
        double get_guess(){return guess;}
        std::string get_name(){return name;}

        void set_min(double x){min = x;}
        void set_max(double x){max = x;}
        void set_guess(double x){guess= x;}
        void set_name(std::string s){name= s;}
};

inline Parameter::Parameter(std::string n, std::string p, double m, double M, double g)
{
    name = n;
    prior = p;
    min = m; 
    max = M; 
    guess = g;
}
#endif
