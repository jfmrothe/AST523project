#ifndef ELLIPSOID_H
#define ELLIPSOID_H

#include "Point.h"

struct Ran{
  unsigned long long int u,v,w;
  Ran(unsigned long long int j) : v(4101842887655102017LL),w(1){
    u = j^v; int64();
    w = v; int64();
  }
  inline unsigned long long int int64() {
    u = u*2862933555777941757LL + 7046029254386353087LL;
    v ^=v >> 17; v^=v<<31; v^=v >>8;
    w = 4294957665U*(w & 0xffffffff) + (w>>32);
    unsigned long long int x = u^(u<<21); x^=x>>35; x^=x<<4;
    return (x +v) ^w;
  }
  inline double doub(){return 5.42101086242752217E-20 * int64();}
  inline unsigned int int32() {return (unsigned int) int64();}
};


//Ran myrand(17);
//float boxmuller();
//float quadr();

void unisphere(float *, int);

class Ellipsoid {
public:
  Ellipsoid(int D, gsl_vector * center, gsl_matrix * C, double f, vector<Point *>&);
  ~Ellipsoid();
  int getD();
  gsl_vector * getCenter();
  gsl_matrix * getCovMat();
  gsl_matrix * getCinv();
  gsl_matrix * getA();
  gsl_matrix * getAinv();
  void SampleEllipsoid();
  void RescaleToCatch();
  bool IsMember(gsl_vector *);
  bool intersect(Ellipsoid& other);
  void fetchPoints(Ellipsoid& other);
  double getEnlFac();
  void setEnlFac(double f);
  double getVol();
  double mdist(Point *);
  void printout();
  gsl_vector * get_newcoor() {return newcoor_;}
  vector<Point *> ell_pts_;
private:
  int D_;
  gsl_vector * center_;
  gsl_matrix * Cinv_;
  gsl_matrix * A_;
  gsl_matrix * Ainv_;
  gsl_matrix * covMat_;
  gsl_vector * newcoor_;
  double f_;
  double vol_;

};
#endif
