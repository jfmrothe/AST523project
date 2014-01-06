#include <ellipsoid.hh>
#include <sampleellipsoid.hh>
#include <findenclosing.hh>

#include <vector>

#include <cluster.hh>       
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>


int main() {
  srand(time(NULL));
  return TestEllipsoidalPartitioningVec();
//  return Test2FindSelectedEnclosingEllipsoid();
}
