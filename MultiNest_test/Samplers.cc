#include "Samplers.h"
#include "Point.h"

gsl_vector * Samplers::get_newcoor()
{
    int RandEll;
    int NumEll = clustering.size();
    int n_e;
    gsl_vector * tmp_coor = gsl_vector_alloc(D);

    do 
    {
        do RandEll= rand() % clustering.size();
        while (clustering[RandEll]->getVol()/Vtot < UNIFORM);

        do clustering[RandEll]->SampleEllipsoid();
        while (!u_in_hypercube(clustering[RandEll]->get_newcoor(), D));

        gsl_vector_memcpy(tmp_coor, clustering[RandEll]->get_newcoor());

        n_e = 0;
        for(int i = 0; i<NumEll; i++)
        {
           if(clustering[i]->IsMember(tmp_coor)) {n_e++;}
        }
    }
    while(1.0/n_e < UNIFORM);
    
    gsl_vector_free(tmp_coor);

    return clustering[RandEll]->get_newcoor();
}

void Samplers::CalcVtot()
{
    Vtot=0.0; 

    for(int i=0; i<clustering.size(); i++)
    {
       Vtot += clustering[i]->getVol();
    }
}

void Samplers::EllipsoidalPartitioning(vector<Point *>& pts, double Xtot) 
{
  // this is the variant currently in development, because segfaults could be
  // avoided. will check if this implementation is "good" in some sense...
  // vector of ellipsoid pointers is returned, these are deleted after use in
  // the function calling EllipsoidalPartitioningVec performs Algorithm I from
  // Feroz, Hobson and Bridges (2009) on N points in D-dimensional
  // [0,1]-hypercube given in coors. Returns resulting number of ellipsoids in
  // Nell, uses new to create array of ellipsoids, first address is returned in
  // clustering

  /////////////////
  // ALGORITHM I //
  /////////////////

  // create mainEllipsoid which may be split further
  // ?? in recursion depth, this first ellipsoid can be passed by parent??

  int N = pts.size();

  Ellipsoid mainEll = FindEnclosingEllipsoid(pts,D);

  //enlarge if neccessary
  if(Xtot>mainEll.getVol()) {
    mainEll.setEnlFac(mainEll.getEnlFac()*pow(Xtot/mainEll.getVol(),1.0/D));
  }
  // initialize splitting using kmeans
  int i;
  int grouping[N];
  int k=2;
  KMeans(pts,D,k,&grouping[0]);
  vector<Point *> pts_group_0;
  vector<Point *> pts_group_1;

  SelectFromGrouping(pts, D, grouping, 0, pts_group_0);
  SelectFromGrouping(pts, D, grouping, 1, pts_group_1);

  bool changed = true;
  double X1;
  double X2;
  double h1;
  double h2;
  double tmp1,tmp2;
  double vol1,vol2;

  while(changed) {
    // return main ellipsoid immediately if partitioning would create singular (flat) ellipsoid
    if(pts_group_0.size()<D+1 or pts_group_1.size()<D+1) {
      clustering.push_back (new Ellipsoid(D, mainEll.getCenter(), mainEll.getCovMat(), mainEll.getEnlFac(), pts) );
      return;
    }
    // find new sub-Ellipsoids
    //"locality" of these variables removes object overwriting trouble
    Ellipsoid subEll1 = FindEnclosingEllipsoid(pts_group_0,D);
    vol1 = subEll1.getVol();
    X1 = ((double)pts_group_0.size()/N)*Xtot;
    Ellipsoid subEll2 = FindEnclosingEllipsoid(pts_group_1,D);
    vol2 = subEll2.getVol();
    X2 = ((double)pts_group_1.size()/N)*Xtot;

    //enlarge if neccessary
    if(X1>vol1) {
      subEll1.setEnlFac(subEll1.getEnlFac()*pow(X1/vol1,(double)(1.0/D)));
      vol1 = subEll1.getVol();
    }
    if(X2>vol2) {
    subEll2.setEnlFac(subEll2.getEnlFac()*pow(X2/vol2,(double)(1.0/D)));
    vol2 = subEll2.getVol();
    }

    // optimize splitting using mahalanobis-h-metric
    changed = false;    
    for(i=0;i<N;i++) 
    {
      tmp1 = subEll1.mdist(pts[i]);
      tmp2 = subEll2.mdist(pts[i]);
   
      h1 = vol1 / X1 * tmp1;
      h2 = vol2 / X2 * tmp2;

      if (h2<h1 and grouping[i]==0) 
      {
          changed = true;
          grouping[i] = 1;
      }
      if (h1<h2 and grouping[i]==1) 
      {
          changed = true;
          grouping[i] = 0;
      }
     }

     pts_group_0.clear();
     pts_group_1.clear();
     SelectFromGrouping(pts, D, grouping, 0, pts_group_0);
     SelectFromGrouping(pts, D, grouping, 1, pts_group_1);
  }
  // judgement if splitting should be continued
  if( vol1+vol2<mainEll.getVol() or mainEll.getVol()>2*Xtot) {

    // recursively start the splittings of subEll1 and subEll2
    EllipsoidalPartitioning(pts_group_0, X1);
    // second one
    EllipsoidalPartitioning(pts_group_1, X2);

  }
  else{
    // allocate memory for mainEll and push it to the vector
    clustering.push_back (new Ellipsoid(D, mainEll.getCenter(), mainEll.getCovMat(), mainEll.getEnlFac(), pts));
  }
  return;
}

void Samplers::mcmc(Point* pt, Data data_obj, double logLmin)
{
    vector<double> new_coords(D);
    double step = 0.1;
    int accept = 0;
    int reject = 0;
    Point *trial;
    trial = new Point(D);
    *trial = *pt;

    for(int j = 20; j > 0; j--)
    {
        for(int i = 0; i<D; i++)
        {
            new_coords[i] = pt->get_u(i) + step*(2.0*UNIFORM - 1.0);
            new_coords[i] -= floor(new_coords[i]);
            trial->set_u(i, new_coords[i]);
        }

        trial->transform_prior();
        data_obj.lighthouse_logL(trial);

        if(trial->get_logL() > logLmin){*pt = *trial; accept++;}
        else reject++;

        if(accept > reject) step *= exp(1.0 / accept);
        else if(accept < reject) step /= exp(1.0 / reject);
    }

    delete trial;
}
