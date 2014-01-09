#include "Point.h"
#include "Data.h"
#include "Samplers.h"

int main(int argc, char *argv[])
{

    if (argc != 2)   // check for command-line arguments
    {
        cout << "usage: " << argv[0] << " runtime_filename"<< endl;
        exit (1);   
    }

    // **** get runtime parameters
    string runtime_filename = argv[1];
    int N, D, num_cols;
    double e;
    ifstream runtime_file(runtime_filename.c_str());
    string datafile_name; 
    runtime_file.ignore(256, ':');
    runtime_file >> datafile_name;
    runtime_file.ignore(256, ':');
    runtime_file >> num_cols;
    runtime_file.ignore(256, ':');
    runtime_file >> D;
    runtime_file.ignore(256, ':');
    runtime_file >> N;
    runtime_file.ignore(256, ':');
    runtime_file >> e;

    vector<string> prior_types(D);
    vector<double> min_vals(D), max_vals(D);
    for(int i = 0; i<D; i++)
    {
        runtime_file.ignore(256, ':');
        runtime_file >> prior_types[i];  
        runtime_file >> min_vals[i];  
        runtime_file >> max_vals[i];  
    }
    runtime_file.close();
    // **************************
    
    // create data object
    Data data_obj(D, num_cols);

    // **** get data if necessary
    if(datafile_name != "None") {data_obj.get_data(datafile_name);}

    // create sampler object
    Samplers sampler(D, e);
   
    // seed random number generator
    srand(time(NULL)); 
    
    // **** create N active points and set params
    vector <Point *> pts(N); 
    for(int j=0; j<N; j++)
    {
        pts[j] = new Point(D); 
        pts[j]->set_params(prior_types, min_vals, max_vals);
        pts[j]->hypercube_prior();
        pts[j]->transform_prior();
        data_obj.eggbox_logL(pts[j]);
    }
    // ******************************************

    // ****************************************** start nested sampling algorithm 
    double logZ = -DBL_MAX;
    double logwidth, logZnew, logZ_err, X_i; 
    double logLmin, H=0.0;
    double logLmax = -999.9;
    double logL_tmp = 0.0;
    
    int j, nest, worst, copy;
    list<Point *> discard_pts; // list of Point objects to sample posterior 

    logwidth = log(1.0 - exp(-1.0/N));

    nest = 0;
 
    do
    {
        X_i = exp(-nest/N);
        worst = 0; 
        for(j=1; j<N; j++)
        {
            logL_tmp = pts[j]->get_logL();
            if(logL_tmp < pts[worst]->get_logL()) worst = j; 
            if(logL_tmp > logLmax) logLmax = logL_tmp;
        }

        pts[worst]->set_logWt(logwidth + pts[worst]->get_logL());

        logZnew = PLUS(logZ, pts[worst]->get_logWt()); 
        H = exp(pts[worst]->get_logWt() - logZnew)*(pts[worst]->get_logL()) 
            + exp(logZ - logZnew)*(H + logZ) - logZnew;
        logZ = logZnew; 
         
        discard_pts.push_back( new Point(*pts[worst]) );
        logLmin = pts[worst]->get_logL();

        // **************** ellipsoidal sampling 
        sampler.ClearCluster();
        sampler.EllipsoidalPartitioning(pts, X_i); 
        sampler.CalcVtot();
        do
        {
            pts[worst]->set_u(sampler.get_newcoor());
            pts[worst]->transform_prior();
            data_obj.eggbox_logL(pts[worst]);
        }
        while(logLmin > pts[worst]->get_logL());
        // **************** 

        // to use MCMC search method, use the 4 lines below (and comment out the ellipsoidal sampling)
        //do copy = (int)(N*UNIFORM) % N; 
        //while (copy == worst);
        // *pts[worst] = *pts[copy];
        //sampler.mcmc(pts[worst], data_obj, logLmin);
    
        logwidth -= 1.0/N;
        nest++;
    }
    while(THRESH < abs(X_i*logLmax));
    // ****************************************** end nested sampling algorithm 

    cout << nest << " iterations used" << endl;

    // ************* output results
    logZ_err = sqrt(H/N);
    data_obj.get_results(discard_pts, logZ);
    cout << "information: H =  " << H/log(2.0) << " bits " << endl;
    cout << "evidence: logZ = " << logZ << " +/- " << logZ_err << endl;

    ofstream outfile;

    outfile.open("posterior_pdfs.dat");
    list<Point *>::iterator s;
    for(s=discard_pts.begin(); s!=discard_pts.end(); s++)
      outfile << (*s)->get_theta(0) << " " << (*s)->get_theta(1) << " " << (*s)->get_logL() << " " << exp((*s)->get_logL() - logZ) << endl;
    // ************* 
    outfile.close();

    for(int j=0; j<N; j++){delete pts[j];} 
    for(list<Point *>::iterator s=discard_pts.begin();s!=discard_pts.end();s++){delete *s;} 

    return 0;
}
