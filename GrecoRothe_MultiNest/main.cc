#include "Point.h"
#include "Data.h"
#include "Samplers.h"

int main(int argc, char *argv[])
{
    if (argc != 2)   // check for command-line arguments
    {
        cout << "usage: " << argv[0] << " numpoints"<< endl;
        exit (1);   
    }

    // **** get runtime parameters
    int N = atoi(argv[1]);
    int D, num_cols;
    ifstream runtime_file("runtime.txt");
    string datafile_name; 
    runtime_file >> datafile_name;
    runtime_file >> num_cols;
    runtime_file >> D;

    vector<string> param_names(D);
    vector<string> prior_types(D);
    double min_vals[D], max_vals[D], guess_vals[D];
    for(int i = 0; i<D; i++)
    {
        runtime_file >> param_names[i];  
        runtime_file >> prior_types[i];  
        runtime_file >> min_vals[i];  
        runtime_file >> max_vals[i];  
        runtime_file >> guess_vals[i];  
    }
    runtime_file.close();
    // **************************
    
    // **** get data
    Data data_obj(D, num_cols);
    data_obj.get_data(datafile_name);
    // ****

    // create sampler object
    Samplers sampler(D);

    // seed random number generator
    srand(time(NULL)); 
   
    // **** create N active points and set params
    Point *pts[N]; 
    for(int j=0; j<N; j++)
    {
        pts[j] = new Point(D); 
        pts[j]->set_params(param_names, prior_types, min_vals, max_vals, guess_vals);
        pts[j]->hypercube_prior();
        pts[j]->transform_prior();
        data_obj.lighthouse_logL(pts[j]);
    }
    // ******************************************

    // ****************************************** start nested sampling algorithm 
    double logZ = -DBL_MAX;
    double logwidth, logZnew, logZ_err, X_i; 
    double logLmin, H=0.0;
    double logLmax = -999.9;
    double logL_tmp = 0.0;
    
    int j, nest, worst, copy;
    list<Point> sample_pts; // list of Point objects to sample posterior 

    logwidth = log(1.0 - exp(-1.0/N));

    nest = 0;
 
    do
    {
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

        sample_pts.push_back(*pts[worst]);
        logLmin = pts[worst]->get_logL();
        
        // **************** ellipsoidal sampling 
        sampler.set_vectors_zero();
        sampler.FindEnclosingEllipsoid(N, pts);
        sampler.set_f_factor(1.06);
        do
        {
            do sampler.SampleEllipsoid();
            while (!sampler.u_in_hypercube());
            pts[worst]->set_u(sampler.get_coor());
            pts[worst]->transform_prior();
            data_obj.lighthouse_logL(pts[worst]);
        }
        while(logLmin > pts[worst]->get_logL());
        // **************** 

        // to use MCMC search method, use the 4 lines below (and comment of the ellipsoidal sampling)
        //do copy = (int)(N*UNIFORM) % N; 
        //while (copy == worst);
        //*pts[worst] = *pts[copy];
        //sampler.mcmc(pts[worst], data_obj, logLmin);
    
        X_i = exp(-nest/N);
        logwidth -= 1.0/N;
        nest++;
    }
    while(THRESH < abs(X_i*logLmax));
    // ****************************************** end nested sampling algorithm 

    cout << logLmax << endl;
    cout << nest << " iterations used" << endl;

    // ************* output results
    logZ_err = sqrt(H/N);
    data_obj.get_results(sample_pts, logZ);
    cout << "information: H =  " << H/log(2.0) << " bits " << endl;
    cout << "evidence: logZ = " << logZ << " +/- " << logZ_err << endl;

    ofstream outfile;

    outfile.open("posterior_pdfs.dat");
    list<Point>::iterator s;
    for(s=sample_pts.begin(); s!=sample_pts.end(); s++)
        outfile << s->get_theta(0) << " " << s->get_theta(1) << " " << exp(s->get_logL() - logZ) << endl;
    // ************* 

    for(int j=0; j<N; j++){delete pts[j];} 
    return 0;
}
