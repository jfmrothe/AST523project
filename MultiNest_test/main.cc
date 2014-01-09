#include "Point.h"
#include "Data.h"
#include "Samplers.h"

int main(int argc, char *argv[])
{

    if (argc != 2)   // check for command-line arguments
    {
        cout << "usage: " << argv[0] << " problem"<< endl;
        exit (1);   
    }

    cout << "getting runtime parameters" << endl;

    // **** get runtime parameters
    string runtime_filename = argv[1];
    int N, D, num_cols;
    double e, RepartitionFactor;
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
    runtime_file.ignore(256, ':');
    runtime_file >> RepartitionFactor;
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
    
    Data data_obj(D, num_cols, datafile_name); // create data object
    if(num_cols != 0) {data_obj.get_data();} // **** get data if necessary
    Samplers sampler(D, e); // create sampler object
    srand(time(NULL));  // seed random number generator
    
    cout << "creating " << N << " active points" << endl;

    // **** create N active points and set params
    vector<Point *> pts(N); 
    for(int j=0; j<N; j++)
    {
        pts[j] = new Point(D); 
        pts[j]->set_params(prior_types, min_vals, max_vals);
        pts[j]->hypercube_prior();
        pts[j]->transform_prior();
        data_obj.logL(pts[j]);
    }
    // ******************************************

    cout << "running MultiNest algorithm... this may take a few minutes" << endl;

    // ****************************************** start nested sampling algorithm 
    double logZ = -DBL_MAX;
    double logwidth, logZnew, logZ_err, X_i; 
    double logLmin, H=0.0;
    double logLmax = -999.9;
    double logL_tmp = 0.0;
    int NumRecluster = 0;
    int j, nest, worst, copy;
    nest = 0;

    list<Point> discard_pts; // list of Point objects to sample posterior 

    logwidth = log(1.0 - exp(-1.0/N)); 
    do
    {
        X_i = exp(-nest/N);
        worst = 0; 
        for(j=1; j<N; j++) // for lowest and highest logL
        {
            logL_tmp = pts[j]->get_logL();
            if(logL_tmp < pts[worst]->get_logL()) worst = j; 
            if(logL_tmp > logLmax) logLmax = logL_tmp;
        }
        logLmin = pts[worst]->get_logL();
        pts[worst]->set_logWt(logwidth + pts[worst]->get_logL()); // find contribution to evidence

        logZnew = PLUS(logZ, pts[worst]->get_logWt()); 
        H = exp(pts[worst]->get_logWt() - logZnew)*(pts[worst]->get_logL()) 
            + exp(logZ - logZnew)*(H + logZ) - logZnew;
        logZ = logZnew; // update global evidence
         
        discard_pts.push_back(*pts[worst]); // save discarded point for posterior sampling

        // **************** ellipsoidal partitioning and sampling 
        if(nest==0 || sampler.ClusteringQuality(X_i) > RepartitionFactor)  // recluster? 
        {
            sampler.ClearCluster();
            sampler.EllipsoidalPartitioning(pts, X_i);
            sampler.CalcVtot();
            NumRecluster++;
        }
        do
        {
            pts[worst]->set_u(sampler.get_newcoor());
            pts[worst]->transform_prior();
            data_obj.logL(pts[worst]);
        }
        while(logLmin > pts[worst]->get_logL());
        // *********************************************
        logwidth -= 1.0/N;
        nest++;
    }
    while(THRESH < abs(X_i*logLmax)); // stoping criterion
    // ****************************************** end nested sampling algorithm 

    // ************* output results
    logZ_err = sqrt(H/N);
    cout << "job complete!" << endl;
    cout << "**** results ****" << endl;
    cout << "number iterations = " << nest << endl;
    cout << "number reclusters = " << NumRecluster <<endl;
    cout << "information: H =  " << H/log(2.0) << " bits " << endl;
    cout << "global evidence: logZ = " << logZ << " +/- " << logZ_err << endl;
    cout << "****************" << endl;
    ofstream outfile;
    outfile.open("posterior_pdfs.dat");
    list<Point>::iterator s;
    for(s=discard_pts.begin(); s!=discard_pts.end(); s++)
        outfile << s->get_theta(0) << " " << s->get_theta(1) << " " << 
                   s->get_logL() << " " << exp(s->get_logL() - logZ) << endl;
    // ************* 
    outfile.close();

    for(int j=0; j<N; j++){delete pts[j];} 

    return 0;
}
