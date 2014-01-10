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
    
    srand(time(NULL));  // seed random number generator

    Data data_obj(D, num_cols, datafile_name); // create data object
    if(num_cols != 0) {data_obj.get_data();} // **** get data if necessary
    Samplers sampler(D, N, data_obj, prior_types, min_vals, max_vals, e); // create sampler object


    cout << "running MultiNest algorithm... this may take a few minutes" << endl;

    // ****************************************** start nested sampling algorithm 
  
    double X_i, logLmax;
    int NumRecluster = 0;
    int j, nest;
    nest = 0;

    do
    {
        X_i = exp(-nest/N);

	// discard and resample, get logLmax for convergence check as byproduct
	logLmax = sampler.ResetWorstPoint(nest, data_obj);

        // **************** ellipsoidal partitioning 
        if(nest==0 || sampler.ClusteringQuality(X_i) > RepartitionFactor)  // recluster? 
        {
	  sampler.ClearCluster();
	  vector <Point *> empty;
	  sampler.EllipsoidalPartitioning(empty, X_i);
	  sampler.EraseFirst();
	  sampler.CalcVtot();
	  NumRecluster++;
        }
        // *********************************************
	nest++;
    }
    while(THRESH < abs(X_i*logLmax)); // stopping criterion
    // ****************************************** end nested sampling algorithm 

    // ************* output results

    cout << "job complete!" << endl;
    cout << "**** results ****" << endl;
    cout << "number iterations = " << nest << endl;
    cout << "number reclusters = " << NumRecluster <<endl;
    sampler.printout();

    return 0;
}
