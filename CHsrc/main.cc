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
    vector<string> tempstring(D);
    string prior_types;
    double min_vals[D], max_vals[D];
    for(int i = 0; i<D; i++)
    {
        runtime_file.ignore(256, ':');
        runtime_file >> tempstring[i];  
        runtime_file >> min_vals[i];  
        runtime_file >> max_vals[i];  
    }
    runtime_file.close();
    prior_types = tempstring[0];
    //cout << min_vals[0] << " " << max_vals[0] << prior_types[0] << endl;
    // **************************
    cout << D << ' '<< N<< endl; 
    srand(time(NULL));  // seed random number generator

    Data data_obj(D, num_cols, datafile_name); // create data object
    cout << "create data object" << endl; 
    if(num_cols != 0) {data_obj.get_data();} // **** get data if necessary
    cout << "finish reading data" << endl; 
    //Samplers sampler(D, N, data_obj, prior_types, min_vals, max_vals, e); // create sampler object
    Samplers sampler( min_vals, D, max_vals,D,e,N,prior_types); // create sampler object
    cout << "running MultiNest algorithm... this may take a few minutes" << endl;

    // ****************************************** start nested sampling algorithm 
  
    double X_i, logLmax;
    int NumRecluster = 0;
    int NumSample = 0;
    int nest;
    nest = 0;
    double *Alltheta = new double [D*N];
    double *logL = new double [N];
    double *theta = new double [D];
    double *templogL = new double [1];
    double logLmin;
    bool FlagSample;
    bool Debug=false;
    //CH: before that, let's ask dataobj for logL
	  sampler.getAlltheta(Alltheta, D, N);  
    data_obj.Get_L(Alltheta,logL,N);
    sampler.SetAllPoint(logL, N);
    //return 0;
    if(Debug){
      for (int i=0; i<N; i++){
        cout<< i << " "<< Alltheta[D*i] << " " << Alltheta[D*i+1] << " " << logL[i] << endl;
      }
    }
    
    do
    {
      X_i = exp((double) -nest/N);
    // discard and resample, get logLmax for convergence check as byproduct
	sampler.DisgardWorstPoint(nest); 
        logLmax = sampler.getlogLmax();
        logLmin = sampler.getlogLmin();
        //cout << nest <<  " " << logLmin << endl;
        //return 0;
        FlagSample = true;
        templogL[0]=0;
        while (FlagSample){
          sampler.ResetWorstPoint(theta,D);
	        NumSample++;
          data_obj.Get_L(theta,templogL,1);
          //cout << theta[0] << " " << theta[1] << " " << templogL[0] <<  " " << logLmin << endl;
          if(logLmin < templogL[0]){
            FlagSample=false;
          }
        }
        double temp = templogL[0];
        sampler.ResetWorstPointLogL(temp);
        // **************** ellipsoidal partitioning 

	//	NumRecluster += sampler.Recluster(X_i,RepartitionFactor);
	NumRecluster++;
	sampler.FullRecluster(X_i);
        // *********************************************
	nest++;
    }
    while(THRESH < fabs(X_i*logLmax)); // stopping criterion
    // ****************************************** end nested sampling algorithm 

    // ************* output results

    cout << "job complete!" << endl;
    cout << "**** results ****" << endl;
    cout << "number iterations = " << nest << endl;
    cout << "number reclusters = " << NumRecluster <<endl;
    cout << "sampling efficiency = " << (double) nest/NumSample << endl;
    double logzinfo[3];
    sampler.getlogZ(logzinfo,3);
    cout << "information: H=" << logzinfo[0] << "bits" <<endl;
    cout << "global evidence: logZ =" << logzinfo[1] << "+/-" << logzinfo[2] << endl;
    
    delete [] Alltheta;
    delete [] logL;
    delete [] theta; 
    delete [] templogL;  
    return 0;
}
