// #include "DTSource.h"
#include "DTSaveError.h"
#include "DTDataFile.h"
#include "DTArguments.h"
#include "DTDictionary.h"
#include "DTDoubleArray.h"
#include "DTMatlabDataFile.h"
#include "DTProgress.h"
#include "DTSeriesArray.h"

//////////////////////////////////////////////////////////////////////////////
//    Main routine
//////////////////////////////////////////////////////////////////////////////


void Computation(const DTDoubleArray &network, double gamma, double omega, string genlouvain, string matlabpath);

int main(int argc,const char *argv[])
{
    DTSetArguments(argc,argv);
    
    DTMatlabDataFile inputFile("Input.mat",DTFile::ReadOnly);
    // Read in the input variables.
    DTDoubleArray network = inputFile.ReadDoubleArray("adjacency");
    double gamma = inputFile.ReadNumber("gamma");
    double omega= inputFile.ReadNumber("omega");
    
    string genlouvain = inputFile.ReadString("genlouvainpath");
    genlouvain = genlouvain.substr(0,genlouvain.size()-13);
    string matlabpath = inputFile.ReadString("matlabpath");
    
    // The computation.
    DTDoubleArray computed;
    clock_t t_before = clock();
    Computation(network,gamma, omega, genlouvain, matlabpath);
    clock_t t_after = clock();
    double exec_time = double(t_after-t_before)/double(CLOCKS_PER_SEC);
    
    // Write the output.
    DTMatlabDataFile outputFile("Output.mat",DTFile::ExistingReadWrite);
    
    // Output from computation
    //outputFile.Save(computed,"Var");
    outputFile.Save("Array","Seq_Var");
    
    // The execution time.
    outputFile.Save(exec_time,"ExecutionTime");
    outputFile.Save("Real Number","Seq_ExecutionTime");
    
    // The errors.
    DTSaveError(outputFile,"ExecutionErrors");
    outputFile.Save("StringList","Seq_ExecutionErrors");
    
    return 0;
}


void Computation(const DTDoubleArray &network,double gamma, double omega, string genlouvain, string matlabpath)
{

    //DTProgress progress;
    
    // generate the matlab file to do the partitions

    //string genlouvain = "/Users/ben/Dropbox/GenLouvain-2.1/";
    //string matlabpath = "/Applications/Matlab_R2018a.app/bin/matlab";

    string script = genlouvain + "/DTWrapper4.m";
    system("touch script.m && rm script.m");
    system("pwd");
    system(("echo \"addpath(genpath('" + genlouvain + "'));\n\" >> script.m").c_str());
    system(("cat \""+script+"\">> script.m").c_str());
    
    system((matlabpath +  " -nojvm -nosplash -r script,quit").c_str());
}
