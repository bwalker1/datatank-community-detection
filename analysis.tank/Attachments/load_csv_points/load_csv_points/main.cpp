// #include "DTSource.h"
#include "DTSaveError.h"

#include "DTArguments.h"
#include "DTDataFile.h"
#include "DTDictionary.h"
#include "DTPointCollection3D.h"
#include "DTProgress.h"
#include "DTSeriesPointCollection3D.h"


#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>

using namespace std;

//////////////////////////////////////////////////////////////////////////////
//    Main routine
//////////////////////////////////////////////////////////////////////////////

void Computation(const string &file,DTSeriesPointCollection3D &computed);

int main(int argc,const char *argv[])
{
    DTSetArguments(argc,argv);
    
    DTDataFile inputFile("Input.dtbin",DTFile::ReadOnly);
    DTDataFile outputFile("Output.dtbin",DTFile::NewReadWrite);
    // Input variables.
    string file = inputFile.ReadString("file");
    
    // Output series.
    DTSeriesPointCollection3D computed(outputFile,"Var");
    if (DTArgumentIncludesFlag("saveInput")) { // Add -saveInput to the argument list to save the input in the output file.
        WriteOne(outputFile,"file",file);
    }
    
    
    // The computation.
    clock_t t_before = clock();
    Computation(file,computed);
    clock_t t_after = clock();
    double exec_time = double(t_after-t_before)/double(CLOCKS_PER_SEC);
    
    // The execution time.
    outputFile.Save(exec_time,"ExecutionTime");
    outputFile.Save("Real Number","Seq_ExecutionTime");
    
    // The errors.
    DTSaveError(outputFile,"ExecutionErrors");
    outputFile.Save("StringList","Seq_ExecutionErrors");
    
    outputFile.SaveIndex();
    
    return 0;
}

//////////////////////////////////////////////////////////////////////////////
//    Computational routine
//////////////////////////////////////////////////////////////////////////////

void Computation(const string &file,DTSeriesPointCollection3D &computed)
{
    DTProgress progress;
    
    // Inside the loop, do
    //     progress.UpdatePercentage(fraction);
    //     computed.Add(returnStructure,time); // Call with time>=0 and strictly increasing.
    
    
    ifstream f(file);
    string line;
    while (getline(f,line))
    {
        stringstream ss(line);
        string item;
        getline(ss,item,',');
        double t = stod(item);
        vector<double> vals;
        while (getline(ss,item,','))
        {
            vals.push_back(stod(item));
        }
        int len = vals.size();
        int N = len/3;
        DTMutableDoubleArray data(3,N);
        for (int i=0;i<len;++i)
        {
            data(i) = vals[i];
        }
        DTPointCollection3D pts(data);
        computed.Add(pts, t);
    }
}

