// #include "DTSource.h"
#include "DTSaveError.h"

#include "DTArguments.h"
#include "DTDataFile.h"
#include "DTDictionary.h"
#include "DTDoubleArray.h"
#include "DTPointCollection3D.h"
#include "DTProgress.h"
#include "DTSeriesArray.h"
#include "DTSeriesPointCollection3D.h"

//////////////////////////////////////////////////////////////////////////////
//    Main routine
//////////////////////////////////////////////////////////////////////////////

void Computation(const DTSeriesPointCollection3D &points,double embed_points,
                 double embed_delay,double time_points,DTSeriesArray &computed);

int main(int argc,const char *argv[])
{
    DTSetArguments(argc,argv);
    
    DTDataFile inputFile("Input.dtbin",DTFile::ReadOnly);
    DTDataFile outputFile("Output.dtbin",DTFile::NewReadWrite);
    // Input variables.
    DTSeriesPointCollection3D points;
    Read(inputFile,"points",points);
    double embed_points = inputFile.ReadNumber("embed_points");
    double embed_delay = inputFile.ReadNumber("embed_delay");
    double time_points = inputFile.ReadNumber("time_points");
    
    // Output series.
    DTSeriesArray computed(outputFile,"Var");
    if (DTArgumentIncludesFlag("saveInput")) { // Add -saveInput to the argument list to save the input in the output file.
        WriteOne(outputFile,"embed_points",embed_points);
        WriteOne(outputFile,"embed_delay",embed_delay);
        WriteOne(outputFile,"time_points",time_points);
    }
    
    
    // The computation.
    clock_t t_before = clock();
    Computation(points,embed_points,embed_delay,time_points,computed);
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

void Computation(const DTSeriesPointCollection3D &points,double embed_points,
                 double embed_delay,double time_points,DTSeriesArray &computed)
{
    // Insert your code here.
    
    DTProgress progress;
    progress.UpdatePercentage(0.0f);
    
    int n = points.Get(points.TimeNumber(0)).NumberOfPoints();
    
    // Inside the loop, do
    //     progress.UpdatePercentage(fraction);
    //     computed.Add(returnStructure,time); // Call with time>=0 and strictly increasing.
    for (int i=0;i<time_points;++i)
    {
        DTDoubleArray cur_points = points.Get(points.TimeNumber(i)).DoubleData();
        DTMutableDoubleArray this_out(3,embed_points,n);
        for (int m=0;m<n;++m)
        {
            for (int j=0;j<embed_points;++j)
            {
                for (int k=0;k<3;++k)
                {
                    this_out(k,j,m) = cur_points(
                }
            }
        }
    }
}

