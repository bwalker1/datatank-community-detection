// Include DTSource.h if you want to include all the headers.

#include "DTArguments.h"
#include "DTSaveError.h"

#include "DTDataFile.h"
#include "DTDoubleArray.h"
#include "DTPointCollection3D.h"
#include "DTSeriesPointCollection3D.h"

// Common utilities
#include "DTDoubleArrayOperators.h"
#include "DTProgress.h"
#include "DTTimer.h"
#include "DTUtilities.h"
#include "DTDictionary.h"
#include <assert.h>
#include <math.h>


DTDoubleArray Computation(double bins,double width,const DTSeriesPointCollection3D &points,
                          double threshold, double scale, double start);

int main(int argc,const char *argv[])
{
    DTSetArguments(argc,argv);
    
    DTDataFile inputFile("Input.dtbin",DTFile::ReadOnly);
    // Read in the input variables.
    double bins = inputFile.ReadNumber("bins");
    double width = inputFile.ReadNumber("width");
    DTSeriesPointCollection3D points;
    Read(inputFile,"points",points);
    double threshold = inputFile.ReadNumber("threshold");
    double scale = inputFile.ReadNumber("scale");
    double start = inputFile.ReadNumber("start");
    
    // The computation.
    DTDoubleArray computed;
    clock_t t_before = clock();
    computed = Computation(bins,width,points,threshold,scale,start);
    clock_t t_after = clock();
    double exec_time = double(t_after-t_before)/double(CLOCKS_PER_SEC);
    
    // Write the output.
    DTDataFile outputFile("Output.dtbin",DTFile::NewReadWrite);
    
    // Output from computation
    outputFile.Save(computed,"Var");
    outputFile.Save("Array","Seq_Var");
    
    // The execution time.
    outputFile.Save(exec_time,"ExecutionTime");
    outputFile.Save("Real Number","Seq_ExecutionTime");
    
    // The errors.
    DTSaveError(outputFile,"ExecutionErrors");
    outputFile.Save("StringList","Seq_ExecutionErrors");
    
    outputFile.SaveIndex();
    
    return 0;
}

DTDoubleArray Computation(double bins,double width,const DTSeriesPointCollection3D &points,
                          double threshold, double scale, double start)
{
    DTProgress progress;
    progress.UpdatePercentage(0);
    
    DTPointCollection3D p0 = points(0);
    int n = p0.NumberOfPoints();
    int nt = points.HowManySaved();
    DTMutableDoubleArray toReturn(n,n,bins);
    
    assert(bins*width <= nt);
    
    for (int bin=0;bin<bins;++bin)
    {
        // track progress in DataTank
        progress.UpdatePercentage(float(bin)/float(bins));
        // initialize memory to 0
        for (int j=0;j<n;++j)
        {
            for (int i=0;i<n;++i)
            {
                toReturn(i,j,bin) = 0;
            }
        }
        
        // do an average distance over all time values inside the current bin
        for (int k=0;k<width;++k)
        {
            int t = start + bin*width + k;
            DTPointCollection3D p = points(t);
            for (int j=0;j<n;++j)
            {
                DTPoint3D p1 = p(j);
                for (int i=0;i<n;++i)
                {
                    if (i==j) continue;
                    DTPoint3D p2 = p(i);
                    
                    toReturn(i,j,bin) += Norm(p2-p1);
                }
            }
        }
        
        // map the average distances into adjacency values
        for (int j=0;j<n;++j)
        {
            for (int i=0;i<n;++i)
            {
                if (i==j) continue;
                toReturn(i,j,bin) /= width;
                if (toReturn(i,j,bin) > threshold)
                {
                    toReturn(i,j,bin) = 0;
                }
                else
                {
                    toReturn(i,j,bin) = exp(-toReturn(i,j,bin)/scale);
                }
            }
        }
    }
    
    return toReturn;
}
