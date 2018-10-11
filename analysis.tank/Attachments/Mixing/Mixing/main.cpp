// Include DTSource.h if you want to include all the headers.

#include "DTArguments.h"
#include "DTSaveError.h"

#include "DTDataFile.h"
#include "DTPointCollection3D.h"
#include "DTSeriesPointCollection3D.h"

// Common utilities
#include "DTDoubleArrayOperators.h"
#include "DTProgress.h"
#include "DTTimer.h"
#include "DTUtilities.h"
#include "DTDictionary.h"

#include <math.h>

double Computation(double dist,const DTSeriesPointCollection3D &points_series);

int main(int argc,const char *argv[])
{
    DTSetArguments(argc,argv);

    DTDataFile inputFile("Input.dtbin",DTFile::ReadOnly);
    // Read in the input variables.
    double dist = inputFile.ReadNumber("dist");
    DTSeriesPointCollection3D points_series;
    Read(inputFile,"points_series",points_series);

    // The computation.
    double computed;
    clock_t t_before = clock();
    computed = Computation(dist,points_series);
    clock_t t_after = clock();
    double exec_time = double(t_after-t_before)/double(CLOCKS_PER_SEC);

    // Write the output.
    DTDataFile outputFile("Output.dtbin",DTFile::NewReadWrite);

    // Output from computation
    outputFile.Save(computed,"Var");
    outputFile.Save("Real Number","Seq_Var");

    // The execution time.
    outputFile.Save(exec_time,"ExecutionTime");
    outputFile.Save("Real Number","Seq_ExecutionTime");

    // The errors.
    DTSaveError(outputFile,"ExecutionErrors");
    outputFile.Save("StringList","Seq_ExecutionErrors");

    outputFile.SaveIndex();

    return 0;
}

double Computation(double dist,const DTSeriesPointCollection3D &points_series)
{
    double toReturn;
    
    DTProgress progress;
    progress.UpdatePercentage(0.0f);
    int nt = points_series.HowManySaved();
    int n = points_series(0).NumberOfPoints();
    
    DTMutableDoubleArray interacted(n*(n-1)/2);
    interacted = 0;
    toReturn = 0;
    
    for (int i=0;i<nt;++i)
    {
        DTPointCollection3D pts = points_series(i);
        int c=0;
        for (int j=0;j<n;++j)
        {
            for (int k=j+1;k<n;++k)
            {
                if (Norm(pts(j)-pts(k)) < dist)
                {
                    interacted(c) = 1;
                }
                ++c;
            }
        }
        progress.UpdatePercentage(float(i)/float(nt));
    }
    
    // compute the average
    for (int i=0;i<interacted.Length();++i)
    {
        toReturn += interacted(i);
    }
    toReturn /= interacted.Length();
    
    return toReturn;
}
