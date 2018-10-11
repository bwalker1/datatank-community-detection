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
    double toReturn = 0;
    double Z = 0;
    
    DTProgress progress;
    
    int nt = points_series.HowManySaved();
    int n = points_series(0).NumberOfPoints();
    
    DTMutableDoubleArray startingTime(n,n);
    startingTime = -2;
    
    for (int t=0;t<nt;++t)
    {
        progress.UpdatePercentage(float(t)/float(nt));
        int c = 0;
        auto pts = points_series(t);
        
        for (int i=0;i<n;++i)
        {
            for (int j=i+1;j<n;++j)
            {
                double d = Norm(pts(i)-pts(j));
                if (d < dist)
                {
                    if (startingTime(i,j) >= 0)
                    {
                        // ending apartness
                        toReturn += points_series.TimeNumber(t)-points_series.TimeNumber(startingTime(i,j));
                        Z += 1;
                    }
                    startingTime(i,j)=-1;
                }
                else
                {
                    if (startingTime(i,j) > 0 || startingTime(i,j) == -2)
                    {
                        // start of apartness
                        startingTime(i,j) = t;
                    }
                }
            }
        }
    }

    toReturn /= (Z);
    
    if (Z==0) toReturn = -1;
    return toReturn;
}
