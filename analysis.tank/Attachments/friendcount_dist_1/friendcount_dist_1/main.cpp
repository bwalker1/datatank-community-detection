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

double Computation(const DTSeriesPointCollection3D &points_series,double dist);

int main(int argc,const char *argv[])
{
    DTSetArguments(argc,argv);

    DTDataFile inputFile("Input.dtbin",DTFile::ReadOnly);
    // Read in the input variables.
    DTSeriesPointCollection3D points_series;
    Read(inputFile,"points_series",points_series);
    double dist = inputFile.ReadNumber("dist");

    // The computation.
    double computed;
    clock_t t_before = clock();
    computed = Computation(points_series,dist);
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

double Computation(const DTSeriesPointCollection3D &points_series,double dist)
{
    DTProgress progress;
    double toReturn = 0;
    
    int nt = points_series.HowManySaved();
    int n = points_series(0).NumberOfPoints();
    
    for (int t=0;t<nt;++t)
    {
        progress.UpdatePercentage(float(t)/float(nt));
        auto pts = points_series(t);
        for (int i=0;i<n;++i)
        {
            for (int j=0;j<n;++j)
            {
                if (i==j) continue;
                if (Norm(pts(i)-pts(j))<dist)
                {
                    ++toReturn;
                }
            }
        }
    }
    
    toReturn /= (nt*n);
    
    return toReturn;
}
