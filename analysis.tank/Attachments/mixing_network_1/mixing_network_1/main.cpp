// Include DTSource.h if you want to include all the headers.

#include "DTArguments.h"
#include "DTSaveError.h"

#include "DTDataFile.h"
#include "DTDoubleArray.h"

// Common utilities
#include "DTDoubleArrayOperators.h"
#include "DTProgress.h"
#include "DTTimer.h"
#include "DTUtilities.h"
#include "DTDictionary.h"

#include <math.h>

double Computation(const DTDoubleArray &clusters);

int main(int argc,const char *argv[])
{
    DTSetArguments(argc,argv);

    DTDataFile inputFile("Input.dtbin",DTFile::ReadOnly);
    // Read in the input variables.
    DTDoubleArray clusters = inputFile.ReadDoubleArray("clusters");

    // The computation.
    double computed;
    clock_t t_before = clock();
    computed = Computation(clusters);
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

double Computation(const DTDoubleArray &clusters)
{
    double toReturn = 0;
    
    int n = clusters.m();
    int nt = clusters.n();
    
    DTMutableDoubleArray interacted(n*(n-1)/2);
    interacted=0;
    
    for (int t=0;t<nt;++t)
    {
        int c = 0;
        for (int i=0;i<n;++i)
        {
            for (int j=i+1;j<n;++j)
            {
                if (clusters(j,t)==clusters(i,t))
                {
                    interacted(c)=1;
                }
                ++c;
            }
        }
    }
    
    for (int i=0;i<n*(n-1)/2;++i)
    {
        toReturn += interacted(i);
    }
    
    toReturn /= (n*(n-1)/2);

    return toReturn;
}
