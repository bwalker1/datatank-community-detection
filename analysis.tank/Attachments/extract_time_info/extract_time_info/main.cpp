// Include DTSource.h if you want to include all the headers.

#include "DTArguments.h"
#include "DTSaveError.h"

#include "DTDataFile.h"
#include "DTDoubleArray.h"
#include "DTPath3D.h"
#include "DTSeriesPath3D.h"

// Common utilities
#include "DTDoubleArrayOperators.h"
#include "DTProgress.h"
#include "DTTimer.h"
#include "DTUtilities.h"
#include "DTDictionary.h"

#include <math.h>

DTDoubleArray Computation(const DTSeriesPath3D &series);

int main(int argc,const char *argv[])
{
    DTSetArguments(argc,argv);

    DTDataFile inputFile("Input.dtbin",DTFile::ReadOnly);
    // Read in the input variables.
    DTSeriesPath3D series;
    Read(inputFile,"series",series);

    // The computation.
    DTDoubleArray computed;
    clock_t t_before = clock();
    computed = Computation(series);
    clock_t t_after = clock();
    double exec_time = double(t_after-t_before)/double(CLOCKS_PER_SEC);

    // Write the output.
    DTDataFile outputFile("Output.dtbin",DTFile::NewReadWrite);

    // Output from computation
    outputFile.Save(computed,"Var");
    outputFile.Save("NumberList","Seq_Var");

    // The execution time.
    outputFile.Save(exec_time,"ExecutionTime");
    outputFile.Save("Real Number","Seq_ExecutionTime");

    // The errors.
    DTSaveError(outputFile,"ExecutionErrors");
    outputFile.Save("StringList","Seq_ExecutionErrors");

    outputFile.SaveIndex();

    return 0;
}

DTDoubleArray Computation(const DTSeriesPath3D &series)
{
    DTMutableDoubleArray toReturn(2);
    
    double start_time = series.TimeNumber(0);
    double dt = series.TimeNumber(1)-start_time;
    
    toReturn(0) = start_time;
    toReturn(1) = dt;

    return toReturn;
}
