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

#include <math.h>

#include "ripser.h"

DTDoubleArray Computation(double embed_points,double embed_delay,const DTSeriesPointCollection3D &points);

int main(int argc,const char *argv[])
{
    DTSetArguments(argc,argv);

    DTDataFile inputFile("Input.dtbin",DTFile::ReadOnly);
    // Read in the input variables.
    double embed_points = inputFile.ReadNumber("embed_points");
    double embed_delay = inputFile.ReadNumber("embed_delay");
    DTSeriesPointCollection3D points;
    Read(inputFile,"points",points);

    // The computation.
    DTDoubleArray computed;
    clock_t t_before = clock();
    computed = Computation(embed_points,embed_delay,points);
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

DTDoubleArray Computation(double embed_points,double embed_delay,const DTSeriesPointCollection3D &points)
{
    const int n = points.Get(points.TimeNumber(0)).NumberOfPoints();
    DTMutableDoubleArray embeddings(3*embed_points,n);
    
    for (int i=0;i<embed_points;i += embed_delay)
    {
        DTFloatArray cur_points = points.Get(points.TimeNumber(i)).FloatData();
        for (int j=0;j<n;++j)
        {
            for (int k=0;k<3;++k)
            {
                embeddings(k+3*i,j) = cur_points(k,j);
            }
        }
    }
    
    std::vector<value_t> distances;
    for (int i=1;i<n;++i)
    {
        for (int j=0;j<i;++j)
        {
            double dist = 0;
            for (int k=0;k<3*embed_points;++k)
            {
                double tmp = embeddings(k,i)-embeddings(k,j);
                dist += tmp*tmp;
            }
            dist = sqrt(dist);
            distances.push_back(dist);
        }
    }
    auto dist = compressed_lower_distance_matrix(std::move(distances));
    ripser<compressed_lower_distance_matrix>(std::move(dist), 1, 1000, 1,2).compute_barcodes();

    return embeddings;
}
