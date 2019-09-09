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

DTDoubleArray Computation(double bins,double width,double threshold,double scale,
                          double start,const DTDoubleArray &partition,const DTSeriesPointCollection3D &points);

int main(int argc,const char *argv[])
{
    DTSetArguments(argc,argv);

    DTDataFile inputFile("Input.dtbin",DTFile::ReadOnly);
    // Read in the input variables.
    double bins = inputFile.ReadNumber("bins");
    double width = inputFile.ReadNumber("width");
    double threshold = inputFile.ReadNumber("threshold");
    double scale = inputFile.ReadNumber("scale");
    double start = inputFile.ReadNumber("start");
    DTDoubleArray partition = inputFile.ReadDoubleArray("partition");
    DTSeriesPointCollection3D points;
    Read(inputFile,"points",points);

    // The computation.
    DTDoubleArray computed;
    clock_t t_before = clock();
    computed = Computation(bins,width,threshold,scale,start,partition,points);
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

DTDoubleArray Computation(double bins,double width,double threshold,double scale,
                          double start,const DTDoubleArray &partition,const DTSeriesPointCollection3D &points)
{
    DTProgress progress;
    progress.UpdatePercentage(0);
    
    DTPointCollection3D p0 = points(0);
    int n = 2803;    // number of total beads
    int n2 = 361;    // number of nucleolus beads
    int nucleolus_offset = 1560;    // index of first nucleolus bead
    
    int nt = points.HowManySaved();
    
    vector<double> dists;
    
    
    for (int bin=0;bin<bins;++bin)
    {
        // track progress in DataTank
        progress.UpdatePercentage(float(bin)/float(bins));
        
        
        // iterate over all time steps within this bin
        for (int k=0;k<width;++k)
        {
            int t = start + bin*width + k;
            DTPointCollection3D p = points(t);
            // figure out cluster COMs in this bin
            map<int,DTPoint3D> com;
            map<int,int> Z;
            for (int i=0;i<n2;++i)
            {
                int community = partition(i,bin);
                com[community] = com[community] + p(i+nucleolus_offset);
                Z[community] += 1;
            }
            
            // figure out distances from COM, iterating over all beads
            for (int i=0;i<n;++i)
            {
                // iterate over all communities and find distance from their COM
                for (auto elem : Z)
                {
                    int community = elem.first;
                    
                    // check if we are working with a nucleolus bead and need to exclude our own community
                    if (i >= nucleolus_offset && i < nucleolus_offset + n2)
                    {
                        // check that we're not on our own community
                        if (community == partition(i-nucleolus_offset,bin))
                        {
                            continue;
                        }
                    }
                    
                    DTPoint3D dist = p(i) - com[community]/Z[community];
                    dists.push_back(Norm(dist));
                }
            }
        }
    }
    
    DTMutableDoubleArray toReturn(dists.size());
    for (int i=0;i<dists.size();++i)
    {
        toReturn(i) = dists[i];
    }
    
    return toReturn;
}
