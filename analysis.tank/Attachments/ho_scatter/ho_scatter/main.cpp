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

#include "KDTree.h"

DTDoubleArray Computation(const DTSeriesPointCollection3D &points,double order,double samples);

int main(int argc,const char *argv[])
{
    DTSetArguments(argc,argv);

    DTDataFile inputFile("Input.dtbin",DTFile::ReadOnly);
    // Read in the input variables.
    DTSeriesPointCollection3D points;
    Read(inputFile,"points",points);
    double order = inputFile.ReadNumber("order");
    double samples = inputFile.ReadNumber("samples");

    // The computation.
    DTDoubleArray computed;
    clock_t t_before = clock();
    computed = Computation(points,order,samples);
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

DTDoubleArray Computation(const DTSeriesPointCollection3D &points,double order,double samples)
{
    
    // To start, need samples of order points that are close together
    
    // Create KDTree
    int tlist[] = {0,1000,2000,3000};
    const int tsteps = 1500;
    DTMutableDoubleArray tracks(4*361,tsteps);
    for (int startset=0;startset<4;++startset)
    {
        int startTime = tlist[startset];
        DTPointCollection3D start_points = points.Get(points.TimeNumber(startTime));
        int n = start_points.NumberOfPoints();
        vector<pt3d> kdpts;
        for (int i=0;i<n;++i)
        {
            pt3d p;
            p.id = i;
            p.data[0] = start_points(i).x;
            p.data[1] = start_points(i).y;
            p.data[1] = start_points(i).z;
            
            kdpts.push_back(p);
        }
        
        KDTree kd(kdpts);
        
        // get some samples
        vector<set<int> >  sample_sets;
        
        for (int i=0;i<n;++i)
        {
            KDTreeNodeAccess p = kd.find(i);
            auto res = kd.KNN(p,order);
            set<int> point_set;
            for (auto v : res)
            {
                point_set.insert(v.id());
            }
            sample_sets.push_back(std::move(point_set));
        }
        
        // create some tracks
        for (int t = 0; t < tsteps; t += 1)
        {
            auto cur_points= points.Get(points.TimeNumber(startTime+t));
            for (int j=0;j<sample_sets.size();++j)
            {
                // compute variance at this point
                DTPoint3D mean(0,0,0);
                for (auto id : sample_sets[j])
                {
                    mean = mean + cur_points(id);
                }
                mean = mean / order;
                double var = 0;
                for (auto id : sample_sets[j])
                {
                    double tmp = Norm(cur_points(id)-mean);
                    var += tmp*tmp;
                }
                tracks(361*startset+j,t) = var;
            }
        }
    }
    
    
    
    // try to normalize them to start around the same value
    double start = 0;
    vector<double> valuesToSort(tracks.m());
    for (int i=0;i<tracks.m();++i)
    {
        valuesToSort[i]=tracks(i,0);
    }
    sort(valuesToSort.begin(),valuesToSort.end());
    start = valuesToSort[int(0.5*valuesToSort.size())];
    
    // first count how many usable tracks we get
    int usable = 0;
    for (int i=0;i<tracks.m();++i)
    {
        if (tracks(i,0) <= start)
        {
            ++usable;
        }
    }
    
    DTMutableDoubleArray finalTracks(tracks.m(),tsteps);
    finalTracks=0;
    int c = 0;
    for (int i=0;i<tracks.m();++i)
    {
        if (tracks(i,0)>start)
        {
            continue;
        }
        // figure out where to start this track
        int startIndex = 0;
        while (startIndex < (tsteps/2) && tracks(i,startIndex)<=(start))
        {
            ++startIndex;
        }
        if (startIndex == tsteps)
        {
            continue;
        }
        for (int j=startIndex;j<tsteps;++j)
        {
            finalTracks(c,j-startIndex) = tracks(i,j);
        }
        ++c;
    }

    auto finalTracksTruncated = Region(finalTracks,DTRange(0,c),DTRange(0,tsteps/2));
    
    // compute average
    DTMutableDoubleArray average(tsteps/2);
    average = 0;
    for (int t=0;t<finalTracksTruncated.n();++t)
    {
        for (int k=0;k<finalTracksTruncated.m();++k)
        {
            average(t) += finalTracksTruncated(k,t);
        }
        average(t) /= finalTracksTruncated.m();
    }
    
    return average;
}
