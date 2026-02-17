#ifndef KRmath
#define KRmath

#include <cmath>
#include <vector>

namespace KR {
    double darc(double x, double y){
        double thetax = std::acos(x);
        double thetay = std::asin(y);
        double theta;
        if (thetay >= 0){
            theta = thetax;
        }
        if (thetay < 0){
            theta = 2*M_PI-thetax;
        }
        return theta;
    }

    void HeapSortModCSwap(std::vector<int>& indexVec, std::vector<double>& MainVec, std::vector<double>& SubVec, int i, int length)
    {
        int swaptag = i;
        int lefti = i*2 + 1;
        int righti = i*2 + 2;
        if (!(lefti < length))
        {
            return;
        }
        if (lefti < length && (MainVec[lefti] > MainVec[swaptag] || (MainVec[lefti] == MainVec[swaptag] && SubVec[lefti] > SubVec[swaptag])))
        {
            swaptag = lefti;
        }
        if (righti < length && (MainVec[righti] > MainVec[swaptag] || (MainVec[righti] == MainVec[swaptag] && SubVec[righti] > SubVec[swaptag])))
        {
            swaptag = righti;
        }
        if (swaptag != i)
        {
            std::swap(indexVec[i],indexVec[swaptag]);
            std::swap(MainVec[i],MainVec[swaptag]);
            std::swap(SubVec[i],SubVec[swaptag]);
            KR::HeapSortModCSwap(indexVec,MainVec,SubVec,swaptag,length);
        }
    }

    std::vector<int> HeapSort(std::vector<double> MainVec, std::vector<double> SubVec, bool udFlag = false)
    {
        int n = static_cast<int>(MainVec.size());

        std::vector<int> indexVec;
        indexVec.clear();
        indexVec.resize(MainVec.size());
        for (int i=0; i<n; i++)
        {
            indexVec[i] = i;
        }

        for(int i=(n/2.0-1.0); i>=0; i--)
        {
            KR::HeapSortModCSwap(indexVec,MainVec,SubVec,i,n);
        }
        while(n > 0)
        {
            std::swap(indexVec[0],indexVec[n-1]);
            std::swap(MainVec[0],MainVec[n-1]);
            std::swap(SubVec[0],SubVec[n-1]);
            n--;
            KR::HeapSortModCSwap(indexVec,MainVec,SubVec,0,n);
        }
        if (udFlag)
        {
            int len = static_cast<int>(indexVec.size());
            for (int i=0; i<len/2; i++)
            {
                std::swap(indexVec[i],indexVec[len-1-i]);
            }
        }
        return indexVec;
    }
}

#endif