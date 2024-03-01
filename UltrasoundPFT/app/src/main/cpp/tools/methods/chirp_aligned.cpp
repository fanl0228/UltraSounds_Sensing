//
// Created by Dell on 2024/2/29.
//
#include <vector>
#include "chirp_aligned.h"
#include "find_peak.h"

int chirpsAligned(const std::vector<double>& inData,
                  const std::vector<int32_t>& peakIndex,
                  const std::vector<int32_t>& xPeakIndex,
                  const int32_t chirpSamples,
                  std::vector<std::vector<double> >& out_alignedData,
                  std::vector<int32_t>& out_xOffsetArray){

    std::vector<double> corrSlice;
    std::vector<int32_t> _peaks1;

    for(int i = 0; i < peakIndex.size(); i++){
        corrSlice.clear();
        corrSlice.assign(inData.begin() + peakIndex[i] - (int32_t)(chirpSamples*0.25),
                         inData.begin() + peakIndex[i] + (int32_t)(chirpSamples*0.75));

        if(corrSlice.size() != chirpSamples){
            int32_t _gapSamples = chirpSamples - corrSlice.size();
            corrSlice.resize(corrSlice.size() + _gapSamples, 0.0);
        }

        int32_t xOffset = 0;
        if(xPeakIndex.empty()){
            _peaks1 = findPeaks(corrSlice,0.8*GetMax(corrSlice),(int32_t)(0.8*chirpSamples));
            xOffset = _peaks1[0];
        }else{
            xOffset = xPeakIndex[i];
        }

        std::vector<double> _corrSliceCopy(corrSlice.begin() + xOffset, corrSlice.end());

        if(_corrSliceCopy.size() < chirpSamples){
            int32_t _gapSamples = chirpSamples - _corrSliceCopy.size();
            _corrSliceCopy.resize(_corrSliceCopy.size() + _gapSamples, 0.0);
        }else if(_corrSliceCopy.size() > chirpSamples){
            int32_t _gapSamples = _corrSliceCopy.size() - chirpSamples;
            _corrSliceCopy.resize(_corrSliceCopy.size() - _gapSamples, 0.0);
        }

        out_alignedData.push_back(_corrSliceCopy);
        out_xOffsetArray.push_back(xOffset);
    }


    return 0;
}


