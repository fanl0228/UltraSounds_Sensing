//
// Created by Dell on 2024/2/29.
//

#ifndef ULTRASOUNDPFT_CHIRP_ALIGNED_H
#define ULTRASOUNDPFT_CHIRP_ALIGNED_H


int chirpsAligned(const std::vector<double>& inData,
                  const std::vector<int32_t>& peakIndex,
                  const std::vector<int32_t>& xPeakIndex,
                  const int32_t chirpSamples,
                  std::vector<std::vector<double> >& out_alignedData,
                  std::vector<int32_t>& out_xOffsetArray);

#endif //ULTRASOUNDPFT_CHIRP_ALIGNED_H
