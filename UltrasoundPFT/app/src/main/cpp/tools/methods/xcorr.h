//
// Created by Dell on 2024/2/26.
//

#ifndef ULTRASOUNDPFT_XCORR_H
#define ULTRASOUNDPFT_XCORR_H



void normalizeSignal(const std::vector<double>& signal,
                     std::vector<double>& normalizedSignal);

int crossCorrelation(const std::vector<double>& x,
                     const std::vector<double>& y,
                     std::vector<double>& res);


//int crossCorrelationFFT(const std::vector<double>& sig1,
//                        const std::vector<double>& sig2,
//                        std::vector<double>& result);


void calculateEnvelope(const std::vector<double>& signal,
                       std::vector<double>& envelope);

#endif //ULTRASOUNDPFT_XCORR_H
