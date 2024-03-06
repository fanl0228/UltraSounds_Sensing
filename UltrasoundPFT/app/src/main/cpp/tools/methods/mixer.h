//
// Created by Dell on 2024/3/1.
//

#ifndef ULTRASOUNDPFT_MIXER_H
#define ULTRASOUNDPFT_MIXER_H


#include "../../WaveSignalProcessing.h"

//fftw_complex * fttInput;
//fftw_complex * fttOutput;

std::vector<double> elementWiseMultiply(const std::vector<double>& vec1,
                                        const std::vector<double>& vec2);

int signalMixer(WaveSignalStruct& waveSignal, int chipNums);

#endif //ULTRASOUNDPFT_MIXER_H
