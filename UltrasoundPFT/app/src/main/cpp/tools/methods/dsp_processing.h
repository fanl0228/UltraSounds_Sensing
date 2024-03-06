//
// Created by Dell on 2024/3/4.
//

#ifndef ULTRASOUNDPFT_DSP_PROCESSING_H
#define ULTRASOUNDPFT_DSP_PROCESSING_H


#include "../../WaveSignalProcessing.h"

#define MODEL_PATH4  0.35
#define MODEL_PATH2  0.15
#define MODEL_RADICAL  0.2
#define VOICE_VELOCITY  346
#define START_FRE 17000

extern const char* LOGTAG;

int
range_fft_for_RxRightChData(WaveSignalStruct& waveSignal, int N,
                            std::vector<std::complex<double>>& fft_result);

int
Range_FFT_for_SignalChirp(std::vector<double>& inData_real,std::vector<double>& inData_imag,
                          int N, bool isFFTshift,
                          std::vector<std::complex<double>>& fft_result);

int
Range_FFT_for_RxRightMixerData(WaveSignalStruct& waveSignal,
                               int N, std::vector<std::complex<double>>& fft_result);

void
fftShift(std::vector<std::complex<double>>& spectrum);

std::vector<std::vector<double>>
reshape(const std::vector<double>& input, size_t cols);

std::vector<std::vector<std::complex<double>>>
reshape(const std::vector<std::complex<double>>& input, size_t cols, bool isFFTshift);

std::vector<double>
rowWiseMean(const std::vector<std::vector<std::complex<double>>>& input);

std::vector<double>
columnWiseMean(const std::vector<std::vector<std::complex<double>>>& input);

std::vector<double>
unwrapPhase(const std::vector<std::complex<double>>& complexVector);

std::vector<double>
subtractVectors(const std::vector<double>& vec1, const std::vector<double>& vec2);

std::vector<double>
vectorDiff(const std::vector<double>& signal);

std::vector<std::complex<double>>
getColumn(const std::vector<std::vector<std::complex<double>>>& matrix, int column);

std::vector<double>
hannWindow(int windowSize);

std::vector<std::vector<std::complex<double>>>
VectorSTFT(const std::vector<double>& signal, int windowSize, int hopSize);

std::vector<double>
ColumnwiseAbsSum(const std::vector<std::vector<std::complex<double>>>& data);

std::vector<double>
rowwiseAbsSum(const std::vector<std::vector<std::complex<double>>>& data);

int
CalAirflowVelocity(std::vector<double>& data);



#endif //ULTRASOUNDPFT_DSP_PROCESSING_H
