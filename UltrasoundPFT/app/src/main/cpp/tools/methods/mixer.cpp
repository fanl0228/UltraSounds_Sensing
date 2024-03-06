//
// Created by Dell on 2024/3/1.
//


#include <string>
#include <stdio.h>
#include <vector>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cstdint>
#include <android/log.h>
#include <fftw3.h>


#include "mixer.h"
#include "../../WaveSignalProcessing.h"
#include "../filter/filter.h"
#include "../common/Utils.h"

//extern const char * LOGTAG;


std::vector<double> elementWiseMultiply(const std::vector<double>& vec1,
                                        const std::vector<double>& vec2) {
    std::vector<double> result;
    result.reserve(std::min(vec1.size(), vec2.size()));

    size_t min_size = std::min(vec1.size(), vec2.size());
    for (size_t i = 0; i < min_size; ++i) {
        result.push_back(vec1[i] * vec2[i]);
    }
    return result;
}


int signalMixer(WaveSignalStruct& waveSignal, int chirpNums){

    if(waveSignal.TxChirpSignal.empty() ||
        waveSignal.RxLeftChData.empty() ||
        waveSignal.RxRightChData.empty()){
        __android_log_print(ANDROID_LOG_ERROR, "UltrasoundPFT: native --->",
                            "signalMixer input inData is empty.");
        return -1;
    }
    std::vector<double> cosTxChirp = waveSignal.TxChirpSignal;
    cosTxChirp = normalizeToMinusOneToOne(cosTxChirp); // Normalization
    ChirpParameters mChirpParams;
    std::vector<double> sinTxChirp = genPCM16MonoToneBytes_Sin(mChirpParams);
    sinTxChirp = normalizeToMinusOneToOne(sinTxChirp); // Normalization

    int TxchirpSize = waveSignal.TxChirpSignal.size();

    std::vector<double> lSlice;
    std::vector<double> lSlice_mixer_I, lSlice_mixer_Q;

    std::vector<double> rSlice;
    std::vector<double> rSlice_mixer_I, rSlice_mixer_Q;

    // Filter initialization
    // Left channel slice
    std::vector<double> out_lSlice_mixer_I, out_lSlice_mixer_Q;
    out_lSlice_mixer_I.resize(TxchirpSize, 0.0);
    out_lSlice_mixer_Q.resize(TxchirpSize, 0.0);

    // right channel slice
    std::vector<double> out_rSlice_mixer_I, out_rSlice_mixer_Q;
    out_rSlice_mixer_I.resize(TxchirpSize, 0.0);
    out_rSlice_mixer_Q.resize(TxchirpSize, 0.0);

    // Chebyshev2 LP-Filter initialization
    Filter mFilter(FilterStructType::Chebyshev2, FilterNameType::Low_Pass);

    for(int i = 0; i < chirpNums; i++){
        lSlice.clear();
        lSlice.assign(waveSignal.RxLeftChData.begin() + i * TxchirpSize,
                      waveSignal.RxLeftChData.begin() + (i+1) * TxchirpSize);

        // mixer I/Q
        lSlice_mixer_I.clear();
        lSlice_mixer_Q.clear();
        lSlice_mixer_I = elementWiseMultiply(lSlice, cosTxChirp);
        lSlice_mixer_Q = elementWiseMultiply(lSlice, sinTxChirp);

        // LPF
        mFilter.applyFilter(lSlice_mixer_I,out_lSlice_mixer_I);
        mFilter.applyFilter(lSlice_mixer_Q,out_lSlice_mixer_Q);

        // push back
        waveSignal.RxLeftMixerData.real.insert(waveSignal.RxLeftMixerData.real.end(),
                                               out_lSlice_mixer_I.begin(),
                                               out_lSlice_mixer_I.end());
        waveSignal.RxLeftMixerData.imag.insert(waveSignal.RxLeftMixerData.imag.end(),
                                               out_lSlice_mixer_Q.begin(),
                                               out_lSlice_mixer_Q.end());
#if 0
        logVector_double("lSlice_mixer_I ===>", lSlice_mixer_I);
        logVector_double("lSlice_mixer_I_Filter ===>", out_lSlice_mixer_I);
#endif

        rSlice.clear();
        rSlice.assign(waveSignal.RxRightChData.begin() + i * TxchirpSize,
                      waveSignal.RxRightChData.begin() + (i+1) * TxchirpSize);
        // mixer I/Q
        rSlice_mixer_I.clear();
        rSlice_mixer_Q.clear();
        rSlice_mixer_I = elementWiseMultiply(rSlice, cosTxChirp);
        rSlice_mixer_Q = elementWiseMultiply(rSlice, sinTxChirp);

        //LPF
        mFilter.applyFilter(rSlice_mixer_I,out_rSlice_mixer_I);
        mFilter.applyFilter(rSlice_mixer_Q,out_rSlice_mixer_Q);

        // push back
        waveSignal.RxRightMixerData.real.insert(waveSignal.RxRightMixerData.real.end(),
                                                out_rSlice_mixer_I.begin(),
                                                out_rSlice_mixer_I.end());
        waveSignal.RxRightMixerData.imag.insert(waveSignal.RxRightMixerData.imag.end(),
                                                out_rSlice_mixer_Q.begin(),
                                                out_rSlice_mixer_Q.end());
#if 0
        logVector_double("rSlice_mixer_I ===>", rSlice_mixer_I);
        logVector_double("rSlice_mixer_Q ===>", rSlice_mixer_Q);
#endif
    }

#if 0
    logVector_double(" waveSignal.RxLeftMixerData.imag ===>", waveSignal.RxLeftMixerData.imag);
    logVector_double(" waveSignal.RxLeftMixerData.real ===>", waveSignal.RxLeftMixerData.real);
#endif


#if 0
//    //LPF  10kHz
    std::vector<double> outData_left_real, outData_left_imag;
    std::vector<double> outData_right_real, outData_right_imag;
    outData_left_real.resize(waveSignal.RxLeftMixerData.real.size(), 0.0);
    outData_left_imag.resize(waveSignal.RxLeftMixerData.imag.size(), 0.0);
    outData_right_real.resize(waveSignal.RxRightMixerData.real.size(), 0.0);
    outData_right_imag.resize(waveSignal.RxRightMixerData.imag.size(), 0.0);

    Filter mFilter(FilterStructType::Chebyshev2, FilterNameType::Low_Pass);

    mFilter.applyFilter(waveSignal.RxLeftMixerData.real,outData_left_real);
    mFilter.applyFilter(waveSignal.RxLeftMixerData.imag,outData_left_imag);
    mFilter.applyFilter(waveSignal.RxRightMixerData.real,outData_right_real);
    mFilter.applyFilter(waveSignal.RxRightMixerData.imag,outData_right_imag);

    waveSignal.RxLeftMixerData.real = outData_left_real;
    waveSignal.RxLeftMixerData.imag = outData_left_imag;
    waveSignal.RxRightMixerData.real = outData_right_real;
    waveSignal.RxRightMixerData.imag = outData_right_imag;

#endif

    // Normalization
    waveSignal.RxLeftMixerData.real = normalizeToMinusOneToOne(waveSignal.RxLeftMixerData.real);
    waveSignal.RxLeftMixerData.imag = normalizeToMinusOneToOne(waveSignal.RxLeftMixerData.imag);
    waveSignal.RxRightMixerData.real = normalizeToMinusOneToOne(waveSignal.RxRightMixerData.real);
    waveSignal.RxRightMixerData.imag = normalizeToMinusOneToOne(waveSignal.RxRightMixerData.imag);

    // Check output
    if(waveSignal.RxLeftMixerData.real.size() != chirpNums * TxchirpSize ||
       waveSignal.RxLeftMixerData.imag.size() != chirpNums * TxchirpSize ||
       waveSignal.RxRightMixerData.real.size() != chirpNums * TxchirpSize ||
       waveSignal.RxRightMixerData.imag.size() != chirpNums * TxchirpSize)
    {
        __android_log_print(ANDROID_LOG_ERROR, "UltrasoundPFT: native --->",
                            "RxRightMixerData size: %d is not equal to RxRightChData size :%d.",
                            waveSignal.RxRightMixerData.real.size(), chirpNums * TxchirpSize);
        return -1;
    }

    return 0;

}







