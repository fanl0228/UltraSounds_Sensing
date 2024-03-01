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


#include "mixer.h"
#include "../../WaveSignalProcessing.h"

extern const char * LOGTAG;

int signalMixer(std::vector<std::vector<double> > &inData,
                WaveSignalStruct& waveSignal){

    if(inData.empty()){
        __android_log_print(ANDROID_LOG_ERROR, LOGTAG,
                            "signalMixer input inData is empty.");
        return -1;
    }
    std::vector<double> cosTxChirp = waveSignal.TxChirpSignal;
    std::vector<double> sinTxChirp = cosToSin(waveSignal.TxChirpSignal);

    for(int i = 0; i < inData.size(); i++){

    }



    return 0;
}







