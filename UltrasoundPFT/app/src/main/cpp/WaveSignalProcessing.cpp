//
// Created by Dell on 2024/2/23.
//

#include <string>
#include <stdio.h>
#include <vector>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <android/log.h>

#include "WaveSignalProcessing.h"

#include "tools/kiss-fft/kiss_fft.h"
#include "tools/kiss-fft/kiss_fftr.h"
#include "Utils.h"

extern const char* LOGTAG;

/** ===========================================================================================
 *  Private function
 * ============================================================================================*/


int read_wav_file(char *fname, WaveSignalStruct& waveSignal) {

    // open WAV file
    FILE *file = fopen(fname, "rb");
    if (file == NULL) {
        // 处理无法打开文件的情况
        __android_log_print(ANDROID_LOG_ERROR, LOGTAG,
                            "Failed to open file: %s", fname);
        return -1;
    }

    // Read WAV file header information
    if (fread(waveSignal.pHeader, sizeof(WavHeader), 1, file) != 1) {
        __android_log_print(ANDROID_LOG_ERROR, LOGTAG,
                            "Read header meaasge failed:%s", fname);
        fclose(file);
        return -1;
    }

    // wav message
    __android_log_print(ANDROID_LOG_DEBUG, LOGTAG,
                        "sampleRate: %d\t "
                        "audioFormat: %d\t "
                        "subchunk1Size: %d\t"
                        "channels: %d\t "
                        "bitDepth: %d\t "
                        "byteRate: %ld\t ",
                        waveSignal.pHeader->sampleRate,
                        waveSignal.pHeader->audioFormat,
                        waveSignal.pHeader->subchunk1Size,
                        waveSignal.pHeader->numChannels,
                        waveSignal.pHeader->bitsPerSample,
                        waveSignal.pHeader->byteRate);


    if (waveSignal.pHeader->audioFormat != 1) {
        __android_log_print(ANDROID_LOG_ERROR, LOGTAG,
                            "Not PCM, audioFormat is: %s", waveSignal.pHeader->audioFormat );
    }

    // Count the number of bytes for each sample
    int bytesPerSample = waveSignal.pHeader->bitsPerSample / 8;

    // Calculate the total number of bytes of sample data
    int dataSize = waveSignal.pHeader->subchunk2Size;
    __android_log_print(ANDROID_LOG_DEBUG, LOGTAG,
                        "wave file size: %d\t ", dataSize);

    // skip header message
    fseek(file, 44, SEEK_SET);
    // Calculate the length of audio data
    fseek(file, 0, SEEK_END);
    long fileDataSize = ftell(file) - 44;
    fseek(file, 44, SEEK_SET);

    __android_log_print(ANDROID_LOG_DEBUG, LOGTAG,
                        "read file data size (header and audio): %ld\t ", fileDataSize);

    // Read audio data
    std::vector<int8_t> leftChannel(fileDataSize / 2); // Assume stereo audio
    std::vector<int8_t> rightChannel(fileDataSize / 2);

    for (size_t i = 0; i < leftChannel.size(); ++i) {
        fread(&leftChannel[i], sizeof(int8_t), 1, file);
        fread(&rightChannel[i], sizeof(int8_t), 1, file);
    }

//    // read audio data
//    std::vector<int8_t> audioData(fileDataSize - 44);
//    size_t audioDataSize = fread(&audioData[0], sizeof(int8_t), fileDataSize - 44, file);
//    // Check if the read is successful
//    if (audioDataSize != (size_t)(fileDataSize - 44)) {
//        __android_log_print(ANDROID_LOG_ERROR, LOGTAG,
//                            "read audio data size is : %d", audioDataSize);
//        return -1;
//    }

    waveSignal.RxLeftChData = deepCopyVector(leftChannel);
    waveSignal.RxRightChData = deepCopyVector(rightChannel);

    // log print vector
    logVectorReverse_int8("RxLeftChData reverse", waveSignal.RxLeftChData);
    logVectorReverse_int8("RxRightChData reverse", waveSignal.RxRightChData);

    // close file
    fclose(file);
    return fileDataSize;
}



int signalPreProcessing(ChirpParameters& params,
                         WaveSignalStruct& waveSignal){
    // parameter check
    if (waveSignal.pHeader == NULL){
        __android_log_print(ANDROID_LOG_ERROR, LOGTAG,
                            "waveSignal.TxChirpSignal is NULL ");
        return -1;
    }

    if (waveSignal.TxChirpSignal.empty()){
        __android_log_print(ANDROID_LOG_ERROR, LOGTAG,
                            "waveSignal.TxChirpSignal is NULL ");
        return -1;
    }

    if (waveSignal.RxLeftChData.empty()){
        __android_log_print(ANDROID_LOG_ERROR, LOGTAG,
                            "waveSignal.RxLeftRightChData is NULL ");
        return -1;
    }

    // Skip Samples
    waveSignal.RxLeftChData.erase(waveSignal.RxLeftChData.begin(),
                                  waveSignal.RxLeftChData.begin() + params.skip_samples);
    waveSignal.RxRightChData.erase(waveSignal.RxRightChData.begin(),
                                   waveSignal.RxRightChData.begin() + params.skip_samples);
    //
    __android_log_print(ANDROID_LOG_DEBUG, LOGTAG,
                        "Left data size: %d, skip data size: %d",
                        waveSignal.RxLeftChData.size(), params.skip_samples);

    int Fs = 48000;
    int nfft = 128;
    int overlap = 64;


    double* cx_in_left = vectorInt8ToDouble(waveSignal.RxLeftChData);
    double* cx_in_right = vectorInt8ToDouble(waveSignal.RxRightChData);

    int stftOutSize = (int) (waveSignal.RxLeftChData.size() - nfft)/(nfft-overlap) * (nfft/2+1);

    waveSignal.RxLeftFFTData = new FFT_CPX[stftOutSize];  // channel 0
    waveSignal.RxRightFFTData = new FFT_CPX[stftOutSize];  // channel 1

    kiss_stftr(cx_in_left, (kiss_fft_cpx*)(waveSignal.RxLeftFFTData), Fs, nfft, overlap);
    kiss_stftr(cx_in_right, (kiss_fft_cpx*)(waveSignal.RxRightFFTData), Fs, nfft, overlap);

    delete[] cx_in_left;
    delete[] cx_in_right;

    logDoubleArray("Left STFT Data  ==>",
                   (kiss_fft_cpx*)(waveSignal.RxLeftFFTData) , stftOutSize);
    logDoubleArray("Right FFT Data ==>",
                   (kiss_fft_cpx*)(waveSignal.RxRightFFTData) , stftOutSize);










    return 0;
}

