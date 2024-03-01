//
// Created by Dell on 2024/2/23.
//

#include <string>
#include <stdio.h>
#include <stddef.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <android/log.h>

#include "WaveSignalProcessing.h"
#include "tools/common/Utils.h"
#include "tools/common//Logger.h"

//#include "tools/kiss-fft/kiss_fft.h"
//#include "tools/kiss-fft/kiss_fftr.h"
#include "tools/filter/filter.h"
#include "tools/methods/xcorr.h"
#include "tools/methods/find_peak.h"
#include "tools/methods/chirp_aligned.h"

extern const char* LOGTAG;

/** ===========================================================================================
 *  Private function
 * ============================================================================================*/


int read_wav_file(char *fname, WaveSignalStruct& waveSignal) {

    // open WAV file
    FILE *file = fopen(fname, "rb");
    if (file == nullptr) {
        // 处理无法打开文件的情况
        Loge("Failed to open file: %s", fname);
        return -1;
    }

    // Read WAV file header information
    if (fread(waveSignal.pHeader, sizeof(WavHeader), 1, file) != 1) {
        Loge("Read header meaasge failed:%s", fname);
        fclose(file);
        return -1;
    }

#if DEBUG
    // wav message
    __android_log_print(ANDROID_LOG_DEBUG, LOGTAG,
                        "sampleRate: %d\t "
                            "audioFormat: %d\t "
                            "subchunk1Size: %d\t"
                            "channels: %d\t "
                            "bitDepth: %d\t "
                            "byteRate: %d\t ",
                            waveSignal.pHeader->sampleRate,
                            waveSignal.pHeader->audioFormat,
                            waveSignal.pHeader->subchunk1Size,
                            waveSignal.pHeader->numChannels,
                            waveSignal.pHeader->bitsPerSample,
                            waveSignal.pHeader->byteRate);
#endif

    if (waveSignal.pHeader->audioFormat != 1) {
        Loge("Not PCM, audioFormat is: %s", waveSignal.pHeader->audioFormat );
    }

    // Calculate the total number of bytes of sample data
    int fileSize = waveSignal.pHeader->subchunk2Size;
    // skip header message
    fseek(file, 44, SEEK_SET);
    // Calculate the length of audio data
    fseek(file, 0, SEEK_END);
    long audioDataSize = ftell(file) - 44;
    fseek(file, 44, SEEK_SET);

    __android_log_print(ANDROID_LOG_DEBUG, LOGTAG,
                        "Wave File Size: %d (bytes) = Audio Data Size: %ld (bytes) + 44 bytes (heads) ",
                            fileSize, audioDataSize);

    // Read audio data: 8bytes -> 16bytes and 2channels, so "audioDataSize / 2 / 2 "
    std::vector<int16_t> leftChannel(audioDataSize / 2 / 2); // Assume stereo audio
    std::vector<int16_t> rightChannel(audioDataSize / 2 / 2);

    for (size_t i = 0; i < leftChannel.size(); ++i) {
        fread(&leftChannel[i], sizeof(int16_t), 1, file);
        fread(&rightChannel[i], sizeof(int16_t), 1, file);
    }

    waveSignal.RxLeftChData = deepCopyVector(leftChannel);
    waveSignal.RxRightChData = deepCopyVector(rightChannel);

    fclose(file);
    return audioDataSize;
}



int signalPreProcessing(ChirpParameters& params,
                         WaveSignalStruct& waveSignal){

    int resFlag = 0;
    /**
     * Input Parameters Check.
     * */
    if (waveSignal.pHeader == nullptr || waveSignal.TxChirpSignal.empty() || waveSignal.RxLeftChData.empty()){
        Loge("TxChirpSignal or RxLeftRightChData is empty.");
        return -1;
    }

    /**
     * Step1: Skip Start Samples.
     * */
    waveSignal.RxLeftChData.erase(waveSignal.RxLeftChData.begin(),
                                  waveSignal.RxLeftChData.begin() + params.skip_samples);
    waveSignal.RxRightChData.erase(waveSignal.RxRightChData.begin(),
                                   waveSignal.RxRightChData.begin() + params.skip_samples);
    waveSignal.RxLeftProcessingData = VInt16ToVDouble(waveSignal.RxLeftChData);
    waveSignal.RxRightProcessingData = VInt16ToVDouble(waveSignal.RxRightChData);

#if DEBUG
    Logd(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Skip Start Samples.")
    __android_log_print(ANDROID_LOG_DEBUG, LOGTAG,
                        "Left data size: %d;\t Right data size: %d;\t skip data size: %d",
                        waveSignal.RxLeftChData.size(), waveSignal.RxRightChData.size(),
                        params.skip_samples);
    logVector_double("After Skip Left Data (RxLeftProcessingData)\t ",
                     waveSignal.RxLeftProcessingData);
    logVector_double("After Skip Right Data (RxLeftProcessingData)\t ",
                     waveSignal.RxLeftProcessingData);
    saveChDataToCSV(waveSignal);    // save left and right data for verify.
    Logd("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<")
#endif


    /**
     * Step 2: Processing crossCorrelation--> Envelope --> Find Peaks--> Aligned
    * */
    std::vector<double> _lProcessingData(waveSignal.RxLeftProcessingData.size() + waveSignal.TxChirpSignal.size(), 0.0);
    std::vector<double> _rProcessingData(waveSignal.RxRightProcessingData.size() + waveSignal.TxChirpSignal.size(), 0.0);

    resFlag += crossCorrelation(waveSignal.RxLeftProcessingData, waveSignal.TxChirpSignal,
                                _lProcessingData);
    resFlag += crossCorrelation(waveSignal.RxRightProcessingData, waveSignal.TxChirpSignal,
                                _rProcessingData);
    if (resFlag < 0){
        __android_log_print(ANDROID_LOG_ERROR, LOGTAG,
                            "crossCorrelation return ERROR: %d", resFlag);
    }
//    waveSignal.RxLeftProcessingData = _lProcessingData;
//    waveSignal.RxRightProcessingData = _rProcessingData;
    /*********************************    Envelope   *****************************************/
    // Envelope
    std::vector<double> _lEnvelope(_lProcessingData.size(), 0.0);
    std::vector<double> _rEnvelope(_rProcessingData.size(), 0.0);

    calculateEnvelope(_lProcessingData, _lEnvelope);
    calculateEnvelope(_rProcessingData, _rEnvelope);

    waveSignal.RxLeftProcessingData = _lEnvelope;
    waveSignal.RxRightProcessingData = _rEnvelope;


#if DEBUG
    Logd(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> xCorr")
    __android_log_print(ANDROID_LOG_DEBUG, LOGTAG,
                        "xCorr ==> TxDataSize:%d, RxLeftChData size: %d, RxRightChData size: %d",
                        waveSignal.TxChirpSignal.size(),
                        waveSignal.RxLeftChData.size(),
                        waveSignal.RxRightChData.size());

    logVector_double("RxLeftProcessingData", waveSignal.RxLeftProcessingData);
    logVector_double("RxRightProcessingData", waveSignal.RxRightProcessingData);

    saveProcessingDataToCSV(waveSignal);
    Logd("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<")
#endif

    /*********************************    Find peaks   ***************************************/
    std::vector<int32_t> refPeakIndexs;
    // bottom mic data as reference data
    const std::vector<double> refPeakSignal = waveSignal.RxRightProcessingData;
    const std::vector<double> otherPeakSignal = waveSignal.RxLeftProcessingData;

    refPeakIndexs = findPeaks(refPeakSignal,
                              0.8*GetMax(refPeakSignal),
                              (int)(0.8*waveSignal.TxChirpSignal.size()));
    // check peak index interval: assert(waveSignal.TxChirpSignal.size() == diff(refPeakIndexs))
    std::vector<int32_t> _diff = diffSignal(refPeakIndexs);
    int32_t _indexMode = modeSignal(_diff);

    if (_indexMode != waveSignal.TxChirpSignal.size()){
        __android_log_print(ANDROID_LOG_ERROR, LOGTAG,
                            "peak index interval: %d", _indexMode);
    }
    logVector_int32("=========> refPeakIndexs", refPeakIndexs);

    /*********************************    Chirp Aligned   ***********************************/
    std::vector<std::vector<double> > out_alignedData_ref;
    std::vector<std::vector<double> > out_alignedData_other;
    std::vector<int32_t> _xPeakIndex={}; // not used
    std::vector<int32_t> out_xOffsetArray_ref;
    std::vector<int32_t> out_xOffsetArray_other;

    resFlag += chirpsAligned(refPeakSignal,
                              refPeakIndexs,
                              _xPeakIndex,
                              waveSignal.TxChirpSignal.size(),
                             out_alignedData_ref,
                             out_xOffsetArray_ref);

    // aligned top mic base on bottom mic
    resFlag += chirpsAligned(otherPeakSignal,
                             refPeakIndexs,
                             out_xOffsetArray_ref,
                             waveSignal.TxChirpSignal.size(),
                             out_alignedData_other,
                             out_xOffsetArray_other);
    if(resFlag < 0){
        __android_log_print(ANDROID_LOG_DEBUG, LOGTAG,
                            "crossCorrelation return ERROR: %d", resFlag);
    }

#if DEBUG
    __android_log_print(ANDROID_LOG_DEBUG, LOGTAG,
                        "_alignedData size %d, subSize:%d, out_xOffsetArray size:%d ",
                        out_alignedData_ref.size(), out_alignedData_ref[0].size(), out_xOffsetArray_ref.size());
    logVector_int32("out_xOffsetArray", out_xOffsetArray_ref);
    // save out_alignedData_ref
    const char* out_alignedData_ref_file = "/storage/emulated/0/Documents/alignedRefData.csv";
    saveVectorDataTo2D(out_alignedData_ref_file, out_alignedData_ref);

    __android_log_print(ANDROID_LOG_DEBUG, LOGTAG,
                        "_alignedData size %d, subSize:%d, out_xOffsetArray size:%d ",
                        out_alignedData_other.size(), out_alignedData_other[0].size(), out_xOffsetArray_other.size());
    logVector_double("out_alignedData_other[0]", out_alignedData_other[0]);

    // save out_alignedData_other
    const char* out_alignedData_other_file = "/storage/emulated/0/Documents/alignedOtherData.csv";
    saveVectorDataTo2D(out_alignedData_other_file, out_alignedData_other);

    logVector_int32("out_xOffsetArray_ref", out_xOffsetArray_ref);
    logVector_int32("out_xOffsetArray_other", out_xOffsetArray_other);
#endif

    /**
     * Step 3: Demodulation chirps
    * */
    // Find received signal start-end and ship samples for aligned chirps
    int _secondSkipSamples = (int32_t)(refPeakIndexs[0]);

    __android_log_print(ANDROID_LOG_DEBUG, LOGTAG,"_secondSkipSamples: %d", _secondSkipSamples);

    waveSignal.RxLeftChData.erase(waveSignal.RxLeftChData.begin(),
                                  waveSignal.RxLeftChData.begin() + _secondSkipSamples);
    waveSignal.RxRightChData.erase(waveSignal.RxRightChData.begin(),
                                   waveSignal.RxRightChData.begin() + _secondSkipSamples);

#if 0
    const char* filename11 = "/storage/emulated/0/Documents/leftChData.csv";
    saveVectorDataToCSV(filename11, waveSignal.RxRightChData);
    const char* filename = "/storage/emulated/0/Documents/rightChData.csv";
    saveVectorDataToCSV(filename, waveSignal.RxRightChData);
#endif
    // demodulate









            /**
     * STFT
     * */
//    int Fs = 48000;
//    int nfft = 128;
//    int overlap = 64;
//
//    double* cx_in_left = vectorInt16ToDouble(waveSignal.RxLeftChData);
//    double* cx_in_right = vectorInt16ToDouble(waveSignal.RxRightChData);
//
//    int stftOutSize = (int) (waveSignal.RxLeftChData.size() - nfft)/(nfft-overlap) * (nfft/2+1);
//
//    waveSignal.RxLeftFFTData = new FFT_CPX[stftOutSize];  // channel 0
//    waveSignal.RxRightFFTData = new FFT_CPX[stftOutSize];  // channel 1
//
//    kiss_stftr(cx_in_left, (kiss_fft_cpx*)(waveSignal.RxLeftFFTData), Fs, nfft, overlap);
//    kiss_stftr(cx_in_right, (kiss_fft_cpx*)(waveSignal.RxRightFFTData), Fs, nfft, overlap);
//
//    delete[] cx_in_left;
//    delete[] cx_in_right;
//
//    logDoubleArray("Left STFT Data  ==>",
//                   (kiss_fft_cpx*)(waveSignal.RxLeftFFTData) , stftOutSize);
//    logDoubleArray("Right FFT Data ==>",
//                   (kiss_fft_cpx*)(waveSignal.RxRightFFTData) , stftOutSize);

    /**
     * Filter
     * */

#if 0
    // Filter
    SOSSection Chebyshev2HPFSOS[] = {
            {{1.0000,   -1.0000,        0   },{  1.0000,   -0.6169,         0}},
            {{1.0000,   -1.9633,    1.0000  },{  1.0000,   -1.2174,    0.3877}},
            {{1.0000,   -1.8595,    1.0000  },{  1.0000,   -1.1698,    0.4084}},
            {{1.0000,   -1.7044,    1.0000  },{  1.0000,   -1.0962,    0.4411}},
            {{1.0000,   -1.5191,    1.0000  },{  1.0000,   -1.0040,    0.4836}},
            {{1.0000,   -1.3239,    1.0000  },{  1.0000,   -0.9014,    0.5334}},
            {{1.0000,   -1.1350,    1.0000  },{  1.0000,   -0.7964,    0.5882}},
            {{1.0000,   -0.9639,    1.0000  },{  1.0000,   -0.6960,    0.6463}},
            {{1.0000,   -0.8174,    1.0000  },{  1.0000,   -0.6060,    0.7065}},
            {{1.0000,   -0.6990,    1.0000  },{  1.0000,   -0.5309,    0.7681}},
            {{1.0000,   -0.6100,    1.0000  },{  1.0000,   -0.4741,    0.8312}},
            {{1.0000,   -0.5508,    1.0000  },{  1.0000,   -0.4387,    0.8963}},
            {{1.0000,   -0.5213,    1.0000  },{  1.0000,   -0.4274,    0.9644}}
    };
    int n_sections = 13;

    double* inDataLeft = vectorInt16ToDouble(waveSignal.RxLeftChData);
    double* inDataRight = vectorInt16ToDouble(waveSignal.RxRightChData);

    res += filter_signal(inDataLeft, waveSignal.RxLeftChData.size(),
                              Chebyshev2HPFSOS, n_sections,  Chebyshev2);
    res += filter_signal(inDataRight, waveSignal.RxLeftChData.size(),
                         Chebyshev2HPFSOS, n_sections,  Chebyshev2);
    if(res < 0){
        __android_log_print(ANDROID_LOG_ERROR, LOGTAG,
                            "crossCorrelation return ERROR: %d", res);
    }
//    waveSignal.RxLeftProcessingData = doubleTovector(inDataLeft, waveSignal.RxLeftChData.size());
//    waveSignal.RxRightProcessingData = doubleTovector(inDataRight, waveSignal.RxRightChData.size());



#if DEBUG
    Logd(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Filter")
    logVector_double("Filter Data RxLeftProcessingData. ",  waveSignal.RxLeftProcessingData );
    logVector_double("Filter Data RxRightProcessingData. ",  waveSignal.RxRightProcessingData );
//
//    const char *lFile = "/storage/emulated/0/Documents/filterLeft.csv";
//    saveVectorDataToCSV(lFile,  waveSignal.RxLeftProcessingData );
//
//    const char *rFile = "/storage/emulated/0/Documents/filterRight.csv";
//    saveVectorDataToCSV(rFile,  waveSignal.RxRightProcessingData );

    const char *txFile = "/storage/emulated/0/Documents/txChirp.csv";
    saveVectorDataToCSV(txFile,  waveSignal.TxChirpSignal );

    std::vector<double> TxChirp = int16_to_double(waveSignal.TxChirpSignal);
    logVector_double("Filter Data TxChirpSignal. ",  TxChirp);
    Logd("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<")
#endif

#endif








    return 0;
}

