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
#include <complex>


#include "WaveSignalProcessing.h"
#include "tools/common/Utils.h"
#include "tools/common//Logger.h"


#include "tools/filter/filter.h"
#include "tools/methods/xcorr.h"
#include "tools/methods/find_peak.h"
#include "tools/methods/chirp_aligned.h"
#include "tools/methods/mixer.h"
#include "tools/methods/dsp_processing.h"

extern const char* LOGTAG;

const char* BaseDir = "/storage/emulated/0/Documents/";

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

    if (waveSignal.pHeader->audioFormat != 1) {
        __android_log_print(ANDROID_LOG_DEBUG, LOGTAG,
                            "Not PCM, audioFormat is: %d", waveSignal.pHeader->audioFormat );
        return -1;
    }

    if (waveSignal.pHeader->sampleRate != 48000) {
        __android_log_print(ANDROID_LOG_DEBUG, LOGTAG,
                            "SampleRate is: %d", waveSignal.pHeader->sampleRate );
        return -1;
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

    __android_log_print(ANDROID_LOG_DEBUG, LOGTAG,
                        "Entry signalPreProcessing+++++++++++++++");
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


    /**
     * Step 2: Processing crossCorrelation--> Envelope --> Find Peaks--> Aligned
    * */
    std::vector<double> temp_lProcessingData(waveSignal.RxLeftProcessingData.size() + waveSignal.TxChirpSignal.size(), 0.0);
    std::vector<double> temp_rProcessingData(waveSignal.RxRightProcessingData.size() + waveSignal.TxChirpSignal.size(), 0.0);

    resFlag += crossCorrelation(waveSignal.RxLeftProcessingData, waveSignal.TxChirpSignal,
                                temp_lProcessingData);
    resFlag += crossCorrelation(waveSignal.RxRightProcessingData, waveSignal.TxChirpSignal,
                                temp_rProcessingData);
    if (resFlag < 0){
        __android_log_print(ANDROID_LOG_ERROR, LOGTAG,
                            "crossCorrelation return ERROR: %d", resFlag);
    }

    /*********************************    Envelope   *****************************************/
    // Envelope
    std::vector<double> temp_lEnvelope(temp_lProcessingData.size(), 0.0);
    std::vector<double> temp_rEnvelope(temp_rProcessingData.size(), 0.0);

    calculateEnvelope(temp_lProcessingData, temp_lEnvelope);
    calculateEnvelope(temp_rProcessingData, temp_rEnvelope);

    waveSignal.RxLeftProcessingData = deepCopyVector(temp_lEnvelope);
    waveSignal.RxRightProcessingData = deepCopyVector(temp_rEnvelope);

    // release temp vector memory
    temp_lEnvelope.clear();
    temp_lEnvelope.shrink_to_fit();
    temp_rEnvelope.clear();
    temp_rEnvelope.shrink_to_fit();
    temp_lProcessingData.clear();
    temp_lProcessingData.shrink_to_fit();
    temp_rProcessingData.clear();
    temp_rProcessingData.shrink_to_fit();

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
    const std::vector<double> refPeakSignal = waveSignal.RxLeftProcessingData;
    const std::vector<double> otherPeakSignal = waveSignal.RxRightProcessingData;

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
#if DEBUG
    logVector_int32("=========> refPeakIndexs", refPeakIndexs);
#endif

#if 0
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

    // aligned bottom mic base on top mic
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
#endif

#if 0
    __android_log_print(ANDROID_LOG_DEBUG, LOGTAG,
                        "_alignedData size %d, subSize:%d, out_xOffsetArray size:%d ",
                        out_alignedData_ref.size(), out_alignedData_ref[0].size(), out_xOffsetArray_ref.size());
    logVector_int32("out_xOffsetArray", out_xOffsetArray_ref);
    // save out_alignedData_ref
    const char* out_alignedData_ref_file = mergeStrings(BaseDir, "alignedRefData.csv");
    saveVectorDataTo2D(out_alignedData_ref_file, out_alignedData_ref);

    __android_log_print(ANDROID_LOG_DEBUG, LOGTAG,
                        "_alignedData size %d, subSize:%d, out_xOffsetArray size:%d ",
                        out_alignedData_other.size(), out_alignedData_other[0].size(), out_xOffsetArray_other.size());
    logVector_double("out_alignedData_other[0]", out_alignedData_other[0]);

    // save out_alignedData_other
    const char* out_alignedData_other_file = mergeStrings(BaseDir, "alignedOtherData.csv");
    saveVectorDataTo2D(out_alignedData_other_file, out_alignedData_other);

    logVector_int32("out_xOffsetArray_ref", out_xOffsetArray_ref);
    logVector_int32("out_xOffsetArray_other", out_xOffsetArray_other);
#endif

    /**
     * Step 3: Demodulation chirps
    * */
    // Find received signal start-end and ship samples for aligned chirps
    int _secondSkipSamples = (int32_t)(refPeakIndexs[0]);
    waveSignal.RxLeftChData.erase(waveSignal.RxLeftChData.begin(),
                                  waveSignal.RxLeftChData.begin() + _secondSkipSamples);
    waveSignal.RxRightChData.erase(waveSignal.RxRightChData.begin(),
                                   waveSignal.RxRightChData.begin() + _secondSkipSamples);

#if DEBUG
    Logd(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Skip Start Samples.")
    __android_log_print(ANDROID_LOG_DEBUG, LOGTAG,
                        "Left data size: %d;\t Right data size: %d;\t Total skip data size: %d",
                        waveSignal.RxLeftChData.size(), waveSignal.RxRightChData.size(),
                        params.skip_samples + _secondSkipSamples);

    logVector_int16("After Skip Left Data (RxLeftChData)\t ", waveSignal.RxLeftChData);
    logVector_int16("After Skip Right Data (RxRightChData)\t ", waveSignal.RxRightChData);

    saveChDataToCSV(waveSignal);    // save left and right data for verify.
    Logd("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<")
#endif

    /**************************** signal mixer and filter*********************************** */
    int _chipNums = (int) (waveSignal.RxLeftChData.size() / waveSignal.TxChirpSignal.size());
    __android_log_print(ANDROID_LOG_DEBUG, LOGTAG, "chirp numbers need mixer: %d", _chipNums);

    resFlag += signalMixer(waveSignal, _chipNums);
    if(resFlag < 0){
        __android_log_print(ANDROID_LOG_DEBUG, LOGTAG,
                            "signalMixer return ERROR: %d", resFlag);
    }

#if DEBUG
    Logd(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> signalMixer and Filter")
    logVector_double("waveSignal.RxLeftMixerData.real",  waveSignal.RxLeftMixerData.real);
    logVector_double("waveSignal.RxLeftMixerData.imag",  waveSignal.RxLeftMixerData.imag);
    logVector_double("waveSignal.RxRightMixerData.real",  waveSignal.RxRightMixerData.real);
    logVector_double("waveSignal.RxRightMixerData.imag",  waveSignal.RxRightMixerData.imag);

    const std::string filename1 = strscpString(waveSignal.pFilename, "_lMixerReal.csv");
    saveVectorDataToCSV(filename1.c_str(),   waveSignal.RxLeftMixerData.real);

    const std::string filename2 = strscpString(waveSignal.pFilename, "_lMixerImag.csv");
    saveVectorDataToCSV(filename2.c_str(),   waveSignal.RxLeftMixerData.imag);

    const std::string filename3 = strscpString(waveSignal.pFilename, "_rMixerReal.csv");
    saveVectorDataToCSV(filename3.c_str(),   waveSignal.RxRightMixerData.real);

    const std::string filename4 = strscpString(waveSignal.pFilename,  "_rMixerImag.csv");
    saveVectorDataToCSV(filename4.c_str(),   waveSignal.RxRightMixerData.imag);
    Logd("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<")
#endif

    std::vector<std::complex<double>> range_fft_result;
    int16_t nFFT = 480;
    resFlag += Range_FFT_for_RxRightMixerData(waveSignal, nFFT, range_fft_result);

    // mean
    std::vector<std::vector<std::complex<double>>> range_fft_result2D;
    range_fft_result2D = reshape(range_fft_result, nFFT, true);
    std::vector<double> range_fft_abs_mean;
    range_fft_abs_mean = columnWiseMean(range_fft_result2D);

#if DEBUG
    __android_log_print(ANDROID_LOG_DEBUG, LOGTAG,
                        "range_fft_result size: %d", range_fft_result.size());

    const std::string filename5 = strscpString(waveSignal.pFilename, "_range_fft_result2D.csv");
    saveComplexMatrixToCSV(filename5, range_fft_result2D);

    const std::string filename6 = strscpString(waveSignal.pFilename, "_range_fft_abs_mean.csv");
    saveVectorDataToCSV(filename6.c_str(), range_fft_abs_mean);
#endif

    // find peaks: target bin peak and references bin peak
    std::vector<int32_t> temp_peak;
    temp_peak = findPeaks(range_fft_abs_mean, 0.1*GetMax(range_fft_abs_mean), 4);
    if(temp_peak.size() < 2){
        __android_log_print(ANDROID_LOG_ERROR, LOGTAG,
                            "range_fft_abs_mean not find target peaks: %d", temp_peak.size());
        return -1;
    }
    // sort peaks
    std::sort(temp_peak.begin(), temp_peak.end(), [&](int32_t a, int32_t b){
        return range_fft_abs_mean[a] > range_fft_abs_mean[b];
    });

    // get doppler signal
    std::vector<std::complex<double>> temp_BottomFFT_bin_ref, temp_BottomFFT_bin_target;

    temp_BottomFFT_bin_ref = getColumn(range_fft_result2D, temp_peak[0]); // ref bin
    temp_BottomFFT_bin_target = getColumn(range_fft_result2D, temp_peak[1]); // target bin

#if DEBUG
    __android_log_print(ANDROID_LOG_DEBUG, LOGTAG,
                        "range_fft_result2D size: %d, range_fft_result2D[0].size(): %d",
                        range_fft_result2D.size(), range_fft_result2D[0].size());
    logVector_int32("peaks:", temp_peak);
#endif

    // release temp
    temp_peak.clear();
    temp_peak.shrink_to_fit();

    std::vector<double> temp_BottomFFT_ref_unwrapPhase = unwrapPhase(temp_BottomFFT_bin_ref);
    std::vector<double> temp_BottomFFT_target_unwrapPhase = unwrapPhase(temp_BottomFFT_bin_target);
    std::vector<double> temp_BottomFFT_ref_unwrapPhase_sub = subtractVectors(temp_BottomFFT_target_unwrapPhase, temp_BottomFFT_ref_unwrapPhase);
    waveSignal.RxRightProcessingData = vectorDiff(temp_BottomFFT_ref_unwrapPhase_sub);

#if DEBUG
    const std::string filename7 = strscpString(waveSignal.pFilename,  "_diff_angle.csv");
    saveVectorDataToCSV(filename7.c_str(),waveSignal.RxRightProcessingData);

    logVector_double("BottomFFT_ref_unwrapPhase", temp_BottomFFT_ref_unwrapPhase);
    logVector_double("temp_BottomFFT_target_unwrapPhase", temp_BottomFFT_target_unwrapPhase);
    logVector_double("temp_BottomFFT_ref_unwrapPhase_diff", temp_BottomFFT_ref_unwrapPhase_sub);
    logVector_double("phase_velocity", waveSignal.RxRightProcessingData);
#endif

    // release temporary variables
    temp_BottomFFT_ref_unwrapPhase.clear();
    temp_BottomFFT_ref_unwrapPhase.shrink_to_fit();
    temp_BottomFFT_target_unwrapPhase.clear();
    temp_BottomFFT_target_unwrapPhase.shrink_to_fit();
    temp_BottomFFT_ref_unwrapPhase_sub.clear();
    temp_BottomFFT_ref_unwrapPhase_sub.shrink_to_fit();

    // STFT for phase velocity
    int windowSize = 16;
    int hopSize = 8;
    std::vector<std::vector<std::complex<double>>> stft_result2D;
    stft_result2D = VectorSTFT(waveSignal.RxRightProcessingData, windowSize, hopSize);
#if DEBUG
    __android_log_print(ANDROID_LOG_DEBUG, LOGTAG,
                        "range_fft_result2D size: %d, range_fft_result2D[0].size(): %d",
                        stft_result2D.size(), stft_result2D[0].size());

    const std::string filename8 = strscpString(waveSignal.pFilename, "_stft_result2D.csv");
    saveComplexMatrixToCSV(filename8, stft_result2D);
#endif

    waveSignal.RxRightProcessingData = rowwiseAbsSum(stft_result2D);
    // calculate the velocity based on model struct
    CalAirflowVelocity(waveSignal.RxRightProcessingData);

    waveSignal.AirflowVelocity = deepCopyVector(waveSignal.RxRightProcessingData);

#if DEBUG
    const std::string filename9 = strscpString(waveSignal.pFilename, "_airflowVelocity.csv");
    saveVectorDataToCSV(filename9.c_str(), waveSignal.RxRightProcessingData);
    logVector_double("CalAirflowVelocity", waveSignal.RxRightProcessingData);
#endif

    __android_log_print(ANDROID_LOG_DEBUG, LOGTAG,
                        "Exit signalPreProcessing------------------");

    return 0;
}



