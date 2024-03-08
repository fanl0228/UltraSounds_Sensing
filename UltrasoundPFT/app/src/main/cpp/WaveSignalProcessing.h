//
// Created by Dell on 2024/2/23.
//
#ifndef ULTRASOUNDPFT_WAVESIGNALPROCESSING_H
#define ULTRASOUNDPFT_WAVESIGNALPROCESSING_H

#include "tools/methods/chirp_generator.h"




// WAV header struct
struct WavHeader{
    char chunkID[4]             = {'R', 'I', 'F', 'F'};          // chunkID
    uint32_t chunkSize          = sizeof(WavHeader) - 8;         // chunkSize，Subtract the size of chunkID and chunkSize itself
    char format[4]              = {'W', 'A', 'V', 'E'};          // format
    char subchunk1ID[4]         = {'f', 'm', 't', ' '};          // subchunk1ID
    uint32_t subchunk1Size      = 16;                            // subchunk1Size, PCM format is 16
    uint16_t audioFormat        = 1;                             // audioFormat, PCM = 1
    uint16_t numChannels        = 2;                             // numChannels, studio = 2
    uint32_t sampleRate         = 44100;                         // sampleRate
    uint32_t byteRate           = 44100 * 2 * sizeof(int16_t);   // byteRate, Sampling rate * number of channels * sample size / 8
    uint16_t blockAlign         = 2 * sizeof(int16_t);           // blockAlign，Number of channels * sample size / 8
    uint16_t bitsPerSample      = 16;                            // bitsPerSample, Sample size
    char subchunk2ID[4]         = {'d', 'a', 't', 'a'};          // subchunk2ID
    uint32_t subchunk2Size      = 0;                             // subchunk2Size, Set it to 0 first, then fill in the actual size

    WavHeader(){}

    // Constructor
    WavHeader(char* fileNa, const char* chunkId, uint32_t chunkSz, const char* fmt, const char* sub1Id,
              uint32_t sub1Sz, uint16_t audioFmt, uint16_t numChnls, uint32_t smpRate,
              uint32_t byteRt, uint16_t blkAlgn, uint16_t btsPrSmp, const char* sub2Id,
              uint32_t sub2Sz)
            : chunkSize(chunkSz), subchunk1Size(sub1Sz), audioFormat(audioFmt), numChannels(numChnls),
              sampleRate(smpRate), byteRate(byteRt), blockAlign(blkAlgn), bitsPerSample(btsPrSmp),
              subchunk2Size(sub2Sz) {
        memcpy(chunkID, chunkId, 4);
        memcpy(format, fmt, 4);
        memcpy(subchunk1ID, sub1Id, 4);
        memcpy(subchunk2ID, sub2Id, 4);
    }
};
typedef struct WavHeader WavHeader;


typedef struct {
    double r;
    double i;
} FFT_CPX;


struct ComplexVector{
    std::vector<double> real;
    std::vector<double> imag;

    ComplexVector(){
        real = {};
        imag = {};
    }

//    ~ComplexVector(){
//        real.clear();
//        real.shrink_to_fit();
//        imag.clear();
//        imag.shrink_to_fit();
//    }

};

typedef struct ComplexVector ComplexVector;

struct WaveSignalStruct {
    char* pFilename = nullptr;

    WavHeader* pHeader = nullptr;

    // transmitter signal
    std::vector<double>& TxChirpSignal;

    // Channel left/right received raw data
    std::vector<int16_t>& RxLeftChData;
    std::vector<int16_t>& RxRightChData;

    // Channel left/right processing signal
    std::vector<double>& RxLeftProcessingData;
    std::vector<double>& RxRightProcessingData;

    // Channel left/right Mixer data
    ComplexVector& RxLeftMixerData;
    ComplexVector& RxRightMixerData;

    //
    std::vector<double>& AirflowVelocity;

    // channel 0 received signal FFT data
    FFT_CPX* RxLeftFFTData = nullptr;
    FFT_CPX* RxRightFFTData = nullptr;


    ~WaveSignalStruct() {
        __android_log_print(ANDROID_LOG_ERROR, "native processing--->",
                            "WaveSignalStruct release is called.");

        delete pFilename;

        delete RxLeftFFTData;
        delete RxRightFFTData;
    }

    // init
    WaveSignalStruct(char* filename,
                     WavHeader* header,
                     std::vector<double>& txSignal,
                     std::vector<int16_t>& leftSignal,
                     std::vector<int16_t>& rightSignal,
                     std::vector<double>& leftProcess,
                     std::vector<double>& rightProcess,
                     ComplexVector& lMixerData,
                     ComplexVector& rMixerData,
                     std::vector<double>& velocity
                     )
                        :pFilename(filename),
                        pHeader(header),
                         TxChirpSignal(txSignal),
                         RxLeftChData(leftSignal),
                         RxRightChData(rightSignal),
                         RxLeftProcessingData(leftProcess),
                         RxRightProcessingData(rightProcess),
                         RxLeftMixerData(lMixerData),
                         RxRightMixerData(rMixerData),
                         AirflowVelocity(velocity)
                         {

        // channel 0 received signal FFT data
        RxLeftFFTData = nullptr;

        // channel 1 received signal FFT data
        RxRightFFTData = nullptr;

    }

};
typedef struct WaveSignalStruct WaveSignalStruct;


int read_wav_file(char *fname, WaveSignalStruct& waveSignal);

int signalPreProcessing(ChirpParameters& params, WaveSignalStruct& waveSignal);

#endif //ULTRASOUNDPFT_WAVESIGNALPROCESSING_H
