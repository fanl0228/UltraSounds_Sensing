//
// Created by Dell on 2024/2/23.
//
#include <vector>
#include <cmath>

#include "chirp_generator.h"


std::vector<int8_t> genPCM16MonoToneBytes(ChirpParameters params) {

    int chirpSampleCnt = std::round(params.chirpDuration * params.sampleRate);
    int frameSampleCnt = std::round(params.chirpDuration * params.sampleRate);

    // 生成窗函数
    std::vector<double> factors;
    factors.reserve(chirpSampleCnt);
    for (int i = 0; i < chirpSampleCnt; i++) {
        if (params.applyHamming) {
            factors.push_back(0.5 * (1 - std::cos((2 * M_PI * i) / (chirpSampleCnt - 1))));
        } else {
            factors.push_back(1.0);
        }
    }

    std::vector<double> sample;
    sample.reserve(frameSampleCnt * params.chirpCircle);

    int t = 0;
    while (t < params.chirpCircle) {
        double instfreq, numerator;
        for (int i = 0; i < chirpSampleCnt; i++) {
            numerator = static_cast<double>(i) / chirpSampleCnt;
            instfreq = params.chirpStartFrq +
                    (0.5 * numerator * (params.chirpStopFrq - params.chirpStartFrq) );

            // apply windows function
            double value = factors[i] * std::cos(2.0 * M_PI * i / (params.sampleRate / instfreq));
            sample.push_back(value);
        }
        for (int i = chirpSampleCnt; i < frameSampleCnt; i++) {
            sample.push_back(0.0);
        }
        t++;
    }

    std::vector<int8_t> data;
    data.reserve(2 * 2 * frameSampleCnt * params.chirpCircle);

    int idx = 0;
    for (const double& dVal : sample) {
        if (params.speakerType == SPEAKER_TYPE_BOTTOM) {
            // left channel
            data.push_back(0);
            data.push_back(0);
            // right channel
            short val = static_cast<short>(dVal * 32767);
            data.push_back(val & 0xFF);
            data.push_back((val & 0xFF00) >> 8);
        } else if (params.speakerType == SPEAKER_TYPE_TOP) {
            // left channel
            short val = static_cast<short>(dVal * 32767);
            data.push_back(val & 0xFF);
            data.push_back((val & 0xFF00) >> 8);
            // right channel
            data.push_back(0);
            data.push_back(0);
        } else if (params.speakerType == SPEAKER_TYPE_BOTH) {
            // left channel
            short val = static_cast<short>(dVal * 32767);
            data.push_back(val & 0xFF);
            data.push_back((val & 0xFF00) >> 8);
            // right channel
            data.push_back(val & 0xFF);
            data.push_back((val & 0xFF00) >> 8);
        } else {
            // others default: 0
            data.push_back(0);
            data.push_back(0);
            data.push_back(0);
            data.push_back(0);
            break;
        }
    }

    return data;
}


