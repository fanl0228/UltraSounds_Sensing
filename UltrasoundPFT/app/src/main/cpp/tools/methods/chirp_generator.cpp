//
// Created by Dell on 2024/2/23.
//
#include <vector>
#include <cmath>

#include "chirp_generator.h"


std::vector<double> genPCM16MonoToneBytes(ChirpParameters params) {

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

    std::vector<double> data;
    data.reserve(frameSampleCnt * params.chirpCircle);

    int idx = 0;
    for (const double& dVal : sample) {
//        int16_t val = static_cast<int16_t>(dVal * 32767);
//        data.push_back(val);
        data.push_back(dVal);
    }

    return data;
}



std::vector<double> genPCM16MonoToneBytes_Sin(ChirpParameters params) {

    int chirpSampleCnt = std::round(params.chirpDuration * params.sampleRate);
    int frameSampleCnt = std::round(params.chirpDuration * params.sampleRate);

    // 生成窗函数
    std::vector<double> factors;
    factors.reserve(chirpSampleCnt);
    for (int i = 0; i < chirpSampleCnt; i++) {
        if (params.applyHamming) {
            factors.push_back(0.5 * (1 - std::sin((2 * M_PI * i) / (chirpSampleCnt - 1))));
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
            double value = factors[i] * std::sin(2.0 * M_PI * i / (params.sampleRate / instfreq));
            sample.push_back(value);
        }
        for (int i = chirpSampleCnt; i < frameSampleCnt; i++) {
            sample.push_back(0.0);
        }
        t++;
    }

    std::vector<double> data;
    data.reserve(frameSampleCnt * params.chirpCircle);

    int idx = 0;
    for (const double& dVal : sample) {
//        int16_t val = static_cast<int16_t>(dVal * 32767);
//        data.push_back(val);
        data.push_back(dVal);
    }

    return data;
}


std::vector<double> cosToSin(const std::vector<double>& cos_sequence) {
    std::vector<double> sin_sequence;
    sin_sequence.reserve(cos_sequence.size());

    for (double cos_val : cos_sequence) {
        double sin_val = std::sin(M_PI / 2 - std::acos(cos_val));
        sin_sequence.push_back(sin_val);
    }

    return sin_sequence;
}


