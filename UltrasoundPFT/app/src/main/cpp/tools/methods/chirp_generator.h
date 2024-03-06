//
// Created by Dell on 2024/2/23.
//

#ifndef ULTRASOUNDPFT_CHIRP_GENERATOR_H
#define ULTRASOUNDPFT_CHIRP_GENERATOR_H

#include <stdio.h>
#include <vector>

enum SpeakerType{
    SPEAKER_TYPE_ZERO   = 0,
    /** Using top speaker only */
    SPEAKER_TYPE_TOP    = 1,
    /** Using bottom speaker only */
    SPEAKER_TYPE_BOTTOM = 2,
    /** Using top and bottom speaker */
    SPEAKER_TYPE_BOTH   = 3
};


typedef struct ChirpParameters{
    /** chirp generator paramaters */
    int sampleRate = 48000;
    double chirpStartFrq = 17000;
    double chirpStopFrq = 23000;
    double chirpDuration = 10 / 1000.0;
    double chirpFrameDuration = 10 / 1000.0;
    int chirpCircle = 1;
    bool applyHamming = false;
    SpeakerType speakerType = SPEAKER_TYPE_BOTTOM;

    /** chirp processing parameters */
    float interval_factor = 0.95;
    float skip_time = 0.25;  // second
    int chirp_samples = (int)(sampleRate * chirpDuration);
    int skip_samples = (int)(skip_time / chirpDuration * sampleRate * chirpDuration);

} ChirpParameters;



std::vector<double> genPCM16MonoToneBytes(ChirpParameters params);

std::vector<double> genPCM16MonoToneBytes_Sin(ChirpParameters params);

std::vector<double> cosToSin(const std::vector<double>& cos_sequence);

#endif //ULTRASOUNDPFT_CHIRP_GENERATOR_H
