# -*- encoding: utf-8 -*-


import numpy as np
import scipy.signal as signal


from Signal_Processing.Common import Xcorr, Aligned_Chirps


def find_reference_peaks(received_signal, chirp_signal, startFreq, stopFreq, sampleRate, interval):
    """

    """
    # get top/bottom mic signal
    TopMicChannel = received_signal[0]
    BottomMicChannel = received_signal[1]

    # Judge the speaker type: top or bottom?
    _, TopMic_corr = Xcorr(TopMicChannel, chirp_signal, startFreq, stopFreq, sampleRate)
    TopMic_corr = np.abs(TopMic_corr)
    TopMic_Peaks, _ = signal.find_peaks(TopMic_corr, distance=interval)
    TopMic_corr_energy = np.mean(TopMic_corr[TopMic_Peaks])

    _, BottomMic_corr = Xcorr(BottomMicChannel, chirp_signal, startFreq, stopFreq, sampleRate)
    BottomMic_corr = np.abs(BottomMic_corr)
    BottomMic_Peaks, _ = signal.find_peaks(BottomMic_corr, distance=interval)
    BottomMic_corr_energy = np.mean(BottomMic_corr[BottomMic_Peaks])

    if TopMic_corr_energy > BottomMic_corr_energy:
        speakerType = "TopSpeaker"
    else:
        speakerType = "BottomSpeaker"

    # 3.Corr:  Channel 0
    if speakerType == "BottomSpeaker":
        # Find Reference Peak
        _, ref_corr = Xcorr(BottomMicChannel, chirp_signal, startFreq, stopFreq, sampleRate)
        ref_Peaks, ref_Peakvalue = signal.find_peaks(abs(ref_corr), distance=interval)
    elif speakerType == "TopSpeaker":
        # Find Reference Peak
        _, ref_corr = Xcorr(TopMicChannel, chirp_signal, startFreq, stopFreq, sampleRate)
        ref_corr = abs(ref_corr)
        ref_Peaks, _ = signal.find_peaks(ref_corr, distance=interval)
    else:
        print("Speaker type is not support......")
        return -1

    return speakerType, ref_Peaks, TopMic_corr, BottomMic_corr
