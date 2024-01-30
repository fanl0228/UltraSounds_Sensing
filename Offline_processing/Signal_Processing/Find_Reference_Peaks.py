# -*- encoding: utf-8 -*-


import numpy as np
import scipy.signal as signal


from Signal_Processing.Common import Xcorr, Aligned_Chirps



def find_reference_peaks_stereo(received_signal, chirp_signal_top, chirp_signal_bottom, startFreq, stopFreq, sampleRate, interval):
    """

    """
    # top speaker
    speaker_type, \
    ref_peaks, other_peaks, \
    TopMic_corr, BottomMic_corr, \
    TopMic_corr_filter, BottomMic_corrfilter = find_reference_peaks(received_signal, chirp_signal_top,
                                                                    startFreq, stopFreq, sampleRate, interval)

    # bottom speaker
    # speaker_type, \
    # ref_peaks, other_peaks, \
    # TopMic_corr, BottomMic_corr, \
    # TopMic_corr_filter, BottomMic_corrfilter = find_reference_peaks(received_signal, chirp_signal_bottom,
    #                                                                 startFreq, stopFreq, sampleRate, interval)
    #

    return speaker_type, \
           ref_peaks, other_peaks, \
           TopMic_corr, BottomMic_corr, \
           TopMic_corr_filter, BottomMic_corrfilter




def find_reference_peaks(received_signals, chirp_signal, startFreq, stopFreq, sampleRate, interval):
    """

    """

    # get top/bottom mic signal
    TopMicChannel = received_signals[0]
    BottomMicChannel = received_signals[1]
    # TopMicChannel = received_signal[0][10000:-10000]
    # BottomMicChannel = received_signal[1][10000:-10000]
    # TopMicChannel = TopMicChannel[196:-207]
    # BottomMicChannel = BottomMicChannel[196:-207]

    ellip_filter = signal.iirfilter(25, Wn=10000, rp=1, rs=60, btype='highpass',
                                    analog=False, ftype='cheby2', fs=48000, output='sos')

    TopMicChannel = signal.sosfilt(ellip_filter, TopMicChannel)
    BottomMicChannel = signal.sosfilt(ellip_filter, BottomMicChannel)

    # Judge the speaker type: top or bottom?
    _, TopMic_corr = Xcorr(TopMicChannel, chirp_signal, startFreq, stopFreq, sampleRate)
    TopMic_corr = np.abs(TopMic_corr)
    TopMic_Peaks, _ = signal.find_peaks(TopMic_corr, distance=0.8*interval, height=0.05*np.max(TopMic_corr))
    TopMic_corr_energy = np.mean(TopMic_corr[TopMic_Peaks])

    _, BottomMic_corr = Xcorr(BottomMicChannel, chirp_signal, startFreq, stopFreq, sampleRate)
    BottomMic_corr = np.abs(BottomMic_corr)
    BottomMic_Peaks, _ = signal.find_peaks(BottomMic_corr, distance=0.8*interval, height=0.05*np.max(BottomMic_corr))
    BottomMic_corr_energy = np.mean(BottomMic_corr[BottomMic_Peaks])

    if TopMic_corr_energy > BottomMic_corr_energy:
        speaker_type = "TopSpeaker"
    else:
        speaker_type = "BottomSpeaker"

    # 3. Find reference peak bath top and bottom
    _, bottom_corr = Xcorr(BottomMicChannel, chirp_signal, startFreq, stopFreq, sampleRate)
    bottom_corr_hilbert = signal.hilbert(bottom_corr)
    bottom_corr = abs(np.sqrt(bottom_corr ** 2 + bottom_corr_hilbert ** 2))

    # LP-Filter
    ellip_filter = signal.iirfilter(15, Wn=6000, rp=1, rs=60, btype='lowpass',
                                        analog=False, ftype='ellip', fs=48000, output='sos')
    bottom_corr = signal.sosfilt(ellip_filter, bottom_corr)

    bottom_corr = abs(bottom_corr)
    bottom_corr_peaks, _ = signal.find_peaks(bottom_corr, distance=0.8*interval, height=0.05*np.max(bottom_corr))

    _, top_corr = Xcorr(TopMicChannel, chirp_signal, startFreq, stopFreq, sampleRate)
    top_corr_hilbert = signal.hilbert(top_corr)
    top_corr = abs(np.sqrt(top_corr ** 2 + top_corr_hilbert ** 2))

    # LP-Filter
    ellip_filter = signal.iirfilter(15, Wn=6000, rp=1, rs=60, btype='lowpass',
                                        analog=False, ftype='ellip', fs=48000, output='sos')
    top_corr = signal.sosfilt(ellip_filter, top_corr)

    top_corr = abs(top_corr)
    top_corr_peaks, _ = signal.find_peaks(top_corr, distance=0.8*interval, height=0.05 * np.max(top_corr))

    # 3.Corr:  Channel 0
    if speaker_type == "BottomSpeaker":
        ref_peaks = bottom_corr_peaks
        other_peaks = top_corr_peaks

    elif speaker_type == "TopSpeaker":
        ref_peaks = top_corr_peaks
        other_peaks = bottom_corr_peaks
    else:
        print("Speaker type is not support......")
        return -1

    #
    TopMic_corr_filter = top_corr
    BottomMic_corrfilter = bottom_corr
    #
    return speaker_type, ref_peaks, other_peaks, TopMic_corr, BottomMic_corr, TopMic_corr_filter, BottomMic_corrfilter
