# -*- encoding: utf-8 -*-

''' processing_chirp_xcorr.py

Description:
    This script processing received ultrasound signal based on mixer.

get datasets:
    * Mi_11: adb connect 192.168.2.186:5555
    * Mi_10: adb connect 192.168.2.154:5555

    adb shell ls /storage/emulated/0/Music/
    adb pull /storage/emulated/0/Music/
'''

import os, sys
import math
import numpy as np
from math import pi
import scipy
import scipy.signal as signal
from scipy.io import wavfile
from pathlib import Path
import matplotlib, copy
import matplotlib as mpl
from matplotlib.axes import Axes
import matplotlib.pyplot as plt
from collections import defaultdict

# import Fmcw, Alsa, Plot, BaseUtils
# from  Audio import AudioClass
#

from Signal_Processing.Audio import AudioClass
from Signal_Processing.Chirp import UpChirp
from Signal_Processing.Common import Xcorr, Aligned_Chirps
from Signal_Processing.Find_Reference_Peaks import find_reference_peaks_stereo

import pdb

DEBUG = 1

if __name__ == "__main__":
    # 1. Generate Chirp signal
    startFreq = 17000
    stopFreq = 20000
    halfFreq = (stopFreq + startFreq) // 2
    sampleRate = 48000
    tSig = 5 / 1000  # duration
    tBlank = 0 * tSig
    interval_factor = 0.95
    interval = int((tSig + tBlank) * sampleRate * interval_factor)

    # 2. Read received ultrasound signal
    basepath = "./data/Music"
    filename = 'ChirpSignal_9_28_26LR.wav'
    filepath = os.path.join(basepath, filename)

    # generate signal chirp signal
    chirp_signal_top = UpChirp(startFreq, halfFreq, sampleRate, tSig, hann=False)
    chirp_signal_bottom = UpChirp(halfFreq, stopFreq, sampleRate, tSig, hann=False)

    if 0:
        # Vis chirp signal
        plt.figure(figsize=(8, 5))
        plt.plot(chirp_signal_top)
        plt.title('Chirp Signal')
        plt.xlabel('Time (s)')
        plt.ylabel('Amplitude')
        plt.grid(False)
        plt.pause(0.05)

    # load received signal
    received_signal = AudioClass.load(filename=filepath)
    if 0:
        # Vis chirp signal
        _, ax = plt.subplots(2, 1, figsize=(8, 5), sharex=False)
        axes = ax.flatten()
        axes[0].plot(received_signal[0])
        axes[0].set_title('received Signal Channel 0')
        axes[0].set_xlabel('Time (s)')
        axes[0].set_ylabel('Amplitude')

        axes[1].plot(received_signal[1])
        axes[1].set_title('received Signal Channel 1')
        axes[1].set_xlabel('Time (s)')
        axes[1].set_ylabel('Amplitude')

        plt.tight_layout()
        plt.pause(0.05)

    # find reference peaks
    speakerType, ref_peaks, other_peaks, \
    TopMic_corr, BottomMic_corr, \
    TopMic_corr_filter, BottomMic_corr_filter = find_reference_peaks_stereo(received_signal,
                                                                     chirp_signal_top, chirp_signal_bottom,
                                                                     startFreq, stopFreq, sampleRate, interval)

    if 1:
        _, ax = plt.subplots(2, 2, figsize=(20, 10), sharex=True)
        axes = ax.flatten()
        plt.suptitle('chirp duration:' + str(tSig) + '   Speaker Type: ' + speakerType)
        axes[0].plot(TopMic_corr)
        axes[0].set_title('received Signal Channel TopMic_corr')
        axes[0].set_xlabel('Range')
        axes[0].set_ylabel('Amplitude')
        # axes[0].scatter(other_peaks, TopMic_corr[other_peaks], c='r')

        axes[2].plot(BottomMic_corr)
        axes[2].set_title('received Signal Channel BottomMic_corr')
        axes[2].set_xlabel('Range')
        axes[2].set_ylabel('Amplitude')
        # axes[2].scatter(ref_peaks, BottomMic_corr[ref_peaks], c='r')

        axes[1].plot(TopMic_corr_filter)
        axes[1].set_title('received Signal Channel TopMic_corr_filter')
        axes[1].set_xlabel('Range')
        axes[1].set_ylabel('Amplitude')
        if speakerType == "TopSpeaker":
            axes[1].scatter(ref_peaks, TopMic_corr_filter[ref_peaks], c='r')
        elif speakerType == "BottomSpeaker":
            axes[1].scatter(other_peaks, TopMic_corr_filter[other_peaks], c='r')
        else:
            pass

        axes[3].plot(BottomMic_corr_filter)
        axes[3].set_title('received Signal Channel BottomMic_corr_filter')
        axes[3].set_xlabel('Range')
        axes[3].set_ylabel('Amplitude')
        if speakerType == "TopSpeaker":
            axes[3].scatter(other_peaks, BottomMic_corr_filter[other_peaks], c='r')
        elif speakerType == "BottomSpeaker":
            axes[3].scatter(ref_peaks, BottomMic_corr_filter[ref_peaks], c='r')
        else:
            pass
        plt.tight_layout()

        plt.pause(0.05)


    # define var
    TopMic_upCorrEnvYList = []
    BottomMic_upCorrEnvYList = []

    if speakerType == "TopSpeaker":
        # processing aligned top mic signals
        TopMic_upCorrEnvYList, xPeaks = Aligned_Chirps(ref_peaks, TopMic_corr_filter, sampleRate=sampleRate,
                                                       chirpSamples=interval, interval_factor=interval_factor)
        TopMic_upCorrEnvYList = np.array(TopMic_upCorrEnvYList, dtype=np.float64)

        # processing aligned bottom mic signals
        BottomMic_upCorrEnvYList, _ = Aligned_Chirps(ref_peaks, BottomMic_corr_filter,  xPeaks=xPeaks, sampleRate=sampleRate,
                                                     chirpSamples=interval, interval_factor=interval_factor)
        BottomMic_upCorrEnvYList = np.array(BottomMic_upCorrEnvYList, dtype=np.float64)

    elif speakerType == "BottomSpeaker":
        # processing aligned bottom mic signals
        BottomMic_upCorrEnvYList, xPeaks = Aligned_Chirps(ref_peaks, BottomMic_corr_filter, sampleRate=sampleRate,
                                                     chirpSamples=interval, interval_factor=interval_factor)
        BottomMic_upCorrEnvYList = np.array(BottomMic_upCorrEnvYList, dtype=np.float64)

        # processing aligned top mic signals
        TopMic_upCorrEnvYList, _ = Aligned_Chirps(ref_peaks, TopMic_corr_filter, xPeaks=xPeaks, sampleRate=sampleRate,
                                                       chirpSamples=interval, interval_factor=interval_factor)
        TopMic_upCorrEnvYList = np.array(TopMic_upCorrEnvYList, dtype=np.float64)
    else:
        print("SpeakerTop is Error.....\n")
        exit(-1)

    # obtain
    if speakerType == "TopSpeaker":
        corr_Envy = abs(np.diff(TopMic_corr_filter[ref_peaks]))
    elif speakerType == "BottomSpeaker":
        corr_Envy = abs(np.diff(BottomMic_corr_filter[ref_peaks]))
    else:
        print("SpeakerTop is Error.....\n")
        exit(-1)

    peaks, _ = signal.find_peaks(corr_Envy, height=5*np.mean(corr_Envy))
    if len(peaks) > 1:
        start_blow = peaks[0]
        stop_blow = peaks[-1]
    else:
        start_blow = 0
        stop_blow = -1
    if 1:
        _, ax = plt.subplots(2, 1, figsize=(8, 8))
        axes = ax.flatten()
        plt.suptitle(speakerType)


        if speakerType == "TopSpeaker":
            top_corr_norm = TopMic_corr_filter[ref_peaks] / np.max(TopMic_corr_filter[ref_peaks])
            bottom_corr_norm = BottomMic_corr_filter[other_peaks] / np.max(BottomMic_corr_filter[other_peaks])
            time_top = np.arange(len(top_corr_norm)) * 5.0 / 1000.0
            time_bottom = np.arange(len(bottom_corr_norm)) * 5.0 / 1000.0
        elif speakerType == "BottomSpeaker":
            top_corr_norm = TopMic_corr_filter[other_peaks] / np.max(TopMic_corr_filter[other_peaks])
            bottom_corr_norm = BottomMic_corr_filter[ref_peaks] / np.max(BottomMic_corr_filter[ref_peaks])
            time_top = np.arange(len(top_corr_norm)) * 5.0 / 1000.0
            time_bottom = np.arange(len(bottom_corr_norm)) * 5.0 / 1000.0

        axes[0].plot(time_top, top_corr_norm, label='Top morm')
        axes[0].plot(time_bottom, bottom_corr_norm, label='Bottom norm')
        axes[0].set_xlabel("Time (s)")
        axes[0].legend()

        axes[1].plot(corr_Envy, label='ref_corr diff')
        axes[1].scatter(peaks, corr_Envy[peaks], c='r', label='Peaks')
        axes[1].legend()

        plt.pause(0.05)

    if 1:
        _, ax = plt.subplots(2, 3, figsize=(15, 6), sharex=True)
        axes = ax.flatten()
        plt.suptitle(speakerType)
        axes[0].plot(TopMic_upCorrEnvYList.T)
        axes[0].set_title('received Signal Channel TopMic')
        axes[0].set_xlabel('Range')
        axes[0].set_ylabel('Amplitude')

        axes[3].plot(BottomMic_upCorrEnvYList.T)
        axes[3].set_title('received Signal Channel BottomMic')
        axes[3].set_xlabel('Range')
        axes[3].set_ylabel('Amplitude')

        # tmp0 = upCorrEnvYList0.reshape(1, -1)
        # f, t, Zxx = signal.stft(tmp0, fs=sampleRate, nperseg=256, noverlap=128, nfft=256)
        # axes[2].pcolor(np.log2(np.abs(Zxx[0])+1e-8), cmap='jet')
        axes[1].pcolormesh(TopMic_upCorrEnvYList, shading='gouraud', cmap='jet', antialiased=True)
        axes[1].set_title('received Signal Channel TopMic All Chirps STFT')
        axes[1].set_xlabel('Range')
        axes[1].set_ylabel('Chirps')

        # tmp1 = upCorrEnvYList1.reshape(1, -1)
        # f, t, Zxx = signal.stft(tmp1, fs=sampleRate, nperseg=256, noverlap=128, nfft=256)
        axes[4].pcolormesh(BottomMic_upCorrEnvYList, shading='gouraud', cmap='jet', antialiased=True)
        axes[4].set_title('received Signal Channel Bottom all Chirps STFT')
        axes[4].set_xlabel('Range')
        axes[4].set_ylabel('Chirps')

        # Zoom [start_blow, stop_blow]
        TopMic_upCorrEnvYList_zoom = TopMic_upCorrEnvYList[start_blow:stop_blow, :]
        BottomMic_upCorrEnvYList_zoom = BottomMic_upCorrEnvYList[start_blow:stop_blow, :]

        axes[2].pcolormesh(TopMic_upCorrEnvYList_zoom, shading='gouraud', cmap='jet', antialiased=True)
        axes[2].set_title('received Signal Channel TopMic All Chirps STFT')
        axes[2].set_xlabel('Range')
        axes[2].set_ylabel('Chirps')

        axes[5].pcolormesh(BottomMic_upCorrEnvYList_zoom, shading='gouraud', cmap='jet', antialiased=True)
        axes[5].set_title('received Signal Channel Bottom all Chirps STFT')
        axes[5].set_xlabel('Range')
        axes[5].set_ylabel('Chirps')

        plt.tight_layout()
        plt.pause(0.05)

    # obtain bin signals
    if speakerType == "TopSpeaker":
        target_peaks = []
        for i in range(int(0.1*BottomMic_upCorrEnvYList.shape[0])):
            _peaks, _ = signal.find_peaks(BottomMic_upCorrEnvYList[i,:], distance=BottomMic_upCorrEnvYList.shape[1])
            target_peaks.append(_peaks)
        bin_num = scipy.stats.mode(target_peaks)[0][0]

        # located the target bin
        phaseSig = np.reshape(BottomMic_upCorrEnvYList[start_blow:stop_blow, bin_num], (1, -1))
        phaseSigfft = np.fft.fftshift(np.fft.fft(phaseSig, n=256))
        if 1:
            _, ax = plt.subplots(2, 1, figsize=(10, 6), sharex=False)
            axes = ax.flatten()
            axes[0].plot(phaseSig.T, label='phase signal')
            axes[0].legend()

            axes[1].plot(np.abs(phaseSigfft.T), label='FFT')
            axes[0].legend()

            plt.tight_layout()
            plt.pause(0.05)

    # obtain bin signals
    if speakerType == "BottomSpeaker":
        target_peaks = []
        for i in range(int(0.1*TopMic_upCorrEnvYList.shape[0])):
            _peaks, _ = signal.find_peaks(TopMic_upCorrEnvYList[i, :], distance=TopMic_upCorrEnvYList.shape[1])
            target_peaks.append(_peaks)
        bin_num = scipy.stats.mode(target_peaks)[0][0]

        phaseSig = np.reshape(TopMic_upCorrEnvYList[start_blow : stop_blow, bin_num], (1, -1)) #- TopMic_upCorrEnvYList[:, 66]
        phaseSigfft = np.fft.fftshift(np.fft.fft(phaseSig, n=128))

        if 1:
            _, ax = plt.subplots(2, 1, figsize=(10, 6), sharex=False)
            axes = ax.flatten()
            axes[0].plot(phaseSig.T, label='phase signal')
            axes[0].legend()

            axes[1].plot(np.abs(phaseSigfft.T), label='FFT')
            axes[0].legend()

            plt.tight_layout()
            plt.pause(0.05)

    # diff,
    if 0:
        TopMic_diff = np.diff(TopMic_upCorrEnvYList, axis=0)
        BottomMic_diff = np.diff(BottomMic_upCorrEnvYList, axis=0)
        fig, ax = plt.subplots(1, 2, figsize=(8, 5), sharex=False)
        fig.suptitle("diff background")
        axes = ax.flatten()
        axes[0].pcolormesh(TopMic_diff, shading='gouraud', cmap='jet', antialiased=True)
        axes[1].pcolormesh(BottomMic_diff, shading='gouraud', cmap='jet', antialiased=True)
        plt.tight_layout()
        plt.pause(0.05)

    # subtract background frame
    if 1:
        # obtaining front 50 frames to estimate background frame.
        TopMic_ref_frame = np.mean(TopMic_upCorrEnvYList[0:50, :], axis=0)
        BottomMic_ref_frame = np.mean(BottomMic_upCorrEnvYList[0:50, :], axis=0)
        TopMic_diff = TopMic_upCorrEnvYList - TopMic_ref_frame
        BottomMic_diff = BottomMic_upCorrEnvYList - BottomMic_ref_frame

        fig, ax = plt.subplots(1, 2, figsize=(8, 5), sharex=False)
        fig.suptitle("subtract background")
        axes = ax.flatten()
        axes[0].pcolormesh(TopMic_diff, shading='gouraud', cmap='jet', antialiased=True)
        axes[0].set_title("Top mic")
        axes[1].pcolormesh(BottomMic_diff, shading='gouraud', cmap='jet', antialiased=True)
        axes[1].set_title("Bottom mic")

        plt.grid()
        plt.tight_layout()
        plt.pause(0.05)

    if speakerType == "TopSpeaker":
        # Processing bottom mic signal
        print("Speaker type is :{}".format(speakerType))
        # frame diff
        BottomMic_diff = np.diff(BottomMic_upCorrEnvYList, axis=0)
        plt.figure()
        plt.pcolormesh(BottomMic_diff, shading='gouraud', cmap='jet', antialiased=True)

    elif speakerType == "BottomSpeaker":
        # Processing top mic signal
        print("Speaker type is :{}".format(speakerType))

        # find target peak

    pass









