# -*- encoding: utf-8 -*-

''' main.py



'''


import os, sys
import math
import numpy as np
from math import pi
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
from Signal_Processing.Find_Reference_Peaks import find_reference_peaks

import pdb

DEBUG = 1


"""
    adb shell ls /storage/emulated/0/Music/
    
    adb pull /storage/emulated/0/Music/ 


"""


if __name__ =="__main__":
    # 1. Generate Chirp signal
    startFreq = 17000
    stopFreq = 23000
    sampleRate = 48000
    tSig = 10 / 1000     # duration
    tBlank = 0 * tSig
    interval_factor = 0.9
    interval = int((tSig + tBlank) * sampleRate * interval_factor)
    # 2. Read received ultrasound signal
    basepath = "./data/Music"
    filename = 'ChirpSignal_14_48_17LR.wav'
    filepath = os.path.join(basepath, filename)

    #
    chirp_signal = UpChirp(startFreq, stopFreq, sampleRate, tSig, hann=False)
    if 1:
        # Vis chirp signal
        plt.figure(figsize=(8, 5))
        plt.plot(chirp_signal)
        plt.title('Chirp Signal')
        plt.xlabel('Time (s)')
        plt.ylabel('Amplitude')
        plt.grid(False)
        plt.pause(0.05)

    received_signal = AudioClass.load(filename=filepath)
    if 1:
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
    speakerType, ref_Peaks, TopMic_corr, BottomMic_corr = find_reference_peaks(received_signal, chirp_signal,
                                                                               startFreq, stopFreq,
                                                                               sampleRate, interval)

    if 1:
        _, ax = plt.subplots(2, 1, figsize=(10, 5), sharex=True)
        axes = ax.flatten()
        plt.suptitle('chirp duration:' + str(tSig) + '   Speaker Type: ' + speakerType)
        axes[0].plot(TopMic_corr)
        axes[0].set_title('received Signal Channel TopMic_corr')
        axes[0].set_xlabel('Range')
        axes[0].set_ylabel('Amplitude')

        axes[1].plot(BottomMic_corr)
        axes[1].set_title('received Signal Channel BottomMic_corr')
        axes[1].set_xlabel('Range')
        axes[1].set_ylabel('Amplitude')

        plt.tight_layout()
        
        plt.pause(0.05)

    # processing top mic
    TopMic_upCorrEnvYList, xPeaks = Aligned_Chirps(ref_Peaks, TopMic_corr, sampleRate=sampleRate, 
                                                   chirpSamples=interval, interval_factor=interval_factor)
    TopMic_upCorrEnvYList = np.array(TopMic_upCorrEnvYList, dtype=np.float64)

    # processing bottom mic
    BottomMic_upCorrEnvYList, _ = Aligned_Chirps(ref_Peaks, BottomMic_corr, xPeaks=xPeaks, sampleRate=sampleRate, 
                                                 chirpSamples=interval, interval_factor=interval_factor)
    BottomMic_upCorrEnvYList = np.array(BottomMic_upCorrEnvYList, dtype=np.float64)

    
    if 1:
        _, ax = plt.subplots(2, 2, figsize=(10, 6), sharex=True)
        axes = ax.flatten()
        plt.suptitle(speakerType)
        axes[0].plot(TopMic_upCorrEnvYList.T)
        axes[0].set_title('received Signal Channel TopMic')
        axes[0].set_xlabel('Range')
        axes[0].set_ylabel('Amplitude')

        axes[1].plot(BottomMic_upCorrEnvYList.T)
        axes[1].set_title('received Signal Channel BottomMic')
        axes[1].set_xlabel('Range')
        axes[1].set_ylabel('Amplitude')

        # tmp0 = upCorrEnvYList0.reshape(1, -1)
        # f, t, Zxx = signal.stft(tmp0, fs=sampleRate, nperseg=256, noverlap=128, nfft=256)
        #axes[2].pcolor(np.log2(np.abs(Zxx[0])+1e-8), cmap='jet')
        axes[2].pcolormesh(TopMic_upCorrEnvYList, shading='gouraud', cmap='jet', antialiased=True)
        axes[2].set_title('received Signal Channel TopMic All Chirps STFT')
        axes[2].set_xlabel('Range')
        axes[2].set_ylabel('Chirps')

        # tmp1 = upCorrEnvYList1.reshape(1, -1)
        # f, t, Zxx = signal.stft(tmp1, fs=sampleRate, nperseg=256, noverlap=128, nfft=256)
        axes[3].pcolormesh(BottomMic_upCorrEnvYList, shading='gouraud', cmap='jet', antialiased=True)
        axes[3].set_title('received Signal Channel Bottom all Chirps STFT')
        axes[3].set_xlabel('Range')
        axes[3].set_ylabel('Chirps')

        plt.tight_layout()
        plt.pause(0.05)

    # obtain bin signals
    if speakerType == "TopSpeaker":
        bin_num = 34
        phaseSig = BottomMic_upCorrEnvYList[:, bin_num]
        _, ax = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
        axes = ax.flatten()
        axes[0].plot(phaseSig)

        phaseSigfft = np.fft.fftshift(np.fft.fft(phaseSig, n=256))
        axes[1].plot(np.abs(phaseSigfft))

        plt.tight_layout()
        plt.pause(0.05)

    # obtain bin signals
    if speakerType == "BottomSpeaker":
        bin_num = 50
        phaseSig = TopMic_upCorrEnvYList[:, bin_num]
        _, ax = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
        axes = ax.flatten()
        axes[0].plot(phaseSig)

        phaseSigfft = np.fft.fftshift(np.fft.fft(phaseSig, n=256))
        axes[1].plot(np.abs(phaseSigfft))

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









