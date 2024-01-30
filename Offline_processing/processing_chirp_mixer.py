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
from Signal_Processing.Chirp import UpChirp, UpChirp_cos, UpChirp_sin
from Signal_Processing.Common import Xcorr, Aligned_Chirps, one_dimensional_kalman_filter
from Signal_Processing.Find_Reference_Peaks import find_reference_peaks

import pdb

DEBUG = 1
DEBUG_VIS = 1

if __name__ == "__main__":
    # 1. Generate Chirp signal
    startFreq = 17000
    stopFreq = 23000
    sampleRate = 48000
    tSig = 10 / 1000  # duration
    tBlank = 0 * tSig
    interval_factor = 0.95
    interval = int((tSig + tBlank) * sampleRate * interval_factor)
    skip_time = 0.25  # 250ms
    chirp_length = int(sampleRate * tSig)

    # 3D模型相关参数
    Model_L4 = 0.35  # m
    Model_L2 = 0.15  # m
    Model_r = 0.2    # m
    Voice_c = 346  # m/s

    # 2. Read received ultrasound signal
    basepath = "./data/Music"
    filename = 'ChirpSignal_10_23_6LR.wav'
    filepath = os.path.join(basepath, filename)

    # generate signal chirp signal
    chirp_signal = UpChirp(startFreq, stopFreq, sampleRate, tSig, hann=False)
    if 0:
        # Vis chirp signal
        plt.figure(figsize=(8, 5))
        plt.plot(chirp_signal)
        plt.title('Chirp Signal')
        plt.xlabel('Time (s)')
        plt.ylabel('Amplitude')
        plt.grid(False)
        plt.pause(0.05)

    # load received signal
    receive = AudioClass.load(filename=filepath)

    # 前后各跳过数据, 1. 计算数据是多少样本点
    skip_samples = int((skip_time / tSig) * sampleRate * tSig)
    received_signals = np.zeros((2, receive.length - skip_samples - 1))
    received_signals[0] = receive[0][skip_samples:-1]
    received_signals[1] = receive[1][skip_samples:-1]

    # get top/bottom mic signal and Filter
    TopMicChannel = received_signals[0]
    BottomMicChannel = received_signals[1]
    if 0:
        ellip_filter = signal.iirfilter(25, Wn=10000, rp=1, rs=60, btype='highpass',
                                        analog=False, ftype='cheby2', fs=48000, output='sos')
        TopMicChannel = signal.sosfilt(ellip_filter, TopMicChannel)
        BottomMicChannel = signal.sosfilt(ellip_filter, BottomMicChannel)

    # 计算mic采集数据的时频谱图
    f_top, t_top, Zxx_top = signal.stft(TopMicChannel, fs=48000, nperseg=128,
                                        noverlap=64, nfft=128, return_onesided=True)
    f_bottom, t_bottom, Zxx_bottom = signal.stft(BottomMicChannel, fs=48000, nperseg=128,
                                                 noverlap=64, nfft=128, return_onesided=True)
    Zxx_top_norm = np.log10(np.abs(Zxx_top)) / np.max(np.log10(np.abs(Zxx_top)))
    Zxx_bottom_norm = np.log10(np.abs(Zxx_bottom)) / np.max(np.log10(np.abs(Zxx_bottom)))

    if DEBUG_VIS:
        # Vis chirp signal
        _, ax = plt.subplots(4, 2, figsize=(8, 5), sharex=False)
        axes = ax.flatten()
        axes[0].plot(TopMicChannel)
        axes[0].set_title('received Signal Channel 0')
        axes[0].set_xlabel('Time (s)')
        # axes[0].set_ylabel('Amplitude')

        axes[1].plot(BottomMicChannel)
        axes[1].set_title('received Signal Channel 1')
        axes[1].set_xlabel('Time (s)')
        # axes[1].set_ylabel('Amplitude')

        axes[2].pcolor(np.log10(np.abs(Zxx_top)), cmap='jet')
        axes[2].set_title('Filtered top mic signal spectrum')

        axes[3].pcolor(np.log10(np.abs(Zxx_bottom)), cmap='jet')
        axes[3].set_title('Filtered bottom mic signal spectrum')


        """
            
        """
        if 0:
            Zxx_bottom_logabs = np.abs(Zxx_bottom) # np.log10(np.abs(Zxx_bottom))
            # Zxx_bottom_logabs = Zxx_bottom_logabs / np.max(Zxx_bottom_logabs)
            Zxx_bottom_logabs_norm = (Zxx_bottom_logabs - np.min(Zxx_bottom_logabs)) / (np.max(Zxx_bottom_logabs) - np.min(Zxx_bottom_logabs))

            Zxx_template_length = 90
            Zxx_template = Zxx_bottom_logabs_norm[:, :int(Zxx_template_length)]
            # Zxx_template = Zxx_template / np.max(Zxx_template)

            cycles = int(np.floor(Zxx_bottom_logabs_norm.shape[1] / Zxx_template.shape[1]))
            ref_noise = np.copy(Zxx_bottom_logabs_norm) #np.zeros(Zxx_bottom_logabs_norm.shape)

            for i in range(cycles):
                ref_noise[:, i*int(Zxx_template_length):(i+1)*int(Zxx_template_length)] \
                    = Zxx_bottom_logabs_norm[:, i*int(Zxx_template_length):(i+1)*int(Zxx_template_length)] - Zxx_template

            axes[4].pcolor(np.log10(ref_noise+1e-8), cmap='jet')
            axes[4].set_title('Airflow noise spectrum of bottom mic')

            Zxx_top_logabs = np.abs(Zxx_top)  # np.log10(np.abs(Zxx_top))
            Zxx_top_logabs_norm = (Zxx_top_logabs - np.min(Zxx_top_logabs)) / (np.max(Zxx_top_logabs) - np.min(Zxx_top_logabs))

            axes[5].pcolor(np.log10(Zxx_top_logabs_norm), cmap='jet')
            axes[5].set_title('Top mic spectrum')

            # ref_noise_norm = (Zxx_bottom_noise - np.min(Zxx_bottom_noise)) / (np.max(Zxx_bottom_noise) - np.min(Zxx_bottom_noise))
            Zxx_top_logabs_denoise1 = np.log10(Zxx_top_logabs_norm) - np.log10(ref_noise)
            axes[6].pcolor(Zxx_top_logabs_denoise1, cmap='jet')
            axes[6].set_title('Top mic spectrum Zxx_top_logabs_denoise1')

            Zxx_top_logabs_denoise2 = Zxx_top_logabs_norm - ref_noise
            axes[7].pcolor(np.log10(Zxx_top_logabs_denoise2+1e-8), cmap='jet')

        if 1:
            Zxx_bottom_logabs = np.log10(np.abs(Zxx_bottom))
            # Zxx_bottom_logabs = Zxx_bottom_logabs / np.max(Zxx_bottom_logabs)
            Zxx_bottom_logabs_norm = (Zxx_bottom_logabs - np.min(Zxx_bottom_logabs)) / (
                        np.max(Zxx_bottom_logabs) - np.min(Zxx_bottom_logabs))

            Zxx_template_length = 90
            Zxx_template = Zxx_bottom_logabs_norm[:, :int(Zxx_template_length)]
            # Zxx_template = Zxx_template / np.max(Zxx_template)

            cycles = int(np.floor(Zxx_bottom_logabs_norm.shape[1] / Zxx_template.shape[1]))
            ref_noise = np.copy(Zxx_bottom_logabs_norm)  # np.zeros(Zxx_bottom_logabs_norm.shape)

            for i in range(cycles):
                ref_noise[:, i * int(Zxx_template_length):(i + 1) * int(Zxx_template_length)] \
                    = Zxx_bottom_logabs_norm[:,
                      i * int(Zxx_template_length):(i + 1) * int(Zxx_template_length)] - Zxx_template

            axes[4].pcolor(ref_noise, cmap='jet')
            axes[4].set_title('Airflow noise spectrum of bottom mic')

            Zxx_top_logabs = np.log10(np.abs(Zxx_top))
            Zxx_top_logabs_norm = (Zxx_top_logabs - np.min(Zxx_top_logabs)) / (
                        np.max(Zxx_top_logabs) - np.min(Zxx_top_logabs))

            axes[5].pcolor(Zxx_top_logabs_norm, cmap='jet')
            axes[5].set_title('Top mic spectrum')

            # ref_noise_norm = (Zxx_bottom_noise - np.min(Zxx_bottom_noise)) / (np.max(Zxx_bottom_noise) - np.min(Zxx_bottom_noise))
            Zxx_top_logabs_denoise1 = np.log10(Zxx_top_logabs_norm) - np.log10(ref_noise)
            axes[6].pcolor(Zxx_top_logabs_denoise1, cmap='jet')
            axes[6].set_title('Top mic spectrum Zxx_top_logabs_denoise1')

            Zxx_top_logabs_denoise2 = Zxx_top_logabs_norm - ref_noise
            axes[7].pcolor(np.log10(Zxx_top_logabs_denoise2 + 1e-8), cmap='jet')
        # 顶发底收 消除气流声音: 谱减法，底部mic的时频图减去顶部mic的时频图
        plt.tight_layout()
        plt.pause(0.05)

    # find reference peaks, 用于对齐信号强度
    speakerType, ref_peaks, other_peaks, \
    TopMic_corr, BottomMic_corr, \
    TopMic_corr_filter, BottomMic_corr_filter = find_reference_peaks(received_signals, chirp_signal,
                                                       startFreq, stopFreq,
                                                       sampleRate, interval)

    if DEBUG_VIS:
        _, ax = plt.subplots(2, 2, figsize=(8, 5),  sharex=True)
        axes = ax.flatten()
        plt.suptitle('chirp duration:' + str(tSig) + '   Speaker Type: ' + speakerType)
        axes[0].plot(TopMic_corr, c='cyan', label='TopMic_corr')
        axes[0].set_title('received Signal Channel TopMic_corr')
        axes[0].set_xlabel('Range Samples')
        axes[0].set_ylabel('Amplitude')
        # axes[0].scatter(other_peaks, TopMic_corr[other_peaks], c='r')
        axes[0].plot(TopMic_corr_filter, c='blue', label='TopMic_corr_filter')
        axes[0].set_title('received Signal Channel TopMic_corr_filter')
        axes[0].set_xlabel('Range Samples')
        axes[0].set_ylabel('Amplitude')
        if speakerType == "TopSpeaker":
            axes[0].scatter(ref_peaks, TopMic_corr_filter[ref_peaks], c='r')
        elif speakerType == "BottomSpeaker":
            axes[0].scatter(other_peaks, TopMic_corr_filter[other_peaks], c='r')
        axes[0].legend()

        axes[2].plot(BottomMic_corr, c='deepskyblue', label='BottomMic_corr')
        axes[2].set_title('received Signal Channel BottomMic_corr')
        axes[2].set_xlabel('Range')
        axes[2].set_ylabel('Amplitude')
        # axes[2].scatter(ref_peaks, BottomMic_corr[ref_peaks], c='r')
        axes[2].plot(BottomMic_corr_filter, c='coral', label='BottomMic_corr_filter')
        axes[2].set_title('received Signal Channel BottomMic_corr_filter')
        axes[2].set_xlabel('Range Samples')
        axes[2].set_ylabel('Amplitude')
        if speakerType == "TopSpeaker":
            axes[2].scatter(other_peaks, BottomMic_corr_filter[other_peaks], c='g')
        elif speakerType == "BottomSpeaker":
            axes[2].scatter(ref_peaks, BottomMic_corr_filter[ref_peaks], c='g')
        axes[2].legend()

        TopMic_corr_filter_norm = TopMic_corr_filter / np.max(TopMic_corr_filter)
        axes[1].plot(TopMic_corr_filter_norm, c='blue', label='Mic of top channel')
        axes[1].set_title('received Signal Channel TopMic_corr_filter Norm.')
        axes[1].set_xlabel('Range Samples')
        axes[1].set_ylabel('Amplitude')
        if speakerType == "TopSpeaker":
            axes[1].scatter(ref_peaks, TopMic_corr_filter_norm[ref_peaks], c='g')
        elif speakerType == "BottomSpeaker":
            axes[1].scatter(other_peaks, TopMic_corr_filter_norm[other_peaks], c='g')

        BottomMic_corr_filter_norm = BottomMic_corr_filter / np.max(BottomMic_corr_filter)
        axes[1].plot(BottomMic_corr_filter_norm, c='red', label='Mic of bottom channel')
        axes[1].set_title('received Signal Channel BottomMic_corr_filter')
        axes[1].set_xlabel('Range Samples')
        axes[1].set_ylabel('Amplitude')
        if speakerType == "TopSpeaker":
            axes[1].scatter(other_peaks, BottomMic_corr_filter_norm[other_peaks], c='g')
        elif speakerType == "BottomSpeaker":
            axes[1].scatter(ref_peaks, BottomMic_corr_filter_norm[ref_peaks], c='g')
        axes[1].legend()

        # axes[3].plot(BottomMic_corr_filter, c='coral')
        # axes[3].set_title('received Signal Channel BottomMic_corr_filter')
        # axes[3].set_xlabel('Range')
        # axes[3].set_ylabel('Amplitude')
        # if speakerType == "TopSpeaker":
        #     axes[3].scatter(other_peaks, BottomMic_corr_filter[other_peaks], c='g')
        # elif speakerType == "BottomSpeaker":
        #     axes[3].scatter(ref_peaks, BottomMic_corr_filter[ref_peaks], c='g')
        axes[3].plot(TopMic_corr, c='cyan', label='TopMic_corr')
        axes[3].plot(BottomMic_corr, c='deepskyblue', label='BottomMic_corr')
        axes[3].set_title('Received Corr-Signal of TopMic and BottomMic')
        axes[3].set_xlabel('Range Samples')
        axes[3].set_ylabel('Amplitude')
        axes[3].legend()

        plt.tight_layout()
        plt.pause(0.05)

    """ 
        检测吹气开始和结束的时刻，并且记录数据
    """
    # define var
    TopMic_upCorrEnvYList = []
    BottomMic_upCorrEnvYList = []
    if speakerType == "TopSpeaker":
        # processing aligned top mic signals
        TopMic_upCorrEnvYList, xPeaks = Aligned_Chirps(ref_peaks, TopMic_corr_filter,
                                                       sampleRate=sampleRate,
                                                       chirpSamples=interval,
                                                       interval_factor=interval_factor)
        TopMic_upCorrEnvYList = np.array(TopMic_upCorrEnvYList, dtype=np.float64)

        # processing aligned bottom mic signals
        BottomMic_upCorrEnvYList, _ = Aligned_Chirps(ref_peaks, BottomMic_corr_filter,
                                                     xPeaks=xPeaks,
                                                     sampleRate=sampleRate,
                                                     chirpSamples=interval,
                                                     interval_factor=interval_factor)
        BottomMic_upCorrEnvYList = np.array(BottomMic_upCorrEnvYList, dtype=np.float64)

        corr_Envy = abs(np.diff(TopMic_corr_filter[ref_peaks]))

        peaks_, _ = signal.find_peaks(corr_Envy, height=5 * np.mean(corr_Envy))
        if len(peaks_) > 1:
            start_airflow = peaks_[0]
            stop_airflow = peaks_[-1]
        else:
            start_airflow = 0
            stop_airflow = -1

    elif speakerType == "BottomSpeaker":
        # processing aligned bottom mic signals
        BottomMic_upCorrEnvYList, xPeaks = Aligned_Chirps(ref_peaks, BottomMic_corr_filter,
                                                          sampleRate=sampleRate,
                                                          chirpSamples=interval,
                                                          interval_factor=interval_factor)
        BottomMic_upCorrEnvYList = np.array(BottomMic_upCorrEnvYList, dtype=np.float64)

        # processing aligned top mic signals
        TopMic_upCorrEnvYList, _ = Aligned_Chirps(ref_peaks, TopMic_corr_filter,
                                                  xPeaks=xPeaks,
                                                  sampleRate=sampleRate,
                                                  chirpSamples=interval,
                                                  interval_factor=interval_factor)
        TopMic_upCorrEnvYList = np.array(TopMic_upCorrEnvYList, dtype=np.float64)

        corr_Envy = abs(np.diff(BottomMic_corr_filter[ref_peaks]))

        peaks_, _ = signal.find_peaks(corr_Envy, height=5 * np.mean(corr_Envy))
        if len(peaks_) > 1:
            start_airflow = peaks_[0]
            stop_airflow = peaks_[-1]
        else:
            start_airflow = 0
            stop_airflow = -1
    else:
        print("SpeakerTop is Error.....\n")
        exit(-1)

    if speakerType == "TopSpeaker":
        top_corr_norm = TopMic_corr_filter[ref_peaks] / np.max(TopMic_corr_filter[ref_peaks])
        bottom_corr_norm = BottomMic_corr_filter[other_peaks] / np.max(BottomMic_corr_filter[other_peaks])
        time_top = np.arange(len(top_corr_norm)) * tSig
        time_bottom = np.arange(len(bottom_corr_norm)) * tSig
    elif speakerType == "BottomSpeaker":
        top_corr_norm = TopMic_corr_filter[other_peaks] / np.max(TopMic_corr_filter[other_peaks])
        bottom_corr_norm = BottomMic_corr_filter[ref_peaks] / np.max(BottomMic_corr_filter[ref_peaks])
        time_top = np.arange(len(top_corr_norm)) * tSig
        time_bottom = np.arange(len(bottom_corr_norm)) * tSig

    if DEBUG_VIS:
        _, ax = plt.subplots(1, 2, figsize=(8, 5))
        axes = ax.flatten()
        plt.suptitle(speakerType)

        if speakerType == "TopSpeaker":
            axes[0].plot(time_top, top_corr_norm, label='Top Norm')
        elif speakerType == "BottomSpeaker":
            axes[0].plot(time_bottom, bottom_corr_norm, label='Bottom Norm')
        axes[0].set_xlabel("Time (s)")
        axes[0].legend()

        if speakerType == "BottomSpeaker":
            corr_bottom = abs(np.diff(bottom_corr_norm))
            peaks_, _ = signal.find_peaks(corr_bottom, height=5 * np.mean(corr_bottom))
            axes[1].plot(corr_bottom, label='Bottom corr diff', c='b')
            axes[1].scatter(peaks_, corr_bottom[peaks_], c='g', label='Peaks')
            axes[1].axvline(x=start_airflow, color='red', linestyle='--', label='start airflow')
            axes[1].axvline(x=stop_airflow, color='red', linestyle='-', label='stop airflow')
            axes[1].legend()

        elif speakerType == "TopSpeaker":
            corr_top = abs(np.diff(top_corr_norm))
            peaks_, _ = signal.find_peaks(corr_top, height=5 * np.mean(corr_top))
            axes[1].plot(corr_top, label='Top corr diff', c='b')
            axes[1].scatter(peaks_, corr_top[peaks_], c='g', label='Peaks')
            axes[1].axvline(x=start_airflow, color='red', linestyle='--', label='start airflow')
            axes[1].axvline(x=stop_airflow, color='red', linestyle='-', label='stop airflow')
            axes[1].legend()
        plt.pause(0.05)

    """ To demodulation
        find the start point and skip, to aligned the chirp data.
    """
    if speakerType == "TopSpeaker":
        TopMic_Peaks, _ = signal.find_peaks(TopMic_corr, distance=0.8 * interval,
                                               height=0.05 * np.max(TopMic_corr))
        second_skip_sample = int( (TopMic_Peaks[0] - chirp_length/2 + chirp_length) % chirp_length)

    elif speakerType == "BottomSpeaker":
        BottomMic_Peaks, _ = signal.find_peaks(BottomMic_corr, distance=0.8 * interval,
                                               height=0.05 * np.max(BottomMic_corr))
        second_skip_sample = int( (BottomMic_Peaks[0] - chirp_length/2 + chirp_length) % chirp_length)

    # Skip the top/bottom mic signal
    TopMicChannel = received_signals[0][second_skip_sample:-1]
    BottomMicChannel = received_signals[1][second_skip_sample:-1]

    if 0:
        # Vis chirp signal
        _, ax = plt.subplots(1, 2, figsize=(8, 5), sharex=True)
        axes = ax.flatten()

        plt.suptitle('Speaker Type: {}'.format(speakerType))
        axes[0].plot(TopMicChannel)
        axes[0].set_title('received Signal Channel 0')
        axes[0].set_xlabel('Time (s)')
        axes[0].set_ylabel('Amplitude')
        axes[0].legend()

        axes[1].plot(BottomMicChannel)
        axes[1].set_title('received Signal Channel 1')
        axes[1].set_xlabel('Time (s)')
        axes[1].set_ylabel('Amplitude')
        axes[1].legend()

        plt.pause(0.05)

    # Test the self-correlation
    if 0:
        _, TopMic_corr = Xcorr(TopMicChannel, TopMicChannel, startFreq, stopFreq, sampleRate)
        _, BottomMic_corr = Xcorr(BottomMicChannel, BottomMicChannel, startFreq, stopFreq, sampleRate)

        # Vis chirp signal
        _, ax = plt.subplots(1, 2, figsize=(8, 5), sharex=False)
        axes = ax.flatten()
        axes[0].plot(TopMic_corr, c='r', label='Topmic self-corr')
        axes[0].legend()
        axes[1].plot(BottomMic_corr, c='b', label='Bottommic self-corr')
        axes[1].legend()

    # 模拟解调发射端的，移位解调
    if speakerType == "TopSpeaker":
        chirp_number = int(np.floor(BottomMicChannel.shape[0] / chirp_length))  # Count the number of chirps
        TopMixSignals = np.zeros((chirp_number, chirp_length), dtype=np.complex64)
        BottomMixSignals = np.zeros((chirp_number, chirp_length), dtype=np.complex64)

    # 模拟解调发射端的，移位解调
    elif speakerType == "BottomSpeaker":
        chirp_number = int(np.floor(TopMicChannel.shape[0] / chirp_length))  # Count the number of chirps
        TopMixSignals = np.zeros((chirp_number, chirp_length), dtype=np.complex64)
        BottomMixSignals = np.zeros((chirp_number, chirp_length), dtype=np.complex64)

    # I/Q demodulate the top and bottom mic signals
    for idx in range(chirp_number):
        # ---------------------- Top signal
        TopRec_chirp = TopMicChannel[idx * chirp_length:(idx + 1) * chirp_length]

        # mixter
        chirp_signal_cos = UpChirp_cos(startFreq, stopFreq, sampleRate, tSig)
        top_signal_I = TopRec_chirp * chirp_signal_cos
        chirp_signal_sin = UpChirp_sin(startFreq, stopFreq, sampleRate, tSig)
        top_signal_Q = TopRec_chirp * chirp_signal_sin

        # filter
        ellip_filter = signal.iirfilter(15, Wn=10000, rp=1, rs=60, btype='lowpass',
                                        analog=False, ftype='ellip', fs=48000, output='sos')
        top_signal_I = signal.sosfilt(ellip_filter, top_signal_I)
        top_signal_Q = signal.sosfilt(ellip_filter, top_signal_Q)
        TopMixSignals[idx, :] = top_signal_I - 1j*top_signal_Q

        # ---------------------- Bottom signal
        BottomRec_chirp = BottomMicChannel[idx * chirp_length:(idx + 1) * chirp_length]
        # mixter
        bottom_signal_I = BottomRec_chirp * chirp_signal_cos
        bottom_signal_Q = BottomRec_chirp * chirp_signal_sin

        # filter
        bottom_signal_I = signal.sosfilt(ellip_filter, bottom_signal_I)
        bottom_signal_Q = signal.sosfilt(ellip_filter, bottom_signal_Q)
        BottomMixSignals[idx, :] = bottom_signal_I - 1j*bottom_signal_Q

    if DEBUG_VIS:
        # Vis chirp signal
        _, ax = plt.subplots(2, 2, figsize=(8, 5), sharex=False)
        axes = ax.flatten()
        plt.suptitle("I/Q vector")

        # TopMixSignals_reshape = TopMixSignals.reshape((TopMixSignals.shape[0]*TopMixSignals.shape[1]))
        # axes[0].plot(TopMixSignals_reshape.real[:1000], c='b', label='top_I')
        # axes[0].plot(TopMixSignals_reshape.imag[:1000], c='r', label='top_Q')
        # axes[0].legend()

        axes[0].plot(TopMixSignals[0, :].real, c='b', label='top_I')
        axes[0].plot(TopMixSignals[0, :].imag, c='r', label='top_Q')
        axes[0].legend()

        axes[1].plot(BottomMixSignals[0, :].real, c='b', label='bottom_I')
        axes[1].plot(BottomMixSignals[0, :].imag, c='r', label='bottom_Q')
        axes[1].legend()

        # 计算比值
        divide_T2B = np.divide(BottomMixSignals, TopMixSignals)
        axes[2].plot(divide_T2B[0, :].real, c='b', label='bottom_I')
        axes[2].plot(divide_T2B[0, :].imag, c='r', label='bottom_Q')
        axes[2].legend()

        subtract_T2B = np.subtract(BottomMixSignals, TopMixSignals)
        axes[3].plot(subtract_T2B[0, :].real, c='b', label='bottom_I')
        axes[3].plot(subtract_T2B[0, :].imag, c='r', label='bottom_Q')
        axes[3].legend()

        # # 进行 IQ 分析
        # axes[2].scatter(TopMixSignals[:, :300].real, TopMixSignals[:, :300].imag, c='b', label='Top')
        # axes[2].scatter(TopMixSignals[:, 300:-1].real, TopMixSignals[:, 300:-1].imag, c='r', label='Top_part2')
        # axes[2].legend()
        #
        # axes[3].scatter(BottomMixSignals[:, :300].real, BottomMixSignals[:, :300].imag, c='b', label='Bottom')
        # axes[3].scatter(BottomMixSignals[:, 300:-1].real, BottomMixSignals[:, 300:-1].imag, c='r', label='Bottom_part2')
        # axes[3].legend()

        plt.tight_layout()
        plt.pause(0.05)

    # 对每个解调后地chirp中频信号计算Range-FFT分析
    # ------------------------Top Signal

    TopFFT = np.fft.fft(subtract_T2B, n=128, axis=-1)
    TopFFT = np.fft.fftshift(TopFFT, axes=-1)
    TopFFT_abs_mean = np.abs(np.mean(TopFFT, axis=0))
    TopFFT_abs_mean = TopFFT_abs_mean / np.max(TopFFT_abs_mean)

    # ------------------------Bottom Signal
    BottomFFT = np.fft.fft(BottomMixSignals, n=128, axis=-1)
    BottomFFT = np.fft.fftshift(BottomFFT, axes=-1)
    BottomFFT_abs_mean = np.abs(np.mean(BottomFFT, axis=0))
    BottomFFT_abs_mean = BottomFFT_abs_mean / np.max(BottomFFT_abs_mean)

    if DEBUG_VIS:
        # Vis chirp signal
        _, ax = plt.subplots(2, 2, figsize=(8, 5), sharex=False)
        axes = ax.flatten()
        plt.suptitle("mixed signal FFT")

        axes[0].plot(TopFFT_abs_mean.T)
        axes[0].set_title("Top mic IF signal mean")

        axes[1].plot(BottomFFT_abs_mean.T)
        axes[1].set_title("Bottom mic IF signal mean")

        axes[2].plot(np.abs(TopFFT).T)
        axes[2].set_title("Top mic fft")

        axes[3].plot(np.abs(BottomFFT).T)
        axes[3].set_title("Bottom mic fft")

        plt.tight_layout()
        plt.pause(0.05)

    """ 
        Find the target bin, and obtain the corresponding phase signal
    """
    # Top signal find target bin
    peaks, _ = signal.find_peaks(TopFFT_abs_mean[:int(len(TopFFT_abs_mean)/2)])
    max_peak_index = np.argmax(TopFFT_abs_mean[peaks])
    max_peak_value = peaks[max_peak_index]
    TopFFT_target_bin = TopFFT[:, max_peak_value]
    TopFFT_phase = np.unwrap(np.angle(TopFFT_target_bin))
    TopFFT_angle = np.diff(TopFFT_phase)

    # Bottom signal find target bin
    peaks, _ = signal.find_peaks(BottomFFT_abs_mean, height=0.1*np.max(BottomFFT_abs_mean))
    # max_peak_index = np.argmax(BottomFFT_abs_mean[peaks])
    peak_indexs = np.argsort(BottomFFT_abs_mean[peaks])[-2:]
    peak_values = peaks[peak_indexs]

    BottomFFT_bin_target = BottomFFT[:, peak_values[0]]
    BottomFFT_bin_ref = BottomFFT[:, peak_values[1]]

    BottomFFT_phase_target = np.unwrap(np.angle(BottomFFT_bin_target))
    BottomFFT_phase_ref = np.unwrap(np.angle(BottomFFT_bin_ref))

    # BottomFFT_angle = np.diff(BottomFFT_phase)
    # path4 - path2
    BottomFFT_phase_diff = BottomFFT_phase_target - BottomFFT_phase_ref
    BottomFFT_angle = np.diff(BottomFFT_phase_diff)

    if DEBUG_VIS:
        # Vis chirp signal
        _, ax = plt.subplots(2, 2, figsize=(8, 5))
        axes = ax.flatten()
        plt.suptitle("extract phase signal, range bin :{}".format(max_peak_value))

        axes[0].plot(BottomFFT_bin_target, c='r', label='target path4')
        axes[0].plot(BottomFFT_bin_ref, c='b', label='ref path2')
        axes[0].set_title("phase")

        axes[1].plot(BottomFFT_phase_target, c='r', label='target path4')
        axes[1].plot(BottomFFT_phase_ref, c='b', label='ref path2')
        axes[1].set_title("angle")

        axes[2].plot(BottomFFT_phase_diff)
        axes[2].set_title("phase subject")
        axes[3].plot(BottomFFT_angle)
        axes[3].set_title("angle after phase subject")

        # axes[4].plot(TopFFT_phase)
        # axes[4].set_title("top mic phase signal")
        # axes[5].plot(TopFFT_angle)
        # axes[5].set_title("top mic phase difference signal")

        plt.tight_layout()
        plt.pause(0.05)

    """
        To explore the velocity of phase displacement, we calculate the diff to obtain the velocity.
        
        Calculating signal power based on the Parseval's law  
    """
    #  phase_velocity = Model_L4 / (BottomFFT_phase_diff / (2*np.pi*startFreq) + Model_L2 / Voice_c + (Model_L4 - Model_L2) / Voice_c) - Voice_c
    phase_velocity = BottomFFT_angle

    if DEBUG_VIS:
        plt.figure(figsize=(8, 5))
        plt.plot(phase_velocity)
        plt.xlabel("Time (s)")
        plt.ylabel("Amplitude")
        plt.title("phase_velocity")
        plt.pause(0.05)

    f, vc_time, Zxx = signal.stft(phase_velocity, fs=100, nperseg=32,
                                  noverlap=16, nfft=32, return_onesided=True)

    sum_phase_velocity = np.sum(np.abs(Zxx), axis=0)  # power
    sum_phase_velocity = - 0.356 * (Model_L4 / (sum_phase_velocity / (2 * np.pi * startFreq)
                                                + Model_L2 / Voice_c +
                                                (Model_L4 - Model_L2) / Voice_c) - Voice_c) \
                          * np.pi * Model_r * Model_r * 1000


    if DEBUG_VIS:
        # Vis chirp signal
        _, ax = plt.subplots(2, 1, figsize=(8, 5))
        axes = ax.flatten()

        axes[0].pcolor(np.log10( np.abs(Zxx) ) )
        axes[0].set_xlabel("Time (s)")
        axes[0].set_ylabel("Frequency (Hz)")
        axes[0].set_title("The STFT spectrum of phase velocity")

        axes[1].plot(np.abs(np.fft.fft(TopFFT_angle)))
        plt.pause(0.05)

        # plt.figure(figsize=(8, 5))
        # plt.plot(vc_time, sum_phase_velocity)
        # plt.xlabel("Time (s)")
        # plt.ylabel("velocity (mL/s)")
        # plt.title("The phase velocity power")
        # plt.pause(0.05)

    """
        Kalman filter
    """
    initial_state_mean = 0
    initial_state_covariance = 1
    process_variance = 0.5
    measurement_variance = 0.3
    filtered_states = one_dimensional_kalman_filter(sum_phase_velocity,
                                                    initial_state_mean,
                                                    initial_state_covariance,
                                                    process_variance,
                                                    measurement_variance)
    if DEBUG_VIS:
        plt.figure(figsize=(8, 5))
        plt.plot(vc_time, sum_phase_velocity, label='Original Signal')
        plt.plot(vc_time, filtered_states, label='Filtered Signal', color='coral')
        plt.xlabel('Time (s)')
        plt.ylabel('Flow measurement (mL/s)')
        plt.legend()
        plt.title('Flow measurement estimation')
        plt.pause(0.05)

    """
        吹气范围检测：利用滑动标准差 计算呼吸活动 边沿
    """

    window_size = 5
    rolling_std = np.std(np.lib.stride_tricks.sliding_window_view(filtered_states,
                                                                  window_shape=(window_size,), axis=0),
                         axis=1, ddof=0)
    peaks_rolling_std, _ = signal.find_peaks(rolling_std)
    sorted_peaks = sorted(peaks_rolling_std, key=lambda i: rolling_std[i], reverse=True)
    highest_peaks = sorted_peaks[:2]
    first_peak_value = highest_peaks[0] if (highest_peaks[0] < highest_peaks[1]) else highest_peaks[1]  # + int(window_size // 2)
    second_peak_value = highest_peaks[1] if (highest_peaks[0] < highest_peaks[1]) else highest_peaks[0]
    second_peak_value += int(window_size // 1)

    if DEBUG_VIS:
        plt.figure(figsize=(8, 5))
        plt.plot(vc_time, sum_phase_velocity, color='k', label='Original Signal')
        plt.plot(vc_time, filtered_states, color='green', label='Kalman Filtered Signal')
        # plt.plot(vc_time[0:len(rolling_std)], rolling_std, label=f'Rolling Std (Window Size = {window_size})', linestyle='--', color='red')
        plt.axvline(x=vc_time[first_peak_value], color='r', linestyle='--')
        plt.axvline(x=vc_time[second_peak_value], color='r', linestyle='--')
        plt.xlabel('Time (s)')
        plt.ylabel('Flow measurement (mL/s)')
        plt.legend()
        plt.title('Rolling Standard Deviation')
        plt.pause(0.05)

    """
        积分计算范围内的面积
    """
    total_delivery_filter = 0.0
    total_delivery_org = 0.0
    for i in range(first_peak_value, second_peak_value):
        total_delivery_filter += filtered_states[i] * (vc_time[i] - vc_time[i-1])
        total_delivery_org += sum_phase_velocity[i] * (vc_time[i] - vc_time[i - 1])
    print("total delivery filter is :{} ml \n org:{} ml ".format(total_delivery_filter, total_delivery_org))



    pass









