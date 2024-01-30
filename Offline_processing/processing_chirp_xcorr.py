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

if __name__ == "__main__":
    # 1. Generate Chirp signal
    startFreq = 17000
    stopFreq = 23000
    sampleRate = 48000
    tSig = 5 / 1000  # duration
    tBlank = 0 * tSig
    interval_factor = 0.95
    interval = int((tSig + tBlank) * sampleRate * interval_factor)
    skip_time = 0.25  # 250ms
    chirp_length = int(sampleRate * tSig)

    # 2. Read received ultrasound signal
    basepath = "./data/Music"
    # filename = 'ChirpSignal_22_15_2LR.wav'    #'ChirpSignal_23_15_14LR.wav' #'ChirpSignal_22_15_2LR.wav'
    # filename = 'ChirpSignal_21_29_5LR.wav'  # 肺活量 3495
    filename = 'ChirpSignal_9_55_42LR.wav'
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

    if 0:
        # get top/bottom mic signal
        TopMicChannel = received_signals[0]
        BottomMicChannel = received_signals[1]

        # Vis chirp signal
        _, ax = plt.subplots(2, 1, figsize=(8, 5), sharex=False)
        axes = ax.flatten()
        axes[0].plot(TopMicChannel)
        axes[0].set_title('received Signal Channel 0')
        axes[0].set_xlabel('Time (s)')
        axes[0].set_ylabel('Amplitude')

        axes[1].plot(BottomMicChannel)
        axes[1].set_title('received Signal Channel 1')
        axes[1].set_xlabel('Time (s)')
        axes[1].set_ylabel('Amplitude')

        plt.tight_layout()
        plt.pause(0.05)


    # find reference peaks
    speakerType, ref_peaks, other_peaks, \
    TopMic_corr, BottomMic_corr, \
    TopMic_corr_filter, BottomMic_corr_filter = find_reference_peaks(received_signals, chirp_signal,
                                                       startFreq, stopFreq,
                                                       sampleRate, interval)

    if 0:
        _, ax = plt.subplots(2, 2, figsize=(20, 10), sharex=True)
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

        axes[1].plot(TopMic_corr_filter, c='blue', label='TopMic_corr_filter')
        axes[1].set_title('received Signal Channel TopMic_corr_filter')
        axes[1].set_xlabel('Range Samples')
        axes[1].set_ylabel('Amplitude')
        if speakerType == "TopSpeaker":
            axes[1].scatter(ref_peaks, TopMic_corr_filter[ref_peaks], c='r')
        elif speakerType == "BottomSpeaker":
            axes[1].scatter(other_peaks, TopMic_corr_filter[other_peaks], c='r')
        axes[1].plot(BottomMic_corr_filter, c='coral', label='BottomMic_corr_filter')
        axes[1].set_title('received Signal Channel BottomMic_corr_filter')
        axes[1].set_xlabel('Range Samples')
        axes[1].set_ylabel('Amplitude')
        if speakerType == "TopSpeaker":
            axes[1].scatter(other_peaks, BottomMic_corr_filter[other_peaks], c='g')
        elif speakerType == "BottomSpeaker":
            axes[1].scatter(ref_peaks, BottomMic_corr_filter[ref_peaks], c='g')
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

    # 寻找起始点
    if speakerType == "TopSpeaker":
        TopMic_Peaks, _ = signal.find_peaks(TopMic_corr, distance=0.8 * interval,
                                               height=0.05 * np.max(TopMic_corr))
        second_skip_sample = int( (TopMic_Peaks[0] - chirp_length/2 + chirp_length) % chirp_length)

    elif speakerType == "BottomSpeaker":
        BottomMic_Peaks, _ = signal.find_peaks(BottomMic_corr, distance=0.8 * interval,
                                               height=0.05 * np.max(BottomMic_corr))

        # if BottomMic_Peaks[0] < (chirp_length / 2):
            # 多次实验测试的固定值，不同设备有所区别
        second_skip_sample = int( (BottomMic_Peaks[0] - chirp_length/2 + chirp_length) % chirp_length)

    if 1:
        # get top/bottom mic signal
        TopMicChannel = received_signals[0][second_skip_sample:-1]
        BottomMicChannel = received_signals[1][second_skip_sample:-1]

        # Vis chirp signal
        _, ax = plt.subplots(2, 3, figsize=(15, 10), sharex=False)
        axes = ax.flatten()
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

        _, TopMic_corr = Xcorr(TopMicChannel, TopMicChannel, startFreq, stopFreq, sampleRate)
        _, BottomMic_corr = Xcorr(BottomMicChannel, BottomMicChannel, startFreq, stopFreq, sampleRate)
        axes[2].plot(TopMic_corr, c='r', label='Topmic self-corr')
        axes[2].legend()

        axes[3].plot(BottomMic_corr, c='b', label='Bottommic self-corr')
        axes[3].legend()
        # 模拟解调发射端的，移位解调

        epoch = int(np.floor(BottomMicChannel.shape[0] / chirp_length))
        TopMixSignals = np.zeros((epoch, chirp_length), dtype=np.complex64)
        BottomMixSignals = np.zeros((epoch, chirp_length), dtype=np.complex64)
        for idx in range(epoch):
            # ---------------------- Top signal
            TopRec_chirp = TopMicChannel[idx * chirp_length:(idx + 1) * chirp_length]
            # mixter
            chirp_signal_cos = UpChirp_cos(startFreq, stopFreq, sampleRate, tSig)
            top_signal_I = TopRec_chirp * chirp_signal_cos
            chirp_signal_sin = UpChirp_sin(startFreq, stopFreq, sampleRate, tSig)
            top_signal_Q = TopRec_chirp * chirp_signal_sin

            # filter
            ellip_filter = signal.iirfilter(15, Wn=6000, rp=1, rs=60, btype='lowpass',
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


        # 进行 IQ 分析
        axes[4].scatter(TopMixSignals.real, TopMixSignals.imag, c='b', label='Top')
        axes[4].scatter(BottomMixSignals.real, BottomMixSignals.imag, c='r', label='Bottom')
        axes[4].legend()
        axes[4].set_title("IQ Space")

        # 进行 FFT 分析
        TopFFT = np.fft.fft(TopMixSignals, n=128, axis=-1)
        TopFFT_abs_mean = np.abs(np.mean(TopFFT[:, 0:int(TopFFT.shape[1]//2)], axis=0))
        # TopFFT_abs_mean = np.abs(TopFFT[:, 0:int(TopFFT.shape[1] // 2)])
        TopFFT_abs_mean = TopFFT_abs_mean / np.max(TopFFT_abs_mean)

        BottomFFT = np.fft.fft(BottomMixSignals, n=128, axis=-1)
        BottomFFT_abs_mean = np.abs(np.mean(BottomFFT[:, 0:int(BottomFFT.shape[1]//2)], axis=0))
        BottomFFT_abs_mean = BottomFFT_abs_mean / np.max(BottomFFT_abs_mean)

        axes[5].plot(TopFFT_abs_mean.T)
        # axes[5].plot(TopFFT_abs_mean.T, c='b', label='Top IF signal mean')
        # axes[5].plot(BottomFFT_abs_mean.T, c='r', label='Bottom IF signal mean')
        axes[5].legend()
        axes[5].set_title("mixed signal FFT")
        plt.tight_layout()
        plt.pause(0.05)

        # Top signal find target bin
        peaks, _ = signal.find_peaks(TopFFT_abs_mean)
        max_peak_index = np.argmax(TopFFT_abs_mean[peaks])
        max_peak_value = peaks[max_peak_index]
        TopFFT_angle = np.unwrap(np.angle(TopFFT[:, max_peak_value]))
        TopFFT_angle = np.diff(TopFFT_angle)

        # Bottom signal find target bin
        peaks, _ = signal.find_peaks(BottomFFT_abs_mean)
        max_peak_index = np.argmax(BottomFFT_abs_mean[peaks])
        max_peak_value = peaks[max_peak_index]
        BottomFFT_angle = np.unwrap(np.angle(BottomFFT[:, max_peak_value]))
        BottomFFT_angle = np.diff(BottomFFT_angle)

        fig, ax = plt.subplots(2, 1, sharex=True)
        axes = ax.flatten()
        axes[0].plot(TopFFT_angle)
        axes[0].set_title("TopFFT_angle")

        axes[1].plot(np.diff(TopFFT_angle))
        axes[1].set_title("BottomFFT_angle")
        plt.tight_layout()
        plt.pause(0.05)

        #
        indata = np.diff(TopFFT_angle)
        f, vc_time, Zxx = signal.stft(indata, fs=200, nperseg=32, noverlap=24, nfft=32, return_onesided=False)

        plt.figure()
        plt.pcolor(np.log10( np.abs(Zxx) ) )
        plt.pause(0.05)

        sum_vec = np.sum(np.abs(Zxx), axis=0)
        # sum_vec = sum_vec / np.max(sum_vec)
        plt.figure()
        plt.plot(vc_time, sum_vec)
        plt.pause(0.05)

        # 卡尔曼滤波, 初始化滤波器参数
        initial_state_mean = 0
        initial_state_covariance = 1
        process_variance = 0.1
        measurement_variance = 0.5

        # 使用卡尔曼滤波
        filtered_states = one_dimensional_kalman_filter(sum_vec,
                                                        initial_state_mean,
                                                        initial_state_covariance,
                                                        process_variance,
                                                        measurement_variance)

        # 绘制结果
        plt.plot(vc_time, sum_vec, label='sum_vec Signal')
        plt.plot(vc_time, filtered_states, label='Filtered Signal', color='green')
        plt.xlabel('Time(s)')
        plt.ylabel('Amplitude')
        plt.legend()
        plt.title('Kalman Filter for Signal Estimation')
        plt.pause(0.05)

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

    peaks_, _ = signal.find_peaks(corr_Envy, height=5*np.mean(corr_Envy))
    if len(peaks_) > 1:
        start_blow = peaks_[0]
        stop_blow = peaks_[-1]
    else:
        start_blow = 0
        stop_blow = -1
    if 1:
        _, ax = plt.subplots(2, 2, figsize=(6, 10))
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

        axes[0].plot(time_top, top_corr_norm, label='Top Norm')

        axes[1].plot(time_bottom, bottom_corr_norm, label='Bottom Norm')
        axes[1].set_xlabel("Time (s)")
        axes[1].legend()

        axes[2].plot(corr_Envy, label='ref_corr diff', c='b')
        axes[2].scatter(peaks_, corr_Envy[peaks_], c='r', label='Peaks')
        axes[2].legend()

        if speakerType == "TopSpeaker":
            corr_other = abs(np.diff(BottomMic_corr_filter[other_peaks]))
            peaks_, _ = signal.find_peaks(corr_other, height=5 * np.mean(corr_other))
        elif speakerType == "BottomSpeaker":
            corr_other = abs(np.diff(TopMic_corr_filter[other_peaks]))
            peaks_, _ = signal.find_peaks(corr_other, height=5 * np.mean(corr_other))

        axes[3].plot(corr_other, label='ref_corr diff', c='b')
        axes[3].scatter(peaks_, corr_other[peaks_], c='r', label='Peaks')
        axes[3].legend()

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

        phaseSig = np.reshape(TopMic_upCorrEnvYList[start_blow:stop_blow, bin_num], (1, -1)) #- TopMic_upCorrEnvYList[:, 66]
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









