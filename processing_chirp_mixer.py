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
    filename = 'exp1228/ChirpSignal_15_37_29LR.wav'
    filepath = os.path.join(basepath, filename)

    # generate signal chirp signal
    chirp_signal = UpChirp(startFreq, stopFreq, sampleRate, tSig, hann=False)
    if DEBUG_VIS:
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

    if DEBUG_VIS:
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

    # find reference peaks, 用于对齐信号强度
    speakerType, ref_peaks, other_peaks, \
    TopMic_corr, BottomMic_corr, \
    TopMic_corr_filter, BottomMic_corr_filter = find_reference_peaks(received_signals, chirp_signal,
                                                       startFreq, stopFreq,
                                                       sampleRate, interval)

    if 0:
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
            start_blow = peaks_[0]
            stop_blow = peaks_[-1]
        else:
            start_blow = 0
            stop_blow = -1

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
            start_blow = peaks_[0]
            stop_blow = peaks_[-1]
        else:
            start_blow = 0
            stop_blow = -1
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
            axes[1].axvline(x=peaks_[0], color='red', linestyle='--', label='start airflow')
            axes[1].axvline(x=peaks_[-1], color='red', linestyle='-', label='stop airflow')
            axes[1].legend()
        elif speakerType == "TopSpeaker":
            corr_top = abs(np.diff(top_corr_norm))
            peaks_, _ = signal.find_peaks(corr_top, height=5 * np.mean(corr_top))
            axes[1].plot(corr_top, label='Top corr diff', c='b')
            axes[1].scatter(peaks_, corr_top[peaks_], c='g', label='Peaks')
            axes[1].axvline(x=peaks_[0], color='red', linestyle='--', label='start airflow')
            axes[1].axvline(x=peaks_[-1], color='red', linestyle='-', label='stop airflow')
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

    if 1:
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

    # Test the self-correlation
    if DEBUG:
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

        if DEBUG_VIS:
            # Vis chirp signal
            _, ax = plt.subplots(2, 1, figsize=(8, 5), sharex=False)
            axes = ax.flatten()
            plt.suptitle("I/Q vector")

            # 进行 IQ 分析
            axes[0].scatter(TopMixSignals[:, 5].real, TopMixSignals[:, 5].imag, c='b', label='Top')
            axes[0].legend()

            axes[1].scatter(BottomMixSignals[:, 5].real, BottomMixSignals[:, 5].imag, c='r', label='Bottom')
            axes[1].legend()

            plt.tight_layout()
            plt.pause(0.05)

        # 对每个解调后地chirp中频信号计算Range-FFT分析
        # ------------------------Top Signal
        TopFFT = np.fft.fft(TopMixSignals, n=128, axis=-1)
        TopFFT_abs_mean = np.abs(np.mean(TopFFT, axis=0))
        TopFFT_abs_mean = TopFFT_abs_mean / np.max(TopFFT_abs_mean)
        # ------------------------Bottom Signal
        BottomFFT = np.fft.fft(BottomMixSignals, n=128, axis=-1)
        BottomFFT_abs_mean = np.abs(np.mean(BottomFFT, axis=0))
        BottomFFT_abs_mean = BottomFFT_abs_mean / np.max(BottomFFT_abs_mean)

        if DEBUG_VIS:
            # Vis chirp signal
            plt.figure(figsize=(8, 5))

            plt.plot(TopFFT_abs_mean.T, c='b', label='Top IF signal mean')
            plt.plot(BottomFFT_abs_mean.T, c='r', label='Bottom IF signal mean')
            plt.legend()
            plt.title("mixed signal FFT")

            plt.tight_layout()
            plt.pause(0.05)

        """ 
            Find the target bin, and obtain the corresponding phase signal
        """
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

        if DEBUG_VIS:
            # Vis chirp signal
            _, ax = plt.subplots(2, 1, figsize=(8, 5), sharex=True)
            axes = ax.flatten()

            axes[0].plot(TopFFT_angle)
            axes[0].set_title("Top mic phase difference signal")

            axes[1].plot(BottomFFT_angle)
            axes[1].set_title("Bottom mic phase difference signal")

            plt.tight_layout()
            plt.pause(0.05)

        """
            To explore the velocity of phase displacement, we calculate the diff to obtain the velocity.
            
            Calculating signal power based on the Parseval's law  
        """
        phase_velocity = np.diff(TopFFT_angle)
        f, vc_time, Zxx = signal.stft(phase_velocity, fs=200, nperseg=32,
                                      noverlap=24, nfft=32, return_onesided=True)

        sum_phase_velocity = np.sum(10*np.log10(np.abs(Zxx)), axis=0)  # power
        
        if DEBUG_VIS:
            plt.figure(figsize=(8, 5))
            plt.pcolor(10*np.log10( np.abs(Zxx) ) )
            plt.colorbar()
            plt.xlabel("Time (s)")
            plt.ylabel("Frequency (Hz)")
            plt.title("The STFT spectrum of phase velocity")
            plt.pause(0.05)

            plt.figure(figsize=(8, 5))
            plt.plot(vc_time, sum_phase_velocity)
            plt.xlabel("Time (s)")
            plt.ylabel("Power (dB)")
            plt.title("The phase velocity power")
            plt.pause(0.05)

        """
            Kalman filter
        """
        initial_state_mean = 0
        initial_state_covariance = 1
        process_variance = 0.01
        measurement_variance = 0.3
        filtered_states = one_dimensional_kalman_filter(sum_phase_velocity,
                                                        initial_state_mean,
                                                        initial_state_covariance,
                                                        process_variance,
                                                        measurement_variance)
        if DEBUG_VIS:
            plt.plot(vc_time, sum_phase_velocity, label='Original Signal')
            plt.plot(vc_time, filtered_states, label='Filtered Signal', color='green')
            plt.xlabel('Time (s)')
            plt.ylabel('Amplitude (dB)')
            plt.legend()
            plt.title('Kalman Filter for Signal Estimation')
            plt.pause(0.05)

        """
            吹气范围检测：利用滑动标准差 计算呼吸活动 边沿
        """
        window_size = 5
        rolling_std = np.std(np.lib.stride_tricks.sliding_window_view(filtered_states,
                                                                      window_shape=(window_size,), axis=0),
                             axis=1, ddof=0)
        peaks_rolling_std, _ = signal.find_peaks(rolling_std, )
        sorted_peaks = sorted(peaks_rolling_std, key=lambda i: rolling_std[i], reverse=True)
        highest_peaks = sorted_peaks[:2]
        first_peak_value = highest_peaks[0]
        second_peak_value = highest_peaks[1]

        if DEBUG_VIS:
            plt.figure(figsize=(8, 5))
            plt.plot(vc_time, filtered_states, label='Kalman Filtered Signal')
            plt.plot(vc_time[0:len(rolling_std)], rolling_std, label=f'Rolling Std (Window Size = {window_size})', linestyle='--', color='red')
            plt.axvline(x=vc_time[first_peak_value], color='red', linestyle='--')
            plt.axvline(x=vc_time[second_peak_value], color='red', linestyle='--')
            plt.xlabel('Time (s)')
            plt.ylabel('Amplitude (dB)')
            plt.legend()
            plt.title('Rolling Standard Deviation')
            plt.pause(0.05)

        """
            积分计算范围内的面积
        """
        tou = 0.00406 # 拟合系数
        tube_r = 20  # 管道半径 mm
        total_delivery = 0.0
        for i in range(first_peak_value, second_peak_value):
            instantaneous_delivery = filtered_states[i] * tou * np.pi * tube_r*tube_r
            total_delivery += instantaneous_delivery

        print("total delivery is :{}".format(total_delivery))


    pass









