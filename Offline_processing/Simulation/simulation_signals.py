import numpy as np
import math
from scipy.signal import chirp, spectrogram, hilbert
from scipy import signal
import matplotlib.pyplot as plt

import sys
sys.path.append("../")
from Signal_Processing import Chirp


import pdb


def plot_spectrogram(title, sig, fs, nperseg=12, nfft=16):
    ff, tt, Sxx = spectrogram(sig, fs=fs, nperseg=nperseg, nfft=nfft)
    fig, ax = plt.subplots()
    ax.pcolormesh(tt, ff, Sxx, cmap='gray_r', shading='gouraud')
    ax.set_title(title)
    ax.set_xlabel('t (sec)')
    ax.set_ylabel('Frequency (Hz)')
    ax.grid(True)
    plt.pause(0.05)


if __name__ == '__main__':
    signalFS = 48000
    startFreq = 17000
    stopFreq = 23000
    chirp_duration = 0.005
    chirpFs = int(signalFS * chirp_duration)

    # Transmitting signal generate
    chirp_signal = Chirp.UpChirp_sin(startFreq, stopFreq, signalFS, chirp_duration)  # chirpfs = 48000*0.005 = 240

    # received signal simulation
    delay_time = 0 # signal delay receive  (seconds)
    delay_sample = int(delay_time * signalFS)
    receive_signal = np.roll(chirp_signal, delay_sample)

    # correlation analysis
    if 0:
        # xcor => received signal corr transmitting signal ==》and obtain the second half of the data.
        corrSig = signal.correlate(receive_signal, chirp_signal, mode='full', method='auto')
        # corrSigAbsHalf = np.abs(corrSig[chirpFs - 1:])  #int(np.floor(corrSig/2))
        corrSigAbsHalf = np.abs(corrSig)  # int(np.floor(corrSig/2))

        fig, ax = plt.subplots(2, 1)
        axes = ax.flatten()
        plt.suptitle("Chirp Signal xcorr Analysis")

        axes[0].plot(chirp_signal)
        axes[0].plot(receive_signal)
        axes[0].set_xlabel('samples')

        axes[1].plot(corrSigAbsHalf)
        axes[1].set_xlabel('samples')
        plt.tight_layout()
        plt.pause(0.05)

    # hilbert analysis
    if 1:
        #
        receive_signal_q = hilbert(receive_signal)
        fig, ax = plt.subplots(2, 1)
        axes = ax.flatten()
        plt.suptitle("hilbert signal analysis")

        receive_signal_real = np.real(receive_signal_q)
        receive_signal_imag = np.imag(receive_signal_q)

        axes[0].scatter(receive_signal_imag, receive_signal_real)

        instantaneous_phase = np.unwrap(np.angle(receive_signal_q))
        instantaneous_freq = np.diff(instantaneous_phase)

        axes[1].plot(instantaneous_phase)
        axes[1].plot(instantaneous_freq)

        plt.pause(0.05)

    # demoludation
    if 1:
        IF_signal = receive_signal * chirp_signal

        fig, ax = plt.subplots(2, 2)
        axes = ax.flatten()
        plt.suptitle("demoludation signal analysis")

        axes[0].plot(IF_signal)
        IF_signal_fft = np.abs(np.fft.fftshift(np.fft.fft(IF_signal) ) )
        axes[1].plot(IF_signal_fft)

        # LP-Filter
        ellip_filter = signal.iirfilter(15, Wn=6000, rp=1, rs=60, btype='lowpass',
                                        analog=False, ftype='ellip', fs=48000, output='sos')
        IF_signal_filter = signal.sosfilt(ellip_filter, IF_signal)

        axes[2].plot(IF_signal_filter)
        IF_signal_filter_fft = np.abs(np.fft.fftshift(np.fft.fft(IF_signal_filter) ) )
        axes[3].plot(IF_signal_filter_fft)

        plt.pause(0.05)

    if 1:
        chirp_signal_I = Chirp.UpChirp_sin(startFreq, stopFreq, signalFS, chirp_duration)
        IF_signal_I = receive_signal * chirp_signal_I

        chirp_signal_Q = Chirp.UpChirp_cos(startFreq, stopFreq, signalFS, chirp_duration)
        IF_signal_Q = receive_signal * chirp_signal_Q

        # IQ channel filter
        ellip_filter = signal.iirfilter(15, Wn=6000, rp=1, rs=60, btype='lowpass',
                                        analog=False, ftype='ellip', fs=48000, output='sos')
        IF_signal_I = signal.sosfilt(ellip_filter, IF_signal_I)
        IF_signal_Q = signal.sosfilt(ellip_filter, IF_signal_Q)

        IF_signal_compex = IF_signal_I + 1j * IF_signal_Q

        fig, ax = plt.subplots(2, 1)
        axes = ax.flatten()
        plt.suptitle("demoludation IQ signal analysis")

        axes[0].scatter(IF_signal_I, IF_signal_Q)

        instantaneous_phase = np.unwrap(np.angle(IF_signal_compex))
        axes[1].plot(instantaneous_phase)


        plt.pause(0.05)



    print("==========================================================\n")

    # ========================================================== ZC-code
    # gendrate zc-code
    u = 1
    N = 479
    signalFS = 48000
    signal_duration = 0.010
    signal_samples = int(signalFS * signal_duration)  # 480
    n = np.arange(signal_samples)
    ZC_signal = np.exp(-1j * np.pi * u * n * (n+1)/(N-1))

    if 1:
        real_signal = np.real(ZC_signal)
        imag_signal = np.imag(ZC_signal)

        fig, ax = plt.subplots(1, 2, figsize=(10, 5), sharex=False)
        axes = ax.flatten()
        axes[0].scatter(real_signal, imag_signal)

        axes[1].plot(real_signal, c='teal')
        axes[1].plot(imag_signal, c='cyan')

        # plt.tight_layout()
        plt.pause(0.05)

    # xcorr
    signalFS = 48000
    delay_time = 0.003  # signal delay receive  (seconds)
    delay_sample = int(delay_time * signalFS)
    ZC_received_signal = np.roll(ZC_signal, delay_sample)

    corrSig = signal.correlate(ZC_received_signal, ZC_signal, mode='full', method='auto')
    corrSigAbsHalf = np.abs(corrSig)  #int(np.floor(corrSig/2))
    if 1:
        fig, ax = plt.subplots(2, 1)
        axes = ax.flatten()
        plt.suptitle("ZC-code Signal Analysis")

        axes[0].plot(ZC_signal)
        axes[0].plot(ZC_received_signal)
        axes[0].set_xlabel('samples')

        axes[1].plot(corrSigAbsHalf)
        axes[1].set_xlabel('samples')
        plt.tight_layout()
        plt.pause(0.05)


    print("==========================================================\n")

    pass


