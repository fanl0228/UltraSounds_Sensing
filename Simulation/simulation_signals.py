import numpy as np
import math
from scipy.signal import chirp, spectrogram
from scipy import signal
import matplotlib.pyplot as plt

import sys
sys.path.append("../")
from Signal_Processing.Chirp import UpChirp


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
    chirp_signal = UpChirp(startFreq, stopFreq, signalFS, chirp_duration, offset=0, phase=0)  # chirpfs = 48000*0.005 = 240

    # received signal simulation
    temp = UpChirp(startFreq, stopFreq, signalFS, chirp_duration, offset=0, phase=0)
    delay_time = 0.0001  # signal delay receive  (seconds)
    delay_sample_for_chirp = int(delay_time * signalFS)
    receive_signal = np.zeros((int(signalFS * chirp_duration),))
    receive_signal[delay_sample_for_chirp:] = temp[:-delay_sample_for_chirp]

    # xcor => received signal corr transmitting signal ==》and obtain the second half of the data.
    corrSig = signal.correlate(receive_signal, chirp_signal, mode='full', method='auto')
    corrSigAbsHalf = np.abs(corrSig[chirpFs-1:])  #int(np.floor(corrSig/2))

    # plot_spectrogram("txchirp", chirp_signal, chirpFs)
    # plot_spectrogram("rxchirp", receive_signal, chirpFs)
    
    #
    fig, ax = plt.subplots(2, 1)
    axes = ax.flatten()
    plt.suptitle("Chirp Signal Analysis")

    axes[0].plot(chirp_signal)
    axes[0].plot(receive_signal)
    axes[0].set_xlabel('t (sec)')

    axes[1].plot(corrSigAbsHalf)
    axes[1].set_xlabel('t (sec)')
    plt.tight_layout()
    plt.pause(0.05)

    pdb.set_trace()


    # ========================================================== ZC-code
    # gendrate zc-code
    u = 1
    N =128
    n = np.arange(N)
    ZC_signal = np.exp(-1j * np.pi * u * n * (n+1)/(N-1))





    pass


