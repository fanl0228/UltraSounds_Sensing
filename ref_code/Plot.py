import numpy as np
from matplotlib.axes import Axes
from numpy import ndarray
from scipy import signal
import matplotlib
import matplotlib.pyplot as plt

from BaseUtils import UpperEnv1
from Audio import AudioClass
from Parameters import zhfont

matplotlib.rcParams["font.sans-serif"] = "Source Han Sans SC, 思源黑体, 微软雅黑, DejaVu Serif, Arial"

def PlotWave(ax: Axes, audio: AudioClass, ch: int = 0):
    ax.plot(audio.timeSequence, audio.data) # audio[ch]
    ax.set_xlabel('Time (sec)')
    ax.set_ylabel('Amplitude (V)')
    ax.set_title(audio.title + ' wave')
    ax.grid(alpha=0.5)


def PlotSpecgram(ax: Axes, audio: AudioClass, bpData: np.ndarray, noverlap=None, nfft=None,ch: int = 0):
    ax.specgram(bpData,
                Fs=audio.sampleRate,
                cmap='viridis', 
                interpolation='bilinear',
                noverlap=noverlap, 
                NFFT=nfft)
    ax.set_xlabel('Time (sec)')
    ax.set_ylabel('Frequency (Hz)')
    ax.set_title(audio.title + ' specgram')


def PlotSpecgram2(ax: Axes, audio: AudioClass, noverlap=None, nfft=None, ch: int = 0):
    f, t, Sxx = signal.spectrogram(audio[ch], 
                            fs=audio.sampleRate, 
                            noverlap=noverlap, 
                            nperseg=nfft)
    ax.pcolormesh(t, f, np.abs(Sxx), shading='gouraud', cmap='viridis', antialiased=True)
    ax.set_xlabel('Time (sec)')
    ax.set_ylabel('Frequency (Hz)')
    ax.set_title(audio.title + ' STFT Magnitude')


def PlotXcorrUpperEnv(ax: Axes, X: AudioClass, y: AudioClass, ch: int = 0):
    corr = np.abs(np.correlate(X[ch], y.data, mode='same'))
    env = UpperEnv1(corr, 'linear')
    ax.plot(X.timeSequence, env)
    ax.set_xlabel('Time (sec)')
    ax.set_ylabel('Corr')
    ax.set_title('Cross-correlation (upper env)')
    ax.grid(alpha=0.5)
    return env


def PlotXcorr(ax: Axes, X: AudioClass, bpData: np.ndarray, y: AudioClass, ch: int = 0):
    corr = np.abs(np.correlate(bpData, y.data, mode='same')) # X[ch]
    ax.plot(corr) #X.timeSequence, 
    ax.set_xlabel('Time (sec)')
    ax.set_ylabel('Corr')
    ax.set_title('Cross-correlation')
    ax.grid(alpha=0.5)
    return corr


def PlotFrequencyResponse(fftFreq:np.ndarray, fftY: np.ndarray):

    plt.figure()
    plt.plot(fftFreq, fftY)
    plt.title('频率响应-chirp',fontproperties=zhfont, fontsize = 14)
    plt.xlabel('频率（Hz）',fontproperties=zhfont, fontsize = 14)
    plt.ylabel('强度',fontproperties=zhfont, fontsize = 14)
    plt.show()
