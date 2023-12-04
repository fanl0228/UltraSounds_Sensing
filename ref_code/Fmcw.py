import math
import numpy as np
from pathlib import Path
from scipy import signal
from scipy.io import wavfile


from Audio import AudioClass


def UpChirp(low: int, high:int, rate: int, T: float, offset:float=0, phase:float =0, hann: bool = False):
    """
        按采样率、频率和持续时间，获取chirp
        return chirp信号（加窗/不加窗），乘/未乘32767
    """
    amplitude = [2.0, 1.6, 1.6, 1.7, 1.2, 1.2, 1.3, 1.2, 1.7, 1.5, 2.9, 6.7]
    length = int(rate*T)
    t = np.zeros(length)
    chirp = np.zeros(length)
    k = (high-low)/T
    maxVal = -1
    for n in range(length):
        index = n // 200
        t[n] = (float(n)-offset)/rate
        chirp[n]= amplitude[index]*(math.sin(2*math.pi*low*t[n]+math.pi*k*t[n]*t[n]))
        maxVal = max(maxVal, abs(chirp[n]))
    
    for i in range(len(chirp)):
        chirp[i] = chirp[i] / maxVal

    return chirp


def GetChirp(fMin: float, fMax: float, tSig: float, tBlank: float = 0,
              rate: int = 48000, repeat: int = 1, windows: bool = True):
    """自动生成线性调频（chirp）音频，单声道
    Args:
        fMin (float): 最低频率
        fMax (float): 最高频率
        tSig (float): 信号长度，单位为秒
        tBlank (float, optional): 信号之后空白长度，单位为秒. Defaults to 0.
        rate (int, optional): 采样率. Defaults to 48000.
        repeat (int, optional): 循环多少次. Defaults to 1.

    Returns:
        信号音频， 全部音频
    """
    nSig = int(tSig * rate)
    nBlank = int(rate * tBlank) 
    
    #ts = np.linspace(0, tsig, n_sig)  # 生成时间序列
    ts =  np.arange(0, tSig, 1.0/rate)

    sig = signal.chirp(ts, f0=fMin, t1=tSig, f1=fMax, method='linear')  # 生成chirp
    #sig = UpChirp(fMin, fMax, rate, tSig)
    # add window or not
    sigWindowed = sig #np.multiply(sig, signal.hanning(nSig))  # 将hanning窗作用于chirp
    if windows:
        sigWindowed = np.multiply(sig, signal.hanning(nSig))

    l = np.concatenate((sigWindowed, np.zeros(nBlank)))
    totalWav = np.tile(l, repeat)  # 重复若干次

    totalTime = totalWav / rate

    chirpAudio =  AudioClass(sigWindowed, rate, f"chirp-{fMin}-{fMax}Hz-{nSig}-{nBlank}")
    totalAudio = AudioClass(totalWav, rate, f"total-{fMin}-{fMax}Hz-{nSig}-{nBlank}-{repeat}")

    return chirpAudio, totalAudio


def GetLRChirp(fMin: float, fMax: float, tSig: float, tBlank: float = 0,
                rate: int = 48000, repeat: int = 1):
    """自动生成线性调频（chirp）音频，双声道

    Args:
        fMin (float): 最低频率
        fMax (float): 最高频率
        tSig (float): 信号长度，单位为秒
        tBlank (float, optional): 信号之后空白长度，单位为秒. Defaults to 0.
        rate (int, optional): 采样率. Defaults to 48000.
        repeat (int, optional): 循环多少次. Defaults to 1.

    Returns:
        信号音频， 全部音频
    """
    nSig = int(tSig * rate)
    nBlank = int(rate * tBlank)
    nTotal = nSig + nBlank 
    
    #ts = np.linspace(0, tsig, n_sig)  # 生成时间序列
    ts =  np.arange(0, tSig, 1.0/rate)

    sig = signal.chirp(ts, f0=fMin, t1=tSig, f1=fMax, method='linear')  # 生成chirp
    sigWindowed = np.multiply(sig, signal.hanning(nSig))  # 将hanning窗作用于chirp

    l = np.concatenate((sigWindowed, np.zeros(nBlank), np.zeros(nTotal)))
    r = np.concatenate((np.zeros(nTotal), sigWindowed, np.zeros(nBlank)))

    totalWavL = np.tile(l, repeat)  # 重复若干次
    totalWavR = np.tile(r, repeat)  # 重复若干次

    blank = int(rate * 0.2)
    totalWavL =  np.concatenate((np.zeros(blank), totalWavL))
    totalWavR =  np.concatenate((np.zeros(blank), totalWavR))

    totalWav = np.column_stack((totalWavL, totalWavR))
    totalTime = totalWav.shape[0]/rate

    chirpAudio =  AudioClass(sig, rate, f"chirp-{fMin}-{fMax}Hz-{nSig}-{nBlank}")
    totalAudio = AudioClass(totalWav, rate, f"total-{fMin}-{fMax}Hz-{nSig}-{nBlank}")

    return chirpAudio, totalAudio


def SaveAudio(fileName, sampleRate, data):
    """
        保存音频
    """
    if fileName == None:
        fileName = "NoneFileName.wav"
    if Path(fileName).suffix != '.wav':
        fileName += ".wav"
    Path(fileName).parent.mkdir(parents=True, exist_ok=True)
    data_s16 = (data * (2**15 - 1)).astype(np.int16) # 默认保存为int16 兼容性好
    wavfile.write(fileName, sampleRate, data_s16)
    return fileName
