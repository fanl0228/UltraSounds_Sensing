# -*- encoding: utf-8 -*-
'''
@File    :   BaseUtil.py
@Time    :   2021/9/2 10:56:00
@Author  :   Dil Duan
@Version :   1.0
@Contact :   1522740702@qq.com
@License :   (C)Copyright 2021

@Update Author github: https://github.com/fanl0228
@Upate Time : 2023-11-22

'''




import math, os
import numpy as np
from math import pi
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
import scipy.signal as signal
from scipy.io import wavfile
from scipy import interpolate, fftpack
from scipy.fftpack import fft, fftfreq


def UpperEnv1(y: np.ndarray, method: str = 'cubic') -> np.ndarray:
    """对给定序列提取upper envelop

    Args:
        y (np.ndarray): 要计算上包络的1维序列y
        method (str, optional): 插值方法：'linear', 'nearest', 'zero', 'slinear', 'quadratic', 'cubic', 'previous', 'next', by default 'linear'

    Returns:
        np.ndarray: upper envelop sequence（长度同x、y）
    """
    x = np.arange(y.shape[0])  # x=0,1,...,n-1
    b = np.full(y.shape[0], False)  # b=0,0,...,0
    for i in x[1:-1]:
        if y[i] > y[i - 1] and y[i] > y[i + 1]:
            b[i] = True
    b[0], b[-1] = True, True  # b=1,...[1 or 0]...,1
    interp = interpolate.interp1d(x[b], y[b], kind=method)  # get interp1d object
    return interp(x)  # do interpolate


def HanningWindow(signal: np.ndarray, pos: int, size: int):
    """ Add Hanning windows
    """
    # signal*=32767
    index = pos
    while (index < pos + size):
        k = index - pos
        signal[index] = (signal[index] * 0.5 * (1.0 - math.cos(2.0 * math.pi * k / size)))
        index += 1
    return signal


def UpChirp(low: int, high: int, rate: int, T: float, offset: float = 0, phase: float = 0, hann: bool = False):
    """
        按采样率、频率和持续时间，获取chirp
        return chirp信号（加窗/不加窗），乘/未乘32767
    """
    length = int(rate * T)
    t = np.zeros(length)
    chirp = np.zeros(length)
    k = (high - low) / T
    for n in range(length):
        t[n] = (float(n) - offset) / rate
        chirp[n] = (math.sin(2 * math.pi * low * t[n] + math.pi * k * t[n] * t[n] + phase))  # *32767

    ## 汉明窗
    # chirp = np.multiply(chirp, signal.hamming(length))

    ## 汉宁窗
    # chirp = np.multiply(chirp, signal.hanning(length))
    if hann:
        chirp = HanningWindow(chirp, 0, length)
    return chirp


def FilterBandpass(wave: np.ndarray, low: int, high: int, rate: int):
    """
        带通滤波
        wave：原始信号
        return: 过滤后的信号，数组
    """
    l = low / (rate / 2.0)
    h = high / (rate / 2.0)
    b, a = signal.butter(10, [l, h], "bandpass")  # 配置滤波器 5/8 表示滤波器的阶数
    return signal.filtfilt(b, a, wave)  # signal.filtfilt(b, a, wave)

    ### 高通滤波 + 低通滤波
    # b, a = signal.butter(5, l, 'highpass')
    # wave = signal.filtfilt(b, a, wave)

    # b1, a1 = signal.butter(5, h, 'lowpass')
    # wave = signal.filtfilt(b1, a1, wave)
    # return wave


def FilterHighpass(wave: np.ndarray, low: int, rate: int):
    """
        高通滤波
        wave：原始信号
        return: 过滤后的信号，数组
    """
    l = low / (rate / 2.0)

    ## 高通滤波 + 低通滤波
    b, a = signal.butter(5, l, 'highpass')
    wave = signal.filtfilt(b, a, wave)

    # b1, a1 = signal.butter(5, h, 'lowpass')
    # wave = signal.filtfilt(b1, a1, wave)
    return wave


def FilterLowpass(wave: np.ndarray, high: int, rate: int):
    """
        低通滤波
        wave：原始信号
        return: 过滤后的信号，数组
    """
    h = high / (rate / 2.0)

    ## 高通滤波 + 低通滤波

    b1, a1 = signal.butter(5, h, 'lowpass')
    wave = signal.filtfilt(b1, a1, wave)
    return wave


def Smooth(data, w=15):  ## w=10
    '''
    smooth raw data
    low-pass filter
    '''
    weight = np.ones(w) / w
    result = np.convolve(weight, data, 'valid')
    return result


def Kalman(zMeasure, xLast: float = 0, pLast: float = 0, Q: float = 0.018, R: float = 0.0542):
    xMid = xLast
    pMid = pLast + Q
    kg = pMid / (pMid + R)
    xNow = xMid + kg * (zMeasure - xMid)
    pNow = (1 - kg) * pMid
    pLast = pNow
    xLast = xNow
    return xNow, pLast, xLast


def FilterGaussian(data, sigma, r):
    gaussTemp = np.ones(2 * r - 1)
    for i in range(0, 2 * r - 1):
        gaussTemp[i] = math.exp(-(i - r) ** 2 / (2 * sigma ** 2)) / (sigma * math.sqrt(2 * math.pi))
    dataFiltered = np.copy(data)
    for i in range(r, len(data) - r):
        dataFiltered[i] = data[i - r:i + r - 1] * np.mat(gaussTemp).T
    return dataFiltered


def Xcorr(receivedSig: np.ndarray, chirp: np.ndarray, low: int, high: int, rate: int):
    """
        与发送的单个chirp做互相关
        received_sig_c：接收到的信号
        dur：chirp信号持续时间
        return：滤波后的信号，互相关结果
    """

    bandpassSig = FilterBandpass(receivedSig, low, high, rate)  # .astype(np.int32)

    # corr1 = np.correlate(bandpassSig, chirp, mode="same")

    corr = signal.correlate(bandpassSig, chirp, mode='full', method='auto')

    # print(np.max(corr))
    # corr = corr / np.max(corr)

    return bandpassSig, corr


def UpInterpolate(data: np.ndarray, times: int):
    """
        上插值 times 倍
        return : 插值后数据
    """
    ## 上插值
    x = np.arange(0, len(data))
    # xNew = np.linspace(x[0],x[-1],len(x)*times) #!!!大小会影响每个区域点的个数
    xNew = np.arange(x[0], x[-1], 1.0 / times)
    # print(len(x), len(xNew))
    func = interpolate.interp1d(x, data, kind='cubic')  # "cubic"
    func1 = interpolate.UnivariateSpline(x, data, s=0)  # 强制通过所有点
    dataNew = func(xNew)
    return dataNew


def UpInterpolate1(x: np.ndarray, data: np.ndarray, times: int):
    """
        上插值 times 倍
        return : 插值后数据
    """
    ## 上插值
    # xNew = np.linspace(x[0],x[-1],len(x)*times) #!!!大小会影响每个区域点的个数
    xNew = np.arange(x[0], x[-1], 1.0 / times)
    # print(len(x), len(xNew))
    func = interpolate.interp1d(x, data, kind=5)  # "cubic"
    func1 = interpolate.UnivariateSpline(x, data, s=0)  # 强制通过所有点
    dataNew = func(xNew)
    return xNew, dataNew


def UpperEnv(data: np.ndarray):
    """
        包络检波, 共两步
        return : smooth
    """
    # 1.1 提取极值点
    x = signal.argrelextrema(data, np.greater)[0]
    x = np.insert(x, 0, 0)
    x = np.append(x, data.size - 1)
    y = data[x]
    if len(x) <= 2: return [], []
    # 1.2 统一坐标轴，插值平滑
    xNew = np.linspace(x[0], x[-1], data.size * 2)  # !!!大小会影响每个区域点的个数
    func = interpolate.interp1d(x, y, kind="cubic")  # "cubic"
    ySmooth = func(xNew)
    return xNew, ySmooth


def HighFrequencyResponse(data: np.ndarray, rate: int = 48000):
    """
        查看高频响应结果
        返回：频率，强度
    """
    fftYAmp = np.abs(fft(data))
    fftFreq = fftfreq(fftYAmp.size, 1.0 / rate)

    index = int(len(fftYAmp) / 2)  # int(np.where(fftFreq>=290)[0][0]) #
    return fftFreq[0:index], fftYAmp[0:index]


def FftAnalyse(mixData: np.ndarray, maxFreq: int, rate: int = 48000):
    """
        进行FFT处理
        返回：频率，强度
    """
    ## 窗可加/可不加：需要加！
    mixData = mixData * signal.hanning(len(mixData))  # 加窗，减少漏频
    fftY = fft(mixData)
    fftYAmp = np.abs(fftY)
    fftFreq = fftfreq(fftY.size, 1.0 / rate)
    index = int(len(fftY) / 2)
    fftYAmp = fftYAmp[0:index]  # / np.max(fftYAmp[0:index])
    phase = np.angle(fftY, deg=True)  # 相位；应该取一半
    return fftFreq[0:index], fftYAmp, fftY[0:index], phase[0:index]


def FftCompute(first: np.ndarray, current: np.ndarray, zeroCounts: int, rate: int = 48000, maxFreq: int = 580):
    """
       先混合信号，再FFT
    """
    first = np.hstack((first, np.zeros(zeroCounts)))
    current = np.hstack((current, np.zeros(zeroCounts)))
    mixData = first * current
    # print(len(mixData))
    fftFreq, fftYAmp, fftY, phase = FftAnalyse(mixData, maxFreq, rate)
    return fftFreq, fftYAmp, fftY, phase


def ForHeatMap(ax: Axes, data: np.ndarray, durTime: float, disAxis: np.ndarray, index: int, vmax: float = 0.2):
    """
        生成热度图
    """
    # plt.figure(figsize=(8,6))
    # while len(data)<=800:
    #    data.append(np.zeros(120))
    ##xticks = np.linspace(0,len(data)/4800,1)
    ##yticks = np.linspace(0,len(data)/48000*340*100/2,)
    # data = np.transpose(data).tolist()
    # sns.heatmap(data, vmax = 0.2)
    # plt.xlabel('Frame Index (Time)',fontproperties=zhfont, fontsize = 14)
    # plt.ylabel('Tap Index (Distance)',fontproperties=zhfont, fontsize = 14)
    # plt.show()

    yAxis = np.round(np.arange(0, len(data), 5), 0)  # 时间
    xAxis = np.round(np.arange(disAxis[0], disAxis[-1], 5), 0)  # 距离
    # 定义画布为1*1个划分，并在第1个位置上进行作图
    # fig = plt.figure(figsize=(10,6))
    # ax = fig.add_subplot(111)
    # 定义横纵坐标的刻度
    ax.set_yticks(yAxis)
    ax.set_xticks(xAxis)
    ax.set_ylabel("Frame Index", fontsize=12)
    ax.set_xlabel("Distance(cm)", fontsize=12)
    # ax.set_title("Ch-"+str(index), fontsize = 12, x = 1)

    # 作图并选择热图的颜色填充风格，这里选择hot
    map = ax.imshow(data, extent=(np.amin(xAxis), np.amax(xAxis), \
                                  min(30, np.amax(yAxis)), np.amin(yAxis)), vmin=0, vmax=vmax)  #

    # 增加右侧的颜色刻度条
    # plt.colorbar(map)
    # ax.axis("scaled")
    # plt.show()
