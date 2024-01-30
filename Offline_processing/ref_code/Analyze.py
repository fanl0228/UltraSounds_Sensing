# -*- encoding: utf-8 -*-
'''
@File    :   Analyze.py
@Time    :   2022/12/04 16:01:00
@MTime   :   2022/12/04 16:01:00
@Author  :   Dil Duan
@Version :   1.0
@Contact :   1522740702@qq.com
@License :   (C)Copyright 2021
'''
"""
    文件说明：声波测风速理论分析
"""


import os, sys
import math
import numpy as np
from math import pi
import scipy.signal as signal
from scipy.io import wavfile
from pathlib import Path
import matplotlib, copy
from matplotlib.axes import Axes
import matplotlib.pyplot as plt
from collections import defaultdict
import matplotlib as mpl

import Fmcw, Alsa, Plot, BaseUtils
from  Audio import AudioClass

mpl.rcParams['font.sans-serif']=['SimHei'] #用来正常显示中文标签
plt.rcParams['font.sans-serif']=['SimHei'] #用来正常显示中文标签
# ['SimSun']宋体；['SimHei']黑体，有很多自己都可以设置
plt.rcParams['axes.unicode_minus'] = False # 正常显示正负号

## global
upTimes = 128  # 1 #
maxDis = 300
mics = [1, 2]


DEBUG = 1

def RangeFft(chirp:np.ndarray, 
               bandpassSig:np.ndarray, 
               peaks: np.ndarray, 
               interval:int, 
               durTime: float,
               bandWidth: int,
               lab: str = "",
               index: int = 10,
               sampleRate:int = 48000):

    #interval # 会影响相位的情况 
    newChirp = np.hstack((chirp, np.zeros(interval-len(chirp))))
    maxFreq = int(maxDis*2.0/100/340*bandWidth/durTime)  # 60cm

    zeroCounts = 2**index - interval # 22500
    startPoint = peaks[0] - 0# 对齐方式，待修改
    if startPoint+int(interval-durTime*sampleRate/2)>len(bandpassSig):
        return

    fCur, yCurAmp, yCur, phase = BaseUtils.FftCompute(newChirp ,bandpassSig[startPoint-int(durTime*sampleRate/2):\
        startPoint+int(interval - durTime*sampleRate/2)], zeroCounts, sampleRate, maxFreq)     # 其他信号
        
    axis = fCur*durTime/bandWidth*340*100/2
    # print(axis[1]-axis[0])
    # plt.plot(axis, yCurAmp, label = lab)
    # plt.show()

    disInterval = 5
    disMin = 25
    disMax = 35
    aimPeaks, _ =signal.find_peaks(yCurAmp, height = 0.5, distance = int(disInterval/(340*100/bandWidth/2)))
    for ap in aimPeaks:
        if axis[ap]>disMin and axis[ap] <disMax and lab!="empty":
            print(lab, axis[ap], phase[ap], axis[0], phase[0])
            break
    
    return axis[ap]


def UpChirp(lowFreq: int, highFreq:int, sampleRate: int, tSig: float, offset:float=0, phase:float =0, hann:bool=True):
    """
        按采样率、频率和持续时间，获取chirp
        return chirp信号（加窗/不加窗），乘/未乘32767
    """
    length = int(sampleRate*tSig)
    t = np.zeros(length)
    chirp = np.zeros(length)
    k = (highFreq-lowFreq)/tSig
    for n in range(length):
        t[n] = (float(n)-offset)/sampleRate # offset为偏移的采样点
        chirp[n]= (np.sin(2*np.pi*lowFreq*t[n]+np.pi*k*t[n]*t[n]+phase)) *32767
    
    if hann: # 需要加窗
        chirp = BaseUtils.HanningWindow(chirp, 0, length)

    return chirp


def Simulation(vel, offset):
    """
    """
    lowFreq = 18000
    highFreq = 23000
    sampleRate = 48000
    tSig = 10 / 1000     # 10  30
    tBlank = 100 / 1000 - tSig
    maxPoints = int(maxDis*2*sampleRate/34000.0)

    nBlank = int(sampleRate * tBlank) 
    interval = int((tSig + tBlank)*sampleRate)
    # chirpWin, SigWin = Fmcw.GetChirp(lowFreq, highFreq, tSig, tBlank, sampleRate, 20)
    # chirpNotWin, SigNotWin = Fmcw.GetChirp(lowFreq, highFreq, tSig, tBlank, sampleRate, 20, False)
    # sig = chirpNotWin.data#

    baseSig = UpChirp(lowFreq, highFreq, sampleRate, tSig)
    sig1 = baseSig
    sig2 = UpChirp(lowFreq, highFreq, sampleRate, tSig, offset=offset) # offset / sampleRate * 340 * 100 / 2, '/2'表示往返

    if DEBUG:
        plt.figure()
        plt.plot(sig1, label ="sig1")
        plt.plot(sig2, label = "sig2")
        plt.legend()
        plt.pause(0.05)

    delay1 = 282#int(30*2/34000*sampleRate)
    simulation = np.hstack((baseSig, np.zeros(nBlank)))
    ob1 = np.hstack((np.zeros(delay1), 0.4*sig1, np.zeros(nBlank-delay1)))
    ob2 = np.hstack((np.zeros(delay1), 0.4*sig2, np.zeros(nBlank-delay1)))

    simulationSig1 = simulation+ob1
    simulationSig2 = simulation+ob2

    bandpassSig1, corr1 = BaseUtils.Xcorr(simulationSig1, baseSig, lowFreq, highFreq, sampleRate)
    bandpassSig2, corr2 = BaseUtils.Xcorr(simulationSig2, baseSig, lowFreq, highFreq, sampleRate)

    peaks1, _ = signal.find_peaks(abs(corr1), height=0.5, distance=interval)
    corrSlice1 = corr1[peaks1[0]:-1]
    peaks2, _ = signal.find_peaks(abs(corr2), height=0.5, distance=interval)
    corrSlice2 = corr2[peaks2[0]:-1]

    # ## FFT 操作
    # dis1 = RangeFft(baseSig, bandpassSig1, peaks1, interval, tSig, highFreq-lowFreq, "sig1", 27)
    # dis2 = RangeFft(baseSig, bandpassSig2, peaks2, interval, tSig, highFreq-lowFreq, "sig2", 27)

    # tdDis = float(offset) / sampleRate * 340 * 100 / 2
    # print("Vel:", vel, " Theoretical dDis:", tdDis, " Calculated dDis: ", abs(dis1-dis2))
    # # plt.xlim(0,50)
    # # plt.show()
    # return

    # 互相关操作
    corrSlice1 = corr1[peaks1[0]-50:peaks1[0]+maxPoints]
    corrSlice2 = corr2[peaks2[0]-50:peaks2[0]+maxPoints]
    ## 上插值看看效果
    upCorrSlice1 = BaseUtils.UpInterpolate(corrSlice1, upTimes)
    upCorrSlice2 = BaseUtils.UpInterpolate(corrSlice2, upTimes)
    #print(np.arange(len(upCorrSlice1))/upTimes)

    if DEBUG:
        plt.figure()
        plt.plot(np.arange(len(upCorrSlice1))/upTimes/sampleRate*340*100/2, upCorrSlice1, label = "ob1")
        plt.plot(np.arange(len(upCorrSlice2))/upTimes/sampleRate*340*100/2, upCorrSlice2,  label = "ob2")
        #plt.xlim(40, 60)
        plt.xlabel('Distance (cm)')
        plt.ylabel('Xcorr')
        plt.tight_layout()
        plt.legend()
        plt.pause(0.05)
        # plt.show()

    corrEnvY1 = np.abs(signal.hilbert(corrSlice1))
    peaks1, _ = signal.find_peaks(corrEnvY1, height=0.5, distance=len(corrEnvY1))
    corrEnvY1 = corrEnvY1[ peaks1[0]: peaks1[0]+(maxPoints-20)-10 ] #-upTimes*30

    corrEnvY2 = np.abs(signal.hilbert(corrSlice2))
    peaks2, _ = signal.find_peaks(corrEnvY2, height=0.5, distance=len(corrEnvY2))
    corrEnvY2 = corrEnvY2[ peaks2[0]: peaks2[0]+(maxPoints-20)-10 ] #-upTimes*30

    corrEnvY1 = corrEnvY1#[260:270]
    corrEnvY2 = corrEnvY2#[260:270]
    upCorrEnvY1 = BaseUtils.UpInterpolate(corrEnvY1, upTimes)
    upCorrEnvY2 = BaseUtils.UpInterpolate(corrEnvY2, upTimes)

    if DEBUG:
        print(len(corrEnvY1), len(upCorrEnvY1))
        plt.figure()
        plt.plot(corrEnvY1, '.', label="1")
        plt.plot(np.arange(len(upCorrEnvY1))/upTimes, upCorrEnvY1,label="2")
        plt.legend()
        plt.pause(0.05)

    corrEnvY1 = corrEnvY1/np.max(corrEnvY1)
    corrEnvY2 = corrEnvY2/np.max(corrEnvY2)
    upCorrEnvY1 = upCorrEnvY1/np.max(upCorrEnvY1)
    upCorrEnvY2 = upCorrEnvY2/np.max(upCorrEnvY2)
    #upTimes=12
    axis = np.arange(len(upCorrEnvY1))/upTimes  #/ sampleRate/ upTime *340*100/2
   
    disInterval = 5
    disMin = 280#25*2/100/340*sampleRate
    disMax = 300#35*2/100/340*sampleRate
    aimPeaks, _ =signal.find_peaks(upCorrEnvY1, height = 0.1, distance = int(disInterval/(340*100)*2*sampleRate*upTimes))
    for ap in aimPeaks:
        if axis[ap]>disMin and axis[ap] <disMax:
            print(axis[ap],  axis[0])
            break
    dis1 = axis[ap]

    aimPeaks, _ =signal.find_peaks(upCorrEnvY2, height = 0.1, distance = int(disInterval/(340*100)*2*sampleRate*upTimes))
    for ap in aimPeaks:
        if axis[ap]>disMin and axis[ap] <disMax:
            print(axis[ap],  axis[0])
            break
    dis2 = axis[ap]
    print(abs(dis1-dis2))
    
    fig, ax = plt.subplots(1, 2, figsize=(16,10), sharex=True, sharey=True)

    curAx = ax.flat[0]
    curAx.plot(corrEnvY1, linewidth=1, label = "仿真信号1")
    curAx.plot(corrEnvY2, linewidth=1, label = "仿真信号2")
    curAx.set_xlabel('采样点', fontsize=35) 
    curAx.set_ylabel('归一化互相关值', fontsize=35)
    curAx.legend(fontsize=20)
    curAx.set_xlim(0, 350)
    #curAx.set_ylim(0,0.45)
    curAx.tick_params(labelsize=20)
    curAx.set_title("插值前", fontsize=30)

    curAx1 = ax.flat[1]
    curAx1.plot(axis, upCorrEnvY1, linewidth=1, label = "仿真信号1")
    curAx1.plot(axis, upCorrEnvY2, linewidth=1, label = "仿真信号2")
    curAx1.set_xlabel('采样点', fontsize=35) 
    #curAx1.set_ylabel('', fontsize=35)
    curAx1.legend(fontsize=20)
    curAx1.set_xlim(0, 350)
    #curAx1.set_ylim(0,0.45)
    curAx1.tick_params(labelsize=20)
    curAx1.set_title("插值后", fontsize=30)
    

    #plt.tight_layout()
    plt.show()
    

def PlotCurve():
    """
    反射式测速: 风速与时延差的关系曲线
    
    """
    L = 0.1
    vs = 340
    sampleRate = 48000
    vw = np.arange(0, 20, 0.1)
    y = []
    yy = []
    prev = 0
    for one in vw:
        #res = 2*L*one*one/(vs*(vs*vs-one*one))*sampleRate
        one = round(one, 2)
        res = (L/(vs+one)+L/(vs-one)-2*L/vs)*sampleRate
        y.append(res)
        yy.append(tuple((one, res)))
        #print((one, res))
        prev = res
    plt.figure()
    plt.plot(vw, y)
    plt.xlabel("Wind speed(m/s)")
    plt.ylabel("Sampling Points")
    # plt.ylabel("Sampling points")
    plt.show()
    return yy


def PlotCurve1():
    """
    直射式测速: 风速与时延差的关系曲线
    
    """
    L = 0.14
    vs = 340
    sampleRate = 48000#48000*2
    vw = np.arange(0, 20, 0.1)
    y = []
    yy = []
    prev = 0
    for one in vw:
        one = round(one, 2)
        res = abs(L/(vs+one)-L/vs)*sampleRate
        y.append(res)
        yy.append(tuple((one, res)))
        print((one, res))
        prev = res
    plt.figure()
    plt.plot(vw, y)
    plt.xlabel("Wind speed(m/s)")
    plt.ylabel("Samping Points")
    plt.show()
    return yy


def PlotIMUWave():
    """
    """
    fileLabel = 'testDirect-v0'
    fileName = f"anemometerData/20221221/note9/10k-23k/hann/48k/{fileLabel}.csv"
    data = [[] for _ in range(6)]
    with open(fileName) as f:
        for line in f.readlines():
            row = line.strip().split(',')
            for i in range(len(row)):
                if row[i]=='':
                    row[i] = '"-1"'
            #print(row)
            try:
                row_data = np.array([row[i].strip('\"') for i in range(len(row))], dtype=np.float32)
                for i in range(len(row_data)):
                    if row_data[i]!=-1:
                        if i == 3:
                            row_data[i]-=8.03
                        if i == 4:
                            row_data[i]+=0.425
                        if i == 5:
                            row_data[i]-=5.675
                        data[i].append(row_data[i])
                        
            except Exception as e:
                print(e)
                continue
    plt.figure()
    plt.plot(data[0], label="Gyro-X")
    plt.plot(data[1], label="Gyro-Y")
    plt.plot(data[2], label="Gyro-Z")
    plt.legend()
    plt.figure()
    plt.plot(data[3], label="Acc-X")
    plt.plot(data[4], label="Acc-Y")
    plt.plot(data[5], label="Acc-Z")
    plt.legend()
    plt.show()
    return np.array(data)


if __name__ =="__main__":
    #PlotIMUWave()
    #sys.exit()

    #velSample = PlotCurve1()
    #sys.exit(0)

    velSample = [tuple((7.9, 0.05))]
    for one in velSample:
        print("-------------------------")
        print(one[0], one[1])
        Simulation(one[0], one[1])
    plt.show()
    sys.exit(0)
