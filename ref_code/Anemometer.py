# -*- encoding: utf-8 -*-
'''
@File    :   Anemometer.py
@Time    :   2021/9/2 10:41:00
@MTime   :   2022/11/27 14:53:00
@Author  :   Dil Duan
@Version :   1.0
@Contact :   1522740702@qq.com
@License :   (C)Copyright 2021
'''
"""
    文件说明：声波测风速
"""

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

import Fmcw, Alsa, Plot, BaseUtils
from  Audio import AudioClass


mpl.rcParams['font.sans-serif']=['SimHei'] #用来正常显示中文标签
plt.rcParams['font.sans-serif']=['SimHei'] #用来正常显示中文标签
# 在本例中，坐标轴刻度和图例均用新罗马字体['TimesNewRoman']来表示
# ['SimSun']宋体；['SimHei']黑体，有很多自己都可以设置
plt.rcParams['axes.unicode_minus'] = False # 正常显示正负号

## 全局变量
upTimes = 1  # 上插值倍数
maxDis = 60
mics=[1,2,]
peakResults = [] # 存储往返两条路径目标峰值对应的时延


DEBUG = 1


def PlotEmpty(chirp: AudioClass, lowFreq:int, highFreq:int, sampleRate:int, tSig:float, interval:int):
    """
        plot background signal
    """
    fileName = "anemometerData/20221204/empty.wav"

    audio = AudioClass.load(filename=fileName)
    fig, ax = plt.subplots(2, 1, figsize=(16,10))

    index=0
    for i in mics:
        aimChannel = audio[i-1]

        _, corr = BaseUtils.Xcorr(aimChannel, chirp.data, lowFreq, highFreq, sampleRate)
        peaks, _ = signal.find_peaks(abs(corr), height=0.5, distance=interval)
        #print(peaks)

        curAx = ax.flat[index]
        index=index+1
        ## 互相关方式
        XcorrCompute(curAx, peaks, corr, i, "empty")  #
    return fig, ax


def XcorrCompute(curAx: Axes, peaks: np.ndarray, corr: np.ndarray, index: int, interval: int,
                 xPeaks: np.ndarray = [],
                 lab: str = "", sampleRate:int = 48000):
    """
        通过互相关方式处理每个麦克风的数据
        return: disAxis, corrEnvData——每一帧数据对应的坐标轴和互相关上包络的集合
    """

    maxPoints = int(maxDis*2*sampleRate/34000.0)

    i = 0
    count = 0
    prevPeak = -1
    aveEnvY = np.zeros(maxPoints+50)
    xOffsets = np.zeros(len(peaks))
    axisList = []
    upCorrEnvYList = []
    while i<len(peaks):

        i+= 1
        if not ((i>=20 and i<=len(peaks)-5) ): # 经验, 去除前后两端一定范围的数据 
            continue

        # if prevPeak!=-1 and peaks[i]-prevPeak>interval*8:
        #     i+=1
        #     prevPeak = peaks[i]
        #     continue

        prevPeak = peaks[i]
        p = peaks[i]

        corrSlice = corr[p - 50: p + maxPoints]
        #peaks1, _ = signal.find_peaks(upCorrSlice, height=0.2, distance=len(upCorrSlice))
        #upCorrSlice1 = upCorrSlice[ peaks1[0]: peaks1[0]+upTimes*(maxPoints-20)-10 ] #-upTimes*30
        #if(len(corrEnvY)!=(upTimes-4)*maxPoints-10): #-10
        #    break

        # plt.figure(7)
        # plt.plot(corrSlice)
        # plt.pause(0.05)

        # 包络
        corrEnvY = np.abs(signal.hilbert(corrSlice))
        #print(len(corrEnvY), len(corrSlice))

        # 对包络信号插值
        upCorrEnvY = BaseUtils.UpInterpolate(corrEnvY, upTimes)

        # 直射信号峰值移动到0坐标位置
        xOffset = -1
        if len(xPeaks)==0:
            peaks1, _ = signal.find_peaks(upCorrEnvY, height=0.2, distance=len(upCorrEnvY))
            xOffset = peaks1[0]
        else:
            xOffset = int(xPeaks[i])

        upCorrEnvY = upCorrEnvY[xOffset: xOffset + maxPoints * upTimes - 30 * upTimes]  # -upTimes*30

        #upCorrEnvY = upCorrEnvY

        # if len(corrEnvY)!=maxPoints-50:
        #     continue
        xOffsets[i] = xOffset
    
        count+=1
        aveEnvY = aveEnvY+corrEnvY
        #break

        aveEnvY = aveEnvY / count
        axis = np.arange(len(upCorrEnvY))/upTimes
        curAx.plot(axis, upCorrEnvY, label = lab+'-Env')

        curAx.set_title("Mic"+str(index+1))
        curAx.set_xlabel("Sampling points")
        curAx.set_ylabel("Amplitude")
        curAx.set_xlim(0, 200)
        axisList.append(axis)
        upCorrEnvYList.append(upCorrEnvY)

    print(count)
    return axisList, upCorrEnvYList, xOffsets


def RangeFft(curAx: Axes, chirp:np.ndarray,
               bandpassSig:np.ndarray,
               peaks: np.ndarray,
               interval:int,
               durTime: float,
               bandWidth: int,
               xPeaks: np.ndarray = [],
               lab: str = "",
               index: int = 10,
               sampleRate:int = 48000):

    #interval=int(interval/2) # 会影响相位的情况
    print(interval)
    newChirp = np.hstack((chirp, np.zeros(interval-len(chirp))))
    maxFreq = int(maxDis*2.0/100/340*bandWidth/durTime)  # 60cm

    zeroCounts = 0 #2**index - interval # 22500

    i = 0
    prevPeak = -1
    #print(len(peaks))
    xOffsets = np.zeros(len(peaks))
    axisList = []
    yCurAmpList = []
    phaseList = []

    while i<len(peaks)-1:

        i += 1
        if not ((i>=2 and i<=len(peaks)-40) ): # 经验, 去除前后两端一定范围的数据  or (i>=len(peaks)-100 and i<=len(peaks)-70)
            continue

        if prevPeak!=-1 and peaks[i]-prevPeak>interval*5:
            #print(len(axisList))
            i+=1
            prevPeak = peaks[i]
            continue
        prevPeak = peaks[i]

        startPoint = peaks[i] - 30# 对齐方式，待修改
        if startPoint+int(interval-durTime*sampleRate/2)>len(bandpassSig):
            return

        fCur, yCurAmp, yCur, phase = BaseUtils.FftCompute(newChirp ,bandpassSig[startPoint-int(durTime*sampleRate/2):\
            startPoint+int(interval - durTime*sampleRate/2)], zeroCounts, sampleRate, maxFreq)     # 其他信号

        xOffset = 0 
        if len(xPeaks) == 0:
            lpeaks1, _ = signal.find_peaks(yCurAmp, height=2*1e-15, distance=len(yCurAmp))
            xOffset = lpeaks1[0]
            xOffsets[i] = xOffset
        else:
            xOffset = int(xPeaks[i])
            pass

        fCur1 = fCur[xOffset:-1]
        yCurAmp = yCurAmp[xOffset:-1]
        yCur = yCur[xOffset:-1]
        phase = phase[xOffset:-1]

        axis = (fCur1-fCur[xOffset])*durTime/bandWidth*sampleRate
        #axis = (fCur)*durTime/bandWidth*sampleRate
        #print(axis[1]-axis[0])
        curAx.plot(axis, yCurAmp, label = lab)
        curAx.set_xlim(0, 20)
        # print(i, axis[0], phase[0])

        axisList.append(axis)
        yCurAmpList.append(yCurAmp)
        phaseList.append(phase)

    return axisList, yCurAmpList, phaseList, xOffsets


def FindPeaks(axisList, upCorrEnvYList, isTop, isFirst):
    """
        寻找空气直射路径下目标峰值对应的时延
    """
    sampleInterval = 4  # 需转换成采样点
    
    if isTop and isFirst:
        xMin = 6#16
        xMax = 16#24
        topXX = np.zeros(len(axisList))
        topZX = np.zeros(len(axisList))
        for i in range(len(axisList)):
            axis = axisList[i]
            upCorrEnvY = upCorrEnvYList[i]
            aimPeaks, _ =signal.find_peaks(upCorrEnvY, height = 0.1, distance = int(sampleInterval*upTimes))
            for ap in aimPeaks:
                if axis[ap]>xMin and axis[ap] <xMax:
                    print("xAx", axis[ap], upCorrEnvY[ap])
                    topXX[i] = axis[ap]
                    break
            if topXX[i]==0:
                topXX[i]=topXX[i-1]
        
        #topXXFiltered = BaseUtils.FilterBandpass(topXX, 2, 10, 50)
        # plt.figure()
        # plt.title("TopSpeaker")
        # plt.plot(topXX)
        # plt.plot(topXXFiltered)


    if not isTop and isFirst:
        xMin = 6#16
        xMax = 16#24
        topXX = np.zeros(len(axisList))
        topZX = np.zeros(len(axisList))
        for i in range(len(axisList)):
            axis = axisList[i]
            upCorrEnvY = upCorrEnvYList[i]
            aimPeaks, _ =signal.find_peaks(upCorrEnvY, height = 0.01, distance = int(sampleInterval*upTimes))
            for ap in aimPeaks:
                if axis[ap]>xMin and axis[ap] <xMax:
                    print("xAx", axis[ap], upCorrEnvY[ap])
                    topXX[i] = axis[ap]
                    break
            if topXX[i]==0:
                topXX[i]=topXX[i-1]
        
        #topXXFiltered = BaseUtils.FilterBandpass(topXX, 2, 10, 50)
            
        # plt.figure()
        # plt.title("BottomSpeaker")
        # plt.plot(topXX)
        # plt.plot(topXXFiltered)
    
    
    if isTop and not isFirst:
        xMin0 = 15
        xMax0 = 20
        topTxBottomRxX0 = np.zeros(len(axisList))
        xMin =  23.8-2
        xMax =  24.5
        topTxBottomRxX = np.zeros(len(axisList))
        for i in range(len(axisList)):
            axis = axisList[i]
            upCorrEnvY = upCorrEnvYList[i]
            aimPeaks, _ =signal.find_peaks(upCorrEnvY, height = 0.01, distance = int(sampleInterval*upTimes))
            for ap in aimPeaks:
                if axis[ap]>xMin0 and axis[ap] <xMax0:
                    print("topSpeakerBottomMic0", axis[ap], upCorrEnvY[ap])
                    topTxBottomRxX0[i] = axis[ap]#-6.4#-4.5
                    break

            for ap in aimPeaks:
                if axis[ap]>xMin and axis[ap] <xMax:
                    print("topSpeakerBottomMic", axis[ap], upCorrEnvY[ap])
                    delta = 0
                    if i <200 and axis[ap] > 22.6: # 修正前序的环境波动（经验值）
                        delta = int(axis[ap] * 100 % 10) /100
                        print(delta)
                    topTxBottomRxX[i] = axis[ap]-delta#-36.7#-18.2
                    break
            if topTxBottomRxX[i]==0:
                topTxBottomRxX[i]=topTxBottomRxX[i-1]

        xLast = topTxBottomRxX[0]
        pLast = topTxBottomRxX[0]
        topTxBottomRxXFiltered = np.zeros(len(topTxBottomRxX))
        for i in range(len(topTxBottomRxXFiltered)):
            pred, pLast, xLast= BaseUtils.Kalman(topTxBottomRxX[i], xLast, pLast, 0.5, 20) #
            topTxBottomRxXFiltered[i]=pred

        topTxBottomRxXFiltered = BaseUtils.Smooth(topTxBottomRxXFiltered, 3)
        #topTxBottomRxXFiltered = BaseUtils.FilterGaussian(topTxBottomRxX, 2, 3)

        peakResults.append(topTxBottomRxXFiltered)
        

    if not isTop and not isFirst:
        xMin0 = 0
        xMax0 = 10
        bottomTxTopRxX0 = np.zeros(len(axisList))
        xMin = 19
        xMax = 21#29
        bottomTxTopRxX = np.zeros(len(axisList))
        for i in range(len(axisList)):
            axis = axisList[i]
            upCorrEnvY = upCorrEnvYList[i]
            aimPeaks, _ =signal.find_peaks(upCorrEnvY, height = 0.01, distance = int(sampleInterval*upTimes))
            for ap in aimPeaks:
                if axis[ap]>xMin0 and axis[ap] <xMax0:
                    print("bottomSpeakerTopMic0", axis[ap], upCorrEnvY[ap])
                    bottomTxTopRxX0[i] = axis[ap]#-7.4#-4.5
                    break
            
            for ap in aimPeaks:
                if axis[ap]>xMin and axis[ap] <xMax:
                    print("bottomSpeakerTopMic", axis[ap], upCorrEnvY[ap])
                    bottomTxTopRxX[i] = axis[ap]#-46.8#-23.5
                    break
            if bottomTxTopRxX[i]==0:
                bottomTxTopRxX[i]=bottomTxTopRxX[i-1]

        xLast = bottomTxTopRxX[0]
        pLast = bottomTxTopRxX[0]
        bottomTxTopRxXFiltered = np.zeros(len(bottomTxTopRxX))
        for i in range(len(bottomTxTopRxXFiltered)):
            pred, pLast, xLast= BaseUtils.Kalman(bottomTxTopRxX[i], xLast, pLast, 1.5, 20) #卡尔曼滤波
            bottomTxTopRxXFiltered[i]=pred

        bottomTxTopRxXFiltered = BaseUtils.Smooth(bottomTxTopRxXFiltered, 3)
        #bottomTxTopRxXFiltered = BaseUtils.FilterGaussian(bottomTxTopRxX, 2, 3)

        peakResults.append(bottomTxTopRxXFiltered)

        plt.figure()
        plt.plot(bottomTxTopRxX, label = "raw data", linewidth=2)
        plt.plot(bottomTxTopRxXFiltered, label = "kalman filter", linewidth=2)
        plt.yticks(fontsize=10)
        plt.xticks(fontsize=10)
        plt.ylabel('peak delay of target', fontsize=10)
        plt.xlabel('发射周期下标', fontsize=10)
        plt.legend(fontsize=10)
        plt.pause(0.05)
        #plt.show()


def ProcessBySpeaker(audio:AudioClass, chirp:np.ndarray, lowFreq:int, highFreq:int,
                    sampleRate:int, interval:int, fileLabel:str, isTop:bool = False):
    """
        距离 Speaker 近的 Mic 强度更大
        小米NotePro: BottomMic 对应 channel[0];
                    TopMic 对应 channel[1]
    """

    fig, ax = plt.subplots(2, 1, figsize=(16,10), sharex=True)
    mics = [1, 0]
    speaker = "BottomSpeaker"
    if isTop:
        mics = [0, 1]
        speaker = "TopSpeaker"

    if 0:
        aimChanneldebug1 = audio[0]
        aimChanneldebug2 = audio[1]
        plt.figure(3)
        plt.subplot(311)
        plt.plot(aimChanneldebug1)
        plt.subplot(312)
        plt.plot(aimChanneldebug2)
        plt.subplot(313)
        plt.plot(chirp)
        plt.pause(0.05)

        plt.figure(5)
        plt.subplot(211)
        f, t, Zxx = signal.stft(aimChanneldebug1, fs=sampleRate, nperseg=256, noverlap=128, nfft=256)
        plt.pcolor(np.log2(np.abs(Zxx)+1e-8), cmap='jet')
        plt.subplot(212)
        f, t, Zxx = signal.stft(aimChanneldebug2, fs=sampleRate, nperseg=256, noverlap=128, nfft=256)
        plt.pcolor(np.log2(np.abs(Zxx)+1e-8), cmap='jet')
        plt.pause(0.05)

    # #print(interval)
    # aimChannel = audio[0]
    aimChannel = audio[mics[0]]

    bandpassSig, corr = BaseUtils.Xcorr(aimChannel, chirp, lowFreq, highFreq, sampleRate)
    peaks, _ = signal.find_peaks(abs(corr), height=1, distance=interval)
    print(peaks)
    curAx = ax.flat[0]
    # corr
    axisList, upCorrEnvYList, xPeaks = XcorrCompute(curAx, peaks, corr, mics[0], interval, lab=fileLabel, sampleRate=sampleRate)
    FindPeaks(axisList, upCorrEnvYList, isTop, True)

    # Channel 2
    aimChannel1 = audio[mics[1]]
    bandpassSig1, corr1 = BaseUtils.Xcorr(aimChannel1, chirp, lowFreq, highFreq, sampleRate)
    peaks1, _ = signal.find_peaks(abs(corr1), height=0.1, distance=interval)
    curAx = ax.flat[1]
    # corr
    axisList, upCorrEnvYList, _ = XcorrCompute(curAx, peaks, corr1, mics[1], interval, xPeaks=xPeaks, lab=fileLabel) #
    FindPeaks(axisList, upCorrEnvYList, isTop, False)

    fig.suptitle(str(fileLabels[0])+"-"+speaker, x=0.2)

    plt.grid()
    plt.pause(0.05)




def CalDistance(delays: np.ndarray):
    """
    """
    delay = 0
    for i in range(50):
        delay+=delays[0]
    
    return delay/50*340


def CalWindSpeed():
    """
    """
    oneSide = peakResults[0]
    otherSide = peakResults[0] # peakResults[1]
    print(len(oneSide), len(otherSide))

    d1 = CalDistance(oneSide)
    d2 = CalDistance(otherSide)

    vsList = []
    wyList = []
    for i in range(len(oneSide)):
        vs = (d1/oneSide[i]+d2/otherSide[i])/2
        wy = abs(d1/oneSide[i]-d2/otherSide[i])/2
        vsList.append(vs)
        wyList.append(wy)

    xLast = vsList[0]
    pLast = vsList[0]
    vsListFiltered = np.zeros(len(vsList))
    for i in range(len(vsListFiltered)):
        pred, pLast, xLast= BaseUtils.Kalman(vsList[i], xLast, pLast, 1, 10) #
        vsListFiltered[i]=pred

    vsListFiltered = BaseUtils.Smooth(vsListFiltered, 3)

    xLast = wyList[0]
    pLast = wyList[0]
    wyListFiltered = np.zeros(len(wyList))
    for i in range(len(wyListFiltered)):
        pred, pLast, xLast= BaseUtils.Kalman(wyList[i], xLast, pLast, 0.5, 10) #
        wyListFiltered[i]=pred

    wyListFiltered = BaseUtils.Smooth(wyListFiltered, 3)

    plt.figure()
    #plt.plot(wyList, label="Speed of Wind")
    plt.plot(wyListFiltered[0:750]/2.5, linewidth = 3)
    plt.ylabel("风速值(m/s)", fontsize=30)
    plt.xlabel("周期下标", fontsize=30)
    plt.yticks(fontsize=25) 
    plt.xticks(fontsize=25)
    # plt.figure()
    # plt.plot(vsList, label="Speed of Sound")
    # plt.plot(vsListFiltered, label="Filtered Speed of Sound")
    # plt.legend()


if __name__ =="__main__":

    lowFreq = 18000       #
    highFreq = 23000     #
    sampleRate = 48000
    tSig = 5 / 1000     # 10  30
    tBlank = 5 * tSig

    interval = int((tSig + tBlank) * sampleRate)
    chirp = BaseUtils.UpChirp(lowFreq, highFreq, sampleRate, tSig, hann=False)

    # fig, ax = PlotEmpty(chirp, lowFreq, highFreq, sampleRate, tSig, interval)

    basepath = "H:\\UltraSound_Dataset"
    fileLabels = ['ChirpSignal_1700652652LR']
    for fileLabel in fileLabels:
        fileName = os.path.join(basepath, f"{fileLabel}.wav")
        #fileName = f"anemometerData/20221223/note9/xk-23k/hann/48k/sushe/{fileLabel}.wav"

        print(fileName)

        audio = AudioClass.load(filename=fileName)
        print(audio.sampleRate)

        # BottomSpeaker
        ProcessBySpeaker(audio, chirp, lowFreq, highFreq, sampleRate, interval, fileLabel)
        # TopSpeaker
        # ProcessBySpeaker(audio, chirp, lowFreq, highFreq, sampleRate, interval, fileLabel, isTop = True)

        print("done")
    CalWindSpeed()
    plt.tight_layout()
    plt.show()


