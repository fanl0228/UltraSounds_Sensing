import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as signal
from scipy import interpolate

def Filter(wave: np.ndarray, low: int, high: int, rate: int, type: str):
    """ Filter
        wave：original signal
        return: filting data
    """
    l = low / (rate / 2.0)
    h = high / (rate / 2.0)

    if type == "bandpass":
        b, a = signal.butter(10, [l, h], "bandpass")
    elif type == "highpass":
        b, a = signal.butter(5, l, 'highpass')
        wave = signal.filtfilt(b, a, wave)
    elif type == "lowpass":
        b, a = signal.butter(5, h, 'lowpass')
        wave = signal.filtfilt(b, a, wave)
    else:
        raise ValueError('type is not support...')
    return signal.filtfilt(b, a, wave)


def Xcorr(receivedSig: np.ndarray, chirp: np.ndarray, low: int, high: int, rate: int):
    """ corrolation
        received_sig_c  :   received signal
        dur             :   chirp duration
        return          :   result
    """

    bandpassSig = Filter(receivedSig, low, high, rate, type='bandpass')
    corr = signal.correlate(bandpassSig, chirp, mode='same', method='auto')
    return bandpassSig, corr



def UpInterpolate(data: np.ndarray, times: int):
    #
    x = np.arange(0, len(data))
    # xNew = np.linspace(x[0],x[-1],len(x)*times) #!!!大小会影响每个区域点的个数
    xNew = np.arange(x[0], x[-1], 1.0 / times)
    # print(len(x), len(xNew))
    func = interpolate.interp1d(x, data, kind='cubic')  # "cubic"
    func1 = interpolate.UnivariateSpline(x, data, s=0)  #
    dataNew = func(xNew)
    return dataNew


def Aligned_Chirps(peaks: np.ndarray, 
                   corr: np.ndarray,
                   xPeaks: np.ndarray = [], 
                   sampleRate:int = 48000, 
                   chirpSamples:int=240, 
                   interval_factor:float=0.9):
    """ Chirp data alignment

    """
    chirpSamples = int(chirpSamples / interval_factor)
    shipFrames = 5 #int(len(corr) / chirpSamples * 0.1)
    

    i = 0
    xOffsetsList = np.zeros(len(peaks))
    upCorrEnvYList = []
    while i<len(peaks):

        i += 1
        if not ((i>=shipFrames and i<=len(peaks)-shipFrames) ): # 经验, 去除前后两端一定范围的数据
            continue

        p = peaks[i]
        # corrSlice = corr[p - 50: p + maxPoints]

        corrSlice = corr[p - int(chirpSamples*0.25): p + int(chirpSamples*0.75)]
        if len(corrSlice) != chirpSamples:
            gap_samples = chirpSamples - len(corrSlice)
            corrSlice = np.append(corrSlice, np.zeros(gap_samples))

        # corrSlice_hilbert = signal.hilbert(corrSlice)
        # upCorrEnvY = abs(np.sqrt(corrSlice ** 2 + corrSlice_hilbert ** 2))
        # # LP-Filter
        # ellip_filter = signal.iirfilter(15, Wn=4000, rp=1, rs=60, btype='lowpass',
        #                                 analog=False, ftype='ellip', fs=48000, output='sos')
        # upCorrEnvY = signal.sosfilt(ellip_filter, upCorrEnvY)
        upCorrEnvY = corrSlice
        if 0:
            plt.figure(100)
            plt.plot(np.abs(corrSlice))
            plt.plot(upCorrEnvY, 'r')
            plt.pause(0.05)

        xOffset = -1
        if len(xPeaks) == 0:
            peaks1, _ = signal.find_peaks(upCorrEnvY, distance=int(chirpSamples*interval_factor))
            xOffset = peaks1[0]
        else:
            xOffset = int(xPeaks[i])

        upCorrEnvY = upCorrEnvY[xOffset: chirpSamples]

        if len(upCorrEnvY) < int(chirpSamples):
            gap_samples = int(chirpSamples - len(upCorrEnvY))
            upCorrEnvY = np.append(upCorrEnvY, np.zeros(gap_samples))
        elif len(upCorrEnvY) > int(chirpSamples):
            upCorrEnvY = upCorrEnvY[:int(chirpSamples)]

        upCorrEnvYList.append(upCorrEnvY)

        xOffsetsList[i] = xOffset

    return upCorrEnvYList, xOffsetsList


def one_dimensional_kalman_filter(measurements,
                                  initial_state_mean,
                                  initial_state_covariance,
                                  process_variance,
                                  measurement_variance):

    state_mean = initial_state_mean
    state_covariance = initial_state_covariance

    filtered_states = []

    for measurement in measurements:
        # Prediction step
        state_mean_prior = state_mean
        state_covariance_prior = state_covariance + process_variance

        # Update step
        kalman_gain = state_covariance_prior / (state_covariance_prior + measurement_variance)
        state_mean = state_mean_prior + kalman_gain * (measurement - state_mean_prior)
        state_covariance = (1 - kalman_gain) * state_covariance_prior

        filtered_states.append(state_mean)

    return np.array(filtered_states)




