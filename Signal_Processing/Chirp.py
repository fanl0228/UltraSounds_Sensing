import math, os
import numpy as np


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
    """ generate upchirp
    """
    length = int(rate * T)
    t = np.zeros(length)
    chirp = np.zeros(length)
    k = (high - low) / T
    for n in range(length):
        t[n] = (float(n) - offset) / rate
        chirp[n] = (math.sin(2 * math.pi * low * t[n] + math.pi * k * t[n] * t[n] + phase))  # *32767

    ## hamming
    # chirp = np.multiply(chirp, signal.hamming(length))

    ## hanning
    # chirp = np.multiply(chirp, signal.hanning(length))
    if hann:
        chirp = HanningWindow(chirp, 0, length)
    return chirp
