# -*- encoding: utf-8 -*-
"""
    py说明：音频类
    xx : public 变量，类对象（即类的实例化）访问
    _xx : protected 变量，类对象访问；
    __xx : private 变量，类对量不可直接访问，但有途径访问——“伪私有”

    包含方法：

"""

from datetime import datetime
from pathlib import Path
from typing import List, Union

import numpy as np
from scipy import signal
from scipy.io import wavfile


class AudioClass:

    ## public && protected:
    data: np.ndarray  # 数据
    channels: int  # 音轨数目
    sampleRate: int  # 采样率
    title: str  # 标题
    length: int  # 采样点数
    duration: float  # 时长
    sampleSequence: np.ndarray  # 采样序列
    timeSequence: np.ndarray  # 时间序列

    ## private:

    @property
    def channels(self) -> int:
        if len(self.data.shape) == 1:
            return 1
        else:
            return self.data.shape[1]

    @property
    def length(self) -> int:
        return self.data.shape[0]

    @property
    def duration(self) -> float:
        return self.length / self.sampleRate

    @property
    def sampleSequence(self) -> np.ndarray:
        return np.arange(0., self.length)

    @property
    def timeSequence(self) -> np.ndarray:
        return np.linspace(0., self.duration, self.data.shape[0])

    @staticmethod
    def load(filename: str,
                  title: Union[str, None] = None):

        file = Path(filename)
        if not file.exists():
            raise FileNotFoundError(filename)

        sampleRate, data = wavfile.read(filename)
        return AudioClass(data, sampleRate, title if title else file.name)

    def __init__(self,
                 data: np.ndarray,
                 sample_rate: int,
                 title: Union[str, None] = None):
        # 必须统一将数据类型转换成float64处理
        if data.dtype == np.int16:
            d = data / (2**15 - 1)
        elif data.dtype == np.int32:
            d = data / (2**31 - 1)
        elif data.dtype == np.float32:
            d = data.astype(np.float64)
        elif data.dtype == np.float64:
            d = np.copy(data) # 将数据copy到专有内存中，防止不可写入的问题
        else:
            raise Exception(f"不支持的数据类型：{data.dtype}")
        self.data = d
        self.sampleRate = sample_rate
        self.title = title if title is not None else datetime.now().strftime("%Y-%m-%d_%H:%M:%S.wav")

    def __str__(self) -> str:
        return f'[{self.title}] {self.channels}ch {self.duration}s({self.length}) {self.sampleRate}Hz {self.data.dtype}'

    def __repr__(self) -> str:
        return self.__str__()

    def __len__(self) -> int:
        return self.length

    def __getitem__(self, key) -> np.ndarray:
        return self.getChannel(key)
    
    def getChannel(self, index) -> np.ndarray:
        if index < 0 or index > self.channels-1:
            raise IndexError("index cannot be < 0 or > channels - 1")
        if self.channels == 1:
            return self.data
        else:
            return self.data[:, index]