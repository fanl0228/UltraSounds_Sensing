import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks

def detect_peaks(signal, threshold):
    # 使用find_peaks函数检测峰值
    peaks, _ = find_peaks(signal, height=threshold)

    # 获取每个峰值的起始和截止点
    peak_ranges = []
    for peak in peaks:
        # 查找峰值左边的起始点
        start_point = np.argmax(signal[:peak] < threshold)

        # 查找峰值右边的截止点
        end_point = np.argmax(signal[peak:] < threshold) + peak

        peak_ranges.append((start_point, end_point))

    return peak_ranges

# 示例数据
time = np.linspace(0, 10, 1000)
signal = np.sin(time) + 0.5 * np.sin(5 * time)  # 一个包含峰值的示例信号

# 设置峰值检测的阈值
threshold = 0.1

# 检测峰值并获取起始和截止点
peak_ranges = detect_peaks(signal, threshold)

# 打印结果
print("峰值检测的起始和截止点:")
for i, (start, end) in enumerate(peak_ranges):
    print(f"峰值{i + 1}: 起始点={start}, 截止点={end}")

plt.figure()
plt.plot(signal)
# plt.scatter(peak_ranges, signal[peak_ranges], c='r')
plt.pause(0.05)

print("done")

