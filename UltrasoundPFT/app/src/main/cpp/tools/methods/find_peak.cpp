#include <stdio.h>
#include <float.h>
#include <android/log.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include <unordered_map>

#include "find_peak.h"
#include "../common/Utils.h"
#include "../../WaveSignalProcessing.h"

extern const char* LOGTAG;

double GetMax(const std::vector<double>& data){
    auto max_it = std::max_element(data.begin(), data.end());

    if (max_it != data.end()) {
        return *max_it;
    } else {
        __android_log_print(ANDROID_LOG_ERROR, LOGTAG, "input data is error");
        return -1;
    }
}

double GetMin(const std::vector<double>& data){
    auto min_it = std::min_element(data.begin(), data.end());

    if (min_it != data.end()) {
        return *min_it;
    } else {
        __android_log_print(ANDROID_LOG_ERROR, LOGTAG, "input data is error");
        return -1;
    }
}

std::vector<int32_t> findPeaks(const std::vector<double>& x,
                                double threshold, int indexDistance) {
    std::vector<int32_t> peaks;
    bool is_peak = false;

    for (int i = 1; i < x.size() - 1; ++i) {
        if (threshold >= 0 && x[i] < threshold) {
            continue;
        }

        if (x[i] > x[i - 1] && x[i] > x[i + 1]) {
            // Check if peak is at least indexDistance apart from other peaks
            bool isFarEnough = true;
            for (size_t j = 0; j < peaks.size(); ++j) {
                if (std::abs(i - peaks[j]) < indexDistance) {
                    isFarEnough = false;
                    break;
                }
            }
            if (isFarEnough) {
                peaks.push_back(i);
            }
        } else {
            is_peak = false;
        }
    }

    return peaks;
}

std::vector<int32_t> findTwoMaxPeaks(const std::vector<double>& x) {
    std::vector<int32_t> maxIndices(2, -1); // 初始化为无效索引值

    // 在信号中查找最大的两个峰值
    double max1 = std::numeric_limits<double>::min();
    double max2 = std::numeric_limits<double>::min();

    for (int i = 0; i < x.size(); ++i) {
        if (x[i] > max1) {
            max2 = max1;
            max1 = x[i];
            maxIndices[1] = maxIndices[0];
            maxIndices[0] = i;
        } else if (x[i] > max2) {
            max2 = x[i];
            maxIndices[1] = i;
        }
    }

    return maxIndices;
}


std::vector<int32_t> diffSignal(const std::vector<int32_t>& signal) {
    std::vector<int32_t> diff_signal;
    for (size_t i = 1; i < signal.size(); ++i) {
        diff_signal.push_back(signal[i] - signal[i - 1]);
    }
    return diff_signal;
}

double meanSignal(const std::vector<double>& signal) {
    double sum = 0.0;
    for (double value : signal) {
        sum += value;
    }
    return sum / signal.size();
}

int32_t meanSignal(const std::vector<int32_t>& signal) {
    int32_t sum = 0;
    for (int32_t value : signal) {
        sum += value;
    }
    return sum / signal.size();
}


int32_t modeSignal(const std::vector<int32_t>& nums) {
        std::unordered_map<int32_t, int> freq_map;
        for (int32_t num : nums) {
            freq_map[num]++;
        }

        int32_t mode_value = 0;
        int max_freq = 0;
        for (const auto& pair : freq_map) {
            if (pair.second > max_freq) {
                max_freq = pair.second;
                mode_value = pair.first;
            }
        }

        return mode_value;
}


