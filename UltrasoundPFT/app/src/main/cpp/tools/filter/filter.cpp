//
// Created by Dell on 2024/2/26.
//

#include <stdio.h>
#include <stdlib.h>
#include <android/log.h>
#include <cstring>

#include "filter.h"

#include "../common/Logger.h"

#define PI 3.14159265358979323846


int filter_signal(double* data, int x_len, SOSSection* sos, int n_sections, FilterType type) {
    // 初始化滤波结果数组
    double* filtered_signal = (double*)malloc(x_len * sizeof(double));
    if (filtered_signal == nullptr) {
        fprintf(stderr, "Memory allocation failed\n");
        return -1;
    }
    // initialization
    memset(filtered_signal, 0, x_len * sizeof(double));

    __android_log_print(ANDROID_LOG_DEBUG, "native processing--->",
                        "n_sections: %d, x_len:%d", n_sections, x_len);

    // 信号滤波
    if (type == Chebyshev2){

        Logd("Chebyshev2 Filter")

        // 遍历输入信号
        for (int i = 0; i < x_len; ++i) {
            // 计算当前输出样本
            double y = data[i];

            for (int j = 0; j < n_sections; ++j) {
                // 对每个 SOS 段进行滤波
                double x1 = (i >= 1) ? filtered_signal[i - 1] : 0.0;
                double x2 = (i >= 2) ? filtered_signal[i - 2] : 0.0;
                double y1 = (i >= 1) ? filtered_signal[i - 1] : 0.0;
                double y2 = (i >= 2) ? filtered_signal[i - 2] : 0.0;

                // 计算当前段的输出
                double b0 = sos[j].b[0];
                double b1 = sos[j].b[1];
                double b2 = sos[j].b[2];
                double a1 = sos[j].a[1];
                double a2 = sos[j].a[2];

                y = (b0 * data[i] + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2) / data[0];
            }

            // 将当前输出样本存入结果数组
            filtered_signal[i] = y;

//            __android_log_print(ANDROID_LOG_DEBUG, "native processing--->",
//                                "filtered_signal[%d] = %f", i, y);
        }
    }
    else if (type == Elliptic){
        ;
    }

    for (int n = 0; n < x_len; ++n) {
        data[n] = filtered_signal[n];
    }

    free(filtered_signal);

    // 返回滤波结果数组
    return 0;
}

