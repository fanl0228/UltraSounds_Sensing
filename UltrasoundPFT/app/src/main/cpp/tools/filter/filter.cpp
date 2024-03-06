//
// Created by Dell on 2024/2/26.
//

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <android/log.h>
#include <cstring>

#include "filter.h"

#include "../common/Logger.h"
#include "../methods/find_peak.h"

#define PI 3.14159265358979323846



Filter::Filter(FilterStructType sType, FilterNameType nType) {

    mStructType = sType;

    mNameType = nType;
}

int Filter::applyFilter(const std::vector<double>& inData,
                        std::vector<double>& outData){
    // Check input
    if(inData.empty()){
        __android_log_print(ANDROID_LOG_DEBUG, "native processing--->",
                            "applyFilter inData is empty.");
        return -1;
    }

    //
    if(mNameType == FilterNameType::High_Pass){
        switch (mStructType) {
            case FilterStructType::Chebyshev2:
                _Chebyshev2Filter(inData, outData, Chebyshev2HPFSOS_Order25_Wn10k_Fs48k, 13);
                break;
            case FilterStructType::Elliptic:
                _EllipticFilter(inData, outData, EllipticHPFSOS_Order25_Wn10k_Fs48k, 13);
            default:
                break;
        }
    } else if(mNameType == FilterNameType::Low_Pass){
        switch (mStructType) {
            case FilterStructType::Chebyshev2:
                _Chebyshev2Filter(inData, outData, Chebyshev2LPFSOS_Order16_Wn6k_Fs48k, 8);
                break;
            case FilterStructType::Elliptic:
                _EllipticFilter(inData, outData, EllipticLPFSOS_Order25_Wn10k_Fs48k, 13);
            default:
                break;
        }
    }else if(mNameType == FilterNameType::Band_Pass){
        ;
    }else{

    }


    return 0;
}


int Filter::_Chebyshev2Filter(const std::vector<double>& inData,
                                 std::vector<double>& outData,
                           const SOSSection* sosParams, int n_sections) {

    std::vector<double> filteredSignal(inData.size());

    for (int i = 0; i < inData.size(); ++i) {
        // 计算当前输出样本
        double x = inData[i];
        double y = x;

        for (int j = 0; j < n_sections; ++j) {
            const double* b = sosParams[j].b;
            const double* a = sosParams[j].a;

            // 计算当前段的输出
            double y1 = b[0] * x + b[1] * sosParams[j].xm1 + b[2] * sosParams[j].xm2 - a[1] * sosParams[j].ym1 - a[2] * sosParams[j].ym2;
            y = y1 / a[0];

            sosParams[j].xm2 = sosParams[j].xm1;
            sosParams[j].xm1 = x;
            sosParams[j].ym2 = sosParams[j].ym1;
            sosParams[j].ym1 = y;

            x = y;
        }
        filteredSignal[i] = y;
//        __android_log_print(ANDROID_LOG_DEBUG, "native processing--->", "filtered_signal[%d] = %f", i, y);

    }

    for (int n = 0; n < filteredSignal.size(); ++n) {
        outData[n] = filteredSignal[n];
    }


    return 0;
}


int Filter::_EllipticFilter(const std::vector<double> &inData,
                         std::vector<double> &outData,
                         const SOSSection* sos, int n_sections) {
    return 0;
}

