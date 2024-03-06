//
// Created by Dell on 2024/2/26.
//

#ifndef ULTRASOUNDPFT_FILTER_H
#define ULTRASOUNDPFT_FILTER_H

#include <stdio.h>
#include <stdlib.h>
#include <vector>

// define Second-Order Section（SOS）struct
typedef struct {
    double b[3]; // Molecular coefficient
    double a[3]; // Denominator coefficient

    mutable double xm1 = 0; // 上一个输入延迟状态
    mutable double xm2 = 0; // 上两个输入延迟状态
    mutable double ym1 = 0; // 上一个输出延迟状态
    mutable double ym2 = 0; // 上两个输出延迟状态

} SOSSection;


typedef enum {
    Butterworth     =   0,
    Chebyshev1      =   1,
    Chebyshev2      =   2,
    Elliptic        =   3,
    Bessel          =   4,
}FilterStructType;

typedef enum {
    Low_Pass       =   0,
    High_Pass      =   1,
    Band_Pass      =   2,
}FilterNameType;

// 定义切比雪夫二型滤波器的 SOS 参数
/**
 * Use Matlab generate sos parameters
 * */
extern const SOSSection Chebyshev2HPFSOS_Order25_Wn10k_Fs48k[13];
extern const SOSSection Chebyshev2LPFSOS_Order25_Wn10k_Fs48k[13];
extern const SOSSection Chebyshev2LPFSOS_Order16_Wn6k_Fs48k[8];

extern const SOSSection EllipticHPFSOS_Order25_Wn10k_Fs48k[13];
extern const SOSSection EllipticLPFSOS_Order25_Wn10k_Fs48k[13];

class Filter{

public:

    FilterStructType mStructType;

    FilterNameType mNameType;

    Filter(FilterStructType sType, FilterNameType nType);

    Filter(){};

    ~Filter(){};

    int applyFilter(const std::vector<double>& inData, std::vector<double>& outData);


private:
    int _Chebyshev2Filter(const std::vector<double>& inData, std::vector<double>& outData,
                       const SOSSection* sos, int n_sections);

    int _EllipticFilter(const std::vector<double>& inData, std::vector<double>& outData,
                     const SOSSection* sos, int n_sections);


};


#endif //ULTRASOUNDPFT_FILTER_H

