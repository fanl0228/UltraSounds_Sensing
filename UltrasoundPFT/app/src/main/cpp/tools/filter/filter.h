//
// Created by Dell on 2024/2/26.
//

#ifndef ULTRASOUNDPFT_FILTER_H
#define ULTRASOUNDPFT_FILTER_H



// define Second-Order Section（SOS）struct
typedef struct {
    double b[3]; // Molecular coefficient
    double a[3]; // Denominator coefficient
} SOSSection;


typedef enum {
    Butterworth     =   0,
    Chebyshev1      =   1,
    Chebyshev2      =   2,
    Elliptic        =   3,
    Bessel          =   4,
}FilterType;

// 定义切比雪夫二型滤波器的 SOS 参数
/**
 * Use Matlab generate sos parameters
 * */





int filter_signal(double* inData, int x_len, SOSSection* sos,
                      int n_sections, FilterType type);


#endif //ULTRASOUNDPFT_FILTER_H
