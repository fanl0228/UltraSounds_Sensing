//
// Created by Dell on 2024/2/26.
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cstring>
#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>
#include <android/log.h>

//#include "../kiss-fft/kiss_fft.h"

#include "xcorr.h"
#include "../common/Utils.h"
#include "../../WaveSignalProcessing.h"

// third lib
#include "../fftw-3.3.10/api/fftw3.h"


void normalizeSignal(const std::vector<double>& signal,
                     std::vector<double>& normalizedSignal) {
    // 找到信号中的最大值和最小值
    auto min_it = std::min_element(signal.begin(), signal.end());
    auto max_it = std::max_element(signal.begin(), signal.end());

    // 计算最大值和最小值之间的范围
    int16_t min_val = *min_it;
    int16_t max_val = *max_it;
    double range = max_val - min_val;

    // 归一化信号
    normalizedSignal.resize(signal.size());
    for (size_t i = 0; i < signal.size(); ++i) {
        // 将每个样本从 [min_val, max_val] 映射到 [-1, 1] 范围内
        normalizedSignal[i] = 2 * (signal[i] - min_val) / static_cast<double>(range) - 1.0;
    }

#if 0
    __android_log_print(ANDROID_LOG_DEBUG, LOGTAG,
                        "min_val: %d, max_val: %d,max_val: %f ",
                        min_val, max_val, range);
    logVector_double("=====> normalizedSignal: ", normalizedSignal);

#endif

}


int crossCorrelation(const std::vector<double>& x,
                     const std::vector<double>& y,
                     std::vector<double>& res) {
    // check input
    assert(res.size() >= x.size());

    size_t len_x = x.size();
    size_t len_y = y.size();

    // Initialization
    std::vector<double> sig1(x.size());
    std::vector<double> sig2(y.size());

    normalizeSignal(x, sig1);
    normalizeSignal(y, sig2);

    // xCorr
    for (size_t i = 0; i < len_x; ++i) {
        double sum = 0.0;
        for (size_t j = 0; j < len_y; ++j) {
            if (i + j < len_x) {
                sum += sig1[i + j] * sig2[j];
            }
        }
        res[i] = sum;
    }

    return 0;
}

//
//int crossCorrelationFFT(const std::vector<double>& sig1,
//                        const std::vector<double>& sig2,
//                        std::vector<double>& result) {
//    // Check input
//    if (sig1.empty() || sig2.empty()) {
//        __android_log_print(ANDROID_LOG_ERROR, LOGTAG,
//                            "input data is empty");
//        return -1;
//    }
//    // Initialization
//    std::vector<double> sig11(sig1.size());
//    std::vector<double> sig22(sig2.size());
//
//    normalizeSignal(sig1, sig11);
//    normalizeSignal(sig2, sig22);
//
//    // calculate FFT length
//    int n = sig11.size() + sig22.size() - 1;
//
//    // allocate memory
//    kiss_fft_cpx* fft_sig11 = new kiss_fft_cpx[n];
//    kiss_fft_cpx* fft_sig22 = new kiss_fft_cpx[n];
//    kiss_fft_cpx* fft_result = new kiss_fft_cpx[n];
//
//    // Initialize FFT input signal
//    for (int i = 0; i < sig11.size(); ++i) {
//        fft_sig11[i].r = sig11[i];
//        fft_sig11[i].i = 0;
//    }
//    for (int i = 0; i < sig22.size(); ++i) {
//        fft_sig22[i].r = sig22[i];
//        fft_sig22[i].i = 0;
//    }
//
//    // 进行 FFT 变换
//    kiss_fft_cfg cfg = kiss_fft_alloc(n, false, 0, 0);
//    kiss_fft(cfg, fft_sig11, fft_sig11);
//    kiss_fft(cfg, fft_sig22, fft_sig22);
//
//    // 计算互相关
//    for (int i = 0; i < n; ++i) {
//        fft_result[i].r = fft_sig11[i].r * fft_sig22[i].r - fft_sig11[i].i * fft_sig22[i].i;
//        fft_result[i].i = fft_sig11[i].r * fft_sig22[i].i + fft_sig11[i].i * fft_sig22[i].r;
//    }
//
//    // 进行逆 FFT 变换
//    kiss_fft(cfg, fft_result, fft_result);
//
//    // 提取实部作为互相关结果
//    result.resize(n);
//    for (int i = 0; i < n; ++i) {
//        result[i] = fft_result[i].r;
//    }
//
//    // 释放内存
//    kiss_fft_free(cfg);
//    delete[] fft_sig11;
//    delete[] fft_sig22;
//    delete[] fft_result;
//
//    return 0;
//}
//



void HilbertTransform(fftw_complex *out, double *x, size_t len_x) {
    if (!out) {
        return;
    }

    size_t N = len_x;
    fftw_complex *in;
    fftw_plan p_fft, p_ifft;

    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    // always create plan before initializing input
    p_fft = fftw_plan_dft_1d(N, in, in, FFTW_FORWARD, FFTW_ESTIMATE);
    p_ifft = fftw_plan_dft_1d(N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);

    // Implementation according to: https://www.mathworks.com/help/signal/ref/hilbert.html
    // 1. calculate FFT(x) in-place
    for (int i = 0; i < N; i++) {
        in[i][0] = x[i];
        in[i][1] = 0.0;
    }
    fftw_execute(p_fft);
    fftw_destroy_plan(p_fft);

    // 2. multiply FFT(x) with H[i], then peform ifft
    for (int i = 1; i < ceil(N/2); i++) {
        in[i][0] *= 2.0;
        in[i][1] *= 2.0;
    }
    for (int i = ceil(N/2) + 1; i < N; i++) {
        in[i][0] = 0.0;
        in[i][1] = 0.0;
    }
    fftw_execute(p_ifft);
    fftw_destroy_plan(p_ifft);
    for (int i = 0; i < N; i++) {
        out[i][0] /= N;
        out[i][1] /= N;
    }

    // 3. release resources
    fftw_free(in);
}

void GetEnvelope(double *out, double *x, size_t len_x) {
    size_t N = len_x;
    fftw_complex *ht_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    HilbertTransform(ht_out, x, N);

    for (int i = 0; i < N; i++) {
        out[i] = sqrt(ht_out[i][0]*ht_out[i][0] + ht_out[i][1]*ht_out[i][1]);
    }
}

void calculateEnvelope(const std::vector<double>& signal,
                       std::vector<double>& envelope) {

    double* in = vectorDoubleToDouble(signal);
    double* out = new double[signal.size()];

    GetEnvelope(out, in, signal.size());

    envelope.assign(out, out + signal.size());

    delete[] in;
    delete[] out;
}

