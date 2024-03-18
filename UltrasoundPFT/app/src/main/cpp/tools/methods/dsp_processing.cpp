//
// Created by Dell on 2024/3/4.
//

#include <stdio.h>
#include <vector>
#include <complex>
#include <android/log.h>
#include "dsp_processing.h"
#include "../../WaveSignalProcessing.h"
#include "../fftw-3.3.10/api/fftw3.h"


int range_fft_for_RxRightChData(WaveSignalStruct& waveSignal, int N,
                                std::vector<std::complex<double>>& fft_result){

    int res = 0;
    int TxchirpSize = waveSignal.TxChirpSignal.size();
    int numChirps = floor(waveSignal.RxLeftChData.size() / TxchirpSize);
    std::vector<double> rSlice;
    std::vector<double> inData_real;
    std::vector<double> inData_imag;
    std::vector<std::complex<double>> range_fft_result;


    for (int i = 0; i < numChirps; ++i){

        rSlice.clear();
        rSlice.assign(waveSignal.RxRightChData.begin() + i * TxchirpSize,
                      waveSignal.RxRightChData.begin() + (i+1) * TxchirpSize);

//        res = Range_FFT_for_SignalChirp(rSlice, inData_imag, 128, true, fft_result);

    }

    std::vector<double> range_fft_real, range_fft_imag;


    return 0;

}


int Range_FFT_for_SignalChirp(std::vector<double>& inData_real,
                              std::vector<double>& inData_imag,
                              int N, bool isFFTshift,
                              std::vector<std::complex<double>>& fft_result){

    assert(inData_real.size() == inData_imag.size());
    // input check
    if(inData_real.empty()||inData_imag.empty()){
        __android_log_print(ANDROID_LOG_ERROR, LOGTAG,
                            "waveSignal.RxRightMixerData is empty.");
        return -1;
    }

    // bottom signal  Range-FFT
    fftw_complex *in, *out;
    fftw_plan p;

    int32_t inDataSize = inData_real.size();

    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    // Create an FFT plan for FFTW3
    p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    //Copy the input signal to the input array of FFTW3
    for (int i = 0; i <= inDataSize-N; i+=N) {
        for (int j =0; j<N;++j){
            in[i][0] = inData_real[i+j]; //real
            in[i][1] = inData_imag[i+j]; //imag
        }

        //Execute FFT
        fftw_execute(p);

        for (int i = 0; i < N; ++i) {
            fft_result.emplace_back(out[i][0], out[i][1]);    // copy FFT-out to waveSignal
        }
    }
    if (isFFTshift){
        // FFT shift
        fftShift(fft_result);
    }

    // release memory
    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);

    return 0;
}

int Range_FFT_for_RxRightMixerData(WaveSignalStruct& waveSignal, int N,
                                   std::vector<std::complex<double>>& fft_result) {
    // bottom signal  Range-FFT
    fftw_complex *in, *out;
    fftw_plan p;

    int32_t RxRightMixerDataSize = waveSignal.RxRightMixerData.real.size();
    // Check if input data is empty
    if (RxRightMixerDataSize == 0) {
        __android_log_print(ANDROID_LOG_ERROR, LOGTAG, "waveSignal.RxRightMixerData is empty.");
        return -1;
    }

    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    // Create an FFT plan for FFTW3
    p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    //Copy the input signal to the input array of FFTW3
    for (int i = 0; i <= RxRightMixerDataSize - N; i += N) {
        for (int j =0; j < N; ++j){
            in[j][0] = waveSignal.RxRightMixerData.real[i + j]; //real
            in[j][1] = waveSignal.RxRightMixerData.imag[i + j]; //imag
        }

        //Execute FFT
        fftw_execute(p);

        for (int j = 0; j < N; ++j) {
            fft_result.emplace_back(out[j][0], out[j][1]);    // copy FFT-out to waveSignal
        }
    }

    // release memory
    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);

    return 0;
}


void fftShift(std::vector<std::complex<double>>& spectrum) {
    int N = spectrum.size();
    int halfSize = N / 2;

    // 创建临时数组以保存前一半频谱值
    std::vector<std::complex<double>> temp(spectrum.begin(), spectrum.begin() + halfSize);

    // 将后一半频谱值移到前半部分
    for (int i = 0; i < halfSize; ++i) {
        spectrum[i] = spectrum[i + halfSize];
    }

    // 将前一半频谱值移到后半部分
    for (int i = 0; i < halfSize; ++i) {
        spectrum[i + halfSize] = temp[i];
    }
}

std::vector<std::vector<double>> reshape(const std::vector<double>& input, size_t cols) {
    std::vector<std::vector<double>> output;
    size_t numRows = input.size() / cols;
    output.reserve(numRows);

    auto inputIter = input.begin();
    for (size_t i = 0; i < numRows; ++i) {
        output.emplace_back(inputIter, inputIter + cols);
        inputIter += cols;
    }

    return output;
}

std::vector<std::vector<std::complex<double>>> reshape(const std::vector<std::complex<double>>& input,
                                                       size_t cols, bool isFFTshift) {
    std::vector<std::vector<std::complex<double>>> output;
    size_t numRows = input.size() / cols;
    output.reserve(numRows);

    auto inputIter = input.begin();
    for (size_t i = 0; i < numRows; ++i) {
        output.emplace_back(inputIter, inputIter + cols);
        inputIter += cols;
    }

    if(isFFTshift){
        for(int i =0; i < output.size(); ++i){
            fftShift(output[i]);
        }
    }

    return output;
}

std::vector<double> rowWiseMean(const std::vector<std::vector<std::complex<double>>>& input) {
    std::vector<double> means;
    for (const auto& row : input) {
        std::complex<double> sum(0.0, 0.0);
        for (const auto& value : row) {
            sum += value;
        }
        means.push_back(std::abs(sum) / row.size());
    }
    return means;
}


std::vector<double> columnWiseMean(const std::vector<std::vector<std::complex<double>>>& input) {
    std::vector<double> means(input[0].size(), 0.0); // init
    for (const auto& row : input) {
        for (size_t i = 0; i < row.size(); ++i) {
            means[i] += std::abs(row[i]);
        }
    }
    for (auto& mean : means) {
        mean /= input.size(); // cal mean
    }
    return means;
}

std::vector<double> unwrapPhase(const std::vector<std::complex<double>>& complexVector) {
    // Calculate phase angles
    std::vector<double> phaseAngles;
    for (const auto& complexNumber : complexVector) {
        phaseAngles.push_back(std::arg(complexNumber));
    }

    // Perform phase unwrapping
    double prevPhase = phaseAngles[0];
    for (size_t i = 1; i < phaseAngles.size(); ++i) {
        double diff = phaseAngles[i] - prevPhase;
        if (diff > M_PI) {
            while (diff > M_PI) {
                phaseAngles[i] -= 2 * M_PI;
                diff = phaseAngles[i] - prevPhase;
            }
        } else if (diff < -M_PI) {
            while (diff < -M_PI) {
                phaseAngles[i] += 2 * M_PI;
                diff = phaseAngles[i] - prevPhase;
            }
        }
        prevPhase = phaseAngles[i];
    }

    return phaseAngles;
}

std::vector<double> subtractVectors(const std::vector<double>& vec1, const std::vector<double>& vec2) {
    // 确保两个向量的大小相同
    if (vec1.size() != vec2.size()) {
        // 如果大小不同，可以选择抛出异常或者返回空向量
        return std::vector<double>();
    }

    // 创建一个新的向量来保存差值
    std::vector<double> diffVector;
    diffVector.reserve(vec1.size()); // 预留内存以提高性能

    // 逐元素计算差值
    for (size_t i = 0; i < vec1.size(); ++i) {
        diffVector.push_back(vec1[i] - vec2[i]);
    }

    return diffVector;
}


std::vector<double> vectorDiff(const std::vector<double>& signal) {
    std::vector<double> diff_signal;
    diff_signal.reserve(signal.size()-1);

    for (size_t i = 1; i < signal.size(); ++i) {
        diff_signal.push_back(signal[i] - signal[i - 1]);
    }
    return diff_signal;
}

std::vector<std::complex<double>> getColumn(const std::vector<std::vector<std::complex<double>>>& matrix, int column) {
    std::vector<std::complex<double>> result;
    for (const auto& row : matrix) {
        if (column < row.size()) {
            result.push_back(row[column]);
        } else {
            // If the column index exceeds the length of a row, insert a default value
            result.push_back(std::complex<double>(0.0, 0.0)); // Or any other default value you think is appropriate
        }
    }
    return result;
}

std::vector<double> hannWindow(int windowSize) {
    std::vector<double> window(windowSize);
    for (int i = 0; i < windowSize; ++i) {
        window[i] = 0.5 * (1 - std::cos(2 * M_PI * i / (windowSize - 1)));
    }
    return window;
}

std::vector<std::vector<std::complex<double>>> VectorSTFT(const std::vector<double>& signal, int windowSize, int hopSize) {
    int signalSize = signal.size();
    int numWindows = (signalSize - windowSize) / hopSize + 1;
    std::vector<std::vector<std::complex<double>>> stftResult(numWindows, std::vector<std::complex<double>>(windowSize / 2 + 1));

    fftw_complex *in, *out;
    fftw_plan p;

    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * windowSize);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * windowSize);

    p = fftw_plan_dft_1d(windowSize, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    std::vector<double> window = hannWindow(windowSize);

    for (int i = 0; i < numWindows; ++i) {
        // Window and load the signal into the input array
        for (int j = 0; j < windowSize; ++j) {
            int index = i * hopSize + j;
            in[j][0] = (index < signalSize) ? signal[index] * window[j] : 0.0; // real
            in[j][1] = 0.0; // imag
        }
        // execute fft
        fftw_execute(p);

        // save result
        for (int j = 0; j <= windowSize / 2; ++j) {
            stftResult[i][j] = std::complex<double>(out[j][0], out[j][1]);
        }
    }

    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);

    return stftResult;
}

std::vector<double> ColumnwiseAbsSum(const std::vector<std::vector<std::complex<double>>>& data) {
    int numRows = data.size();
    int numCols = data[0].size();
    std::vector<double> absSums(numCols, 0.0);

    // 对每一列求绝对值并求和
    for (int j = 0; j < numCols; ++j) {
        for (int i = 0; i < numRows; ++i) {
            absSums[j] += std::abs(data[i][j]);
        }
    }

    return absSums;
}

std::vector<double> rowwiseAbsSum(const std::vector<std::vector<std::complex<double>>>& data) {
    std::vector<double> absSums;
    for (const auto& row : data) {
        double sum = 0.0;
        for (const auto& complexNumber : row) {
            sum += std::abs(complexNumber);
        }
        absSums.push_back(sum);
    }
    return absSums;
}


int CalAirflowVelocity(std::vector<double>& data){
    for (int i = 0; i < data.size(); ++i){
        data[i] = -0.41 * (MODEL_PATH4 / ( data[i] / (2 * M_PI * START_FRE)
                                                + MODEL_PATH2 / VOICE_VELOCITY
                                                + (MODEL_PATH4-MODEL_PATH2)/VOICE_VELOCITY)
                              - VOICE_VELOCITY) * M_PI * MODEL_RADICAL * MODEL_RADICAL * 1000 * 1000;
    }

    return 0;
}

