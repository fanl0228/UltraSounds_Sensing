//
// Created by Dell on 2024/2/24.
//
#include <string>
#include <stdio.h>
#include <vector>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cstdint>
#include <android/log.h>

//#include <SkCanvas.h>
//#include <SkSurface.h>
//#include <SkImageEncoder.h>
//#include <SkStream.h>


#include "Utils.h"

extern const char* LOGTAG = "UltrasoundPFT: native --->";


const char* mergeStrings(const char* msg, const char* additional) {
    std::string mergedString = std::string(msg) + std::string(additional);
    return mergedString.c_str();
}


void logVector_uint8(const char * msg, const std::vector<uint8_t>& vec) {
    std::ostringstream oss;
    oss << "Vector content: [";
    for (const auto& elem : vec) {
        oss << static_cast<int>(elem) << ", ";
    }
    oss << "]";
    __android_log_print(ANDROID_LOG_DEBUG, LOGTAG,
                        "%s size: %d , %s",
                        msg, vec.size(), oss.str().c_str());
}


void logVector_int8(const char * msg, const std::vector<int8_t>& vec) {
    std::ostringstream oss;
    oss << "Vector content: [";
    for (const auto& elem : vec) {
        oss << static_cast<int>(elem) << ", ";
    }
    oss << "]";
    __android_log_print(ANDROID_LOG_DEBUG, LOGTAG,
                        "%s size: %d , %s",
                        msg, vec.size(), oss.str().c_str());
}


void logVectorReverse_int8(const char * msg, const std::vector<int8_t>& vec) {
    std::ostringstream oss;
    oss << "Vector content: [";
    for (auto it = vec.rbegin(); it != vec.rend(); ++it) {
        oss << static_cast<int>(*it) << ", ";
    }
    oss << "]";
    __android_log_print(ANDROID_LOG_DEBUG, LOGTAG,
                        "%s size: %d , %s",
                        msg, vec.size(), oss.str().c_str());
}


void logDoubleArray(const char * msg, kiss_fft_cpx* array, size_t size) {
    std::ostringstream oss;
    oss << "Vector content: [";
    for (size_t i = 0; i < size; ++i) {
        oss << "(" << array[i].r << ", " << array[i].i << ")" << ";";
    }
    oss << "]";
    __android_log_print(ANDROID_LOG_DEBUG, LOGTAG,
                        "%s, size: %d , %s",
                        msg, size, oss.str().c_str());
}


// deep copy
std::vector<int8_t> deepCopyVector(const std::vector<int8_t>& original) {
    // Create a new std::vector using the copy constructor and
    // copy the contents of the original vector into the new vector
    std::vector<int8_t> copied(original);
    return copied;
}


//  std::vector<int8_t> to double*
double* vectorInt8ToDouble(const std::vector<int8_t>& inputVector) {
    size_t size = inputVector.size();
    double* outputArray = new double[size];

    for (size_t i = 0; i < size; ++i) {
        outputArray[i] = static_cast<double>(inputVector[i]);
    }

    return outputArray;
}





void drawAndSaveCurve(const std::vector<double>& x, const std::vector<double>& y, const char* filename) {
//    // 创建 Skia 画布
//    sk_sp<SkSurface> surface = SkSurface::MakeRasterN32Premul(800, 600);
//    SkCanvas* canvas = surface->getCanvas();
//
//    // 清空画布
//    canvas->clear(SK_ColorWHITE);
//
//    // 绘制曲线
//    SkPaint paint;
//    paint.setColor(SK_ColorBLACK);
//    paint.setStrokeWidth(2);
//    for (size_t i = 0; i < x.size() - 1; ++i) {
//        canvas->drawLine(x[i], y[i], x[i + 1], y[i + 1], paint);
//    }
//
//    // 保存图像到文件
//    sk_sp<SkImage> image = surface->makeImageSnapshot();
//    sk_sp<SkData> png = image->encodeToData(SkEncodedImageFormat::kPNG, 100);
//    SkFILEWStream file(filename);
//    file.write(png->data(), png->size());

    ;
}


