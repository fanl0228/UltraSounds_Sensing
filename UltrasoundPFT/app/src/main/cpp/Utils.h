//
// Created by Dell on 2024/2/24.
//

#ifndef ULTRASOUNDPFT_UTILS_H
#define ULTRASOUNDPFT_UTILS_H

#include "tools/kiss-fft/kiss_fft.h"
/**
 * Log print functions
 * */
const char* mergeStrings(const char* msg, const char* additional);

void logVector_uint8(const char * msg, const std::vector<uint8_t>& vec);

void logVector_int8(const char * msg, const std::vector<int8_t>& vec);

void logVectorReverse_int8(const char * msg, const std::vector<int8_t>& vec);

void logDoubleArray(const char * msg, kiss_fft_cpx* array, size_t size);


/**
 * tool functions
 * */
std::vector<int8_t> deepCopyVector(const std::vector<int8_t>& original);

double* vectorInt8ToDouble(const std::vector<int8_t>& inputVector);



void drawAndSaveCurve(const std::vector<double>& x, const std::vector<double>& y, const char* filename);


#endif //ULTRASOUNDPFT_UTILS_H
