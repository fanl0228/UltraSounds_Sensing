//
// Created by Dell on 2024/2/24.
//

#ifndef ULTRASOUNDPFT_UTILS_H
#define ULTRASOUNDPFT_UTILS_H


#include <complex>
#include "../../WaveSignalProcessing.h"


#define DEBUG   1

/**
 * Log print functions
 * */
char* mergeStrings(const char* msg, const char* additional);

std::string replaceSubstring(const std::string& str, const std::string& from, const std::string& to);

std::string strscpString(const char* filePath, const char* str);

void logVector_uint8(const char * msg, const std::vector<uint8_t>& vec);

void logVector_int8(const char * msg, const std::vector<int8_t>& vec);

void logVector_int16(const char * msg, const std::vector<int16_t>& vec);

void logVector_int32(const char * msg, const std::vector<int32_t>& vec);

void logVectorReverse_int8(const char * msg, const std::vector<int8_t>& vec);

void logVectorReverse_int16(const char * msg, const std::vector<int16_t>& vec);

//void logDoubleArray(const char * msg, kiss_fft_cpx* array, size_t size);

void logVector_double(const char * msg, const std::vector<double>& vec);


/**
 * convert functions
 * */
int isDirectoryExists(const char* path);

char* extractDirectory(const char* path);

std::vector<int8_t> deepCopyVector(const std::vector<int8_t>& original);

std::vector<int16_t> deepCopyVector(const std::vector<int16_t>& original);

std::vector<double> deepCopyVector(const std::vector<double>& original);

double* vectorInt8ToDouble(const std::vector<int8_t>& inputVector);

double* vectorInt16ToDouble(const std::vector<int16_t>& inputVector);

double* vectorDoubleToDouble(const std::vector<double>& inputVector);

std::vector<double> doubleTovector(const double* input, int size);

std::vector<double> VInt16ToVDouble(const std::vector<int16_t>& input);

std::vector<double> normalizeToMinusOneToOne(const std::vector<double>& inputVector);

std::vector<double> normalizeToMinusOneToOne(const std::vector<int16_t>& inputVector);
/**
 * Save data functions
 * */
int saveChDataToCSV(WaveSignalStruct& waveSignal);

int saveProcessingDataToCSV(WaveSignalStruct& waveSignal);

int saveVectorDataToCSV(const char *filename, std::vector<double>& inData);

int saveVectorDataToCSV(const char *filename, std::vector<int16_t>& inData);

int saveVectorDataTo2D(const char* filename,
                       const std::vector<std::vector<double> >& inData);


void saveComplexMatrixToCSV(const std::string& filename,
                            const std::vector<std::vector<std::complex<double>>>& matrix);





#endif //ULTRASOUNDPFT_UTILS_H
