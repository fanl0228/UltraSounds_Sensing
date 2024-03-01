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

void logVector_int16(const char * msg, const std::vector<int16_t>& vec) {
    std::ostringstream oss;
    oss << "Vector content: [";
    for (const auto& elem : vec) {
        oss << static_cast<int16_t>(elem) << ", ";
    }
    oss << "]";
    __android_log_print(ANDROID_LOG_DEBUG, LOGTAG,
                        "%s size: %d , %s",
                        msg, vec.size(), oss.str().c_str());
}

void logVector_int32(const char * msg, const std::vector<int32_t>& vec) {
    std::ostringstream oss;
    oss << "Vector content: [";
    for (const auto& elem : vec) {
        oss << static_cast<int32_t>(elem) << ", ";
    }
    oss << "]";
    __android_log_print(ANDROID_LOG_DEBUG, LOGTAG,
                        "%s size: %d , %s",
                        msg, vec.size(), oss.str().c_str());
}

void logVector_double(const char * msg, const std::vector<double>& vec) {
    std::ostringstream oss;
    oss << "Vector content: [";
    for (const auto& elem : vec) {
        oss << static_cast<double>(elem) << ", ";
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

void logVectorReverse_int16(const char * msg, const std::vector<int16_t>& vec) {
    std::ostringstream oss;
    oss << "Vector content: [";
    for (auto it = vec.rbegin(); it != vec.rend(); ++it) {
        oss << static_cast<int16_t>(*it) << ", ";
    }
    oss << "]";
    __android_log_print(ANDROID_LOG_DEBUG, LOGTAG,
                        "%s size: %d , %s",
                        msg, vec.size(), oss.str().c_str());
}

//void logDoubleArray(const char * msg, kiss_fft_cpx* array, size_t size) {
//    std::ostringstream oss;
//    oss << "Vector content: [";
//    for (size_t i = 0; i < size; ++i) {
//        oss << "(" << array[i].r << ", " << array[i].i << ")" << ";";
//    }
//    oss << "]";
//    __android_log_print(ANDROID_LOG_DEBUG, LOGTAG,
//                        "%s, size: %d , %s",
//                        msg, size, oss.str().c_str());
//}


// deep copy
std::vector<int8_t> deepCopyVector(const std::vector<int8_t>& original) {
    // Create a new std::vector using the copy constructor and
    // copy the contents of the original vector into the new vector
    std::vector<int8_t> copied(original);
    return copied;
}

// deep copy
std::vector<int16_t> deepCopyVector(const std::vector<int16_t>& original) {
    // Create a new std::vector using the copy constructor and
    // copy the contents of the original vector into the new vector
    std::vector<int16_t> copied(original);
    return copied;
}

// deep copy
std::vector<double> deepCopyVector(const std::vector<double>& original) {
    // Create a new std::vector using the copy constructor and
    // copy the contents of the original vector into the new vector
    std::vector<double> copied(original);
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

//  std::vector<int8_t> to double*
double* vectorInt16ToDouble(const std::vector<int16_t>& inputVector) {
    size_t size = inputVector.size();
    double* outputArray = new double[size];

    for (size_t i = 0; i < size; ++i) {
        outputArray[i] = static_cast<double>(inputVector[i]);
    }

    return outputArray;
}

//  std::vector<int8_t> to double*
double* vectorDoubleToDouble(const std::vector<double>& inputVector) {
    size_t size = inputVector.size();
    double* outputArray = new double[size];

    for (size_t i = 0; i < size; ++i) {
        outputArray[i] = static_cast<double>(inputVector[i]);
    }

    return outputArray;
}

std::vector<double> doubleTovector(const double* input, int size) {
    std::vector<double> outputVector(size,0.0);
    for (int i = 0; i < size; ++i) {
        outputVector[i] = input[i];
    }
    return outputVector;
}

std::vector<double> VInt16ToVDouble(const std::vector<int16_t>& input) {
    std::vector<double> output;
    output.reserve(input.size());
    for (int i = 0; i < input.size(); ++i) {
        output.push_back(static_cast<double>(input[i]));
    }
    return output;
}


int saveChDataToCSV(WaveSignalStruct& waveSignal){
    // Debug, 保存数据文件
    const char *filename = mergeStrings("/storage/emulated/0/Documents",
                                        strrchr(waveSignal.pFilename, '/'));
    std::string lStrFile(filename);
    std::string target = "LR.wav";
    size_t pos = lStrFile.find(target);
    if (pos != std::string::npos) {
        lStrFile.replace(pos, target.length(), "L.csv");
    }
    FILE *lFilename;
    lFilename = fopen(lStrFile.c_str(), "w");
    if (lFilename == nullptr) {
        __android_log_print(ANDROID_LOG_ERROR, LOGTAG,
                            "Unable to open file: %s", lStrFile.c_str());
        return -1;
    }

    if(waveSignal.RxLeftChData.empty()){
        __android_log_print(ANDROID_LOG_ERROR, LOGTAG,
                            "Save input waveSignal.RxLeftChData is empty: %s", filename);
        return -1;
    }

    for (int i = 0; i < waveSignal.RxLeftChData.size(); i++) {
        fprintf(lFilename, "%d\n", waveSignal.RxLeftChData[i]);
    }
    __android_log_print(ANDROID_LOG_DEBUG, LOGTAG,
                        "Write file success: %s", lStrFile.c_str());


    std::string rStrFile = lStrFile;
    std::string targetR = "L.csv";
    size_t posR = rStrFile.find(targetR);
    if (posR != std::string::npos) {
        rStrFile.replace(posR, targetR.length(), "R.csv");
    }
    FILE *rFilename;
    rFilename = fopen(rStrFile.c_str(), "w");
    if (rFilename == nullptr) {
        __android_log_print(ANDROID_LOG_ERROR, LOGTAG,
                            "Unable to open file: %s", rStrFile.c_str());
        return -1;
    }

    if(waveSignal.RxRightChData.empty()){
        __android_log_print(ANDROID_LOG_ERROR, LOGTAG,
                            "Save input waveSignal.RxRightChData is empty: %s", filename);
        return -1;
    }

    for (int i = 0; i < waveSignal.RxRightChData.size(); i++) {
        fprintf(rFilename, "%d\n", waveSignal.RxRightChData[i]);
    }
    __android_log_print(ANDROID_LOG_DEBUG, LOGTAG,
                        "Write file success: %s", rStrFile.c_str());

    std::string TxStrFile = lStrFile;
    std::string targetTx = "L.csv";
    size_t posTx = TxStrFile.find(targetTx);
    if (posTx != std::string::npos) {
        TxStrFile.replace(posTx, targetTx.length(), "Tx.csv");
    }
    FILE *pTxFilename;
    pTxFilename = fopen(TxStrFile.c_str(), "w");
    if (pTxFilename == nullptr) {
        __android_log_print(ANDROID_LOG_ERROR, LOGTAG,
                            "Unable to open file: %s", TxStrFile.c_str());
        return -1;
    }

    if(waveSignal.TxChirpSignal.empty()){
        __android_log_print(ANDROID_LOG_ERROR, LOGTAG,
                            "Save input waveSignal.TxChirpSignal is empty: %s", filename);
        return -1;
    }

    for (int i = 0; i < waveSignal.TxChirpSignal.size(); i++) {
        fprintf(pTxFilename, "%lf\n", waveSignal.TxChirpSignal[i]);
    }
    __android_log_print(ANDROID_LOG_DEBUG, LOGTAG,
                        "Write file success: %s", TxStrFile.c_str());

    fclose(lFilename);
    fclose(rFilename);
    fclose(pTxFilename);

    return 0;
}



int saveProcessingDataToCSV(WaveSignalStruct& waveSignal){
    // Debug, 保存数据文件
    const char *filename = mergeStrings("/storage/emulated/0/Documents",
                                        strrchr(waveSignal.pFilename, '/'));
    std::string lStrFile(filename);
    std::string target = "LR.wav";
    size_t pos = lStrFile.find(target);
    if (pos != std::string::npos) {
        lStrFile.replace(pos, target.length(), "_ProcessingDataL.csv");
    }
    FILE *lFilename;
    lFilename = fopen(lStrFile.c_str(), "w");
    if (lFilename == nullptr) {
        __android_log_print(ANDROID_LOG_ERROR, LOGTAG,
                            "Unable to open file: %s", lStrFile.c_str());
        return -1;
    }

    if(waveSignal.RxLeftProcessingData.empty()){
        __android_log_print(ANDROID_LOG_ERROR, LOGTAG,
                            "Save input waveSignal.RxLeftProcessingData is empty: %s", filename);
        return -1;
    }

    for (int i = 0; i < waveSignal.RxLeftProcessingData.size(); i++) {
        fprintf(lFilename, "%f\n", waveSignal.RxLeftProcessingData[i]);
    }
    __android_log_print(ANDROID_LOG_DEBUG, LOGTAG,
                        "Write file success: %s", lStrFile.c_str());


    std::string rStrFile = lStrFile;
    std::string targetR = "_ProcessingDataL.csv";
    size_t posR = rStrFile.find(targetR);
    if (posR != std::string::npos) {
        rStrFile.replace(posR, targetR.length(), "_ProcessingDataR.csv");
    }
    FILE *rFilename;
    rFilename = fopen(rStrFile.c_str(), "w");
    if (rFilename == nullptr) {
        __android_log_print(ANDROID_LOG_ERROR, LOGTAG,
                            "Unable to open file: %s", rStrFile.c_str());
        return -1;
    }

    if(waveSignal.RxRightProcessingData.empty()){
        __android_log_print(ANDROID_LOG_ERROR, LOGTAG,
                            "Save input waveSignal.RxRightProcessingData is empty: %s", filename);
        return -1;
    }

    for (int i = 0; i < waveSignal.RxRightProcessingData.size(); i++) {
        fprintf(rFilename, "%f\n", waveSignal.RxRightProcessingData[i]);
    }
    __android_log_print(ANDROID_LOG_DEBUG, LOGTAG,
                        "Write file success: %s", rStrFile.c_str());

    fclose(lFilename);
    fclose(rFilename);

    return 0;
}



int saveVectorDataToCSV(const char *filename,
                        std::vector<double>& inData){
    FILE *pFile = fopen(filename, "w+");
    if (pFile == nullptr) {
        __android_log_print(ANDROID_LOG_ERROR, LOGTAG,
                            "Unable to open file: %s", filename);
        return -1;
    }

    if(inData.empty()){
        __android_log_print(ANDROID_LOG_ERROR, LOGTAG,
                            "Save input inData is empty: %s", filename);
        return -1;
    }

    for (int i = 0; i < inData.size(); i++) {
        fprintf(pFile, "%lf\n", inData[i]);
    }
    __android_log_print(ANDROID_LOG_DEBUG, LOGTAG,
                        "Write file success: %s", filename);

    fclose(pFile);
    return 0;
}


int saveVectorDataToCSV(const char *filename,
                        std::vector<int16_t>& inData){
    FILE *pFile = fopen(filename, "w+");
    if (pFile == nullptr) {
        __android_log_print(ANDROID_LOG_ERROR, LOGTAG,
                            "Unable to open file: %s", filename);
        return -1;
    }

    if(inData.empty()){
        __android_log_print(ANDROID_LOG_ERROR, LOGTAG,
                            "Save input inData is empty: %s", filename);
        return -1;
    }

    for (int i = 0; i < inData.size(); i++) {
        fprintf(pFile, "%d\n", inData[i]);
    }
    __android_log_print(ANDROID_LOG_DEBUG, LOGTAG,
                        "Write file success: %s", filename);

    fclose(pFile);
    return 0;
}


int saveVectorDataTo2D(const char* filename,
                       const std::vector<std::vector<double> >& inData){
    FILE *pFile = fopen(filename, "w+");
    if (pFile == nullptr) {
        __android_log_print(ANDROID_LOG_ERROR, LOGTAG,
                            "Unable to open file: %s", filename);
        return -1;
    }

    if(inData.empty()){
        __android_log_print(ANDROID_LOG_ERROR, LOGTAG,
                            "Save input inData is empty.");
        return -1;
    }

    for (int r = 0; r < inData.size(); r++) {
        for (int c = 0; c < inData[r].size() -1; c++) {
            fprintf(pFile, "%lf,", inData[r][c]);
        }
        fprintf(pFile, "%lf\n", inData[r][inData[r].size() -1]);
    }

    __android_log_print(ANDROID_LOG_DEBUG, LOGTAG,
                        "Save filw success: %s.", filename);
    fclose(pFile);
    return 0;
}

