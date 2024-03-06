#include <jni.h>
#include <string>
#include <stdio.h>
#include <vector>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <android/log.h>

#include "tools/common/Logger.h"
#include "tools/methods/chirp_generator.h"
#include "WaveSignalProcessing.h"
#include "tools/common/Utils.h"

#ifdef __cplusplus
extern "C" {
#endif

//extern const char* LOGTAG = "native processing--->";
//extern const char* LOGTAG;

const char* BASED_DIR = "/storage/emulated/0/Music/";
const char* SUFFIX = "LR.wav";
extern const char* LOGTAG;

JNIEXPORT jstring JNICALL
Java_com_example_ultrasoundpft_MainActivity_stringFromJNI(
        JNIEnv* env,
        jobject jobj, /* this */
        jstring str) {

    Logd("Native--> Processing");

    const char *c_str = env->GetStringUTFChars(str, NULL);
    if (c_str != NULL) {
        Logd(c_str);
        env->ReleaseStringUTFChars(str, c_str);
    }

    const char* hello = "Hello from C++";
    return env->NewStringUTF(hello);
}

JNIEXPORT jfloatArray JNICALL
Java_com_example_ultrasoundpft_MainActivity_ProcessingFromJNI(
        JNIEnv* env,
        jobject jobj, /* this */
        jstring filename) {

    jsize size= env->GetStringUTFLength(filename);

    __android_log_print(ANDROID_LOG_DEBUG, LOGTAG,
                        "Call native processing, filename length: %d", (int)size);

    // 将Java字符串转换为C字符串
    const char *c_filename = env->GetStringUTFChars(filename, NULL);
    if (c_filename == nullptr) {
        Loge("Path sting is NULL");
        env->ReleaseStringUTFChars(filename, c_filename);
        exit(-1);
    }
    Logd(c_filename);

    // read *.wav file
    /** Read Audio header message and Chirp to mixer */
    // init WaveSignal struct
    WavHeader mHeader;
    ChirpParameters mChirpParams;
    std::vector<double> TxChirpSignal = genPCM16MonoToneBytes(mChirpParams);
    std::vector<int16_t> RxLeftChSignal = {};
    std::vector<int16_t> RxRightChSignal = {};
    std::vector<double> RxLeftProcessingData = {};
    std::vector<double> RxRightProcessingData = {};
    ComplexVector RxLeftMixerData;
    ComplexVector RxRightMixerData;
    std::vector<double> airflowVelocity = {};
    WaveSignalStruct mWaveSignal((char*)c_filename, &mHeader, TxChirpSignal,
                                 RxLeftChSignal, RxRightChSignal,
                                 RxLeftProcessingData, RxRightProcessingData,
                                 RxLeftMixerData, RxRightMixerData, airflowVelocity);

    // Read audio data
    int datasize = read_wav_file((char*)c_filename, mWaveSignal);
    if (datasize == -1){
        __android_log_print(ANDROID_LOG_ERROR, LOGTAG, "Read wave file failed...");
        env->ReleaseStringUTFChars(filename, c_filename);
        exit(-1);
    }

#if DEBUG
    // Log vis
    Logd(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Original Data")
    logVector_double("Tx Chirp Signal Data ",mWaveSignal.TxChirpSignal);
    logVector_int16("Rx Left Channel Data ", mWaveSignal.RxLeftChData);
    logVector_int16("Rx Right Channel Data ", mWaveSignal.RxRightChData);
    Logd("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<")
#endif

    // input CH data check
    if(mWaveSignal.RxLeftChData.empty() || mWaveSignal.RxRightChData.empty()){
        __android_log_print(ANDROID_LOG_ERROR, LOGTAG, "Read wave data is empty...");
        exit(-1);
    }

    // Processing
    int res = signalPreProcessing(mChirpParams, mWaveSignal);
    if (res == -1){
        __android_log_print(ANDROID_LOG_ERROR, LOGTAG, "signalPreProcessing failed...");
    }


    airflowVelocity = mWaveSignal.AirflowVelocity;
    jfloatArray jAirflowVelocity = env->NewFloatArray(airflowVelocity.size());

    if (jAirflowVelocity == nullptr) {
        __android_log_print(ANDROID_LOG_ERROR, LOGTAG, "jAirflowVelocity Memory allocation failed");
        return nullptr; /* out of memory error thrown */
    }

    //Convert each element in std::vector<double> to jfloat and copy into jfloatArray
    std::vector<float> floatVec;
    floatVec.reserve(airflowVelocity.size());
    for (double value : airflowVelocity) {
        floatVec.push_back(static_cast<float>(value));
    }
    env->SetFloatArrayRegion(jAirflowVelocity, 0, airflowVelocity.size(), floatVec.data());


    return jAirflowVelocity;
}



#ifdef __cplusplus
}
#endif


