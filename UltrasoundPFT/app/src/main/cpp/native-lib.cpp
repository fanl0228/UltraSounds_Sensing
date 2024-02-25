#include <jni.h>
#include <string>
#include <stdio.h>
#include <vector>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <android/log.h>

#include "Logger.h"
#include "chirp_generator.h"
#include "WaveSignalProcessing.h"
#include "Utils.h"

#ifdef __cplusplus
extern "C" {
#endif

//extern const char* LOGTAG = "native processing--->";
//extern const char* LOGTAG;

const char* BASED_DIR = "/storage/emulated/0/Music/";
const char* SUFFIX = "LR.wav";


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

JNIEXPORT jintArray JNICALL
Java_com_example_ultrasoundpft_MainActivity_ProcessingFromJNI(
        JNIEnv* env,
        jobject jobj, /* this */
        jstring filename) {

    jsize size= env->GetStringUTFLength(filename);
    Logd("str length: %d", (int)size);

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
    std::vector<int8_t> TxChirpSignal = genPCM16MonoToneBytes(mChirpParams);
    std::vector<int8_t> RxLeftChSignal = {};
    std::vector<int8_t> RxRightChSignal = {};
    WaveSignalStruct mWaveSignal(&mHeader, TxChirpSignal,
                                 RxLeftChSignal, RxRightChSignal);

    // Read audio data
    int datasize = read_wav_file((char*)c_filename, mWaveSignal);
    if (datasize == -1){
        __android_log_print(ANDROID_LOG_ERROR, LOGTAG, "Read wave file failed...");
        env->ReleaseStringUTFChars(filename, c_filename);
        exit(-1);
    }

    // Processing
    int res = signalPreProcessing(mChirpParams, mWaveSignal);
    if (res == -1){
        __android_log_print(ANDROID_LOG_ERROR, LOGTAG, "signalPreProcessing failed...");
    }

    // Log vis
    // _logVector_int8(TxChirpSignal);
    logVector_int8("RxLeftChData", mWaveSignal.RxLeftChData);
    logVector_int8("RxRightChData", mWaveSignal.RxRightChData);















//    // 将C++向量中的数据复制到Java字节数组中
//    jclass cls_ArrayList = env->GetObjectClass(result_buffer);
//    jmethodID mID_size = env->GetMethodID(cls_ArrayList, "size", "()I");
//    jmethodID mID_get = env->GetMethodID(cls_ArrayList, "get", "(I)Ljava/lang/Object;");
//    jbyteArray array = env->NewByteArray(audioData_buffer.size());
//    env->SetByteArrayRegion(array, 0, audioData_buffer.size(), (jbyte *)&audioData_buffer[0]);


    jintArray result;
    result = env->NewIntArray(10);
    if (result == NULL) {
        return NULL; /* out of memory error thrown */
    }

    /* fill a temp structure to use to populate the java int array */
    int fill[10];
    for (int i = 0; i < 10; i++) {
        fill[i] = i * 2;
    }

    /* move from the temp structure to the java structure */
    env->SetIntArrayRegion(result, 0, 10, fill);

    // release path
    env->ReleaseStringUTFChars(filename, c_filename);
    return result;
}



#ifdef __cplusplus
}
#endif


