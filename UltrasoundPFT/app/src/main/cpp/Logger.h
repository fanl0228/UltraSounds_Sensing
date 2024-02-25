//
// Created by Administrator on 2021/12/28.
//

#ifndef EMLIVESDK_LOGGER_H
#define EMLIVESDK_LOGGER_H
#include <stdarg.h>
#define Logv( msg , ...) log(Log_Level_V , msg , ##__VA_ARGS__);
#define Logd( msg , ...) log(Log_Level_D , msg , ##__VA_ARGS__);
#define Logi( msg , ...) log(Log_Level_I , msg , ##__VA_ARGS__);
#define Logw( msg , ...) log(Log_Level_W , msg , ##__VA_ARGS__);
#define Loge( msg , ...) log(Log_Level_E , msg , ##__VA_ARGS__);
typedef void (*androidLog)(int level , char * msg) ;

enum LogLevel{
    Log_Level_V = 1,
    Log_Level_D,
    Log_Level_I,
    Log_Level_W,
    Log_Level_E,
};

class Logger{

public:
    /**
     * 日志回调
     *
     * 非线程安全的，所以一般在jni load的时候初始化
     *
     * @param level 日志等级
     * @param msg   日志内容
     * @param args  日志内容的占位参数值
     */
    static void log(int level , const char * msg , va_list args);
    /**
    * 初始化log回调
    * @param logCb 全局日志的回调函数
    * @param level 全局的日志等级
    */
    static void initLog(androidLog logCb , LogLevel level);

private :
    static androidLog g_LogCb;
    static LogLevel g_LogLevel;
};

void log(int level , const char * msg , ...);
#endif //EMLIVESDK_LOGGER_H
