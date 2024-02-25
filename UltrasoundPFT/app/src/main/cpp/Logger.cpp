//
// Created by Dell on 2024/2/22.
//

#include "Logger.h"
#include <android/log.h>
#include <stdio.h>

LogLevel Logger::g_LogLevel = Log_Level_D;
androidLog Logger::g_LogCb = 0;

void Logger::initLog(androidLog logCb , LogLevel level){
    g_LogCb = logCb;
    g_LogLevel = level;
}
void Logger::log(int level , const char * msg , va_list args){
    if(level < g_LogLevel){
        return;
    }

    if(g_LogCb){
        char msgBuffer[800] = {0};
        vsnprintf(msgBuffer, sizeof(msgBuffer) - 1 , msg, args);
        g_LogCb(level , msgBuffer);
    }else{
        int androidLogLevel;
        switch (level) {
            case Log_Level_V:
                androidLogLevel = ANDROID_LOG_VERBOSE;
                break;
            case Log_Level_D:
                androidLogLevel = ANDROID_LOG_DEBUG;
                break;
            case Log_Level_I:
                androidLogLevel = ANDROID_LOG_INFO;
                break;
            case Log_Level_W:
                androidLogLevel = ANDROID_LOG_WARN;
                break;
            case Log_Level_E:
                androidLogLevel = ANDROID_LOG_ERROR;
                break;
        }
        __android_log_print(androidLogLevel, "native processing--->", msg , args);
    }
}

void log(int level, const char *msg, ...) {
    va_list args;
    va_start(args , msg);
    Logger::log(level , msg , args);
    va_end(args);
}
