package com.example.ultrasoundpft.utils;

import android.Manifest;
import android.content.pm.PackageManager;
import android.media.AudioFormat;
import android.widget.Toast;

import androidx.appcompat.app.AppCompatActivity;
import androidx.core.app.ActivityCompat;
import androidx.core.content.ContextCompat;

import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.text.SimpleDateFormat;
import java.util.Date;

public class Common {

    public static void requestPermissons (AppCompatActivity activity) {
        String[]  permissions= new String[]{
                Manifest.permission.WRITE_EXTERNAL_STORAGE,
                Manifest.permission.READ_EXTERNAL_STORAGE,
                Manifest.permission.RECORD_AUDIO,
                Manifest.permission.ACCESS_COARSE_LOCATION,
                Manifest.permission.ACCESS_FINE_LOCATION,
        };
        if (ContextCompat.checkSelfPermission(activity, Manifest.permission.WRITE_EXTERNAL_STORAGE)!= PackageManager.PERMISSION_GRANTED){
            if (ActivityCompat.shouldShowRequestPermissionRationale(activity, Manifest.permission.WRITE_EXTERNAL_STORAGE)){
                Toast.makeText(activity, "User rejected permissons before xxxx", Toast.LENGTH_SHORT).show();
            }else {
                ActivityCompat.requestPermissions(activity, permissions, 1000);
                return;
            }
        }

    }

    public static void createDirIfNotExist(String path) {
        File dir = new File(path);
        System.err.println(dir.isDirectory() + " " + dir.isFile());
        System.err.println(dir.exists());
        if (dir.isDirectory() && !dir.exists()) {
            boolean res = dir.mkdirs();
            System.err.println(path + " is created (status: " + res + ")");
        }
    }

    public static int getChannelNumFromConfig(int channelConfig) {
        int channel = 2;
        switch (channelConfig) {
            case AudioFormat.CHANNEL_IN_MONO:   channel = 1; break;
            case AudioFormat.CHANNEL_IN_STEREO: channel = 2; break;
            default:
                System.err.println("(DemoRecorder: getChannel) UNSUPPORTED CHANNEL CONFIG");
                System.exit(1);
        }
        return channel;
    }

    public static void writeWavFileHeader(
            RandomAccessFile raf,
            long fileLength,
            long sampleRate,
            int channelCnt,
            long byteRate) throws IOException
    {
        long totalDataLen = fileLength + 36;
        byte[] header = new byte[44];
        header[0] = 'R'; // RIFF/WAVE header
        header[1] = 'I';
        header[2] = 'F';
        header[3] = 'F';
        header[4] = (byte) (totalDataLen & 0xff);
        header[5] = (byte) ((totalDataLen >> 8) & 0xff);
        header[6] = (byte) ((totalDataLen >> 16) & 0xff);
        header[7] = (byte) ((totalDataLen >> 24) & 0xff);
        header[8] = 'W';
        header[9] = 'A';
        header[10] = 'V';
        header[11] = 'E';
        header[12] = 'f'; // 'fmt ' chunk
        header[13] = 'm';
        header[14] = 't';
        header[15] = ' ';
        header[16] = 16; // 4 bytes: size of 'fmt ' chunk
        header[17] = 0;
        header[18] = 0;
        header[19] = 0;
        header[20] = 1; // format = 1
        header[21] = 0;
        header[22] = (byte) channelCnt;
        header[23] = 0;
        header[24] = (byte) (sampleRate & 0xff);
        header[25] = (byte) ((sampleRate >> 8) & 0xff);
        header[26] = (byte) ((sampleRate >> 16) & 0xff);
        header[27] = (byte) ((sampleRate >> 24) & 0xff);
        header[28] = (byte) (byteRate & 0xff);
        header[29] = (byte) ((byteRate >> 8) & 0xff);
        header[30] = (byte) ((byteRate >> 16) & 0xff);
        header[31] = (byte) ((byteRate >> 24) & 0xff);
        header[32] = (byte) (2 * 16 / 8); // block align
        header[33] = 0;
        header[34] = 16; // bits per sample
        header[35] = 0;
        header[36] = 'd';
        header[37] = 'a';
        header[38] = 't';
        header[39] = 'a';
        header[40] = (byte) (fileLength & 0xff);
        header[41] = (byte) ((fileLength >> 8) & 0xff);
        header[42] = (byte) ((fileLength >> 16) & 0xff);
        header[43] = (byte) ((fileLength >> 24) & 0xff);
        raf.seek(0);
        raf.write(header, 0, 44);
    }

}


