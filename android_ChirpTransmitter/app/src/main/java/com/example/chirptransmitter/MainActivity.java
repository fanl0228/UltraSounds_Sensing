package com.example.chirptransmitter;

import androidx.appcompat.app.AppCompatActivity;
import android.os.Bundle;
import android.view.View;
import android.widget.Button;
import android.widget.TextView;
import android.util.Log;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

import com.example.chirptransmitter.utils.Common;

import com.orhanobut.logger.Logger;
import com.orhanobut.logger.AndroidLogAdapter;
import com.orhanobut.logger.DiskLogAdapter;
import com.orhanobut.logger.FormatStrategy;
import com.orhanobut.logger.Logger;
import com.orhanobut.logger.PrettyFormatStrategy;


/**
 * WIFI-ADB:
 * 1. Connect Android phone and host machine to same WiFi network.
 * 2. Connect Android phone to host machine using USB cable (to start with).
 * 3. Run <adb> adb tcpip 5555 </adb> from a command prompt
 * 4. Run <adb> adb shell "ip addr show wlan0 | grep -e wlan0$ | cut -d\  -f 6 | cut -d/ -f 1" </adb>  to obtain the phone's IP address
 * 5. Disconnect USB cable and run <adb> adb connect <ip_address>:5555 </adb>
 *
 * Mi_11: adb connect 192.168.2.186:5555
 * Mi_10: adb connect 192.168.2.154:5555
 */


public class MainActivity extends AppCompatActivity {

    public StereoRecorder stereoRecorder;
    public ChirpPlayer chirpPlayer;
    public Thread recorderThread;
    public Thread chirpPlayerThread;

    public Button ctrlBtn;
    public TextView logView;

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);

        // Logger adapter init
        FormatStrategy formatStrategy = PrettyFormatStrategy.newBuilder()
                .showThreadInfo(false)  // (Optional) Whether to show thread info or not. Default true
                .methodCount(0)         // (Optional) How many method line to show. Default 2
                .methodOffset(3)        // (Optional) Skips some method invokes in stack trace. Default 5
//        .logStrategy(customLog) // (Optional) Changes the log strategy to print out. Default LogCat
                .tag("ChirpTransmitter")   // (Optional) Custom tag for each log. Default PRETTY_LOGGER
                .build();
        Logger.addLogAdapter(new AndroidLogAdapter(formatStrategy));



        stereoRecorder = new StereoRecorder();
        chirpPlayer = new ChirpPlayer(stereoRecorder);
        recorderThread = new Thread(stereoRecorder);
        chirpPlayerThread = new Thread(chirpPlayer);


        setContentView(R.layout.activity_main);
        Common.requestPermissons(this);
        ctrlBtn = (Button) findViewById(R.id.button_ctrl);
        logView = (TextView) findViewById(R.id.text_log);
    }

    @Override
    protected void onStart() {
        super.onStart();

        ctrlBtn.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View v) {
                try {
                    if (((String)v.getTag()).equals("0")) {
                        chirpPlayerThread.join();
                        recorderThread.join();

                        //current time
                        long totalMilliSeconds = System.currentTimeMillis();
                        long totalSeconds = totalMilliSeconds / 1000;

                        long currentSecond = totalSeconds % 60;

                        long totalMinutes = totalSeconds / 60;
                        long currentMinute = totalMinutes % 60;

                        long totalHour = totalMinutes / 60;
                        long currentHour = totalHour % 24 + 8;

                        String currTime = currentHour + "_" + currentMinute + "_" + currentSecond;

                        Logger.d(currTime);

                        stereoRecorder.setRecFileName("ChirpSignal_" + currTime);

                        recorderThread = new Thread(stereoRecorder);
                        chirpPlayerThread = new Thread(chirpPlayer);

                        recorderThread.start();
                        chirpPlayerThread.start();

                        v.setTag("1");
                        ctrlBtn.setText("END");
                    } else {
                        // recorder is stopped from player
                        chirpPlayer.setPlaying(false);
                        while (stereoRecorder.getRecording()) {Thread.yield();}
                        v.setTag("0");
                        ctrlBtn.setText("START");
                    }
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        });
    }


    @Override
    protected void onResume() {
        super.onResume();

    }

    @Override
    protected void onPause() {
        super.onPause();

    }

    @Override
    protected void onStop() {
        super.onStop();

    }

    @Override
    protected void onDestroy() {
        super.onDestroy();

        Logger.clearLogAdapters();
    }
}