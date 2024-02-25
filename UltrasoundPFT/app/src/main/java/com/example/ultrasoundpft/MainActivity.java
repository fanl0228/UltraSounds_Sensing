package com.example.ultrasoundpft;

import androidx.appcompat.app.AppCompatActivity;

import android.os.Bundle;
import android.view.View;
import android.widget.Button;
import android.widget.EditText;
import android.widget.TextView;

import com.example.ultrasoundpft.databinding.ActivityMainBinding;


import com.example.ultrasoundpft.utils.Common;
import com.example.ultrasoundpft.utils.Timestamp;
import com.orhanobut.logger.AndroidLogAdapter;
import com.orhanobut.logger.FormatStrategy;
import com.orhanobut.logger.Logger;
import com.orhanobut.logger.PrettyFormatStrategy;

import java.sql.Array;

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
    // Used to load the 'native-lib' library on application startup.
    static {
        System.loadLibrary("native-lib");
    }

    private ActivityMainBinding binding;

    /**
     * A native method that is implemented by the 'native-lib' native library,
     * which is packaged with this application.
     */
    public native String stringFromJNI(String str);
    public native int[] ProcessingFromJNI(String filename);


    private StereoRecorder stereoRecorder;
    private ChirpPlayer chirpPlayer;
    private Thread recorderThread;
    private Thread chirpPlayerThread;

    private EditText mUsername;
    private EditText mPassword;
    private Button mCtrlBtn;
    private Button mProcessBtn;
    private TextView mLogView;

    private String mBasePath = "/storage/emulated/0/Music/"; //const audio storage path
    private String mCurrTimeStamp;
    private String mInputUsername;
    private String mDataFilename = null;  // need init
    private int[] mProcessing_results;

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);

        // Logger adapter init
        FormatStrategy formatStrategy = PrettyFormatStrategy.newBuilder()
                .showThreadInfo(false)  // (Optional) Whether to show thread info or not. Default true
                .methodCount(0)         // (Optional) How many method line to show. Default 2
                .methodOffset(3)        // (Optional) Skips some method invokes in stack trace. Default 5
                .tag("ChirpTransmitter")   // (Optional) Custom tag for each log. Default PRETTY_LOGGER
                .build();
        Logger.addLogAdapter(new AndroidLogAdapter(formatStrategy));

        // FMCW signal transmitter
        stereoRecorder = new StereoRecorder();
        chirpPlayer = new ChirpPlayer(stereoRecorder);
        recorderThread = new Thread(stereoRecorder);
        chirpPlayerThread = new Thread(chirpPlayer);

        // UI handle
        setContentView(R.layout.activity_main);
        Common.requestPermissons(this);
        mUsername = (EditText) findViewById(R.id.username);
        mPassword = (EditText) findViewById(R.id.password);
        mCtrlBtn = (Button) findViewById(R.id.button_ctrl);
        mProcessBtn = (Button) findViewById(R.id.button_process);
        mLogView = (TextView) findViewById(R.id.text_log);

    }

    @Override
    protected void onStart() {
        super.onStart();
        Logger.d("onStart");
    }


    @Override
    protected void onResume() {
        super.onResume();
        Logger.d("onResume");

        mCtrlBtn.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View v) {
                try {
                    if (((String)v.getTag()).equals("0")) {
                        chirpPlayerThread.join();
                        recorderThread.join();

                        mCurrTimeStamp = Timestamp.getCurrtimeStamp();
                        mInputUsername = mUsername.getText().toString();
                        mDataFilename = mBasePath + mInputUsername + mCurrTimeStamp + "LR.wav";
                        Logger.d(mDataFilename);

                        stereoRecorder.setRecFileName(mDataFilename);

                        recorderThread = new Thread(stereoRecorder);
                        chirpPlayerThread = new Thread(chirpPlayer);

                        recorderThread.start();
                        chirpPlayerThread.start();

                        v.setTag("1");
                        mCtrlBtn.setText("STOP");
                    } else {
                        // recorder is stopped from player
                        chirpPlayer.setPlaying(false);
                        while (stereoRecorder.getRecording()) {Thread.yield();}
                        v.setTag("0");
                        mCtrlBtn.setText("START");
                    }
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        });

        mProcessBtn.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View v) {
                try {
                    Logger.d(mDataFilename);
                    if (mDataFilename != null){
                        mProcessing_results = ProcessingFromJNI(mDataFilename);
                        Logger.d(mProcessing_results);

                    }else {
                        Logger.d("Wait....");
                    }


                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        });

    }

    @Override
    protected void onPause() {
        super.onPause();
        Logger.d("onPause");
    }

    @Override
    protected void onStop() {
        super.onStop();
        Logger.d("onStop");
    }

    @Override
    protected void onDestroy() {
        super.onDestroy();
        Logger.d("onDestroy");
        Logger.clearLogAdapters();
    }

}