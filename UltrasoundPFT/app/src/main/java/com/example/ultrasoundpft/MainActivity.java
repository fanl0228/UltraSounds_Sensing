package com.example.ultrasoundpft;

import androidx.appcompat.app.AppCompatActivity;
import androidx.core.content.ContextCompat;

import android.annotation.SuppressLint;
import android.app.AlertDialog;
import android.content.DialogInterface;
import android.graphics.Color;
import android.graphics.DashPathEffect;
import android.graphics.Typeface;
import android.graphics.drawable.Drawable;
import android.os.AsyncTask;
import android.os.Bundle;
import android.os.Handler;
import android.view.View;
import android.widget.Button;
import android.widget.EditText;
import android.widget.TextView;


import com.example.ultrasoundpft.utils.Common;
import com.example.ultrasoundpft.utils.Timestamp;
import com.github.mikephil.charting.charts.LineChart;
import com.github.mikephil.charting.components.AxisBase;
import com.github.mikephil.charting.components.Legend;
import com.github.mikephil.charting.components.LimitLine;
import com.github.mikephil.charting.components.XAxis;
import com.github.mikephil.charting.components.YAxis;
import com.github.mikephil.charting.data.Entry;
import com.github.mikephil.charting.data.LineData;
import com.github.mikephil.charting.data.LineDataSet;
import com.github.mikephil.charting.formatter.IAxisValueFormatter;
import com.github.mikephil.charting.formatter.IFillFormatter;
import com.github.mikephil.charting.highlight.Highlight;
import com.github.mikephil.charting.interfaces.dataprovider.LineDataProvider;
import com.github.mikephil.charting.interfaces.datasets.ILineDataSet;
import com.github.mikephil.charting.listener.OnChartValueSelectedListener;
import com.github.mikephil.charting.utils.Utils;
import com.orhanobut.logger.AndroidLogAdapter;
import com.orhanobut.logger.FormatStrategy;
import com.orhanobut.logger.Logger;
import com.orhanobut.logger.PrettyFormatStrategy;

import java.sql.Array;

import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;

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
 *
 * Android Studio WIFIADB (Windows)
 * 1. cd C:\Users\Dell\AppData\Local\Android\Sdk\platform-tools
 * 2. adb tcpip 5555
 * 3. adb connect 192.168.2.186:5555
 *
 * adb pull /storage/emulated/0/Music/
 * adb pull /storage/emulated/0/Download/
 *
 * adb shell rm -rf /storage/emulated/0/Music/
 * adb shell rm -rf /storage/emulated/0/Download/
 *
 */


public class MainActivity extends AppCompatActivity implements OnChartValueSelectedListener {
    // Used to load the 'native-lib' library on application startup.
    static {
        System.loadLibrary("native-lib");
    }

//    private ActivityMainBinding binding;

    /**
     * A native method that is implemented by the 'native-lib' native library,
     * which is packaged with this application.
     */
    public native String stringFromJNI(String str);
    public native float[] ProcessingFromJNI(String filename);

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
    private float[] mProcessing_results;
    private boolean mIsProcessing = false;
    private boolean mIsRecordData = false;
    private boolean mIsProcessingSuccess = false;
    private AlertDialog mProgressDialog;

    // Signal chart
    public LineChart mChart;
    private TextView mtvDelete;
    private ArrayList<Entry> mValues = new ArrayList<>();    //plot values
    private final Handler mHandler = new Handler();

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);

        // Logger adapter init
        FormatStrategy formatStrategy = PrettyFormatStrategy.newBuilder()
                .showThreadInfo(false)  // (Optional) Whether to show thread info or not. Default true
                .methodCount(0)         // (Optional) How many method line to show. Default 2
                .methodOffset(3)        // (Optional) Skips some method invokes in stack trace. Default 5
                .tag("UltrasoundPFT")   // (Optional) Custom tag for each log. Default PRETTY_LOGGER
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

        // init Linechart
        {
            mChart = findViewById(R.id.chart);
            //set background
            mChart.setBackgroundColor(Color.WHITE);
            mChart.setDrawBorders(true);

            //Whether to display coordinate data
            mChart.getDescription().setEnabled(false);
            //Whether to support double-click
            mChart.setTouchEnabled(true);
            //Value selected callback listening
            mChart.setOnChartValueSelectedListener(this);
            //Whether to draw the grid background
            mChart.setDrawGridBackground(false);

            // box displaying coordinate data
            MyMarkerView mv = new MyMarkerView(this, R.layout.custom_marker_view);
            // Set the marker to the chart
            mv.setChartView(mChart);
            mChart.setMarker(mv);
            // enable scaling and dragging
            mChart.setDragEnabled(true);
            mChart.setScaleEnabled(true);
            // mChart.setScaleXEnabled(true);
            // mChart.setScaleYEnabled(true);
            // force pinch zoom along both axis
            mChart.setPinchZoom(true);
        }


    }


    @Override
    protected void onResume() {
        super.onResume();
        Logger.d("onResume");

        mProcessBtn.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View v) {
                try {
//                    Logger.d(mDataFilename);
                    if (mDataFilename != null && !mIsProcessing && !mIsRecordData){

                        showProgressDialog();

                        mIsProcessing = true;
                        Logger.d("Call executeNativeFunctionInBackground");

                        // clear chart value
                        mValues.clear();
                        //Execute the native function after clicking the button
                        executeNativeFunctionInBackground();

                    }else if(mIsProcessing) {
                        Logger.e("Signal processing, please wait....");
                    }
                    else if(mIsRecordData) {
                        Logger.e("Signal collecting, please wait....");
                    }
                    else {
                        Logger.e("Wait....");
                    }

                } catch (Exception e) {
                    e.printStackTrace();
                    Logger.e("onClick failed");
                }
            }
        });
        Logger.d("onResume11111");

    }

    @Override
    protected void onStart(){
        super.onStart();
        Logger.d("onStart");

        mCtrlBtn.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View v) {
                try {
//                    v.setSelected(!v.isSelected());
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
//                        mCtrlBtn.setBackgroundColor(Color.parseColor("#006400"));
                        mIsRecordData = true;
                    } else {
                        // recorder is stopped from player
                        chirpPlayer.setPlaying(false);
                        while (stereoRecorder.getRecording()) {Thread.yield();}
                        v.setTag("0");
                        mCtrlBtn.setText("START");
                        mIsRecordData = false;
                    }
                } catch (Exception e) {

                    mIsRecordData = false;
                    e.printStackTrace();
                }
            }
        });

        startPlotChartThread();

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
        mHandler.removeCallbacksAndMessages(null);

        Logger.d("onDestroy");
        Logger.clearLogAdapters();
    }

    @Override
    public void onValueSelected(Entry e, Highlight h) {
//        Logger.i("Entry selected", e.toString());
//        Logger.i("LOW and HIGH value==> low: " + mChart.getLowestVisibleX() + ", high: " + mChart.getHighestVisibleX());
//        Logger.i("MIN and MAX ==> xMin: " + mChart.getXChartMin() + ", xMax: " + mChart.getXChartMax() + ", yMin: " + mChart.getYChartMin() + ", yMax: " + mChart.getYChartMax());
    }

    @Override
    public void onNothingSelected() {
        Logger.i("Nothing selected.");
    }


    private void startPlotChartThread() {
        // Create a Runnable object and execute your logic in it
        Runnable runnable = new Runnable() {
            @Override
            public void run() {
                // If the processing is successful, perform the corresponding operations
                if (mIsProcessingSuccess) {
                    mIsProcessingSuccess = false;

                    dismissProgressDialog();

                    // Set the thickness of the coordinate axis axis
                    Typeface tfRegular = Typeface.create("sans-serif", Typeface.NORMAL);
                    XAxis xAxis;
                    {   // X-Axis Style //
                        xAxis = mChart.getXAxis();
                        xAxis.setTypeface(tfRegular);
                        xAxis.setValueFormatter(new IAxisValueFormatter() {
                            @Override
                            public String getFormattedValue(float value, AxisBase axis) {
                                String formattedValue = String.format("%.2f", value * 0.08);
                                return formattedValue + "s";
                            }

                        });

                        xAxis.setPosition(XAxis.XAxisPosition.BOTTOM);
                        xAxis.setDrawLabels(true);
                        xAxis.setDrawAxisLine(true);
                        // vertical grid lines
                        xAxis.enableGridDashedLine(10f, 10f, 0f);
                    }
                    YAxis yAxis;
                    {   // // Y-Axis Style // //
                        yAxis = mChart.getAxisLeft();
                        xAxis.setTypeface(tfRegular);
                        yAxis.setValueFormatter(new IAxisValueFormatter() {
                            @Override
                            public String getFormattedValue(float value, AxisBase axis) {
                                return value + "ml";
                            }

                        });
                        // disable dual axis (only use LEFT axis)
                        mChart.getAxisRight().setEnabled(false);
                        // horizontal grid lines
                        yAxis.enableGridDashedLine(10f, 10f, 0f);
                        // axis range
                        yAxis.setAxisMaximum(3000f);
                        yAxis.setAxisMinimum(0f);
                    }

//                    {
//                        // // Create Limit Lines // //
//                        LimitLine llXAxis = new LimitLine(9f, "Index 10");
//                        llXAxis.setLineWidth(4f);
//                        llXAxis.enableDashedLine(10f, 10f, 0f);
//                        llXAxis.setLabelPosition(LimitLine.LimitLabelPosition.RIGHT_BOTTOM);
//                        llXAxis.setTextSize(10f);
//                        llXAxis.setTypeface(tfRegular);
//
//                        LimitLine ll1 = new LimitLine(2000f, "Upper Limit");
//                        ll1.setLineWidth(4f);
//                        ll1.enableDashedLine(10f, 10f, 0f);
//                        ll1.setLabelPosition(LimitLine.LimitLabelPosition.RIGHT_TOP);
//                        ll1.setTextSize(10f);
//                        ll1.setTypeface(tfRegular);
//
//                        LimitLine ll2 = new LimitLine(10f, "Lower Limit");
//                        ll2.setLineWidth(4f);
//                        ll2.enableDashedLine(10f, 10f, 0f);
//                        ll2.setLabelPosition(LimitLine.LimitLabelPosition.RIGHT_BOTTOM);
//                        ll2.setTextSize(10f);
//                        ll2.setTypeface(tfRegular);
//
//                        // draw limit lines behind data instead of on top
//                        yAxis.setDrawLimitLinesBehindData(true);
//                        xAxis.setDrawLimitLinesBehindData(true);
//
//                        // add limit lines
//                        yAxis.addLimitLine(ll1);
//                        yAxis.addLimitLine(ll2);
//                    }

                    // add data
                    setData(mProcessing_results);

                    // draw points over time
                    mChart.animateX(1500);

                    // get the legend (only possible after setting data)
                    Legend l = mChart.getLegend();

                    // draw legend entries as lines
                    l.setForm(Legend.LegendForm.LINE);
                    l.setPosition(Legend.LegendPosition.ABOVE_CHART_RIGHT);

                    mChart.setOnChartValueSelectedListener(MainActivity.this);
                }

                // Execute this thread again after a delay
                mHandler.postDelayed(this, 1000); // 这里的 1000 表示延迟的毫秒数，可以根据需要调整
            }
        };

        // start thread
        mHandler.post(runnable);
    }

    private void setData(float[] velocityData) {

        mValues.clear();

        if(velocityData.length == 0){
            Logger.e("Input velocity data is empty.");
            return;
        }
        for (int i = 0; i < velocityData.length; i++) {
            float val = velocityData[i];
            mValues.add(new Entry(i, val, getResources().getDrawable(R.drawable.star)));
        }

        LineDataSet set1 = new LineDataSet(mValues, "velocity");
        set1.setMode(LineDataSet.Mode.CUBIC_BEZIER);

        if (mChart.getData() != null &&
                mChart.getData().getDataSetCount() > 0) {
            set1 = (LineDataSet) mChart.getData().getDataSetByIndex(0);
            set1.setValues(mValues);
            set1.notifyDataSetChanged();
            mChart.getData().notifyDataChanged();
            mChart.notifyDataSetChanged();
        } else {
            // create a dataset and give it a type
            set1 = new LineDataSet(mValues, "Airflow velocity");

            set1.setDrawIcons(false);

            // draw dashed line
            set1.enableDashedLine(10f, 5f, 0f);

            // black lines and points
            set1.setColor(Color.BLACK);
            set1.setCircleColor(Color.BLACK);

            // line thickness and point size
            set1.setLineWidth(1f);
            set1.setCircleRadius(3f);

            // draw points as solid circles
            set1.setDrawCircleHole(false);

            // customize legend entry
            set1.setFormLineWidth(1f);
            set1.setFormLineDashEffect(new DashPathEffect(new float[]{10f, 5f}, 0f));
            set1.setFormSize(15.f);

            // text size of values
            set1.setValueTextSize(9f);

            // draw selection line as dashed
            set1.enableDashedHighlightLine(10f, 5f, 0f);

            // set the filled area
            set1.setDrawFilled(true);
            set1.setFillFormatter(new IFillFormatter() {
                @Override
                public float getFillLinePosition(ILineDataSet dataSet, LineDataProvider dataProvider) {
                    return mChart.getAxisLeft().getAxisMinimum();
                }
            });

            // set color of filled area
            if (Utils.getSDKInt() >= 18) {
                // drawables only supported on api level 18 and above

                Drawable drawable = ContextCompat.getDrawable(this, R.drawable.fade_red);
                set1.setFillDrawable(drawable);
            } else {
                set1.setFillColor(Color.BLACK);
            }

            ArrayList<ILineDataSet> dataSets = new ArrayList<>();
            dataSets.add(set1); // add the data sets

            // create a data object with the data sets
            LineData data = new LineData(dataSets);

            // set data
            mChart.setData(data);


        }

        Logger.i("setData success.");
    }

    //Execute native function in background thread
    private void executeNativeFunctionInBackground() {
        new Thread(new Runnable() {
            @Override
            public void run() {
                // Call native function
                mProcessing_results = ProcessingFromJNI(mDataFilename);
                mIsProcessing = false;

                if(mProcessing_results.length > 0){
                    mIsProcessingSuccess = true; // for plot chart
                    Logger.d("Signal processing success.");
                    Logger.d(mProcessing_results);
                } else{
                    mIsProcessingSuccess = false;
                    dismissProgressDialog();

                    //
                    processingSignalFailed();
                    Logger.e("Signal processing failed.");
                }

            }
        }).start();
    }

    private void showProgressDialog() {
        AlertDialog.Builder builder = new AlertDialog.Builder(this);
        builder.setMessage("Processing, please wait...")
                .setCancelable(false);
        mProgressDialog = builder.create();
        mProgressDialog.show();
    }

    private void dismissProgressDialog() {
        if (mProgressDialog != null && mProgressDialog.isShowing()) {
            mProgressDialog.dismiss();
        }
    }

    private void processingSignalFailed(){
        runOnUiThread(new Runnable() {
            @Override
            public void run() {
                new AlertDialog.Builder(MainActivity.this)
                        .setTitle("Signal processing failed.")
                        .setMessage("Please collect the data again and process it!")
                        .setPositiveButton("Yes", new DialogInterface.OnClickListener() {
                            @Override
                            public void onClick(DialogInterface dialog, int which) {
                                // 处理点击"Yes"按钮的操作，如果不需要执行任何操作，可以传递null
                            }
                        })
                        .setNegativeButton("No", new DialogInterface.OnClickListener() {
                            @Override
                            public void onClick(DialogInterface dialog, int which) {
                                // 处理点击"No"按钮的操作，如果不需要执行任何操作，可以传递null
                            }
                        })
                        .show();
            }
        });
    }

}