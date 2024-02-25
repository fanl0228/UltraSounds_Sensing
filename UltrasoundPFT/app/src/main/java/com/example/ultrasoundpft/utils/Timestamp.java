package com.example.ultrasoundpft.utils;

import java.text.SimpleDateFormat;
import java.util.Date;

public class Timestamp {

    public static String getCurrtimeStamp(){
        Date currentTime=new Date();
        SimpleDateFormat TimeFmt=new SimpleDateFormat("yyyy_MM_dd_HH_mm_ss");
        return TimeFmt.format(currentTime);
    }

}


