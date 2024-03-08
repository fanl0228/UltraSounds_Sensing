package com.example.ultrasoundpft;

import android.content.Context;
import android.widget.TextView;

import com.github.mikephil.charting.components.MarkerView;
import com.github.mikephil.charting.data.Entry;
import com.github.mikephil.charting.highlight.Highlight;

import static com.example.ultrasoundpft.R.id.tvContent;

public class MyMarkerView extends MarkerView {

    private TextView tvContent;

    public MyMarkerView(Context context, int layoutResource) {
        super(context, layoutResource);

        // Initialize UI components
        tvContent = findViewById(R.id.tvContent);
    }

    @Override
    public void refreshContent(Entry e, Highlight highlight) {
        // Display the value of the Entry
        tvContent.setText("Value: " + e.getY());

        // This will force the marker-view to redraw
        super.refreshContent(e, highlight);
    }

//    @Override
//    public int getXOffset(float xpos) {
//        // Adjust the offset of the MarkerView
//        return -(getWidth() / 2);
//    }
//
//    @Override
//    public int getYOffset(float ypos) {
//        // Adjust the offset of the MarkerView
//        return -getHeight();
//    }
}
