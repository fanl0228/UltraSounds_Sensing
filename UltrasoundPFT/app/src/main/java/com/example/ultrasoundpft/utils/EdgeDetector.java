package com.example.ultrasoundpft.utils;

import java.util.ArrayDeque;
import java.util.Queue;
public class EdgeDetector {
    private Queue<Double> dataQueue = new ArrayDeque<>();
    private int windowSize;
    private double threshold;

    public EdgeDetector(int windowSize, double threshold) {
        this.windowSize = windowSize;
        this.threshold = threshold;
    }

    // 添加新的数据点
    public void addDataPoint(double newDataPoint) {
        dataQueue.add(newDataPoint);
        if (dataQueue.size() > windowSize) {
            dataQueue.poll(); // 移除队列头部的数据点，保持窗口大小不变
        }
    }

    // 检测上升沿
    public boolean detectRisingEdge() {
        // 计算当前窗口内的滑动标准差
        double movingStdDev = getMovingStandardDeviation();
        // 判断是否发生上升沿
        return movingStdDev >= threshold;
    }

    // 检测下降沿
    public boolean detectFallingEdge() {
        // 计算当前窗口内的滑动标准差
        double movingStdDev = getMovingStandardDeviation();
        // 判断是否发生下降沿
        return movingStdDev >= threshold;
    }

    // 返回当前窗口内的滑动标准差
    private double getMovingStandardDeviation() {
        double sum = 0.0;
        double sumOfSquares = 0.0;
        double mean;

        for (double dataPoint : dataQueue) {
            sum += dataPoint;
            sumOfSquares += dataPoint * dataPoint;
        }

        mean = sum / dataQueue.size();
        double variance = (sumOfSquares / dataQueue.size()) - (mean * mean);

        return Math.sqrt(variance);
    }

    // 返回当前窗口内的数据点均值
    private double getMovingAverage() {
        double sum = 0.0;

        for (double dataPoint : dataQueue) {
            sum += dataPoint;
        }

        return sum / dataQueue.size();
    }

    // 返回当前窗口内最大上升沿的横坐标值
    public double getMaxRisingEdgeXValue() {
        double maxRisingEdge = Double.MIN_VALUE;
        double maxRisingEdgeXValue = 0.0;

        for (double dataPoint : dataQueue) {
            addDataPoint(dataPoint); // 更新数据队列
            double movingAverage = getMovingAverage();
            if (dataPoint - movingAverage >= threshold && dataPoint > maxRisingEdge) {
                maxRisingEdge = dataPoint;
                maxRisingEdgeXValue = dataQueue.size(); // 横坐标值即为数据队列的大小
            }
        }

        return maxRisingEdgeXValue;
    }

    // 返回当前窗口内最大下降沿的横坐标值
    public double getMaxFallingEdgeXValue() {
        double maxFallingEdge = Double.MAX_VALUE;
        double maxFallingEdgeXValue = 0.0;

        for (double dataPoint : dataQueue) {
            addDataPoint(dataPoint); // 更新数据队列
            double movingAverage = getMovingAverage();
            if (movingAverage - dataPoint >= threshold && dataPoint < maxFallingEdge) {
                maxFallingEdge = dataPoint;
                maxFallingEdgeXValue = dataQueue.size(); // 横坐标值即为数据队列的大小
            }
        }

        return maxFallingEdgeXValue;
    }


}