package com.example.chirptransmitter.utils;

/**
 * // only play sound on left
 * for(int i = 0; i < count; i += 2){
 *     short sample = (short)(Math.sin(2 * Math.PI * i / (44100.0 / freqHz)) * 0x7FFF);
 *     samples[i + 0] = sample;
 *     samples[i + 1] = 0;
 * }
 * // only play sound on right
 * for(int i = 0; i < count; i += 2){
 *     short sample = (short)(Math.sin(2 * Math.PI * i / (44100.0 / freqHz)) * 0x7FFF);
 *     samples[i + 0] = 0;
 *     samples[i + 1] = sample;
 * }
 *
 * */




public class ChirpGenerator {

    public final static int SPEAKER_TYPE_ZERO   = 0;
    /** Using top speaker only */
    public final static int SPEAKER_TYPE_TOP    = 1;
    /** Using bottom speaker only */
    public final static int SPEAKER_TYPE_BOTTOM = 2;
    /** Using top and bottom speaker */
    public final static int SPEAKER_TYPE_BOTH   = 3;


    static public byte[] genPCM16SingleToneBytes(
            int sampleRate,
            double freq,
            double tSweep,
            double tFrame,
            int times,
            boolean applyHamming) {


        int chirpSampleCnt = (int) Math.round(tSweep * sampleRate);
        int frameSampleCnt = (int) Math.round(tFrame * sampleRate);

        // generate Window
        int windowSize = chirpSampleCnt;

        double[] sample_top = new double[frameSampleCnt * times];
        double[] sample_bottom = new double[frameSampleCnt * times];

        double[] factors = new double[windowSize];
        for(int i = 0; i < windowSize; i++) {
            if (applyHamming) {
                factors[i] = 0.5 * (1 - Math.cos((2 * Math.PI * i) / (windowSize-1)));
            } else {
                factors[i] = 1.0d;
            }
        }

        int t = 0;
        while (t < times) {
            double instfreq_top, instfreq_bottom, numerator;
            for (int i = 0; i < chirpSampleCnt; i++) {
                numerator = (double)(i) / (double) chirpSampleCnt;
                instfreq_top = numerator * freq; //fStart + (0.5 * numerator * (fHalf - fStart));
                sample_top[i + t * frameSampleCnt] = factors[i] * Math.cos(2.0 * Math.PI * i / (sampleRate / instfreq_top));
            }
            for (int i = chirpSampleCnt; i < frameSampleCnt; i++) {
                sample_top[i + t * frameSampleCnt] = 0.0;
            }

            for (int i = 0; i < chirpSampleCnt; i++) {
                numerator = (double)(i) / (double) chirpSampleCnt;
                instfreq_bottom = numerator * freq; //fHalf + (0.5 * numerator * (fEnd - fHalf));
                sample_bottom[i + t * frameSampleCnt] = factors[i] * Math.cos(2.0 * Math.PI * i / (sampleRate / instfreq_bottom ));
            }
            for (int i = chirpSampleCnt; i < frameSampleCnt; i++) {
                sample_bottom[i + t * frameSampleCnt] = 0.0;
            }

            t++;
        }

        int idx = 0;
        byte[] data = new byte[2 * 2 * frameSampleCnt * times];
        for (final double dVal : sample_top) {
            // left channel data
            final short val = (short) ((dVal * 32767));
            data[idx++] = (byte) (val & 0x00ff);
            data[idx++] = (byte) ((val & 0xff00) >>> 8);

            // right channel data
            data[idx++] = (byte) (0);
            data[idx++] = (byte) (0);
        }

//        idx = 0; // reset idx
//        for (final double dVal : sample_bottom) {
//            // left channel data
//            idx++;
//            idx++;
//
//            // right channel data
//            final short val = (short) ((dVal * 32767));
//            data[idx++] = (byte) (val & 0x00ff);
//            data[idx++] = (byte) ((val & 0xff00) >>> 8);
//        }

        return data;
    }


    /** Generate the speaker chirp signal data */
    static public byte[] genPCM16StereoBytes(
            int sampleRate,
            double fStart,
            double fEnd,
            double tSweep,
            double tFrame,
            int times,
            boolean applyHamming,
            int speakerType) {

        double fHalf = (fEnd + fStart) / 2.0;

        int chirpSampleCnt = (int) Math.round(tSweep * sampleRate);
        int frameSampleCnt = (int) Math.round(tFrame * sampleRate);

        // generate Window
        int windowSize = chirpSampleCnt;

        double[] sample_top = new double[frameSampleCnt * times];
        double[] sample_bottom = new double[frameSampleCnt * times];

        double[] factors = new double[windowSize];
        for(int i = 0; i < windowSize; i++) {
            if (applyHamming) {
                factors[i] = 0.5 * (1 - Math.cos((2 * Math.PI * i) / (windowSize-1)));
            } else {
                factors[i] = 1.0d;
            }
        }

        int t = 0;
        while (t < times) {
            double instfreq_top, instfreq_bottom, numerator;
            for (int i = 0; i < chirpSampleCnt; i++) {
                numerator = (double)(i) / (double) chirpSampleCnt;
                instfreq_top = fStart + (0.5 * numerator * (fHalf - fStart));
                sample_top[i + t * frameSampleCnt] = factors[i] * Math.cos(2.0 * Math.PI * i / (sampleRate / instfreq_top));
            }
            for (int i = chirpSampleCnt; i < frameSampleCnt; i++) {
                sample_top[i + t * frameSampleCnt] = 0.0;
            }

            for (int i = 0; i < chirpSampleCnt; i++) {
                numerator = (double)(i) / (double) chirpSampleCnt;
                instfreq_bottom = fHalf + (0.5 * numerator * (fEnd - fHalf));
                sample_bottom[i + t * frameSampleCnt] = factors[i] * Math.cos(2.0 * Math.PI * i / (sampleRate / instfreq_bottom));
            }
            for (int i = chirpSampleCnt; i < frameSampleCnt; i++) {
                sample_bottom[i + t * frameSampleCnt] = 0.0;
            }

            t++;
        }

        int idx = 0;
        byte[] data = new byte[2 * 2 * frameSampleCnt * times];
        for (final double dVal : sample_top) {
            // left channel data
            final short val = (short) ((dVal * 32767));
            data[idx++] = (byte) (val & 0x00ff);
            data[idx++] = (byte) ((val & 0xff00) >>> 8);

            // right channel data
            data[idx++] = (byte) (0);
            data[idx++] = (byte) (0);
        }

        idx = 0; // reset idx
        for (final double dVal : sample_bottom) {
            // left channel data
            idx++;
            idx++;

            // right channel data
            final short val = (short) ((dVal * 32767));
            data[idx++] = (byte) (val & 0x00ff);
            data[idx++] = (byte) ((val & 0xff00) >>> 8);
        }

        return data;
    }


    /** Generate the speaker chirp signal data */
    static public byte[] genPCM16MonoToneBytes(
            int sampleRate,
            double fStart,
            double fEnd,
            double tSweep,
            double tFrame,
            int times,
            boolean applyHamming,
            int speakerType) {

        int chirpSampleCnt = (int) Math.round(tSweep * sampleRate);
        int frameSampleCnt = (int) Math.round(tFrame * sampleRate);

        // generate Window
        int windowSize = chirpSampleCnt;
        double[] sample = new double[frameSampleCnt * times];
        double[] factors = new double[windowSize];
        for(int i = 0; i < windowSize; i++) {
            if (applyHamming) {
                factors[i] = 0.5 * (1 - Math.cos((2 * Math.PI * i) / (windowSize-1)));
            } else {
                factors[i] = 1.0d;
            }
        }

        int t = 0;
        while (t < times) {
            double instfreq, numerator;
            for (int i = 0; i < chirpSampleCnt; i++)
            {
                numerator = (double)(i) / (double) chirpSampleCnt;
                instfreq = fStart + (0.5 * numerator * (fEnd - fStart));
                // apply window
                sample[i + t * frameSampleCnt] = factors[i] * Math.cos(2.0 * Math.PI * i / (sampleRate / instfreq));
            }
            for (int i = chirpSampleCnt; i < frameSampleCnt; i++) {
                sample[i + t * frameSampleCnt] = 0.0;
            }
            t++;
        }

        int idx = 0;
        byte[] data = new byte[2 * 2 * frameSampleCnt * times];
        for (final double dVal : sample) {

            if (speakerType == SPEAKER_TYPE_BOTTOM){
                // left channel data
                data[idx++] = (byte) (0);
                data[idx++] = (byte) (0);

                // right channel data
                final short val = (short) ((dVal * 32767)); //scale to maximum amplitude
                data[idx++] = (byte) (val & 0x00ff);
                data[idx++] = (byte) ((val & 0xff00) >>> 8);
            }else if (speakerType == SPEAKER_TYPE_TOP){
                // left channel data
                final short val = (short) ((dVal * 32767));
                data[idx++] = (byte) (val & 0x00ff);
                data[idx++] = (byte) ((val & 0xff00) >>> 8);

                // right channel data
                data[idx++] = (byte) (0);
                data[idx++] = (byte) (0);
            }
            else if (speakerType == SPEAKER_TYPE_BOTH) {
                // left channel data
                final short val = (short) ((dVal * 32767));
                data[idx++] = (byte) (val & 0x00ff);
                data[idx++] = (byte) ((val & 0xff00) >>> 8);

                // right channel data
                data[idx++] = (byte) (val & 0x00ff);
                data[idx++] = (byte) ((val & 0xff00) >>> 8);

            }
            else if (speakerType == SPEAKER_TYPE_ZERO) {
                data[idx++] = (byte) (0);
                data[idx++] = (byte) (0);

                data[idx++] = (byte) (0);
                data[idx++] = (byte) (0);

            }
            else {
                data[idx++] = (byte) (0);
                data[idx++] = (byte) (0);

                data[idx++] = (byte) (0);
                data[idx++] = (byte) (0);
                break;
            }

        }

        return data;
    }

}
