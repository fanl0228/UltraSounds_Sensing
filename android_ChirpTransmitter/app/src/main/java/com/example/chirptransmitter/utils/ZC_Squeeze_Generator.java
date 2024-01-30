package com.example.chirptransmitter.utils;

/**
 *
 * */



public class ZC_Squeeze_Generator {

    /** Using top speaker only */
    public final static int SPEAKER_TYPE_TOP    = 1;
    /** Using bottom speaker only */
    public final static int SPEAKER_TYPE_BOTTOM = 2;
    /** Using top and bottom speaker */
    public final static int SPEAKER_TYPE_BOTH   = 3;

    /** Generate the speaker chirp signal data */
    static public byte[] genPCM16MonoToneBytes(
            int sampleRate,
            int root_index,
            int ZC_length,
            int duration_time,
            boolean applyHamming,
            int speakerType) {

        int Zc_signal_sampleCnt = (int) Math.round(duration_time * sampleRate);
        double[] signal_samples = new double[Zc_signal_sampleCnt];

        // generate Window
        int windowSize = Zc_signal_sampleCnt;
        double[] window_function = new double[windowSize];
        for(int i = 0; i < windowSize; i++) {
            if (applyHamming) {
                window_function[i] = 0.5 * (1 - Math.cos((2 * Math.PI * i) / (windowSize-1)));
            } else {
                window_function[i] = 1.0d;
            }
        }

        // generate signal
        for (int n = 0; n < Zc_signal_sampleCnt; n++) {
            signal_samples[n] = Math.exp(-1 * Math.PI * root_index * n * (n+1) / (ZC_length));
        }


        // sample signal -> left or right channels
        int idx = 0;
        byte[] data = new byte[2 * 2 * Zc_signal_sampleCnt];
        for (final double dVal : signal_samples) {

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
