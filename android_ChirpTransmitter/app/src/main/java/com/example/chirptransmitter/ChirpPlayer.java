package com.example.chirptransmitter;

import android.media.AudioAttributes;
import android.media.AudioFormat;
import android.media.AudioTrack;


import com.orhanobut.logger.Logger;

import com.example.chirptransmitter.utils.ChirpGenerator;

public class ChirpPlayer implements Runnable {

//    private MonoRecorder recorder = null;
    private com.example.chirptransmitter.StereoRecorder recorder = null;
    public boolean playing = false;
    public int sampleRate = 48000;
    public int sessionID = 0;

    public double chirpStartFrq = 17000;
    public double chirpStopFrq = 23000;
    public double chirpDuration = 10 / 1000.0;
    public double chirpFrameDuration = 10 / 1000.0;
    public int chirpCircle = 1;
    public boolean applyHamming = false;

    AudioTrack track = new AudioTrack.Builder()
            .setAudioAttributes(new AudioAttributes.Builder()
                    .setUsage(AudioAttributes.USAGE_MEDIA)  // USAGE_VOICE_COMMUNICATION   USAGE_MEDIA  USAGE_VOICE_COMMUNICATION_SIGNALLING
                    .setContentType(AudioAttributes.CONTENT_TYPE_MUSIC) // CONTENT_TYPE_SPEECH CONTENT_TYPE_MUSIC
                    .build())
            .setAudioFormat(new AudioFormat.Builder()
                    .setEncoding(AudioFormat.ENCODING_PCM_16BIT)
                    .setSampleRate(this.sampleRate)
                    .setChannelMask(AudioFormat.CHANNEL_OUT_STEREO)  //CHANNEL_OUT_STEREO   CHANNEL_OUT_MONO
                    .build())
            .setTransferMode(AudioTrack.MODE_STREAM)
            .build();


    ChirpPlayer(com.example.chirptransmitter.StereoRecorder recorder) {
        this.recorder = recorder;
    }

    @Override
    public void run() {

        byte[] chirpFrame = ChirpGenerator.genPCM16MonoToneBytes(
                this.sampleRate, chirpStartFrq, chirpStopFrq,
                chirpDuration, chirpFrameDuration, chirpCircle,
                applyHamming, ChirpGenerator.SPEAKER_TYPE_BOTTOM);

//        byte[] chirpFrame = ChirpGenerator.genPCM16SingleToneBytes(
//                this.sampleRate, chirpStartFrq, chirpDuration,
//                chirpFrameDuration, chirpCircle, applyHamming);

        Logger.d(chirpFrame);

        // wait until recorder start recording
        while (!recorder.getRecording()) { Thread.yield(); }

        track.play();
        setPlaying(true);
        while (getPlaying()) {
            track.write(chirpFrame, 0, chirpFrame.length);
        }
        track.pause();
        track.flush();

        // stop recorder from player
        recorder.setRecording(false);

    }

    public int getSessionID() {
        return this.sessionID;
    }
    public void setSessionID(int id) {
        this.sessionID = id;
    }

    synchronized public boolean getPlaying() {
        return playing;
    }
    synchronized public void setPlaying(boolean v) {
        playing = v;
    }

}
