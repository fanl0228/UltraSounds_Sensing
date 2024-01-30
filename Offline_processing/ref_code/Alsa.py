import subprocess
from datetime import datetime
from pathlib import Path
import numpy as np


# Linux Only
def AplayArecord(playFile:str,
                 recordFile: str,
                 duration: float,
                 nChannels: int = 2,
                 sampleRate: int = 48000,
                 playHw: str = None,
                 recordHw: str = None) -> bool:
    
    playHw = f' -D "{play_hw}" ' if playHw else ''
    recordHw = f' -D "{record_hw}" ' if recordHw else ''
    samples = int(duration*sampleRate)
    cmd = f'(aplay "{playFile}" {playHw} &) && arecord {recordHw} -r {sampleRate} -c {nChannels} -f S16_LE -s {samples} "{recordFile}"'
    #cmd = f'arecord {recordHw} -r {sampleRate} -c {nChannels} -f S16_LE -s {samples} "{recordFile}"'
    print(cmd)
    p = subprocess.run(cmd, shell=True)
    return p.returncode==0
