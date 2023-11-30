import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.fft import fft, ifft, fftfreq

def plot_waveform(filename):
    y = np.loadtxt(filename)
    instrument = filename[:-4]
    #plot the waveform
    fig, ax = plt.subplots()
    ax.plot(range(len(y)), y, 'k-', alpha = .5, linewidth = .1, label=instrument)
    ax.set_xlabel('time [arbs]')
    ax.set_ylabel('amplitude [arbs]')
    plt.tight_layout()
    plt.savefig(f'waveform_{instrument}.png', dpi=200)

def FFT_waveform(filename):
    y = np.loadtxt(filename)
    instrument = filename[:-4]
    #do FFT
    f0 = 44100 #from part b)
    y_fft = fft(np.array(y))
    f = fftfreq(len(y), 1/f0) #number of points is first arg, second arg is sample spacing in time
    psf = np.sqrt(np.real(y_fft)**2 + np.imag(y_fft)**2)
    fig, ax = plt.subplots()
    ax.plot(f[:10000], psf[:10000], 'b-', alpha = 1., linewidth = .5, label = instrument)
    ax.set_xlabel('f [Hz]')
    ax.set_ylabel('PSF [arbs]')
    plt.savefig(f'fft_{instrument}.png')

    #print the note being played
    max_index = np.argmax(psf[:10000])
    print(f'for the {instrument}, the note being played was at {f[max_index]}Hz')
    #should be 524 and 1043 Hz. https://nickfever.com/music/note-frequencies <- look at middle C, octave 5 and octave 6. Same note, just different octaves.
    
if __name__ == "__main__":
    plot_waveform('piano.txt')
    plot_waveform('trumpet.txt')
    FFT_waveform('piano.txt')
    FFT_waveform('trumpet.txt')
