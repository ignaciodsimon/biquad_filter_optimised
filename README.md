# Biquad filter libary for Python (implemented in C++)

This is a simple library for Python, which can be used to
create objects of type BiquadFilter that can filter signals
and also calculate coefficients from the standard filter
parameters (type, frequency, gain and Q).

## How to install it

From your terminal, run:

    pip3 install https://github.com/ignaciodsimon/biquad_filter_optimised/archive/refs/heads/master.zip

## Use example

    # Import this library (install it first)
    from biquad_filter_optimised import *
    
    # Create a filter
    sampleRate = 48000.0
    filt = BiquadFilter()
    
    # Tell it to generate its coefficients to be a second-order low-pass filter at 1 kHz
    filt.generateBiQuadCoefficients(
        BiquadFilterParameters(
            sampleRate   = sampleRate,
            filterf0     = 1000,
            filterQ      = 2.0,
            filterGaindB = 0.0,
            filterType   = BiquadFilterType.LPF))
    
    # Measure the impulse response
    IR = [filt.processSample(1.0 if i == 0 else 0.0) for i in range(int(0.5 * sampleRate))]
    timeBins = [i / sampleRate for i in range(len(IR))]
    
    # Display the impulse response
    plot.plot(timeBins, IR)
    plot.grid(1)
    plot.show()

