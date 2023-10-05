"""
    This script is used to demonstrate the usage of the biquad
    filter library.

    https://github.com/ignaciodsimon/biquad_filter_optimised

    Joe Simon 2018 - 2023.
"""

from biquad_filter_optimised import *
import numpy
import matplotlib.pyplot as plot


if __name__ == '__main__':

    sampleRate = 48000.0
    filt = BiquadFilter()

    types = [
        BiquadFilterType.BPF_NORMALIZED,
        BiquadFilterType.HPF,
        BiquadFilterType.BPF,
        BiquadFilterType.LPF,
        BiquadFilterType.NOTCH,
        BiquadFilterType.APF,
        BiquadFilterType.LOW_SHELVING,
        BiquadFilterType.HIGH_SHELVING,
        BiquadFilterType.PEAK,
        BiquadFilterType.LPF_FIRST_ORDER,
        BiquadFilterType.HPF_FIRST_ORDER,
        BiquadFilterType.APF_FIRST_ORDER,
        BiquadFilterType.LOW_SHELVING_FIRST_ORDER,
        BiquadFilterType.HIGH_SHELVING_FIRST_ORDER]

    for index, filterType in enumerate(types):

        # Tell the filter to calculate new coefficients following the standard parameters
        filt.generateBiQuadCoefficients(
            BiquadFilterParameters(
                sampleRate   = sampleRate,
                filterf0     = 1000,
                filterQ      = 2.0,
                filterGaindB = -20.0,
                filterType   = filterType))

        # Measure the impulse response
        IR = [filt.processSample(1.0 if i == 0 else 0.0) for i in range(int(0.5 * sampleRate))]

        # Compute the spectrum module
        filterResponse = 20.0 * numpy.log10(numpy.abs(numpy.fft.rfft(IR)))
        freqAxis = [i * sampleRate / len(filterResponse) / 2.0 for i in range(len(filterResponse))]

        # Plot it
        plot.semilogx(freqAxis, filterResponse, label=filterType)

    plot.xlabel("Frequency [Hz]")
    plot.ylabel("Magnitude [dB]")
    plot.legend()
    plot.grid(1)
    plot.xlim([20, 20e3])
    plot.ylim([-60, 10])
    plot.show()
