"""
    This script is used to test both implementations of 
    the biquad filter and measure their time performance.

    Joe Simon 2018.
"""

import biquad_filter_optimised
import biquad_filter_original
import numpy
import matplotlib.pyplot as plot
import time


if __name__ == '__main__':

    _sampleRate = 48000.0
    print("Instantiating filters ...")
    _filter1 = biquad_filter_original.BiquadFilter(sampleRate=_sampleRate)
    _filter2 = biquad_filter_optimised.BiquadFilter()

    print("Generating filter coefficients ...")
    _filter1.generateBiQuadCoefficients(filterf0=1000, filterQ=0.707, filterType=biquad_filter_original.BiquadFilterType.LPF)
    _filter2.generateBiQuadCoefficients(biquad_filter_optimised.BiquadFilterParameters(filterf0=1000, filterQ=0.707, filterType=biquad_filter_optimised.BiquadFilterType.LPF))

    # Example of how to use the method 'computeBiquadFilterResponse()' without
    # having to 'measure' it passing a test signal and alysing the output.
    #
    # _filter2Response = _filter2.computeBiquadFilterResponse([i/_sampleRate for i in range(int(_sampleRate/2))])
    # plot.subplot(2,1,1)
    # plot.semilogx(20.0*numpy.log10(numpy.abs(_filter2Response)))
    # plot.ylim([-60, 10])
    # plot.grid(1)
    # plot.subplot(2,1,2)
    # plot.semilogx(180.0 / numpy.pi * numpy.angle(_filter2Response))
    # plot.grid(1)
    # plot.show()

    _impulse = [0.0]*int(_sampleRate * 10.0)
    _impulse[0] = 1.0
    print("Processing impulse signal with filters ...")
    _amountOfTests = 10
    _timesFilter1 = []
    _timesFilter2 = []
    _testStartTime = time.time()
    for i in range(_amountOfTests):
        _initialTime1 = time.time()
        _filter1Output = [_filter1.processSample(inputSample=_sample) for _sample in _impulse]
        _finalTime1 = time.time()
        _initialTime2 = time.time()
        _filter2Output = [_filter2.processSample(_sample) for _sample in _impulse]
        _finalTime2 = time.time()
        _timesFilter1.append(_finalTime1 - _initialTime1)
        _timesFilter2.append(_finalTime2 - _initialTime2)

    _improvementRatios = [_timesFilter1[i] / _timesFilter2[i] for i in range(_amountOfTests)]
    print("Time improvement ratio: " + str(["%.2f" % _improvementRatios[i] for i in range(_amountOfTests)]))
    print("Mean improvement ratio: %.2f" % numpy.mean(_improvementRatios))
    print("Total elapsed time: %.2f [s]" % (time.time() - _testStartTime))

    print("Plotting filters response ...")
    _filter1Response = 20.0 * numpy.log10(numpy.abs(numpy.fft.rfft(_filter1Output[:48000])) + 10.0**-200/20.0)
    _filter2Response = 20.0 * numpy.log10(numpy.abs(numpy.fft.rfft(_filter2Output[:48000])) + 10.0**-200/20.0)
    _freqAxis = [i * _sampleRate / len(_filter1Response) / 2.0 for i in range(len(_filter1Response))]
    plot.semilogx(_freqAxis, _filter1Response)
    plot.semilogx(_freqAxis, _filter2Response)
    plot.legend(['Original %.1f [s]'  % numpy.mean(_timesFilter1), 'Optimised %.1f [s]' % numpy.mean(_timesFilter2)])
    plot.grid(1)
    plot.title("Time improvement ratio: %.2f [.]" % numpy.mean(_improvementRatios))
    plot.ylim([-90, 10])
    plot.show()
