'''
    Standard implementation of a biquad filter
    Joe Simon 2018.
'''


import math
from enum import Enum


class BiquadFilterCoefficients():
    def __init__(self, b0=1.0, b1=0, b2=0, a0=0, a1=0, a2=0):
        self.b0 = b0
        self.b1 = b1
        self.b2 = b2
        self.a0 = a0
        self.a1 = a1
        self.a2 = a2


class BiquadFilterType(Enum):
    BPF_NORMALIZED = 1
    HPF            = 2
    BPF            = 3
    LPF            = 4
    NOTCH          = 5
    APF            = 6
    LOW_SHELVING   = 7
    HIGH_SHELVING  = 8
    PEAK           = 9


class BiquadFilterParameters:

    def __init__(self, sampleRate=48000.0, filterType=BiquadFilterType.APF, filterf0=1000.0, filterQ=1.0, filterGain=0.0):
        self.filterf0 = filterf0
        self.filterQ = filterQ
        self.filterGain = filterGain
        self.filterType = filterType
        self.sampleRate = sampleRate


class BiquadFilter:

    def __init__(self, filterParameters=None, sampleRate=None):
        # Filter parameters (initialize to a default value if not provided)
        if filterParameters is None:
            filterParameters = BiquadFilterParameters()
        self.filterf0 = filterParameters.filterf0
        self.filterQ = filterParameters.filterQ
        self.filterGain = filterParameters.filterGain
        self.filterType = filterParameters.filterType
        self._sampleRate = filterParameters.sampleRate
        if not sampleRate is None:
            self._sampleRate = sampleRate
        # Initilize filter states
        self.x1 = 0.0
        self.x2 = 0.0
        self.y1 = 0.0
        self.y2 = 0.0
        # Generate filter coefficients
        self.filterCoefficients = BiquadFilterCoefficients()
        self.generateBiQuadCoefficients(filterType=self.filterType, filterf0=self.filterf0, filterQ=self.filterQ, filterGain=self.filterGain)


    def __str__(self):
        return "f0:%10.3f Hz, Q: %.3f, Gain: %.3f (%5.1f dB), Type: %s" % (self.filterf0, self.filterQ,
                                                                          10**(self.filterGain/20.0), self.filterGain,
                                                                          str(self.filterType))


    def processSample(self, inputSample, outputSample=None):
        # return inputSample
        # # Calculate new output sample
        # _output = (inputSample * self.filterCoefficients.b0) + \
        #           (self.x1 * self.filterCoefficients.b1) + \
        #           (self.x2 * self.filterCoefficients.b2) - \
        #           (self.y1 * self.filterCoefficients.a1) - \
        #           (self.y2 * self.filterCoefficients.a2)
        # _output = _output / self.filterCoefficients.a0
        # return _output

        if outputSample is None:
            inputSample   = float(inputSample)
            self._output  = (inputSample * self.filterCoefficients.b0)
            self._output += (self.x1 * self.filterCoefficients.b1)
            self._output += (self.x2 * self.filterCoefficients.b2)
            self._output -= (self.y1 * self.filterCoefficients.a1)
            self._output -= (self.y2 * self.filterCoefficients.a2)
            self._output /= self.filterCoefficients.a0

            # Rotate states
            self.x2 = self.x1
            self.x1 = inputSample
            self.y2 = self.y1
            self.y1 = self._output
            # Return new output
            return self._output
        else:
            # Calculate new output sample
            # inputSample[0]   = float(inputSample[0])
            outputSample[0]  = (float(inputSample[0]) * self.filterCoefficients.b0)
            outputSample[0] += (self.x1 * self.filterCoefficients.b1)
            outputSample[0] += (self.x2 * self.filterCoefficients.b2)
            outputSample[0] -= (self.y1 * self.filterCoefficients.a1)
            outputSample[0] -= (self.y2 * self.filterCoefficients.a2)
            outputSample[0] /= self.filterCoefficients.a0
            # Rotate states
            self.x2 = self.x1
            self.x1 = inputSample[0]
            self.y2 = self.y1
            self.y1 = outputSample[0]


    def generateBiQuadCoefficients(self, filterType, filterf0, filterQ, filterGain=0):
        '''
            This example generates the coefficients for a bandpass filter
            with the peak gain at 0 dB

            http://shepazu.github.io/Audio-EQ-Cookbook/audio-eq-cookbook.html
        '''

        # print("Generating filter coefficients for parameters: Type=%s, f0=%.1f [Hz], Q=%.2f, Gain=%.1f[dB]" % (str(filterType), filterf0, filterQ, filterGain))

        if filterf0 > self._sampleRate / 2.0:
            _limitFreq = 0.99 * self._sampleRate / 2.0
            print("Warning: Filter's f0 was set to %.1f [Hz], limiting it f0 to: %.1f" % (filterf0, _limitFreq))
            filterf0 = _limitFreq

        _omega0 = 2 * math.pi * filterf0 / self._sampleRate
        _alpha = math.sin(_omega0) / (2 * filterQ)

        if filterType == BiquadFilterType.BPF_NORMALIZED:
            # BPF:
            _b0 = _alpha
            _b1 = 0
            _b2 = -_alpha
            _a0 = 1 + _alpha
            _a1 = -2 * math.cos(_omega0)
            _a2 = 1 - _alpha

        elif filterType == BiquadFilterType.HPF:
            # HPF:
            _b0 = (1 + math.cos(_omega0)) / 2.0
            _b1 = - (1 + math.cos(_omega0))
            _b2 = _b0
            _a0 = 1 + _alpha
            _a1 = -2 * math.cos(_omega0)
            _a2 = 1 - _alpha

        elif filterType == BiquadFilterType.BPF:
            # BPF-2:
            _b0 = filterQ * _alpha
            _b1 = 0
            _b2 = -filterQ * _alpha
            _a0 = 1 + _alpha
            _a1 = -2 * math.cos(_omega0)
            _a2 = 1 - _alpha

        elif filterType == BiquadFilterType.LPF:
            _b0 = (1 - math.cos(_omega0)) / 2.0
            _b1 = 1 - math.cos(_omega0)
            _b2 = _b0
            _a0 = 1 + _alpha
            _a1 = -2 * math.cos(_omega0)
            _a2 = 1 - _alpha

        elif filterType == BiquadFilterType.NOTCH:
            _b0 = 1
            _b1 = -2 * math.cos(_omega0)
            _b2 = 1
            _a0 = 1 + _alpha
            _a1 = -2 * math.cos(_omega0)
            _a2 = 1 - _alpha
        
        elif filterType == BiquadFilterType.APF:
            _b0 = 1 - _alpha
            _b1 = -2 * math.cos(_omega0)
            _b2 = 1 + _alpha
            _a0 = 1 + _alpha
            _a1 = -2 * math.cos(_omega0)
            _a2 = 1 - _alpha

        elif filterType == BiquadFilterType.PEAK:
            A = 10**(filterGain / 40.0)
            _b0 = 1 + (_alpha * A)
            _b1 = -2 * math.cos(_omega0)
            _b2 = 1 - (_alpha * A)
            _a0 = 1 + (_alpha / A)
            _a1 = -2 * math.cos(_omega0)
            _a2 = 1 - (_alpha / A)

        elif filterType == BiquadFilterType.LOW_SHELVING:
            A = 10**(filterGain / 40.0)
            _b0 = A * ((A + 1) - (A - 1)*math.cos(_omega0) + 2*math.sqrt(A)*_alpha)
            _b1 = 2 * A * ((A - 1) - (A + 1)*math.cos(_omega0))
            _b2 = A * ((A + 1) - (A - 1)*math.cos(_omega0) - 2*math.sqrt(A)*_alpha)
            _a0 = (A + 1) + (A - 1)*math.cos(_omega0) + 2*math.sqrt(A)*_alpha
            _a1 = -2 * ((A - 1) + (A + 1)*math.cos(_omega0))
            _a2 = (A + 1) + (A - 1)*math.cos(_omega0) - 2*math.sqrt(A)*_alpha

        elif filterType == BiquadFilterType.HIGH_SHELVING:
            A = 10**(filterGain / 40.0)
            _b0 = A * ((A + 1) + (A - 1)*math.cos(_omega0) + 2*math.sqrt(A)*_alpha)
            _b1 = -2 * A * ((A - 1) + (A + 1)*math.cos(_omega0))
            _b2 = A * ((A + 1) + (A - 1)*math.cos(_omega0) - 2*math.sqrt(A)*_alpha)
            _a0 = (A + 1) - (A - 1)*math.cos(_omega0) + 2*math.sqrt(A)*_alpha
            _a1 = 2 * ((A - 1) - (A + 1)*math.cos(_omega0))
            _a2 = (A + 1) - (A - 1)*math.cos(_omega0) - 2*math.sqrt(A)*_alpha

        else:
            # Unknown type of filter:
            _b0 = 1
            _b1 = 0
            _b2 = 0
            _a0 = 0
            _a1 = 0
            _a2 = 0

        # Using numpy types is slower, so for real-time applications, better use floats
        self.filterCoefficients.b0 = float(_b0)
        self.filterCoefficients.b1 = float(_b1)
        self.filterCoefficients.b2 = float(_b2)
        self.filterCoefficients.a0 = float(_a0)
        self.filterCoefficients.a1 = float(_a1)
        self.filterCoefficients.a2 = float(_a2)


def computeBiquadFilterResponse(filterCoefficients, normalizedFreqBins):
    if filterCoefficients is None:
        return

    if not len(filterCoefficients) == 6:
        return

    b0 = filterCoefficients[0]
    b1 = filterCoefficients[1]
    b2 = filterCoefficients[2]
    a0 = filterCoefficients[3]
    a1 = filterCoefficients[4]
    a2 = filterCoefficients[5]

    alpha = b0
    import math
    import numpy as np
    w0 = np.abs(math.acos(-a1 / 2.0))
    Q = np.sin(w0) / (2.0 * alpha)

    e = np.exp
    pi = np.pi
    _freqResponse = [(b0 + (b1*e(-1j*2*pi*f)) + (b2 * e(-2j*2*pi*f))) / (a0 + (a1*e(-1j*2*pi*f)) + (a2*e(-2j*2*pi*f))) for f in normalizedFreqBins]

    return _freqResponse


def computeBiquadFilterIR(filterCoefficients, IRLength):
    _freqResponse = computeBiquadFilterResponse(filterCoefficients, [i / IRLength for i in range(int(IRLength))])
    return np.real(np.fft.ifft(_freqResponse))


def getFilterSpectrumModule(filter):
    _coeffs = [filter.filterCoefficients.b0, filter.filterCoefficients.b1, filter.filterCoefficients.b2, \
               filter.filterCoefficients.a0, filter.filterCoefficients.a1, filter.filterCoefficients.a2]
    # Compute complex frequency response
    sampleRate = 48000
    _normalizedFreqBins = [i/sampleRate for i in range(int(sampleRate/2.0))]
    _freqResponse = computeBiquadFilterResponse(_coeffs, _normalizedFreqBins)
    _freqResponseMod = 20.0 * np.log10(np.abs(_freqResponse))
    _freqResponsePha = np.angle(_freqResponse) / np.pi * 180.0

    _ir = computeBiquadFilterIR(_coeffs, len(_normalizedFreqBins))
    return np.multiply(_normalizedFreqBins, sampleRate), _freqResponseMod, _freqResponsePha, _ir


if __name__ == "__main__":

    import numpy as np
    import matplotlib.pyplot as plot

    # Real-time filtering of an MLS signal and IR calculation
    _filter = BiquadFilter()
    _filter.generateBiQuadCoefficients(filterf0=1000, filterQ=0.125, filterType=BiquadFilterType.BPF_NORMALIZED)
    _inputSignal = [0.0 for i in range(48000)]
    _inputSignal[0] = 1.0
    _outputSignal = [_filter.processSample(_inputSignal[i]) for i in range(len(_inputSignal))]
    _ir = _outputSignal
    _spec = 20 * np.log10(np.abs(np.fft.fft(_ir)))
    _spec = _spec[0 : round(len(_spec)/2.0)]
    _freqBins = [i / len(_spec) * 48000.0 / 2.0 for i in range(len(_spec))]
    plot.semilogx(_freqBins, _spec)
    plot.grid(1)
    plot.ylim([-60, 10])
    plot.xlim([2, 20000])
    plot.show()
    # quit()

    # Real-time filtering of an input signal
    _filter = BiquadFilter()
    _filter.generateBiQuadCoefficients(filterf0=1000, filterQ=10.0, filterType=BiquadFilterType.NOTCH)
    _inputSignal = [np.sin(2 * np.pi * 1000 * i / 48000) for i in range(1000)]
    _outputSignal = [_filter.processSample(_inputSignal[i]) for i in range(len(_inputSignal))]
    plot.plot(_inputSignal)
    plot.plot(_outputSignal)
    plot.legend(['Input', 'Output'])
    plot.grid(1)
    plot.show()
    # quit()

    # Spectrum plotting of a couple of filters
    _filter = BiquadFilter()
    _filterSpecs = [[1000, 0.3, BiquadFilterType.NOTCH], 
                    [1000, 30, BiquadFilterType.NOTCH]]
    for _specs in _filterSpecs:
        _filter.generateBiQuadCoefficients(filterf0=_specs[0], filterQ=_specs[1], filterType=_specs[2])
        _bins, _mod, _pha, _ir = getFilterSpectrumModule(_filter)
        plot.subplot(3,1,1)
        plot.semilogx(_bins, _mod)
        plot.grid(1)
        plot.ylim([-60, 10])
        plot.subplot(3,1,2)
        plot.semilogx(_bins, _pha)
        plot.ylim([-180, 180])
        plot.grid(1)
        plot.subplot(3,1,3)
        plot.plot(_ir)
        plot.grid(1)
    plot.grid(1)
    plot.show()
