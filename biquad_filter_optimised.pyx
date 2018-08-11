'''
    C++ optimised implementation of a biquad filter,
    (approximately 8.5 times faster).

    Joe Simon 2018.
'''


import math


# ----------------- Additional classes -----------------


cpdef enum BiquadFilterType:
    BPF_NORMALIZED
    HPF
    BPF
    LPF
    NOTCH
    APF
    LOW_SHELVING
    HIGH_SHELVING
    PEAK


cdef class BiquadFilterCoefficients:
    cdef public double b0, b1, b2, a0, a1, a2
    def __init__(self, double b0=1.0, double b1=0, double b2=0,
                 double a0=0, double a1=0, double a2=0):
        self.b0 = b0
        self.b1 = b1
        self.b2 = b2
        self.a0 = a0
        self.a1 = a1
        self.a2 = a2


cdef class BiquadFilterParameters:
    cdef public double filterf0, filterQ, filterGaindB, sampleRate
    cdef public BiquadFilterType filterType
    def __init__(self,
                 double sampleRate=48000.0,
                 BiquadFilterType filterType=BiquadFilterType.APF,
                 double filterf0=1000.0,
                 double filterQ=0.707,
                 double filterGaindB=0.0):
        self.filterf0     = filterf0
        self.filterQ      = filterQ
        self.filterGaindB = filterGaindB
        self.filterType   = filterType
        self.sampleRate   = sampleRate


# ------------------ Main filter part ------------------


cdef class BiquadFilter:

    cdef double _output
    cdef double _x1, _x2, _y1, _y2
    cdef BiquadFilterType _filterType
    cdef BiquadFilterCoefficients _filterCoefficients
    cdef BiquadFilterParameters _filterParameters

    def __init__(self, filterParameters=None):
        # Initialize filter parameters to default values if none provided
        self._filterParameters = BiquadFilterParameters() \
            if filterParameters is None else filterParameters
        # Initilize filter states
        self._x1 = 0.0
        self._x2 = 0.0
        self._y1 = 0.0
        self._y2 = 0.0
        self._output = 0.0
        # Generate filter coefficients
        self._filterCoefficients = BiquadFilterCoefficients()
        self.generateBiQuadCoefficients(self._filterParameters)

    def __str__(self):
        return "f0:%10.3f Hz, Q: %.3f, Gain: %.3f (%5.1f dB), Type: %s" % \
                (self._filterParameters.filterf0,
                 self._filterParameters.filterQ,
                 10**(self._filterParameters.filterGaindB/20.0),
                 self._filterParameters.filterGaindB,
                 str(self._filterParameters.filterType))

    cpdef double processSample(self, double inputSample):
        self._output  = (inputSample * self._filterCoefficients.b0)
        self._output +=    (self._x1 * self._filterCoefficients.b1)
        self._output +=    (self._x2 * self._filterCoefficients.b2)
        self._output -=    (self._y1 * self._filterCoefficients.a1)
        self._output -=    (self._y2 * self._filterCoefficients.a2)
        self._output /=                self._filterCoefficients.a0
        # Rotate states
        self._x2 = self._x1
        self._x1 = inputSample
        self._y2 = self._y1
        self._y1 = self._output
        # Return new output
        return self._output


    cpdef void generateBiQuadCoefficients(self, BiquadFilterParameters _parameters):
        # print("Generating filter coefficients for parameters: Type=%s, f0=%.1f [Hz], \
        #       Q=%.2f, Gain=%.1f[dB]" % (str(filterType), filterf0, filterQ, filterGain))
        cdef double filterf0   = _parameters.filterf0
        cdef double filterQ    = _parameters.filterQ
        cdef double filterGain = _parameters.filterGaindB
        cdef BiquadFilterType filterType = _parameters.filterType
        cdef double _b0, _b1, _b2, _a0, _a1, _a2, _alpha, _omega0, _limitFreq, A

        if filterf0 > _parameters.sampleRate / 2.0:
            _limitFreq = 0.99 * _parameters.sampleRate / 2.0
            print("Warning: Filter's f0 was set to %f [Hz], limiting it f0 to: %f due to the sample rate." % \
                  (filterf0, _limitFreq))
            filterf0 = _limitFreq

        # Pre-calculations
        _omega0 = 2 * math.pi * filterf0 / _parameters.sampleRate
        _alpha = math.sin(_omega0) / (2 * filterQ)

        # Coefficients calculation
        if filterType == BiquadFilterType.BPF_NORMALIZED:
            # BPF-Normalized:
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
            # BPF:
            _b0 = filterQ * _alpha
            _b1 = 0
            _b2 = -filterQ * _alpha
            _a0 = 1 + _alpha
            _a1 = -2 * math.cos(_omega0)
            _a2 = 1 - _alpha

        elif filterType == BiquadFilterType.LPF:
            # LPF:
            _b0 = (1 - math.cos(_omega0)) / 2.0
            _b1 = 1 - math.cos(_omega0)
            _b2 = _b0
            _a0 = 1 + _alpha
            _a1 = -2 * math.cos(_omega0)
            _a2 = 1 - _alpha

        elif filterType == BiquadFilterType.NOTCH:
            # NOTCH:
            _b0 = 1
            _b1 = -2 * math.cos(_omega0)
            _b2 = 1
            _a0 = 1 + _alpha
            _a1 = -2 * math.cos(_omega0)
            _a2 = 1 - _alpha
        
        elif filterType == BiquadFilterType.APF:
            # APF:
            _b0 = 1 - _alpha
            _b1 = -2 * math.cos(_omega0)
            _b2 = 1 + _alpha
            _a0 = 1 + _alpha
            _a1 = -2 * math.cos(_omega0)
            _a2 = 1 - _alpha

        elif filterType == BiquadFilterType.PEAK:
            # PEAK:
            A = 10**(filterGain / 40.0)
            _b0 = 1 + (_alpha * A)
            _b1 = -2 * math.cos(_omega0)
            _b2 = 1 - (_alpha * A)
            _a0 = 1 + (_alpha / A)
            _a1 = -2 * math.cos(_omega0)
            _a2 = 1 - (_alpha / A)

        elif filterType == BiquadFilterType.LOW_SHELVING:
            # LO-SHELVING:
            A = 10**(filterGain / 40.0)
            _b0 = A * ((A + 1) - (A - 1)*math.cos(_omega0) + 2*math.sqrt(A)*_alpha)
            _b1 = 2 * A * ((A - 1) - (A + 1)*math.cos(_omega0))
            _b2 = A * ((A + 1) - (A - 1)*math.cos(_omega0) - 2*math.sqrt(A)*_alpha)
            _a0 = (A + 1) + (A - 1)*math.cos(_omega0) + 2*math.sqrt(A)*_alpha
            _a1 = -2 * ((A - 1) + (A + 1)*math.cos(_omega0))
            _a2 = (A + 1) + (A - 1)*math.cos(_omega0) - 2*math.sqrt(A)*_alpha

        elif filterType == BiquadFilterType.HIGH_SHELVING:
            # HI-SHELVING:
            A = 10**(filterGain / 40.0)
            _b0 = A * ((A + 1) + (A - 1)*math.cos(_omega0) + 2*math.sqrt(A)*_alpha)
            _b1 = -2 * A * ((A - 1) + (A + 1)*math.cos(_omega0))
            _b2 = A * ((A + 1) + (A - 1)*math.cos(_omega0) - 2*math.sqrt(A)*_alpha)
            _a0 = (A + 1) - (A - 1)*math.cos(_omega0) + 2*math.sqrt(A)*_alpha
            _a1 = 2 * ((A - 1) - (A + 1)*math.cos(_omega0))
            _a2 = (A + 1) - (A - 1)*math.cos(_omega0) - 2*math.sqrt(A)*_alpha

        else:
            # Unknown type of filter:
            _b0 = 1.0
            _b1 = 0.0
            _b2 = 0.0
            _a0 = 0.0
            _a1 = 0.0
            _a2 = 0.0

        self._filterCoefficients.b0 = _b0
        self._filterCoefficients.b1 = _b1
        self._filterCoefficients.b2 = _b2
        self._filterCoefficients.a0 = _a0
        self._filterCoefficients.a1 = _a1
        self._filterCoefficients.a2 = _a2


# --------------- Additional filter tools --------------


    cpdef computeComplexBiquadFilterResponse(self, normalizedFreqBins):
        cdef complex b0 = self._filterCoefficients.b0
        cdef complex b1 = self._filterCoefficients.b1
        cdef complex b2 = self._filterCoefficients.b2
        cdef complex a0 = self._filterCoefficients.a0
        cdef complex a1 = self._filterCoefficients.a1
        cdef complex a2 = self._filterCoefficients.a2
        cdef double w0  = abs(math.acos(-self._filterCoefficients.a1 / 2.0))
        cdef complex Q  = math.sin(w0) / (2.0 * self._filterCoefficients.b0)
        cdef complex pi = math.pi
        cdef complex f
        import cmath
        return [(b0 + (b1*cmath.exp(-1j*2*pi*f)) + (b2 * cmath.exp(-2j*2*pi*f))) / (a0 + (a1*cmath.exp(-1j*2*pi*f)) + (a2*cmath.exp(-2j*2*pi*f))) for f in normalizedFreqBins]
