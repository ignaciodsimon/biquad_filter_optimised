'''
    C++ implementation of a biquad filter, with a
    useful calculator function.

    https://github.com/ignaciodsimon/biquad_filter_optimised

    Joe Simon 2018 - 2023.
'''


import math


# ----------------- Additional classes -----------------

cpdef enum BiquadFilterType:
    HPF
    LPF
    HPF_FIRST_ORDER
    LPF_FIRST_ORDER
    BPF
    BPF_NORMALIZED
    PEAK
    NOTCH
    APF
    APF_FIRST_ORDER
    LOW_SHELVING
    LOW_SHELVING_FIRST_ORDER
    HIGH_SHELVING
    HIGH_SHELVING_FIRST_ORDER


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
    cdef public double filterf0, filterQ, filterGaindB, sampleRate, filterMasterGaindB
    cdef public BiquadFilterType filterType
    def __init__(self,
                 double           sampleRate   = 48000.0,
                 BiquadFilterType filterType   = BiquadFilterType.APF,
                 double           filterf0     = 1000.0,
                 double           filterQ      = 0.707,
                 double           filterGaindB = 0.0,
                 double           filterMasterGaindB = 0.0):
        self.filterf0     = filterf0
        self.filterQ      = filterQ
        self.filterGaindB = filterGaindB
        self.filterType   = filterType
        self.sampleRate   = sampleRate
        self.filterMasterGaindB = filterMasterGaindB


# ------------------ Main filter part ------------------


cdef class BiquadFilter:

    cdef double _output
    cdef double _s1, _s2
    cdef BiquadFilterType _filterType
    cdef BiquadFilterCoefficients _filterCoefficients
    cdef BiquadFilterParameters _filterParameters
    cdef double _masterGain
    cdef int _isSecondOrder

    def version(self):
        return "1.0.2"

    def __init__(self, filterParameters=None):
        # Initialize filter parameters to default values if none provided
        self._filterParameters = BiquadFilterParameters() \
            if filterParameters is None else filterParameters
        # Initilize filter states
        self._s1 = 0.0
        self._s2 = 0.0
        self._output = 0.0
        self._masterGain = 1.0
        # Generate filter coefficients
        self._filterCoefficients = BiquadFilterCoefficients()
        self.generateBiQuadCoefficients(self._filterParameters)

    def resetStates(self):
        self._s1 = 0.0
        self._s2 = 0.0

    def setCoefficients(self, filterParameters, isSecondOrder):
        if not isinstance(filterParameters, BiquadFilterCoefficients):
            raise "The input parameter 'filterParameters' must be of type <BiquadFilterCoefficients>!"
        if not isinstance(isSecondOrder, bool):
            raise "The input parameter 'isSecondOrder' must be of type <bool>!"

        if not filterParameters.a0 == 0.0:
            self._filterCoefficients.b0 = filterParameters.b0 / filterParameters.a0
            self._filterCoefficients.b1 = filterParameters.b1 / filterParameters.a0
            self._filterCoefficients.b2 = filterParameters.b2 / filterParameters.a0
            self._filterCoefficients.a0 = 1.0
            self._filterCoefficients.a1 = filterParameters.a1 / filterParameters.a0
            self._filterCoefficients.a2 = filterParameters.a2 / filterParameters.a0
            self._isSecondOrder = isSecondOrder
        else:
            self._filterCoefficients.b0 = 0.0
            self._filterCoefficients.b1 = 0.0
            self._filterCoefficients.b2 = 0.0
            self._filterCoefficients.a0 = 0.0
            self._filterCoefficients.a1 = 0.0
            self._filterCoefficients.a2 = 0.0
            self._isSecondOrder = False


    def getCoefficients(self):
        return self._filterCoefficients

    def _filterTypeToString(self, filterType):
        if filterType == BiquadFilterType.BPF_NORMALIZED:
            return "BPF_NORMALIZED"
        if filterType == BiquadFilterType.HPF:
            return "HPF"
        if filterType == BiquadFilterType.BPF:
            return "BPF"
        if filterType == BiquadFilterType.LPF:
            return "LPF"
        if filterType == BiquadFilterType.NOTCH:
            return "NOTCH"
        if filterType == BiquadFilterType.APF:
            return "APF"
        if filterType == BiquadFilterType.LOW_SHELVING:
            return "LOW_SHELVING"
        if filterType == BiquadFilterType.HIGH_SHELVING:
            return "HIGH_SHELVING"
        if filterType == BiquadFilterType.PEAK:
            return "PEAK"
        if filterType == BiquadFilterType.LPF_FIRST_ORDER:
            return "LPF_FIRST_ORDER"
        if filterType == BiquadFilterType.HPF_FIRST_ORDER:
            return "HPF_FIRST_ORDER"
        if filterType == BiquadFilterType.APF_FIRST_ORDER:
            return "APF_FIRST_ORDER"
        if filterType == BiquadFilterType.LOW_SHELVING_FIRST_ORDER:
            return "LOW_SHELVING_FIRST_ORDER"
        if filterType == BiquadFilterType.HIGH_SHELVING_FIRST_ORDER:
            return "HIGH_SHELVING_FIRST_ORDER"
        return "Unknown type"

    def __str__(self):
        return "f0:%.3f Hz, Q: %.3f, Gain: %.3f (%.1f dB), Type: %s" % \
                (self._filterParameters.filterf0,
                 self._filterParameters.filterQ,
                 10**(self._filterParameters.filterGaindB/20.0),
                 self._filterParameters.filterGaindB,
                 self._filterTypeToString(self._filterParameters.filterType))

    cpdef double processSample(self, double inputSample):
        cdef double x = inputSample
        cdef double y = self._filterCoefficients.b0 * x + self._s1
        self._s1 = self._filterCoefficients.b1 * x - self._filterCoefficients.a1 * y + self._s2
        self._s2 = self._filterCoefficients.b2 * x - self._filterCoefficients.a2 * y
        y *= self._masterGain
        self._output = y
        return y

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
        if filterType == BiquadFilterType.APF_FIRST_ORDER:
            K = math.tan(math.pi * filterf0 / _parameters.sampleRate)
            _b0 = (1 - K) / (1 + K)
            _b1 = -1
            _b2 = 0
            _a0 =  1.0
            _a1 = -_b0
            _a2 = 0
            self._isSecondOrder = 0

        elif filterType == BiquadFilterType.HPF_FIRST_ORDER:
            # First-order high-pass
            _a0 =  1.0 + (_omega0 / 2.0)
            _a1 = -2.0 + _a0
            _a2 =  0.0
            _b0 =  1.0
            _b1 = -1.0
            _b2 =  0.0
            self._isSecondOrder = 0
            # print("f0: %f\nSample rate: %f\nB0: %f\nB1: %f\nB2: %f\nA0: %f\nA1: %f\nA2: %f" % (filterf0, _parameters.sampleRate, _b0, _b1, _b2, _a0, _a1, _a2))

        elif filterType == BiquadFilterType.LPF_FIRST_ORDER:
            # First-order low-pass
            _a0 = 1.0
            _a1 = -math.exp(-2.0 * math.pi * (filterf0 / _parameters.sampleRate))
            _a2 = 0.0
            _b0 = 1.0 + _a1;
            _b1 = 0.0
            _b2 = 0.0
            self._isSecondOrder = 0
            # print("f0: %f\nSample rate: %f\nB0: %f\nB1: %f\nB2: %f\nA0: %f\nA1: %f\nA2: %f" % (filterf0, _parameters.sampleRate, _b0, _b1, _b2, _a0, _a1, _a2))

        elif filterType == BiquadFilterType.BPF_NORMALIZED:
            # BPF-Normalized:
            _b0 = _alpha
            _b1 = 0
            _b2 = -_alpha
            _a0 = 1 + _alpha
            _a1 = -2 * math.cos(_omega0)
            _a2 = 1 - _alpha
            self._isSecondOrder = 1

        elif filterType == BiquadFilterType.HPF:
            # HPF:
            _b0 = (1 + math.cos(_omega0)) / 2.0
            _b1 = - (1 + math.cos(_omega0))
            _b2 = _b0
            _a0 = 1 + _alpha
            _a1 = -2 * math.cos(_omega0)
            _a2 = 1 - _alpha
            self._isSecondOrder = 1

        elif filterType == BiquadFilterType.BPF:
            # BPF:
            _b0 = filterQ * _alpha
            _b1 = 0
            _b2 = -filterQ * _alpha
            _a0 = 1 + _alpha
            _a1 = -2 * math.cos(_omega0)
            _a2 = 1 - _alpha
            self._isSecondOrder = 1

        elif filterType == BiquadFilterType.LPF:
            # LPF:
            _b0 = (1 - math.cos(_omega0)) / 2.0
            _b1 = 1 - math.cos(_omega0)
            _b2 = _b0
            _a0 = 1 + _alpha
            _a1 = -2 * math.cos(_omega0)
            _a2 = 1 - _alpha
            self._isSecondOrder = 1

        elif filterType == BiquadFilterType.NOTCH:
            # NOTCH:
            _b0 = 1
            _b1 = -2 * math.cos(_omega0)
            _b2 = 1
            _a0 = 1 + _alpha
            _a1 = -2 * math.cos(_omega0)
            _a2 = 1 - _alpha
            self._isSecondOrder = 1
        
        elif filterType == BiquadFilterType.APF:
            # APF:
            _b0 = 1 - _alpha
            _b1 = -2 * math.cos(_omega0)
            _b2 = 1 + _alpha
            _a0 = 1 + _alpha
            _a1 = -2 * math.cos(_omega0)
            _a2 = 1 - _alpha
            self._isSecondOrder = 1

        elif filterType == BiquadFilterType.PEAK:
            # PEAK:
            A = 10**(filterGain / 40.0)
            _b0 = 1 + (_alpha * A)
            _b1 = -2 * math.cos(_omega0)
            _b2 = 1 - (_alpha * A)
            _a0 = 1 + (_alpha / A)
            _a1 = -2 * math.cos(_omega0)
            _a2 = 1 - (_alpha / A)
            self._isSecondOrder = 1

        elif filterType == BiquadFilterType.LOW_SHELVING:
            # LO-SHELVING:
            A = 10**(filterGain / 40.0)
            _b0 = A * ((A + 1) - (A - 1)*math.cos(_omega0) + 2*math.sqrt(A)*_alpha)
            _b1 = 2 * A * ((A - 1) - (A + 1)*math.cos(_omega0))
            _b2 = A * ((A + 1) - (A - 1)*math.cos(_omega0) - 2*math.sqrt(A)*_alpha)
            _a0 = (A + 1) + (A - 1)*math.cos(_omega0) + 2*math.sqrt(A)*_alpha
            _a1 = -2 * ((A - 1) + (A + 1)*math.cos(_omega0))
            _a2 = (A + 1) + (A - 1)*math.cos(_omega0) - 2*math.sqrt(A)*_alpha
            self._isSecondOrder = 1

        elif filterType == BiquadFilterType.HIGH_SHELVING:
            # HI-SHELVING:
            A = 10**(filterGain / 40.0)
            _b0 = A * ((A + 1) + (A - 1)*math.cos(_omega0) + 2*math.sqrt(A)*_alpha)
            _b1 = -2 * A * ((A - 1) + (A + 1)*math.cos(_omega0))
            _b2 = A * ((A + 1) + (A - 1)*math.cos(_omega0) - 2*math.sqrt(A)*_alpha)
            _a0 = (A + 1) - (A - 1)*math.cos(_omega0) + 2*math.sqrt(A)*_alpha
            _a1 = 2 * ((A - 1) - (A + 1)*math.cos(_omega0))
            _a2 = (A + 1) - (A - 1)*math.cos(_omega0) - 2*math.sqrt(A)*_alpha
            self._isSecondOrder = 1

        elif filterType == BiquadFilterType.HIGH_SHELVING_FIRST_ORDER:
            g2 = 10.0**(filterGain/ 20.0)
            k2 = math.sqrt(g2)
            B0 = g2
            B1 = g2 * 2.0 * math.pi * filterf0 / k2
            B2 = 0.0
            A0 = 1.0
            A1 = 2 * math.pi * filterf0 * k2
            A2 = 0.0
            C = 2.0 * _parameters.sampleRate # no pre-warping
            CC = C * C
            _a0 = (A2 + A1 * C + A0 * CC) / (4 * CC)
            _a1 = (A2 - A0 * CC) / (2 * CC * _a0)
            _a2 = (A2 - A1 * C + A0 * CC) / (4 * CC * _a0)
            _b0 = (B2 + B1 * C + B0 * CC) / (4 * CC * _a0)
            _b1 = (B2 - B0 * CC) / (2 * CC * _a0)
            _b2 = (B2 - B1 * C + B0 * CC) / (4 * CC * _a0)
            _a0 = 1.0
            self._isSecondOrder = 1

        elif filterType == BiquadFilterType.LOW_SHELVING_FIRST_ORDER:
            g1 = 10.0**(filterGain / 20.0)
            k1 = math.sqrt(g1)
            B0 = 1.0
            B1 = 2 * math.pi * filterf0 * k1
            B2 = 0.0
            A0 = 1.0
            A1 = 2.0 * math.pi * filterf0 / k1
            A2 = 0.0
            C = 2.0 * _parameters.sampleRate # no pre-warping
            CC = C * C
            _a0 = (A2 + A1 * C + A0 * CC) / (4 * CC)
            _a1 = (A2 - A0 * CC) / (2 * CC * _a0)
            _a2 = (A2 - A1 * C + A0 * CC) / (4 * CC * _a0)
            _b0 = (B2 + B1 * C + B0 * CC) / (4 * CC * _a0)
            _b1 = (B2 - B0 * CC) / (2 * CC * _a0)
            _b2 = (B2 - B1 * C + B0 * CC) / (4 * CC * _a0)
            _a0 = 1.0
            self._isSecondOrder = 1
        else:
            # Unknown type of filter
            _b0 = 1.0
            _b1 = 0.0
            _b2 = 0.0
            _a0 = 1.0
            _a1 = 0.0
            _a2 = 0.0
            self._isSecondOrder = 0

        # Just a sanity check
        if not _a0 == 0.0:
            self._filterCoefficients.b0 = _b0 / _a0
            self._filterCoefficients.b1 = _b1 / _a0
            self._filterCoefficients.b2 = _b2 / _a0
            self._filterCoefficients.a0 = 1.0
            self._filterCoefficients.a1 = _a1 / _a0
            self._filterCoefficients.a2 = _a2 / _a0
        else:
            self._filterCoefficients.b0 = 1.0
            self._filterCoefficients.b1 = 0.0
            self._filterCoefficients.b2 = 0.0
            self._filterCoefficients.a0 = 1.0
            self._filterCoefficients.a1 = 0.0
            self._filterCoefficients.a2 = 0.0

        # The master gain is not applied to the coefficients so that getting / setting them works independently
        # from the "output volume setting"
        self._masterGain = 10**(_parameters.filterMasterGaindB / 20.0)

    cpdef reachStationaryState(self, double stationaryInput):
        cdef double den = (1.0 + self._filterCoefficients.a1 + self._filterCoefficients.a2)
        if den == 0.0:
            return
        cdef double DC_gain = stationaryInput * (self._filterCoefficients.b0 + self._filterCoefficients.b1 + self._filterCoefficients.b2) / den
        self._s1 = DC_gain - self._filterCoefficients.b0 * stationaryInput
        self._s2 = self._s1 - self._filterCoefficients.b1 * stationaryInput + self._filterCoefficients.a1 * DC_gain

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

    cpdef calculateCascadedQ(self, combinedOrder):
        '''
            Calculates the Q factor of cascaded filters
            to produce a Butterworth type of response.
        '''
        _pairs = combinedOrder >> 1
        _poleIncrement = math.pi / combinedOrder
        _firstAngle = _poleIncrement

        _qValues = []
        if not combinedOrder & 1: # If a first-order filter is needed
            _firstAngle /= 2
        else:
            _qValues.append([1.0, 0.5])

        for idx in range(_pairs):
            _qValues.append([2.0, 1.0 / (2.0 * math.cos(_firstAngle + idx * _poleIncrement))])

        return _qValues
