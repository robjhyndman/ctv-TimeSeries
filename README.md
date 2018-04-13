CRAN Task View: Time Series Analysis
------------------------------------


|-----------------|----------------------------------------------  
| **Maintainer:** | Rob J Hyndman                                  
| **Contact:**    | Rob.Hyndman at monash.edu                      
| **Version:**    | 2018-04-03                                     
| **URL:**        | <https://CRAN.R-project.org/view=TimeSeries>   

Base R ships with a lot of functionality useful for time series, in
particular in the stats package. This is complemented by many packages
on CRAN, which are briefly summarized below. There is also a
considerable overlap between the tools for time series and those in the
[Econometrics](https://cran.r-project.org/web/views/Econometrics.html)
and [Finance](https://cran.r-project.org/web/views/Finance.html) task
views. The packages in this view can be roughly structured into the
following topics. If you think that some package is missing from the
list, please let us know.

**Basics**

-   *Infrastructure* : Base R contains substantial infrastructure for
    representing and analyzing time series data. The fundamental class
    is `"ts"` that can represent regularly spaced time series (using
    numeric time stamps). Hence, it is particularly well-suited for
    annual, monthly, quarterly data, etc.
-   *Rolling statistics* : Moving averages are computed by `ma` from
    [forecast](https://cran.r-project.org/package=forecast), and `rollmean` from
    [zoo](https://cran.r-project.org/package=zoo). The latter also provides a
    general function `rollapply`, along with other specific rolling
    statistics functions. [roll](https://cran.r-project.org/package=roll) provides
    parallel functions for computing rolling statistics.
-   *Graphics* : Time series plots are obtained with `plot()` applied to
    `ts` objects. (Partial) autocorrelation functions plots are
    implemented in `acf()` and `pacf()`. Alternative versions are
    provided by `Acf()` and `Pacf()` in
    [forecast](https://cran.r-project.org/package=forecast), along with a
    combination display using `tsdisplay()`.
    [SDD](https://cran.r-project.org/package=SDD) provides more general serial
    dependence diagrams, while [dCovTS](https://cran.r-project.org/package=dCovTS)
    computes and plots the distance covariance and correlation functions
    of time series. Seasonal displays are obtained using `monthplot()`
    in stats and `seasonplot` in
    [forecast](https://cran.r-project.org/package=forecast).
    [Wats](https://cran.r-project.org/package=Wats) implements wrap-around time
    series graphics. [ggseas](https://cran.r-project.org/package=ggseas) provides
    ggplot2 graphics for seasonally adjusted series and rolling
    statistics. [dygraphs](https://cran.r-project.org/package=dygraphs) provides an
    interface to the Dygraphs interactive time series charting library.
    [TSstudio](https://cran.r-project.org/package=TSstudio) provides some
    interactive visualization tools for time series.
    [ZRA](https://cran.r-project.org/package=ZRA) plots forecast objects from the
    [forecast](https://cran.r-project.org/package=forecast) package using dygraphs.
    Basic fan plots of forecast distributions are provided by
    [forecast](https://cran.r-project.org/package=forecast) and
    [vars](https://cran.r-project.org/package=vars). More flexible fan plots of any
    sequential distributions are implemented in
    [fanplot](https://cran.r-project.org/package=fanplot).

**Times and Dates**

-   Class `"ts"` can only deal with numeric time stamps, but many more
    classes are available for storing time/date information and
    computing with it. For an overview see *R Help Desk: Date and Time
    Classes in R* by Gabor Grothendieck and Thomas Petzoldt in [R News
    4(1)](http://CRAN.R-project.org/doc/Rnews/Rnews_2004-1.pdf), 29-32.
-   Classes `"yearmon"` and `"yearqtr"` from
    [zoo](https://cran.r-project.org/package=zoo) allow for more convenient
    computation with monthly and quarterly observations, respectively.
-   Class `"Date"` from the base package is the basic class for dealing
    with dates in daily data. The dates are internally stored as the
    number of days since 1970-01-01.
-   The [chron](https://cran.r-project.org/package=chron) package provides classes
    for `dates()`, `hours()` and date/time (intra-day) in `chron()`.
    There is no support for time zones and daylight savings time.
    Internally, `"chron"` objects are (fractional) days since
    1970-01-01.
-   Classes `"POSIXct"` and `"POSIXlt"` implement the POSIX standard for
    date/time (intra-day) information and also support time zones and
    daylight savings time. However, the time zone computations require
    some care and might be system-dependent. Internally, `"POSIXct"`
    objects are the number of seconds since 1970-01-01 00:00:00 GMT.
    Package [lubridate](https://cran.r-project.org/package=lubridate) provides
    functions that facilitate certain POSIX-based computations.
    [wktmo](https://cran.r-project.org/package=wktmo) converts weekly data to
    monthly data in several different ways.
-   Several packages aim to handle time-based tibbles:
    [tsibble](https://cran.r-project.org/package=tsibble) provides tidy temporal
    data frames and associated tools;
    [tibbletime](https://cran.r-project.org/package=tibbletime) handles time aware
    tibbles; [timetk](https://cran.r-project.org/package=timetk) contains tools for
    working with and coercing between time-based tibbles, xts, zoo and
    ts objects.
-   Class `"timeDate"` is provided in the
    [timeDate](https://cran.r-project.org/package=timeDate) package (previously:
    fCalendar). It is aimed at financial time/date information and deals
    with time zones and daylight savings times via a new concept of
    "financial centers". Internally, it stores all information in
    `"POSIXct"` and does all computations in GMT only. Calendar
    functionality, e.g., including information about weekends and
    holidays for various stock exchanges, is also included.
-   The [tis](https://cran.r-project.org/package=tis) package provides the `"ti"`
    class for time/date information.
-   The `"mondate"` class from the
    [mondate](https://cran.r-project.org/package=mondate) package facilitates
    computing with dates in terms of months.
-   The [tempdisagg](https://cran.r-project.org/package=tempdisagg) package includes
    methods for temporal disaggregation and interpolation of a low
    frequency time series to a higher frequency series.
-   Time series disaggregation is also provided by
    [tsdisagg2](https://cran.r-project.org/package=tsdisagg2).
-   [TimeProjection](https://cran.r-project.org/package=TimeProjection) extracts
    useful time components of a date object, such as day of week,
    weekend, holiday, day of month, etc, and put it in a data frame.

**Time Series Classes**

-   As mentioned above, `"ts"` is the basic class for regularly spaced
    time series using numeric time stamps.
-   The [zoo](https://cran.r-project.org/package=zoo) package provides
    infrastructure for regularly and irregularly spaced time series
    using arbitrary classes for the time stamps (i.e., allowing all
    classes from the previous section). It is designed to be as
    consistent as possible with `"ts"`. Coercion from and to `"zoo"` is
    available for all other classes mentioned in this section.
-   The package [xts](https://cran.r-project.org/package=xts) is based on
    [zoo](https://cran.r-project.org/package=zoo) and provides uniform handling of
    R's different time-based data classes.
-   Various packages implement irregular time series based on
    `"POSIXct"` time stamps, intended especially for financial
    applications. These include `"irts"` from
    [tseries](https://cran.r-project.org/package=tseries), and `"fts"` from
    [fts](https://cran.r-project.org/package=fts).
-   The class `"timeSeries"` in
    [timeSeries](https://cran.r-project.org/package=timeSeries) (previously:
    fSeries) implements time series with `"timeDate"` time stamps.
-   The class `"tis"` in [tis](https://cran.r-project.org/package=tis) implements
    time series with `"ti"` time stamps.
-   The package [tframe](https://cran.r-project.org/package=tframe) contains
    infrastructure for setting time frames in different formats.

**Forecasting and Univariate Modeling**

-   The [forecast](https://cran.r-project.org/package=forecast) package provides a
    class and methods for univariate time series forecasts, and provides
    many functions implementing different forecasting models including
    all those in the stats package.
-   *Exponential smoothing* : `HoltWinters()` in stats provides some
    basic models with partial optimization, `ets()` from the
    [forecast](https://cran.r-project.org/package=forecast) package provides a
    larger set of models and facilities with full optimization.
    [robets](https://cran.r-project.org/package=robets) provides a robust
    alternative to the `ets()` function.
    [smooth](https://cran.r-project.org/package=smooth) implements some
    generalizations of exponential smoothing. The
    [MAPA](https://cran.r-project.org/package=MAPA) package combines exponential
    smoothing models at different levels of temporal aggregation to
    improve forecast accuracy.
-   [prophet](https://cran.r-project.org/package=prophet) forecasts time series
    based on an additive model where non-linear trends are fit with
    yearly and weekly seasonality, plus holidays. It works best with
    daily data.
-   The theta method is implemented in the `thetaf` function from the
    [forecast](https://cran.r-project.org/package=forecast) package. An alternative
    and extended implementation is provided in
    [forecTheta](https://cran.r-project.org/package=forecTheta).
-   *Autoregressive models* : `ar()` in stats (with model selection) and
    [FitAR](https://cran.r-project.org/package=FitAR) for subset AR models.
-   *ARIMA models* : `arima()` in stats is the basic function for ARIMA,
    SARIMA, ARIMAX, and subset ARIMA models. It is enhanced in the
    [forecast](https://cran.r-project.org/package=forecast) package via the function
    `Arima()` along with `auto.arima()` for automatic order selection.
    `arma()` in the [tseries](https://cran.r-project.org/package=tseries) package
    provides different algorithms for ARMA and subset ARMA models. Other
    estimation methods including the innovations algorithm are provided
    by [itsmr](https://cran.r-project.org/package=itsmr).
    [FitARMA](https://cran.r-project.org/package=FitARMA) implements a fast MLE
    algorithm for ARMA models. Package
    [gsarima](https://cran.r-project.org/package=gsarima) contains functionality for
    Generalized SARIMA time series simulation. Robust ARIMA modeling is
    provided in the [robustarima](https://cran.r-project.org/package=robustarima)
    package. The [mar1s](https://cran.r-project.org/package=mar1s) package handles
    multiplicative AR(1) with seasonal processes.
    [TSTutorial](https://cran.r-project.org/package=TSTutorial) provides an
    interactive tutorial for Box-Jenkins modelling. Improved prediction
    intervals for ARIMA and structural time series models are provided
    by [tsPI](https://cran.r-project.org/package=tsPI).
-   *Periodic ARMA models* : [pear](https://cran.r-project.org/package=pear) and
    [partsm](https://cran.r-project.org/package=partsm) for periodic autoregressive
    time series models, and [perARMA](https://cran.r-project.org/package=perARMA)
    for periodic ARMA modelling and other procedures for periodic time
    series analysis.
-   *ARFIMA models* : Some facilities for fractional differenced ARFIMA
    models are provided in the
    [fracdiff](https://cran.r-project.org/package=fracdiff) package. The
    [arfima](https://cran.r-project.org/package=arfima) package has more advanced
    and general facilities for ARFIMA and ARIMA models, including
    dynamic regression (transfer function) models. Fractional Gaussian
    noise and simple models for hyperbolic decay time series are handled
    in the [FGN](https://cran.r-project.org/package=FGN) package.
-   *Transfer function* models are provided by the `arimax` function in
    the [TSA](https://cran.r-project.org/package=TSA) package, and the `arfima`
    function in the [arfima](https://cran.r-project.org/package=arfima) package.
-   Outlier detection following the Chen-Liu approach is provided by
    [tsoutliers](https://cran.r-project.org/package=tsoutliers).
-   *Structural models* are implemented in `StructTS()` in stats, and in
    [stsm](https://cran.r-project.org/package=stsm) and
    [stsm.class](https://cran.r-project.org/package=stsm.class).
    [KFKSDS](https://cran.r-project.org/package=KFKSDS) provides a naive
    implementation of the Kalman filter and smoothers for univariate
    state space models. Bayesian structural time series models are
    implemented in [bsts](https://cran.r-project.org/package=bsts)
-   Non-Gaussian time series can be handled with GLARMA state space
    models via [glarma](https://cran.r-project.org/package=glarma), and using
    Generalized Autoregressive Score models in the
    [GAS](https://cran.r-project.org/package=GAS) package. Conditional
    Auto-Regression models using Monte Carlo Likelihood methods are
    implemented in [mclcar](https://cran.r-project.org/package=mclcar).
-   *GARCH models* : `garch()` from
    [tseries](https://cran.r-project.org/package=tseries) fits basic GARCH models.
    Many variations on GARCH models are provided by
    [rugarch](https://cran.r-project.org/package=rugarch). Other univariate GARCH
    packages include [fGarch](https://cran.r-project.org/package=fGarch) which
    implements ARIMA models with a wide class of GARCH innovations.
    There are many more GARCH packages described in the
    [Finance](https://cran.r-project.org/web/views/Finance.html) task
    view.
-   *Stochastic volatility* models are handled by
    [stochvol](https://cran.r-project.org/package=stochvol) in a Bayesian framework.
-   *Count time series* models are handled in the
    [tscount](https://cran.r-project.org/package=tscount) and
    [acp](https://cran.r-project.org/package=acp) packages.
    [ZIM](https://cran.r-project.org/package=ZIM) provides for Zero-Inflated Models
    for count time series.
    [tsintermittent](https://cran.r-project.org/package=tsintermittent) implements
    various models for analysing and forecasting intermittent demand
    time series.
-   *Censored time series* can be modelled using
    [cents](https://cran.r-project.org/package=cents) and
    [carx](https://cran.r-project.org/package=carx).
-   *Portmanteau tests* are provided via `Box.test()` in the stats
    package. Additional tests are given by
    [portes](https://cran.r-project.org/package=portes) and
    [WeightedPortTest](https://cran.r-project.org/package=WeightedPortTest).
-   *Change point detection* is provided in
    [strucchange](https://cran.r-project.org/package=strucchange) (using linear
    regression models), and in [trend](https://cran.r-project.org/package=trend)
    (using nonparametric tests). The
    [changepoint](https://cran.r-project.org/package=changepoint) package provides
    many popular changepoint methods, and
    [ecp](https://cran.r-project.org/package=ecp) does nonparametric changepoint
    detection for univariate and multivariate series.
    [InspectChangepoint](https://cran.r-project.org/package=InspectChangepoint) uses
    sparse projection to estimate changepoints in high-dimensional time
    series.
-   Tests for possibly non-monotonic trends are provided by
    [funtimes](https://cran.r-project.org/package=funtimes).
-   *Time series imputation* is provided by the
    [imputeTS](https://cran.r-project.org/package=imputeTS) package. Some more
    limited facilities are available using `na.interp()` from the
    [forecast](https://cran.r-project.org/package=forecast) package.
    [mtsdi](https://cran.r-project.org/package=mtsdi) implements an EM algorithm for
    imputing missing values in multivariate normal time series,
    accounting for spatial and temporal correlations.
-   Forecasts can be combined using
    [ForecastComb](https://cran.r-project.org/package=ForecastComb) which supports
    many forecast combination methods including simple, geometric and
    regression-based combinations.
    [forecastHybrid](https://cran.r-project.org/package=forecastHybrid) provides
    functions for ensemble forecasts, combining approaches from the
    [forecast](https://cran.r-project.org/package=forecast) package.
    [opera](https://cran.r-project.org/package=opera) has facilities for online
    predictions based on combinations of forecasts provided by the user.
    [mafs](https://cran.r-project.org/package=mafs) fits several forecast models and
    selects the best one according to an error metric.
-   Forecast evaluation is provided in the `accuracy()` function from
    [forecast](https://cran.r-project.org/package=forecast). Distributional forecast
    evaluation using scoring rules is available in
    [scoringRules](https://cran.r-project.org/package=scoringRules)
-   Tidy tools for forecasting are provided by
    [sweep](https://cran.r-project.org/package=sweep), converting objects produced
    in [forecast](https://cran.r-project.org/package=forecast) to "tidy" data
    frames.
-   *Miscellaneous* : [ltsa](https://cran.r-project.org/package=ltsa) contains
    methods for linear time series analysis,
    [timsac](https://cran.r-project.org/package=timsac) for time series analysis and
    control, and [tsbugs](https://cran.r-project.org/package=tsbugs) for time series
    BUGS models.

**Frequency analysis**

-   *Spectral density estimation* is provided by `spectrum()` in the
    stats package, including the periodogram, smoothed periodogram and
    AR estimates. Bayesian spectral inference is provided by
    [bspec](https://cran.r-project.org/package=bspec).
    [quantspec](https://cran.r-project.org/package=quantspec) includes methods to
    compute and plot Laplace periodograms for univariate time series.
    The Lomb-Scargle periodogram for unevenly sampled time series is
    computed by [lomb](https://cran.r-project.org/package=lomb).
    [spectral](https://cran.r-project.org/package=spectral) uses Fourier and Hilbert
    transforms for spectral filtering. [psd](https://cran.r-project.org/package=psd)
    produces adaptive, sine-multitaper spectral density estimates.
    [kza](https://cran.r-project.org/package=kza) provides Kolmogorov-Zurbenko
    Adaptive Filters including break detection, spectral analysis,
    wavelets and KZ Fourier Transforms.
    [multitaper](https://cran.r-project.org/package=multitaper) also provides some
    multitaper spectral analysis tools.
-   *Wavelet methods* : The [wavelets](https://cran.r-project.org/package=wavelets)
    package includes computing wavelet filters, wavelet transforms and
    multiresolution analyses. Wavelet methods for time series analysis
    based on Percival and Walden (2000) are given in
    [wmtsa](https://cran.r-project.org/package=wmtsa).
    [WaveletComp](https://cran.r-project.org/package=WaveletComp) provides some
    tools for wavelet-based analysis of univariate and bivariate time
    series including cross-wavelets, phase-difference and significance
    tests. Tests of white noise using wavelets are provided by
    [hwwntest](https://cran.r-project.org/package=hwwntest). Further wavelet methods
    can be found in the packages
    [brainwaver](https://cran.r-project.org/package=brainwaver),
    [rwt](https://cran.r-project.org/package=rwt),
    [waveslim](https://cran.r-project.org/package=waveslim),
    [wavethresh](https://cran.r-project.org/package=wavethresh) and
    [mvcwt](https://cran.r-project.org/package=mvcwt).
-   *Harmonic regression* using Fourier terms is implemented in
    [HarmonicRegression](https://cran.r-project.org/package=HarmonicRegression). The
    [forecast](https://cran.r-project.org/package=forecast) package also provides
    some simple harmonic regression facilities via the `fourier`
    function.

**Decomposition and Filtering**

-   *Filters and smoothing* : `filter()` in stats provides
    autoregressive and moving average linear filtering of multiple
    univariate time series. The
    [robfilter](https://cran.r-project.org/package=robfilter) package provides
    several robust time series filters, while
    [mFilter](https://cran.r-project.org/package=mFilter) includes miscellaneous
    time series filters useful for smoothing and extracting trend and
    cyclical components. `smooth()` from the stats package computes
    Tukey's running median smoothers, 3RS3R, 3RSS, 3R, etc.
    [sleekts](https://cran.r-project.org/package=sleekts) computes the 4253H twice
    smoothing method.
-   *Decomposition* : Seasonal decomposition is discussed below.
    Autoregressive-based decomposition is provided by
    [ArDec](https://cran.r-project.org/package=ArDec).
    [tsdecomp](https://cran.r-project.org/package=tsdecomp) implements ARIMA-based
    decomposition of quarterly and monthly data.
    [rmaf](https://cran.r-project.org/package=rmaf) uses a refined moving average
    filter for decomposition.
-   *Singular Spectrum Analysis* is implemented in
    [Rssa](https://cran.r-project.org/package=Rssa) and
    [spectral.methods](https://cran.r-project.org/package=spectral.methods).
-   *Empirical Mode Decomposition* (EMD) and Hilbert spectral analysis
    is provided by [EMD](https://cran.r-project.org/package=EMD). Additional tools,
    including ensemble EMD, are available in
    [hht](https://cran.r-project.org/package=hht). An alternative implementation of
    ensemble EMD and its complete variant are available in
    [Rlibeemd](https://cran.r-project.org/package=Rlibeemd).

**Seasonality**

-   *Seasonal decomposition* : the stats package provides classical
    decomposition in `decompose()`, and STL decomposition in `stl()`.
    Enhanced STL decomposition is available in
    [stlplus](https://cran.r-project.org/package=stlplus).
    [stR](https://cran.r-project.org/package=stR) provides Seasonal-Trend
    decomposition based on Regression.
-   X-13-ARIMA-SEATS binaries are provided in the
    [x13binary](https://cran.r-project.org/package=x13binary) package, with
    [seasonal](https://cran.r-project.org/package=seasonal) providing an R interface
    and [seasonalview](https://cran.r-project.org/package=seasonalview) providing a
    GUI. An alternative interface is provided by
    [x12](https://cran.r-project.org/package=x12), with an associated alternative
    GUI provided by [x12GUI](https://cran.r-project.org/package=x12GUI).
-   *Analysis of seasonality* : the
    [bfast](https://cran.r-project.org/package=bfast) package provides methods for
    detecting and characterizing abrupt changes within the trend and
    seasonal components obtained from a decomposition.
    [npst](https://cran.r-project.org/package=npst) provides a generalization of
    Hewitt's seasonality test.
-   [season](https://cran.r-project.org/package=season): Seasonal analysis of health
    data including regression models, time-stratified case-crossover,
    plotting functions and residual checks.
-   [seas](https://cran.r-project.org/package=seas): Seasonal analysis and graphics,
    especially for climatology.
-   [deseasonalize](https://cran.r-project.org/package=deseasonalize): Optimal
    deseasonalization for geophysical time series using AR fitting.

**Stationarity, Unit Roots, and Cointegration**

-   *Stationarity and unit roots* :
    [tseries](https://cran.r-project.org/package=tseries) provides various
    stationarity and unit root tests including Augmented Dickey-Fuller,
    Phillips-Perron, and KPSS. Alternative implementations of the ADF
    and KPSS tests are in the [urca](https://cran.r-project.org/package=urca)
    package, which also includes further methods such as
    Elliott-Rothenberg-Stock, Schmidt-Phillips and Zivot-Andrews tests.
    [uroot](https://cran.r-project.org/package=uroot) provides seasonal unit root
    tests. [CADFtest](https://cran.r-project.org/package=CADFtest) provides
    implementations of both the standard ADF and a covariate-augmented
    ADF (CADF) test.
    [MultipleBubbles](https://cran.r-project.org/package=MultipleBubbles) tests for
    the existence of bubbles based on Phillips-Shi-Yu (2015).
-   *Local stationarity* : [locits](https://cran.r-project.org/package=locits)
    provides a test of local stationarity and computes the localized
    autocovariance. Time series costationarity determination is provided
    by [costat](https://cran.r-project.org/package=costat).
    [LSTS](https://cran.r-project.org/package=LSTS) has functions for locally
    stationary time series analysis. Locally stationary wavelet models
    for nonstationary time series are implemented in
    [wavethresh](https://cran.r-project.org/package=wavethresh) (including
    estimation, plotting, and simulation functionality for time-varying
    spectra).
-   *Cointegration* : The Engle-Granger two-step method with the
    Phillips-Ouliaris cointegration test is implemented in
    [tseries](https://cran.r-project.org/package=tseries) and
    [urca](https://cran.r-project.org/package=urca). The latter additionally
    contains functionality for the Johansen trace and lambda-max tests.
    [tsDyn](https://cran.r-project.org/package=tsDyn) provides Johansen's test and
    AIC/BIC simultaneous rank-lag selection.
    [CommonTrend](https://cran.r-project.org/package=CommonTrend) provides tools to
    extract and plot common trends from a cointegration system.
    Parameter estimation and inference in a cointegrating regression are
    implemented in [cointReg](https://cran.r-project.org/package=cointReg).
    [nardl](https://cran.r-project.org/package=nardl) estimates nonlinear
    cointegrating autoregressive distributed lag models.

**Nonlinear Time Series Analysis**

-   *Nonlinear autoregression* : Various forms of nonlinear
    autoregression are available in
    [tsDyn](https://cran.r-project.org/package=tsDyn) including additive AR, neural
    nets, SETAR and LSTAR models, threshold VAR and VECM. Neural network
    autoregression is also provided in
    [GMDH](https://cran.r-project.org/package=GMDH).
    [nnfor](https://cran.r-project.org/package=nnfor) provides time series
    forecasting with neural networks.
    [bentcableAR](https://cran.r-project.org/package=bentcableAR) implements
    Bent-Cable autoregression. [BAYSTAR](https://cran.r-project.org/package=BAYSTAR)
    provides Bayesian analysis of threshold autoregressive models.
-   [tseriesChaos](https://cran.r-project.org/package=tseriesChaos) provides an R
    implementation of the algorithms from the
    *[TISEAN](http://www.mpipks-dresden.mpg.de/~tisean/) project* .
-   Autoregression Markov switching models are provided in
    [MSwM](https://cran.r-project.org/package=MSwM), while dependent mixtures of
    latent Markov models are given in
    [depmix](https://cran.r-project.org/package=depmix) and
    [depmixS4](https://cran.r-project.org/package=depmixS4) for categorical and
    continuous time series.
-   *Tests* : Various tests for nonlinearity are provided in
    [fNonlinear](https://cran.r-project.org/package=fNonlinear).
    [tseriesEntropy](https://cran.r-project.org/package=tseriesEntropy) tests for
    nonlinear serial dependence based on entropy metrics.
-   Additional functions for nonlinear time series are available in
    [nlts](https://cran.r-project.org/package=nlts) and
    [nonlinearTseries](https://cran.r-project.org/package=nonlinearTseries).
-   Fractal time series modeling and analysis is provided by
    [fractal](https://cran.r-project.org/package=fractal).
    [fractalrock](https://cran.r-project.org/package=fractalrock) generates fractal
    time series with non-normal returns distributions.

**Dynamic Regression Models**

-   *Dynamic linear models* : A convenient interface for fitting dynamic
    regression models via OLS is available in
    [dynlm](https://cran.r-project.org/package=dynlm); an enhanced approach that
    also works with other regression functions and more time series
    classes is implemented in [dyn](https://cran.r-project.org/package=dyn). More
    advanced dynamic system equations can be fitted using
    [dse](https://cran.r-project.org/package=dse). Gaussian linear state space
    models can be fitted using [dlm](https://cran.r-project.org/package=dlm) (via
    maximum likelihood, Kalman filtering/smoothing and Bayesian
    methods), or using [bsts](https://cran.r-project.org/package=bsts) which uses
    MCMC. [dLagM](https://cran.r-project.org/package=dLagM) provides time series
    regression with distributed lags. Functions for distributed lag
    non-linear modelling are provided in
    [dlnm](https://cran.r-project.org/package=dlnm).
-   *Time-varying parameter* models can be fitted using the
    [tpr](https://cran.r-project.org/package=tpr) package.
-   [orderedLasso](https://cran.r-project.org/package=orderedLasso) fits a sparse
    linear model with an order constraint on the coefficients in order
    to handle lagged regressors where the coefficients decay as the lag
    increases.
-   Dynamic modeling of various kinds is available in
    [dynr](https://cran.r-project.org/package=dynr) including discrete and
    continuous time, linear and nonlinear models, and different types of
    latent variables.

**Multivariate Time Series Models**

-   *Vector autoregressive (VAR) models* are provided via `ar()` in the
    basic stats package including order selection via the AIC. These
    models are restricted to be stationary.
    [MTS](https://cran.r-project.org/package=MTS) is an all-purpose toolkit for
    analyzing multivariate time series including VAR, VARMA, seasonal
    VARMA, VAR models with exogenous variables, multivariate regression
    with time series errors, and much more. Possibly non-stationary VAR
    models are fitted in the [mAr](https://cran.r-project.org/package=mAr) package,
    which also allows VAR models in principal component space.
    [sparsevar](https://cran.r-project.org/package=sparsevar) allows estimation of
    sparse VAR and VECM models,
    [bigtime](https://cran.r-project.org/package=bigtime) estimates large sparse
    VAR, VARX and VARMA models, while
    [BigVAR](https://cran.r-project.org/package=BigVAR) estimates VAR and VARX
    models with structured lasso penalties and
    [svars](https://cran.r-project.org/package=svars) implements data-driven
    structural VARs. Automated VAR models and networks are available in
    [autovarCore](https://cran.r-project.org/package=autovarCore). More elaborate
    models are provided in package [vars](https://cran.r-project.org/package=vars),
    [tsDyn](https://cran.r-project.org/package=tsDyn), `estVARXls()` in
    [dse](https://cran.r-project.org/package=dse), and a Bayesian approach is
    available in [MSBVAR](https://cran.r-project.org/package=MSBVAR). Another
    implementation with bootstrapped prediction intervals is given in
    [VAR.etp](https://cran.r-project.org/package=VAR.etp).
    [mlVAR](https://cran.r-project.org/package=mlVAR) provides multi-level vector
    autoregression. [VARsignR](https://cran.r-project.org/package=VARsignR) provides
    routines for identifying structural shocks in VAR models using sign
    restrictions. [gdpc](https://cran.r-project.org/package=gdpc) implements
    generalized dynamic principal components.
    [pcdpca](https://cran.r-project.org/package=pcdpca) extends dynamic principal
    components to periodically correlated multivariate time series.
-   *VARIMA models* and *state space models* are provided in the
    [dse](https://cran.r-project.org/package=dse) package.
    [EvalEst](https://cran.r-project.org/package=EvalEst) facilitates Monte Carlo
    experiments to evaluate the associated estimation methods.
-   *Vector error correction models* are available via the
    [urca](https://cran.r-project.org/package=urca),
    [ecm](https://cran.r-project.org/package=ecm),
    [vars](https://cran.r-project.org/package=vars),
    [tsDyn](https://cran.r-project.org/package=tsDyn) packages, including versions
    with structural constraints and thresholding.
-   *Time series component analysis* : Time series factor analysis is
    provided in [tsfa](https://cran.r-project.org/package=tsfa).
    [ForeCA](https://cran.r-project.org/package=ForeCA) implements forecastable
    component analysis by searching for the best linear transformations
    that make a multivariate time series as forecastable as possible.
    [PCA4TS](https://cran.r-project.org/package=PCA4TS) finds a linear
    transformation of a multivariate time series giving
    lower-dimensional subseries that are uncorrelated with each other.
    One-sided dynamic principal components are computed in
    [odpc](https://cran.r-project.org/package=odpc). Frequency-domain-based dynamic
    PCA is implemented in [freqdom](https://cran.r-project.org/package=freqdom).
-   *Multivariate state space models* An implementation is provided by
    the [KFAS](https://cran.r-project.org/package=KFAS) package which provides a
    fast multivariate Kalman filter, smoother, simulation smoother and
    forecasting. Yet another implementation is given in the
    [dlm](https://cran.r-project.org/package=dlm) package which also contains tools
    for converting other multivariate models into state space form.
    [dlmodeler](https://cran.r-project.org/package=dlmodeler) provides a unified
    interface for [dlm](https://cran.r-project.org/package=dlm) and
    [KFAS](https://cran.r-project.org/package=KFAS).
    [MARSS](https://cran.r-project.org/package=MARSS) fits constrained and
    unconstrained multivariate autoregressive state-space models using
    an EM algorithm. All of these packages assume the observational and
    state error terms are uncorrelated.
-   *Partially-observed Markov processes* are a generalization of the
    usual linear multivariate state space models, allowing non-Gaussian
    and nonlinear models. These are implemented in the
    [pomp](https://cran.r-project.org/package=pomp) package.
-   Multivariate stochastic volatility models (using latent factors) are
    provided by [factorstochvol](https://cran.r-project.org/package=factorstochvol).

**Analysis of large groups of time series**

-   *Time series clustering* is implemented in
    [TSclust](https://cran.r-project.org/package=TSclust),
    [dtwclust](https://cran.r-project.org/package=dtwclust),
    [BNPTSclust](https://cran.r-project.org/package=BNPTSclust) and
    [pdc](https://cran.r-project.org/package=pdc).
-   [TSdist](https://cran.r-project.org/package=TSdist) provides distance measures
    for time series data.
-   [TSrepr](https://cran.r-project.org/package=TSrepr) includes methods for
    representing time series using dimension reduction and feature
    extraction.
-   [jmotif](https://cran.r-project.org/package=jmotif) implements tools based on
    time series symbolic discretization for finding motifs in time
    series and facilitates interpretable time series classification.
-   [rucrdtw](https://cran.r-project.org/package=rucrdtw) provides R bindings for
    functions from the UCR Suite to enable ultrafast subsequence search
    for a best match under Dynamic Time Warping and Euclidean Distance.
-   Methods for plotting and forecasting collections of hierarchical and
    grouped time series are provided by
    [hts](https://cran.r-project.org/package=hts).
    [thief](https://cran.r-project.org/package=thief) uses hierarchical methods to
    reconcile forecasts of temporally aggregated time series. An
    alternative approach to reconciling forecasts of hierarchical time
    series is provided by [gtop](https://cran.r-project.org/package=gtop).
    [thief](https://cran.r-project.org/package=thief)

**Functional time series**

-   Tools for visualizing, modeling, forecasting and analysis of
    functional time series are implemented in
    [ftsa](https://cran.r-project.org/package=ftsa).
-   [freqdom.fda](https://cran.r-project.org/package=freqdom.fda) provides
    implements of dynamical functional principal components for
    functional time series.

**Continuous time models**

-   *Continuous time autoregressive modelling* is provided in
    [cts](https://cran.r-project.org/package=cts), while
    [carfima](https://cran.r-project.org/package=carfima) allows for continuous-time
    ARFIMA models.
-   [Sim.DiffProc](https://cran.r-project.org/package=Sim.DiffProc) simulates and
    models stochastic differential equations.
-   Simulation and inference for stochastic differential equations is
    provided by [sde](https://cran.r-project.org/package=sde) and
    [yuima](https://cran.r-project.org/package=yuima).

**Resampling**

-   *Bootstrapping* : The [boot](https://cran.r-project.org/package=boot) package
    provides function `tsboot()` for time series bootstrapping,
    including block bootstrap with several variants. `tsbootstrap()`
    from [tseries](https://cran.r-project.org/package=tseries) provides fast
    stationary and block bootstrapping. Maximum entropy bootstrap for
    time series is available in [meboot](https://cran.r-project.org/package=meboot).
    [timesboot](https://cran.r-project.org/package=timesboot) computes the bootstrap
    CI for the sample ACF and periodogram.
    [BootPR](https://cran.r-project.org/package=BootPR) computes bias-corrected
    forecasting and bootstrap prediction intervals for autoregressive
    time series.

**Time Series Data**

-   Data from Makridakis, Wheelwright and Hyndman (1998) *Forecasting:
    methods and applications* are provided in the
    [fma](https://cran.r-project.org/package=fma) package.
-   Data from Hyndman, Koehler, Ord and Snyder (2008) *Forecasting with
    exponential smoothing* are in the
    [expsmooth](https://cran.r-project.org/package=expsmooth) package.
-   Data from Hyndman and Athanasopoulos (2013) *Forecasting: principles
    and practice* are in the [fpp](https://cran.r-project.org/package=fpp) package.
-   Data from Hyndman and Athanasopoulos (2017) *Forecasting: principles
    and practice* (2nd ed) are in the
    [fpp2](https://cran.r-project.org/package=fpp2) package.
-   Data from the M-competition and M3-competition are provided in the
    [Mcomp](https://cran.r-project.org/package=Mcomp) package.
    [Tcomp](https://cran.r-project.org/package=Tcomp) provides data from the 2010
    IJF Tourism Forecasting Competition.
-   [pdfetch](https://cran.r-project.org/package=pdfetch) provides facilities for
    downloading economic and financial time series from public sources.
-   Data from the [Quandl](http://www.quandl.com) online portal to
    financial, economical and social datasets can be queried
    interactively using the [Quandl](https://cran.r-project.org/package=Quandl)
    package.
-   Data from the [Datamarket](http://datamarket.com/data/) online
    portal can be fetched using the
    [rdatamarket](https://cran.r-project.org/package=rdatamarket) package.
-   Data from Switzerland via [dataseries.org](http://dataseries.org)
    can be downloaded and imported using
    [dataseries](https://cran.r-project.org/package=dataseries).
-   [BETS](https://cran.r-project.org/package=BETS) provides access to the most
    important economic time series in Brazil.
-   [rbcb](https://cran.r-project.org/package=rbcb) R interface to Brazilian Central Bank
    web services.
-   Data from Cryer and Chan (2010) are in the
    [TSA](https://cran.r-project.org/package=TSA) package.
-   Data from Shumway and Stoffer (2011) are in the
    [astsa](https://cran.r-project.org/package=astsa) package.
-   [tswge](https://cran.r-project.org/package=tswge) accompanies the text *Applied
    Time Series Analysis with R* , 2nd edition by Woodward, Gray, and
    Elliott.
-   [TSdbi](https://cran.r-project.org/package=TSdbi) provides a common interface to
    time series databases.
-   [fame](https://cran.r-project.org/package=fame) provides an interface for FAME
    time series databases
-   [influxdbr](https://cran.r-project.org/package=influxdbr) provides an interface
    to the InfluxDB time series database.
-   [AER](https://cran.r-project.org/package=AER) and
    [Ecdat](https://cran.r-project.org/package=Ecdat) both contain many data sets
    (including time series data) from many econometrics text books

**Miscellaneous**

-   [dtw](https://cran.r-project.org/package=dtw): Dynamic time warping algorithms
    for computing and plotting pairwise alignments between time series.
-   [ensembleBMA](https://cran.r-project.org/package=ensembleBMA): Bayesian Model
    Averaging to create probabilistic forecasts from ensemble forecasts
    and weather observations.
-   [earlywarnings](https://cran.r-project.org/package=earlywarnings): Early
    warnings signals toolbox for detecting critical transitions in time
    series
-   [events](https://cran.r-project.org/package=events): turns machine-extracted
    event data into regular aggregated multivariate time series.
-   [FeedbackTS](https://cran.r-project.org/package=FeedbackTS): Analysis of
    fragmented time directionality to investigate feedback in time
    series.
-   [LPStimeSeries](https://cran.r-project.org/package=LPStimeSeries) aims to find
    "learned pattern similarity" for time series.
-   [MAR1](https://cran.r-project.org/package=MAR1) provides tools for preparing
    ecological community time series data for multivariate AR modeling.
-   [nets](https://cran.r-project.org/package=nets): routines for the estimation of
    sparse long run partial correlation networks for time series data.
-   [paleoTS](https://cran.r-project.org/package=paleoTS): Modeling evolution in
    paleontological time series.
-   [pastecs](https://cran.r-project.org/package=pastecs): Regulation, decomposition
    and analysis of space-time series.
-   [PSF](https://cran.r-project.org/package=PSF): Forecasting univariate time
    series using pattern-sequences.
-   [ptw](https://cran.r-project.org/package=ptw): Parametric time warping.
-   [RGENERATE](https://cran.r-project.org/package=RGENERATE) provides tools to
    generate vector time series.
-   [RMAWGEN](https://cran.r-project.org/package=RMAWGEN) is set of S3 and S4
    functions for spatial multi-site stochastic generation of daily
    time-series of temperature and precipitation making use of VAR
    models. The package can be used in climatology and statistical
    hydrology.
-   [RSEIS](https://cran.r-project.org/package=RSEIS): Seismic time series analysis
    tools.
-   [rts](https://cran.r-project.org/package=rts): Raster time series analysis
    (e.g., time series of satellite images).
-   [sae2](https://cran.r-project.org/package=sae2): Time series models for small
    area estimation.
-   [spTimer](https://cran.r-project.org/package=spTimer): Spatio-temporal Bayesian
    modelling.
-   [surveillance](https://cran.r-project.org/package=surveillance): Temporal and
    spatio-temporal modeling and monitoring of epidemic phenomena.
-   [TED](https://cran.r-project.org/package=TED): Turbulence time series Event
    Detection and classification.
-   [Tides](https://cran.r-project.org/package=Tides): Functions to calculate
    characteristics of quasi periodic time series, e.g. observed
    estuarine water levels.
-   [tiger](https://cran.r-project.org/package=tiger): Temporally resolved groups of
    typical differences (errors) between two time series are determined
    and visualized.
-   [tsfknn](https://cran.r-project.org/package=tsfknn): Time series forecasting
    with k-nearest-neighbours.
-   [TSMining](https://cran.r-project.org/package=TSMining): Mining Univariate and
    Multivariate Motifs in Time-Series Data.
-   [tsModel](https://cran.r-project.org/package=tsModel): Time series modeling for
    air pollution and health.

### CRAN packages:

-   [acp](https://cran.r-project.org/package=acp)
-   [AER](https://cran.r-project.org/package=AER)
-   [ArDec](https://cran.r-project.org/package=ArDec)
-   [arfima](https://cran.r-project.org/package=arfima)
-   [astsa](https://cran.r-project.org/package=astsa)
-   [autovarCore](https://cran.r-project.org/package=autovarCore)
-   [BAYSTAR](https://cran.r-project.org/package=BAYSTAR)
-   [bentcableAR](https://cran.r-project.org/package=bentcableAR)
-   [BETS](https://cran.r-project.org/package=BETS)
-   [bfast](https://cran.r-project.org/package=bfast)
-   [bigtime](https://cran.r-project.org/package=bigtime)
-   [BigVAR](https://cran.r-project.org/package=BigVAR)
-   [BNPTSclust](https://cran.r-project.org/package=BNPTSclust)
-   [boot](https://cran.r-project.org/package=boot)
-   [BootPR](https://cran.r-project.org/package=BootPR)
-   [brainwaver](https://cran.r-project.org/package=brainwaver)
-   [bspec](https://cran.r-project.org/package=bspec)
-   [bsts](https://cran.r-project.org/package=bsts)
-   [CADFtest](https://cran.r-project.org/package=CADFtest)
-   [carfima](https://cran.r-project.org/package=carfima)
-   [carx](https://cran.r-project.org/package=carx)
-   [cents](https://cran.r-project.org/package=cents)
-   [changepoint](https://cran.r-project.org/package=changepoint)
-   [chron](https://cran.r-project.org/package=chron)
-   [cointReg](https://cran.r-project.org/package=cointReg)
-   [CommonTrend](https://cran.r-project.org/package=CommonTrend)
-   [costat](https://cran.r-project.org/package=costat)
-   [cts](https://cran.r-project.org/package=cts)
-   [dataseries](https://cran.r-project.org/package=dataseries)
-   [dCovTS](https://cran.r-project.org/package=dCovTS)
-   [depmix](https://cran.r-project.org/package=depmix)
-   [depmixS4](https://cran.r-project.org/package=depmixS4)
-   [deseasonalize](https://cran.r-project.org/package=deseasonalize)
-   [dLagM](https://cran.r-project.org/package=dLagM)
-   [dlm](https://cran.r-project.org/package=dlm)
-   [dlmodeler](https://cran.r-project.org/package=dlmodeler)
-   [dlnm](https://cran.r-project.org/package=dlnm)
-   [dse](https://cran.r-project.org/package=dse)
-   [dtw](https://cran.r-project.org/package=dtw)
-   [dtwclust](https://cran.r-project.org/package=dtwclust)
-   [dygraphs](https://cran.r-project.org/package=dygraphs)
-   [dyn](https://cran.r-project.org/package=dyn)
-   [dynlm](https://cran.r-project.org/package=dynlm)
-   [dynr](https://cran.r-project.org/package=dynr)
-   [earlywarnings](https://cran.r-project.org/package=earlywarnings)
-   [Ecdat](https://cran.r-project.org/package=Ecdat)
-   [ecm](https://cran.r-project.org/package=ecm)
-   [ecp](https://cran.r-project.org/package=ecp)
-   [EMD](https://cran.r-project.org/package=EMD)
-   [ensembleBMA](https://cran.r-project.org/package=ensembleBMA)
-   [EvalEst](https://cran.r-project.org/package=EvalEst)
-   [events](https://cran.r-project.org/package=events)
-   [expsmooth](https://cran.r-project.org/package=expsmooth)
-   [factorstochvol](https://cran.r-project.org/package=factorstochvol)
-   [fame](https://cran.r-project.org/package=fame)
-   [fanplot](https://cran.r-project.org/package=fanplot)
-   [FeedbackTS](https://cran.r-project.org/package=FeedbackTS)
-   [fGarch](https://cran.r-project.org/package=fGarch)
-   [FGN](https://cran.r-project.org/package=FGN)
-   [FitAR](https://cran.r-project.org/package=FitAR)
-   [FitARMA](https://cran.r-project.org/package=FitARMA)
-   [fma](https://cran.r-project.org/package=fma)
-   [fNonlinear](https://cran.r-project.org/package=fNonlinear)
-   [ForeCA](https://cran.r-project.org/package=ForeCA)
-   [forecast](https://cran.r-project.org/package=forecast) (core)
-   [ForecastComb](https://cran.r-project.org/package=ForecastComb)
-   [forecastHybrid](https://cran.r-project.org/package=forecastHybrid)
-   [forecTheta](https://cran.r-project.org/package=forecTheta)
-   [fpp](https://cran.r-project.org/package=fpp)
-   [fpp2](https://cran.r-project.org/package=fpp2)
-   [fracdiff](https://cran.r-project.org/package=fracdiff)
-   [fractal](https://cran.r-project.org/package=fractal)
-   [fractalrock](https://cran.r-project.org/package=fractalrock)
-   [freqdom](https://cran.r-project.org/package=freqdom)
-   [freqdom.fda](https://cran.r-project.org/package=freqdom.fda)
-   [fts](https://cran.r-project.org/package=fts)
-   [ftsa](https://cran.r-project.org/package=ftsa)
-   [funtimes](https://cran.r-project.org/package=funtimes)
-   [GAS](https://cran.r-project.org/package=GAS)
-   [gdpc](https://cran.r-project.org/package=gdpc)
-   [ggseas](https://cran.r-project.org/package=ggseas)
-   [glarma](https://cran.r-project.org/package=glarma)
-   [GMDH](https://cran.r-project.org/package=GMDH)
-   [gsarima](https://cran.r-project.org/package=gsarima)
-   [gtop](https://cran.r-project.org/package=gtop)
-   [HarmonicRegression](https://cran.r-project.org/package=HarmonicRegression)
-   [hht](https://cran.r-project.org/package=hht)
-   [hts](https://cran.r-project.org/package=hts)
-   [hwwntest](https://cran.r-project.org/package=hwwntest)
-   [imputeTS](https://cran.r-project.org/package=imputeTS)
-   [influxdbr](https://cran.r-project.org/package=influxdbr)
-   [InspectChangepoint](https://cran.r-project.org/package=InspectChangepoint)
-   [itsmr](https://cran.r-project.org/package=itsmr)
-   [jmotif](https://cran.r-project.org/package=jmotif)
-   [KFAS](https://cran.r-project.org/package=KFAS)
-   [KFKSDS](https://cran.r-project.org/package=KFKSDS)
-   [kza](https://cran.r-project.org/package=kza)
-   [locits](https://cran.r-project.org/package=locits)
-   [lomb](https://cran.r-project.org/package=lomb)
-   [LPStimeSeries](https://cran.r-project.org/package=LPStimeSeries)
-   [LSTS](https://cran.r-project.org/package=LSTS)
-   [ltsa](https://cran.r-project.org/package=ltsa)
-   [lubridate](https://cran.r-project.org/package=lubridate)
-   [mafs](https://cran.r-project.org/package=mafs)
-   [MAPA](https://cran.r-project.org/package=MAPA)
-   [mAr](https://cran.r-project.org/package=mAr)
-   [MAR1](https://cran.r-project.org/package=MAR1)
-   [mar1s](https://cran.r-project.org/package=mar1s)
-   [MARSS](https://cran.r-project.org/package=MARSS)
-   [mclcar](https://cran.r-project.org/package=mclcar)
-   [Mcomp](https://cran.r-project.org/package=Mcomp)
-   [meboot](https://cran.r-project.org/package=meboot)
-   [mFilter](https://cran.r-project.org/package=mFilter)
-   [mlVAR](https://cran.r-project.org/package=mlVAR)
-   [mondate](https://cran.r-project.org/package=mondate)
-   [MSBVAR](https://cran.r-project.org/package=MSBVAR)
-   [MSwM](https://cran.r-project.org/package=MSwM)
-   [MTS](https://cran.r-project.org/package=MTS)
-   [mtsdi](https://cran.r-project.org/package=mtsdi)
-   [MultipleBubbles](https://cran.r-project.org/package=MultipleBubbles)
-   [multitaper](https://cran.r-project.org/package=multitaper)
-   [mvcwt](https://cran.r-project.org/package=mvcwt)
-   [nardl](https://cran.r-project.org/package=nardl)
-   [nets](https://cran.r-project.org/package=nets)
-   [nlts](https://cran.r-project.org/package=nlts)
-   [nnfor](https://cran.r-project.org/package=nnfor)
-   [nonlinearTseries](https://cran.r-project.org/package=nonlinearTseries)
-   [npst](https://cran.r-project.org/package=npst)
-   [odpc](https://cran.r-project.org/package=odpc)
-   [opera](https://cran.r-project.org/package=opera)
-   [orderedLasso](https://cran.r-project.org/package=orderedLasso)
-   [paleoTS](https://cran.r-project.org/package=paleoTS)
-   [partsm](https://cran.r-project.org/package=partsm)
-   [pastecs](https://cran.r-project.org/package=pastecs)
-   [PCA4TS](https://cran.r-project.org/package=PCA4TS)
-   [pcdpca](https://cran.r-project.org/package=pcdpca)
-   [pdc](https://cran.r-project.org/package=pdc)
-   [pdfetch](https://cran.r-project.org/package=pdfetch)
-   [pear](https://cran.r-project.org/package=pear)
-   [perARMA](https://cran.r-project.org/package=perARMA)
-   [pomp](https://cran.r-project.org/package=pomp)
-   [portes](https://cran.r-project.org/package=portes)
-   [prophet](https://cran.r-project.org/package=prophet)
-   [psd](https://cran.r-project.org/package=psd)
-   [PSF](https://cran.r-project.org/package=PSF)
-   [ptw](https://cran.r-project.org/package=ptw)
-   [Quandl](https://cran.r-project.org/package=Quandl)
-   [quantspec](https://cran.r-project.org/package=quantspec)
-   [rdatamarket](https://cran.r-project.org/package=rdatamarket)
-   [RGENERATE](https://cran.r-project.org/package=RGENERATE)
-   [Rlibeemd](https://cran.r-project.org/package=Rlibeemd)
-   [rmaf](https://cran.r-project.org/package=rmaf)
-   [RMAWGEN](https://cran.r-project.org/package=RMAWGEN)
-   [robets](https://cran.r-project.org/package=robets)
-   [robfilter](https://cran.r-project.org/package=robfilter)
-   [robustarima](https://cran.r-project.org/package=robustarima)
-   [roll](https://cran.r-project.org/package=roll)
-   [RSEIS](https://cran.r-project.org/package=RSEIS)
-   [Rssa](https://cran.r-project.org/package=Rssa)
-   [rts](https://cran.r-project.org/package=rts)
-   [rucrdtw](https://cran.r-project.org/package=rucrdtw)
-   [rugarch](https://cran.r-project.org/package=rugarch)
-   [rwt](https://cran.r-project.org/package=rwt)
-   [sae2](https://cran.r-project.org/package=sae2)
-   [scoringRules](https://cran.r-project.org/package=scoringRules)
-   [SDD](https://cran.r-project.org/package=SDD)
-   [sde](https://cran.r-project.org/package=sde)
-   [seas](https://cran.r-project.org/package=seas)
-   [season](https://cran.r-project.org/package=season)
-   [seasonal](https://cran.r-project.org/package=seasonal)
-   [seasonalview](https://cran.r-project.org/package=seasonalview)
-   [Sim.DiffProc](https://cran.r-project.org/package=Sim.DiffProc)
-   [sleekts](https://cran.r-project.org/package=sleekts)
-   [smooth](https://cran.r-project.org/package=smooth)
-   [sparsevar](https://cran.r-project.org/package=sparsevar)
-   [spectral](https://cran.r-project.org/package=spectral)
-   [spectral.methods](https://cran.r-project.org/package=spectral.methods)
-   [spTimer](https://cran.r-project.org/package=spTimer)
-   [stlplus](https://cran.r-project.org/package=stlplus)
-   [stochvol](https://cran.r-project.org/package=stochvol)
-   [stR](https://cran.r-project.org/package=stR)
-   [strucchange](https://cran.r-project.org/package=strucchange)
-   [stsm](https://cran.r-project.org/package=stsm)
-   [stsm.class](https://cran.r-project.org/package=stsm.class)
-   [surveillance](https://cran.r-project.org/package=surveillance)
-   [svars](https://cran.r-project.org/package=svars)
-   [sweep](https://cran.r-project.org/package=sweep)
-   [Tcomp](https://cran.r-project.org/package=Tcomp)
-   [TED](https://cran.r-project.org/package=TED)
-   [tempdisagg](https://cran.r-project.org/package=tempdisagg)
-   [tframe](https://cran.r-project.org/package=tframe)
-   [thief](https://cran.r-project.org/package=thief)
-   [tibbletime](https://cran.r-project.org/package=tibbletime)
-   [Tides](https://cran.r-project.org/package=Tides)
-   [tiger](https://cran.r-project.org/package=tiger)
-   [timeDate](https://cran.r-project.org/package=timeDate)
-   [TimeProjection](https://cran.r-project.org/package=TimeProjection)
-   [timesboot](https://cran.r-project.org/package=timesboot)
-   [timeSeries](https://cran.r-project.org/package=timeSeries)
-   [timetk](https://cran.r-project.org/package=timetk)
-   [timsac](https://cran.r-project.org/package=timsac)
-   [tis](https://cran.r-project.org/package=tis)
-   [tpr](https://cran.r-project.org/package=tpr)
-   [trend](https://cran.r-project.org/package=trend)
-   [TSA](https://cran.r-project.org/package=TSA)
-   [tsbugs](https://cran.r-project.org/package=tsbugs)
-   [TSclust](https://cran.r-project.org/package=TSclust)
-   [tscount](https://cran.r-project.org/package=tscount)
-   [TSdbi](https://cran.r-project.org/package=TSdbi)
-   [tsdecomp](https://cran.r-project.org/package=tsdecomp)
-   [tsdisagg2](https://cran.r-project.org/package=tsdisagg2)
-   [TSdist](https://cran.r-project.org/package=TSdist)
-   [tsDyn](https://cran.r-project.org/package=tsDyn)
-   [tseries](https://cran.r-project.org/package=tseries) (core)
-   [tseriesChaos](https://cran.r-project.org/package=tseriesChaos)
-   [tseriesEntropy](https://cran.r-project.org/package=tseriesEntropy)
-   [tsfa](https://cran.r-project.org/package=tsfa)
-   [tsfknn](https://cran.r-project.org/package=tsfknn)
-   [tsibble](https://cran.r-project.org/package=tsibble)
-   [tsintermittent](https://cran.r-project.org/package=tsintermittent)
-   [TSMining](https://cran.r-project.org/package=TSMining)
-   [tsModel](https://cran.r-project.org/package=tsModel)
-   [tsoutliers](https://cran.r-project.org/package=tsoutliers)
-   [tsPI](https://cran.r-project.org/package=tsPI)
-   [TSrepr](https://cran.r-project.org/package=TSrepr)
-   [TSstudio](https://cran.r-project.org/package=TSstudio)
-   [TSTutorial](https://cran.r-project.org/package=TSTutorial)
-   [tswge](https://cran.r-project.org/package=tswge)
-   [urca](https://cran.r-project.org/package=urca)
-   [uroot](https://cran.r-project.org/package=uroot)
-   [VAR.etp](https://cran.r-project.org/package=VAR.etp)
-   [vars](https://cran.r-project.org/package=vars)
-   [VARsignR](https://cran.r-project.org/package=VARsignR)
-   [Wats](https://cran.r-project.org/package=Wats)
-   [WaveletComp](https://cran.r-project.org/package=WaveletComp)
-   [wavelets](https://cran.r-project.org/package=wavelets)
-   [waveslim](https://cran.r-project.org/package=waveslim)
-   [wavethresh](https://cran.r-project.org/package=wavethresh)
-   [WeightedPortTest](https://cran.r-project.org/package=WeightedPortTest)
-   [wktmo](https://cran.r-project.org/package=wktmo)
-   [wmtsa](https://cran.r-project.org/package=wmtsa)
-   [x12](https://cran.r-project.org/package=x12)
-   [x12GUI](https://cran.r-project.org/package=x12GUI)
-   [x13binary](https://cran.r-project.org/package=x13binary)
-   [xts](https://cran.r-project.org/package=xts)
-   [yuima](https://cran.r-project.org/package=yuima)
-   [ZIM](https://cran.r-project.org/package=ZIM)
-   [zoo](https://cran.r-project.org/package=zoo) (core)
-   [ZRA](https://cran.r-project.org/package=ZRA)

### Related links:

-   CRAN Task View: [Finance](Finance.html)
-   CRAN Task View: [Econometrics](Econometrics.html)
-   CRAN Task View: [Environmetrics](Environmetrics.html)
-   [Time Series Data Library](http://data.is/TSDLdemo)
-   [TISEAN Project](http://www.mpipks-dresden.mpg.de/~tisean/)
-   [Quandl](http://www.quandl.com/)
-   [GitHub repository for this Task
    View](https://github.com/robjhyndman/ctv-TimeSeries)
