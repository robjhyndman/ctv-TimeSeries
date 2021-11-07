## CRAN Task View: Time Series Analysis

                                                               
**Maintainer:** Rob J Hyndman                                  
**Contact:**    Rob.Hyndman at monash.edu                      
**Version:**    2021-11-08                                     
**URL:**        <https://CRAN.R-project.org/view=TimeSeries>   

<div>

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
    statistics functions. [slider](https://cran.r-project.org/package=slider)
    calculates a diverse and comprehensive set of type-stable running
    functions for any R data types.
    [tsibble](https://cran.r-project.org/package=tsibble) provides `slide()` for
    rolling statistics, `tile()` for non-overlapping sliding windows,
    and `stretch()` for expanding windows.
    [tbrf](https://cran.r-project.org/package=tbrf) provides rolling functions based
    on date and time windows instead of n-lagged observations.
    [roll](https://cran.r-project.org/package=roll) provides parallel functions for
    computing rolling statistics.
    [runner](https://cran.r-project.org/package=runner) provides tools for running
    any R function in rolling windows or date windows.
    [runstats](https://cran.r-project.org/package=runstats) provides fast
    computational methods for some running sample statistics. For
    [data.table](https://cran.r-project.org/package=data.table), `froll()` can be
    used for high-performance rolling statistics.
-   *Graphics* : Time series plots are obtained with `plot()` applied to
    `ts` objects. (Partial) autocorrelation functions plots are
    implemented in `acf()` and `pacf()`. Alternative versions are
    provided by `Acf()` and `Pacf()` in
    [forecast](https://cran.r-project.org/package=forecast), along with a
    combination display using `tsdisplay()`. Seasonal displays are
    obtained using `monthplot()` in stats, `seasonplot` in
    [forecast](https://cran.r-project.org/package=forecast), and `seasplot` in
    [tsutils](https://cran.r-project.org/package=tsutils).
    [feasts](https://cran.r-project.org/package=feasts) and
    [brolgar](https://cran.r-project.org/package=brolgar) provide various time
    series graphics for tsibble objects including time plots, season
    plots, subseries plots, ACF and PACF plots, and some combination
    displays. Interactive graphics for tsibbles using htmlwidgets are
    provided by [tsibbletalk](https://cran.r-project.org/package=tsibbletalk).
    [SDD](https://cran.r-project.org/package=SDD) provides more general serial
    dependence diagrams, while [dCovTS](https://cran.r-project.org/package=dCovTS)
    computes and plots the distance covariance and correlation functions
    of time series. [ggseas](https://cran.r-project.org/package=ggseas) provides
    additional ggplot2 graphics for seasonally adjusted series and
    rolling statistics. Calendar plots are implemented in
    [sugrrants](https://cran.r-project.org/package=sugrrants).
    [gravitas](https://cran.r-project.org/package=gravitas) allows for visualizing
    probability distributions conditional on bivariate temporal
    granularities. [dygraphs](https://cran.r-project.org/package=dygraphs) provides
    an interface to the Dygraphs interactive time series charting
    library. [TSstudio](https://cran.r-project.org/package=TSstudio) provides some
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
    4(1)](http://CRAN.R-project.org/doc/Rnews/Rnews_2004-1.pdf) , 29-32.
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
    functions that facilitate certain POSIX-based computations, while
    [clock](https://cran.r-project.org/package=clock) provides a comprehensive
    library for date-time manipulations using a new family of orthogonal
    date-time classes (durations, time points, zoned-times, and
    calendars). [timechange](https://cran.r-project.org/package=timechange) allows
    for efficient manipulation of date-times accounting for time zones
    and daylight saving times. [wktmo](https://cran.r-project.org/package=wktmo)
    converts weekly data to monthly data in several different ways.
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
    [tsdisagg2](https://cran.r-project.org/package=tsdisagg2) and
    [disaggR](https://cran.r-project.org/package=disaggR).
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
    consistent as possible with `"ts"`.
-   The package [xts](https://cran.r-project.org/package=xts) is based on
    [zoo](https://cran.r-project.org/package=zoo) and provides uniform handling of
    R's different time-based data classes.
-   Several packages aim to handle time-based tibbles:
    [tsibble](https://cran.r-project.org/package=tsibble) provides tidy temporal
    data frames and associated tools;
    [tsbox](https://cran.r-project.org/package=tsbox) contains tools for working
    with and coercing between many time series classes including
    tsibble, ts, xts, zoo and more.
    [timetk](https://cran.r-project.org/package=timetk) is another toolkit for
    converting between various time series data classes.
-   Some manipulation tools for time series are available in
    [data.table](https://cran.r-project.org/package=data.table) including `shift()`
    for lead/lag operations. Further basic time series functionalities
    are offered by [DTSg](https://cran.r-project.org/package=DTSg) which is based on
    [data.table](https://cran.r-project.org/package=data.table).
-   [collapse](https://cran.r-project.org/package=collapse) provides fast
    computation of several time series functions such as lead/lag
    operations, (quasi-, log-) differences and growth rates on
    time-series and panel data, and ACF/PACF/CCF estimation for panel
    data.
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
-   [timeseriesdb](https://cran.r-project.org/package=timeseriesdb) manages time
    series for official statistics by mapping `ts` objects to PostgreSQL
    relations.

**Forecasting and Univariate Modeling**

-   The [fable](https://cran.r-project.org/package=fable) package provides tools for
    fitting univariate time series models to many series simultaneously
    including ETS, ARIMA, TSLM and other models. It also provides many
    functions for computing and analysing forecasts. The time series
    must be in the `tsibble` format.
    [fabletools](https://cran.r-project.org/package=fabletools) provides tools for
    extending the [fable](https://cran.r-project.org/package=fable) framework.
-   The [forecast](https://cran.r-project.org/package=forecast) package provides
    similar tools for `ts` objects, while
    [modeltime](https://cran.r-project.org/package=modeltime) and
    [modeltime.ensemble](https://cran.r-project.org/package=modeltime.ensemble)
    provides time series forecasting tools for use with the 'tidymodels'
    ecosystem.
-   *Exponential smoothing* : `HoltWinters()` in stats provides some
    basic models with partial optimization, `ETS()` from
    [fable](https://cran.r-project.org/package=fable) and `ets()` from
    [forecast](https://cran.r-project.org/package=forecast) provide a larger set of
    models and facilities with full optimization.
    [robets](https://cran.r-project.org/package=robets) provides a robust
    alternative to the `ets()` function.
    [smooth](https://cran.r-project.org/package=smooth) implements some
    generalizations of exponential smoothing.
    [legion](https://cran.r-project.org/package=legion) implements multivariate
    versions of exponential smoothing. The
    [MAPA](https://cran.r-project.org/package=MAPA) package combines exponential
    smoothing models at different levels of temporal aggregation to
    improve forecast accuracy. Some Bayesian extensions of exponential
    smoothing are contained in [Rlgt](https://cran.r-project.org/package=Rlgt).
-   [prophet](https://cran.r-project.org/package=prophet) forecasts time series
    based on an additive model where nonlinear trends are fit with
    yearly and weekly seasonality, plus holidays. It works best with
    daily data. [fable.prophet](https://cran.r-project.org/package=fable.prophet)
    allows prophet models to be used in the
    [fable](https://cran.r-project.org/package=fable) framework.
-   The theta method is implemented in the `THETA()` function from
    [fable](https://cran.r-project.org/package=fable), `thetaf()` function from
    [forecast](https://cran.r-project.org/package=forecast), and `theta()` from
    [tsutils](https://cran.r-project.org/package=tsutils). An alternative and
    extended implementation is provided in
    [forecTheta](https://cran.r-project.org/package=forecTheta).
-   *Autoregressive models* : `ar()` in stats (with model selection) and
    [FitAR](https://cran.r-project.org/package=FitAR) for subset AR models.
-   *ARIMA models* : `arima()` in stats is the basic function for ARIMA,
    SARIMA, RegARIMA, and subset ARIMA models. It is enhanced in the
    [fable](https://cran.r-project.org/package=fable) package via the `ARIMA()`
    function which allows for automatic modelling. Similar functionality
    is provided in the [forecast](https://cran.r-project.org/package=forecast)
    package via the `auto.arima()` function. `arma()` in the
    [tseries](https://cran.r-project.org/package=tseries) package provides different
    algorithms for ARMA and subset ARMA models. Other estimation methods
    including the innovations algorithm are provided by
    [itsmr](https://cran.r-project.org/package=itsmr).
    [FitARMA](https://cran.r-project.org/package=FitARMA) implements a fast MLE
    algorithm for ARMA models. Package
    [gsarima](https://cran.r-project.org/package=gsarima) contains functionality for
    Generalized SARIMA time series simulation. Robust ARIMA modeling is
    provided in the [robustarima](https://cran.r-project.org/package=robustarima)
    package. [bayesforecast](https://cran.r-project.org/package=bayesforecast) fits
    Bayesian time series models including seasonal ARIMA and ARIMAX
    models. [BayesARIMAX](https://cran.r-project.org/package=BayesARIMAX) implements
    Bayesian estimation of ARIMAX models. The
    [mar1s](https://cran.r-project.org/package=mar1s) package handles multiplicative
    AR(1) with seasonal processes.
    [TSTutorial](https://cran.r-project.org/package=TSTutorial) provides an
    interactive tutorial for Box-Jenkins modelling. Improved prediction
    intervals for ARIMA and structural time series models are provided
    by [tsPI](https://cran.r-project.org/package=tsPI).
-   ARIMA models with multiple seasonal periods can be handled with
    [tfarima](https://cran.r-project.org/package=tfarima) and
    [smooth](https://cran.r-project.org/package=smooth).
-   *Periodic ARMA models* : [partsm](https://cran.r-project.org/package=partsm) for
    periodic autoregressive time series models, and
    [perARMA](https://cran.r-project.org/package=perARMA) and
    [pcts](https://cran.r-project.org/package=pcts) for periodic ARMA modelling and
    other procedures for periodic time series analysis.
-   *Long memory models* : Some facilities for fractional differenced
    ARFIMA models are provided in the
    [fracdiff](https://cran.r-project.org/package=fracdiff) package. The
    [arfima](https://cran.r-project.org/package=arfima) package has more advanced
    and general facilities for ARFIMA and ARIMA models, including
    dynamic regression (transfer function) models. Additional methods
    for fitting and simulating non-stationary ARFIMA models are in
    [nsarfima](https://cran.r-project.org/package=nsarfima).
    [LongMemoryTS](https://cran.r-project.org/package=LongMemoryTS) provides a
    collection of functions for analysing long memory time series.
    Fractionally differenced Gegenbaur ARMA processes are handled by
    [garma](https://cran.r-project.org/package=garma).
    [esemifar](https://cran.r-project.org/package=esemifar) provides tools for
    nonparametric smoothing of long-memory time series.
-   *Transfer function* models are provided by the `arfima` function in
    the [arfima](https://cran.r-project.org/package=arfima) and the
    [tfarima](https://cran.r-project.org/package=tfarima) packages.
-   *Structural (or unobserved component) models* are implemented in
    `StructTS()` in stats and in [stsm](https://cran.r-project.org/package=stsm),
    while automatic modelling and forecasting are provided by
    [UComp](https://cran.r-project.org/package=UComp) and
    [autostsm](https://cran.r-project.org/package=autostsm).
    [KFKSDS](https://cran.r-project.org/package=KFKSDS) provides a naive
    implementation of the Kalman filter and smoothers for univariate
    state space models.
    [statespacer](https://cran.r-project.org/package=statespacer) implements
    univariate state space models including structural and SARIMA
    models. Bayesian structural time series models are implemented in
    [bsts](https://cran.r-project.org/package=bsts) Robust Kalman filtering is
    provided by [RobKF](https://cran.r-project.org/package=RobKF).
-   Non-Gaussian time series can be handled with GLARMA state space
    models via [glarma](https://cran.r-project.org/package=glarma), and using
    Generalized Autoregressive Score models in the
    [GAS](https://cran.r-project.org/package=GAS) package.
    [GlarmaVarSel](https://cran.r-project.org/package=GlarmaVarSel) provides
    variable selection in high-dimensional sparse GLARMA models.
    Conditional Auto-Regression models using Monte Carlo Likelihood
    methods are implemented in [mclcar](https://cran.r-project.org/package=mclcar).
    Efficient Bayesian inference for nonlinear and non-Gaussian state
    space models is provided in [bssm](https://cran.r-project.org/package=bssm).
    Non-Gaussian state space models with exact marginal likelihood are
    given by [NGSSEML](https://cran.r-project.org/package=NGSSEML).
-   *GARCH models* : `garch()` from
    [tseries](https://cran.r-project.org/package=tseries) fits basic GARCH models.
    Many variations on GARCH models are provided by
    [rugarch](https://cran.r-project.org/package=rugarch). Other univariate GARCH
    packages include [fGarch](https://cran.r-project.org/package=fGarch) which
    implements ARIMA models with a wide class of GARCH innovations.
    [bayesforecast](https://cran.r-project.org/package=bayesforecast) fits Bayesian
    time series models including several variations of GARCH models.
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
    [carx](https://cran.r-project.org/package=carx).
    [ARCensReg](https://cran.r-project.org/package=ARCensReg) fits univariate
    censored regression models with autoregressive errors.
-   *Diffusion models* such as Bass and Gompertz curves are provided by
    [diffusion](https://cran.r-project.org/package=diffusion).
-   *Portmanteau tests* are provided via `Box.test()` in the stats
    package. Additional tests are given by
    [portes](https://cran.r-project.org/package=portes),
    [WeightedPortTest](https://cran.r-project.org/package=WeightedPortTest) and
    [testcorr](https://cran.r-project.org/package=testcorr).
-   Outlier detection following the Chen-Liu approach is provided by
    [tsoutliers](https://cran.r-project.org/package=tsoutliers). The `tsoutliers`
    and `tsclean` functions in the
    [forecast](https://cran.r-project.org/package=forecast) package provide some
    simple heuristic methods for identifying and correcting outliers.
    [otsad](https://cran.r-project.org/package=otsad) implements a set of online
    anomaly detectors for time series.
    [tsrobprep](https://cran.r-project.org/package=tsrobprep) provides methods for
    replacing missing values and outliers using a model-based approach.
-   *Change point detection* is provided in
    [strucchange](https://cran.r-project.org/package=strucchange) and
    [strucchangeRcpp](https://cran.r-project.org/package=strucchangeRcpp) (using
    linear regression models) and in
    [trend](https://cran.r-project.org/package=trend) (using nonparametric tests).
    The [changepoint](https://cran.r-project.org/package=changepoint) package
    provides many popular changepoint methods, and
    [ecp](https://cran.r-project.org/package=ecp) does nonparametric changepoint
    detection for univariate and multivariate series.
    [changepoint.np](https://cran.r-project.org/package=changepoint.np) implements
    the nonparametric PELT algorithm,
    [changepoint.mv](https://cran.r-project.org/package=changepoint.mv) detects
    changepoints in multivariate time series, while
    [changepoint.geo](https://cran.r-project.org/package=changepoint.geo) implements
    the high-dimensional changepoint detection method GeomCP.
    [VARDetect](https://cran.r-project.org/package=VARDetect) implements algorithms
    for detecting multiple changes in structural VAR models.
    [InspectChangepoint](https://cran.r-project.org/package=InspectChangepoint) uses
    sparse projection to estimate changepoints in high-dimensional time
    series. [Rbeast](https://cran.r-project.org/package=Rbeast) provides Bayesian
    change-point detection and time series decomposition.
    [breakfast](https://cran.r-project.org/package=breakfast) includes methods for
    fast multiple change-point detection and estimation.
-   Tests for possibly non-monotonic trends are provided by
    [funtimes](https://cran.r-project.org/package=funtimes).
-   *Time series imputation* is provided by the
    [imputeTS](https://cran.r-project.org/package=imputeTS) package. Some more
    limited facilities are available using `na.interp()` from the
    [forecast](https://cran.r-project.org/package=forecast) package.
    [imputeTestbench](https://cran.r-project.org/package=imputeTestbench) provides
    tools for testing and comparing imputation methods.
    [mtsdi](https://cran.r-project.org/package=mtsdi) implements an EM algorithm for
    imputing missing values in multivariate normal time series,
    accounting for spatial and temporal correlations.
-   The [seer](https://cran.r-project.org/package=seer) package implements a
    framework for feature-based forecast model selection.
-   Forecasts can be combined in the
    [fable](https://cran.r-project.org/package=fable) package using simple linear
    expressions. [ForecastComb](https://cran.r-project.org/package=ForecastComb)
    supports many forecast combination methods including simple,
    geometric and regression-based combinations.
    [forecastHybrid](https://cran.r-project.org/package=forecastHybrid) provides
    functions for ensemble forecasts, combining approaches from the
    [forecast](https://cran.r-project.org/package=forecast) package.
    [opera](https://cran.r-project.org/package=opera) has facilities for online
    predictions based on combinations of forecasts provided by the user.
    [profoc](https://cran.r-project.org/package=profoc) combines probabilistic
    forecasts using CRPS learning.
-   Point forecast evaluation is provided in the `accuracy()` function
    from the [fable](https://cran.r-project.org/package=fable) and
    [forecast](https://cran.r-project.org/package=forecast) packages. Distributional
    forecast evaluation using scoring rules is available in
    [fable](https://cran.r-project.org/package=fable),
    [scoringRules](https://cran.r-project.org/package=scoringRules) and
    [scoringutils](https://cran.r-project.org/package=scoringutils). The
    Diebold-Mariano test for comparing the forecast accuracy of two
    models is implemented in the `dm.test()` function in
    [forecast](https://cran.r-project.org/package=forecast). A multivariate version
    of the Diebold-Mariano test is provided by
    [multDM](https://cran.r-project.org/package=multDM).
    [tsutils](https://cran.r-project.org/package=tsutils) implements the Nemenyi
    test for comparing forecasts.
    [greybox](https://cran.r-project.org/package=greybox) provides `ro()` for
    general rolling origin evaluation of forecasts.
-   Tidy tools for forecasting are provided by
    [sweep](https://cran.r-project.org/package=sweep), converting objects produced
    in [forecast](https://cran.r-project.org/package=forecast) to "tidy" data
    frames.
-   Multi-step-ahead direct forecasting with several machine learning
    approaches are provided in
    [forecastML](https://cran.r-project.org/package=forecastML).
-   *Miscellaneous* : [ltsa](https://cran.r-project.org/package=ltsa) contains
    methods for linear time series analysis,
    [timsac](https://cran.r-project.org/package=timsac) for time series analysis and
    control.

**Frequency analysis**

-   *Spectral density estimation* is provided by `spectrum()` in the
    stats package, including the periodogram, smoothed periodogram and
    AR estimates. Bayesian spectral inference is provided by
    [bspec](https://cran.r-project.org/package=bspec) and
    [regspec](https://cran.r-project.org/package=regspec).
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
    multitaper spectral analysis tools. Higher-order spectral analysis
    is implemented in [rhosa](https://cran.r-project.org/package=rhosa), including
    bispectrum, bicoherence, cross-bispectrum and cross-bicoherence.
-   *Wavelet methods* : The [wavelets](https://cran.r-project.org/package=wavelets)
    package includes computing wavelet filters, wavelet transforms and
    multiresolution analyses. Multiresolution forecasting using wavelets
    is also implemented in [mrf](https://cran.r-project.org/package=mrf).
    [WaveletComp](https://cran.r-project.org/package=WaveletComp) provides some
    tools for wavelet-based analysis of univariate and bivariate time
    series including cross-wavelets, phase-difference and significance
    tests. [biwavelet](https://cran.r-project.org/package=biwavelet) is a port of
    the WTC Matlab package for univariate and bivariate wavelet
    analyses. [mvLSW](https://cran.r-project.org/package=mvLSW) provides tools for
    multivariate locally stationary wavelet processes. Tests of white
    noise using wavelets are provided by
    [hwwntest](https://cran.r-project.org/package=hwwntest). Wavelet scalogram tools
    are contained in
    [wavScalogram](https://cran.r-project.org/package=wavScalogram). Further wavelet
    methods can be found in the packages
    [rwt](https://cran.r-project.org/package=rwt),
    [waveslim](https://cran.r-project.org/package=waveslim),
    [wavethresh](https://cran.r-project.org/package=wavethresh).
-   *Harmonic regression* using Fourier terms is implemented in
    [HarmonicRegression](https://cran.r-project.org/package=HarmonicRegression). The
    [fable](https://cran.r-project.org/package=fable) and
    [forecast](https://cran.r-project.org/package=forecast) packages also provide
    some simple harmonic regression facilities via the `fourier`
    function.

**Decomposition and Filtering**

-   *Filters and smoothing* : `filter()` in stats provides
    autoregressive and moving average linear filtering of multiple
    univariate time series. The
    [robfilter](https://cran.r-project.org/package=robfilter) package provides
    several robust time series filters. `smooth()` from the stats
    package computes Tukey's running median smoothers, 3RS3R, 3RSS, 3R,
    etc. [sleekts](https://cran.r-project.org/package=sleekts) computes the 4253H
    twice smoothing method. [mFilter](https://cran.r-project.org/package=mFilter)
    implements several filters for smoothing and extracting trend and
    cyclical components including Hodrick-Prescott and Butterworth
    filters. [smoots](https://cran.r-project.org/package=smoots) provides
    nonparametric estimation of the time trend and its derivatives.
-   *Decomposition* : Seasonal decomposition is discussed below.
    Autoregressive-based decomposition is provided by
    [ArDec](https://cran.r-project.org/package=ArDec).
    [tsdecomp](https://cran.r-project.org/package=tsdecomp) implements ARIMA-based
    decomposition of quarterly and monthly data.
    [rmaf](https://cran.r-project.org/package=rmaf) uses a refined moving average
    filter for decomposition.
-   *Singular Spectrum Analysis* is implemented in
    [Rssa](https://cran.r-project.org/package=Rssa) and
    [ASSA](https://cran.r-project.org/package=ASSA).
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
    [smooth](https://cran.r-project.org/package=smooth) and
    [tsutils](https://cran.r-project.org/package=tsutils) implement extended
    versions of classical decomposition.
-   X-13-ARIMA-SEATS binaries are provided in the
    [x13binary](https://cran.r-project.org/package=x13binary) package, with
    [seasonal](https://cran.r-project.org/package=seasonal) providing an R interface
    and [seasonalview](https://cran.r-project.org/package=seasonalview) providing a
    GUI. An alternative interface is provided by
    [x12](https://cran.r-project.org/package=x12).
-   An interface to the JDemetra+ seasonal adjustment software is
    provided by [RJDemetra](https://cran.r-project.org/package=RJDemetra).
    [ggdemetra](https://cran.r-project.org/package=ggdemetra) provides associated
    ggplot2 functions.
-   Seasonal adjustment of daily time series, allowing for day-of-week,
    time-of-month, time-of-year and holiday effects is provided by
    [dsa](https://cran.r-project.org/package=dsa).
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
-   [sazedR](https://cran.r-project.org/package=sazedR): Method to estimate the
    period of a seasonal time series.

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
    AIC/BIC simultaneous rank-lag selection. Parameter estimation and
    inference in a cointegrating regression are implemented in
    [cointReg](https://cran.r-project.org/package=cointReg).
    [nardl](https://cran.r-project.org/package=nardl) estimates nonlinear
    cointegrating autoregressive distributed lag models.

**Nonlinear Time Series Analysis**

-   *Nonlinear autoregression* : Tools for nonlinear time series
    analysis are provided in [NTS](https://cran.r-project.org/package=NTS) including
    threshold autoregressive models, Markov-switching models,
    convolutional functional autoregressive models, and nonlinearity
    tests. Various forms of nonlinear autoregression are available in
    [tsDyn](https://cran.r-project.org/package=tsDyn) including additive AR, neural
    nets, SETAR and LSTAR models, threshold VAR and VECM. Neural network
    autoregression is also provided in
    [GMDH](https://cran.r-project.org/package=GMDH).
    [nnfor](https://cran.r-project.org/package=nnfor) provides time series
    forecasting with neural networks.
    [NlinTS](https://cran.r-project.org/package=NlinTS) includes neural network VAR,
    and a nonlinear version of the Granger causality test based on
    feedforward neural networks.
    [bentcableAR](https://cran.r-project.org/package=bentcableAR) implements
    Bent-Cable autoregression. [BAYSTAR](https://cran.r-project.org/package=BAYSTAR)
    provides Bayesian analysis of threshold autoregressive models.
    Mixture AR models are implemented in
    [mixAR](https://cran.r-project.org/package=mixAR) and
    [uGMAR](https://cran.r-project.org/package=uGMAR).
-   [tseriesChaos](https://cran.r-project.org/package=tseriesChaos) provides an R
    implementation of the algorithms from the
    *[TISEAN](http://www.mpipks-dresden.mpg.de/~tisean/) project* .
    [DChaos](https://cran.r-project.org/package=DChaos) provides several algorithms
    for detecting chaotic signals inside univariate time series.
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

**Entropy**

-   [RTransferEntropy](https://cran.r-project.org/package=RTransferEntropy) measures
    information flow between time series with Shannon and Renyi transfer
    entropy.
-   An entropy measure based on the Bhattacharya-Hellinger-Matusita
    distance is implemented in
    [tseriesEntropy](https://cran.r-project.org/package=tseriesEntropy).
-   Various approximate and sample entropies are computed using
    [TSEntropies](https://cran.r-project.org/package=TSEntropies).

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
    nonlinear modelling are provided in
    [dlnm](https://cran.r-project.org/package=dlnm).
    [sym.arma](https://cran.r-project.org/package=sym.arma) will fit ARMA models
    with regressors where the observations follow a conditional
    symmetric distribution.
-   *Time-varying parameter* models can be fitted using the
    [tpr](https://cran.r-project.org/package=tpr) package.
-   [greybox](https://cran.r-project.org/package=greybox) provides several tools for
    modelling and forecasting with dynamic regression models.
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
    structural VARs. Shrinkage estimation methods for VARs are
    implemented in [VARshrink](https://cran.r-project.org/package=VARshrink). More
    elaborate models are provided in package
    [vars](https://cran.r-project.org/package=vars),
    [tsDyn](https://cran.r-project.org/package=tsDyn), `estVARXls()` in
    [dse](https://cran.r-project.org/package=dse). Another implementation with
    bootstrapped prediction intervals is given in
    [VAR.etp](https://cran.r-project.org/package=VAR.etp).
    [bvartools](https://cran.r-project.org/package=bvartools) assists in the set-up
    of Bayesian VAR models, while [BMTAR](https://cran.r-project.org/package=BMTAR)
    implements Baysian Multivariate Threshold AR models with missing
    data. [mfbvar](https://cran.r-project.org/package=mfbvar) includes tools for
    estimating mixed-frequency Bayesian VAR models.
    [BVAR](https://cran.r-project.org/package=BVAR) provides a toolkit for
    hierarchical Bayesian VAR models.
    [BGVAR](https://cran.r-project.org/package=BGVAR) implements Bayesian Global VAR
    models. [mlVAR](https://cran.r-project.org/package=mlVAR) provides multi-level
    vector autoregression. [VARsignR](https://cran.r-project.org/package=VARsignR)
    provides routines for identifying structural shocks in VAR models
    using sign restrictions. [gmvarkit](https://cran.r-project.org/package=gmvarkit)
    estimates Gaussian mixture VAR models.
    [GNAR](https://cran.r-project.org/package=GNAR) provides methods for fitting
    network AR models, while
    [graphicalVAR](https://cran.r-project.org/package=graphicalVAR) estimates
    graphical VAR models. [gdpc](https://cran.r-project.org/package=gdpc) implements
    generalized dynamic principal components.
    [pcdpca](https://cran.r-project.org/package=pcdpca) extends dynamic principal
    components to periodically correlated multivariate time series.
    [onlineVAR](https://cran.r-project.org/package=onlineVAR) implements online
    fitting of time-adaptive lasso VARs.
    [mgm](https://cran.r-project.org/package=mgm) estimates time-varying mixed
    graphical models and mixed VAR models via regularized regression.
-   *VARIMA models* and *state space models* are provided in the
    [dse](https://cran.r-project.org/package=dse) package.
-   *Vector error correction models* are available via the
    [urca](https://cran.r-project.org/package=urca),
    [ecm](https://cran.r-project.org/package=ecm),
    [vars](https://cran.r-project.org/package=vars),
    [tsDyn](https://cran.r-project.org/package=tsDyn) packages, including versions
    with structural constraints and thresholding.
-   *Vector exponential smoothing* is provided by
    [smooth](https://cran.r-project.org/package=smooth).
-   *Time series component analysis* :
    [ForeCA](https://cran.r-project.org/package=ForeCA) implements forecastable
    component analysis by searching for the best linear transformations
    that make a multivariate time series as forecastable as possible.
    [PCA4TS](https://cran.r-project.org/package=PCA4TS) finds a linear
    transformation of a multivariate time series giving
    lower-dimensional subseries that are uncorrelated with each other.
    [HDTSA](https://cran.r-project.org/package=HDTSA) provides procedures for
    several high-dimensional time series analysis tools. One-sided
    dynamic principal components are computed in
    [odpc](https://cran.r-project.org/package=odpc). Frequency-domain-based dynamic
    PCA is implemented in [freqdom](https://cran.r-project.org/package=freqdom).
    [tsBSS](https://cran.r-project.org/package=tsBSS) provides blind source
    separation and supervised dimension reduction for time series.
-   *Multivariate state space models* An implementation is provided by
    the [KFAS](https://cran.r-project.org/package=KFAS) package which provides a
    fast multivariate Kalman filter, smoother, simulation smoother and
    forecasting. [FKF](https://cran.r-project.org/package=FKF) provides a fast and
    flexible implementation of the Kalman filter, which can deal with
    missing values. [FKF.SP](https://cran.r-project.org/package=FKF.SP) implements
    fast Kalman filtering through sequential processing. Another
    implementation is given in the [dlm](https://cran.r-project.org/package=dlm)
    package which also contains tools for converting other multivariate
    models into state space form. [mssm](https://cran.r-project.org/package=mssm)
    also provides methods for multivariate state space models.
    [MARSS](https://cran.r-project.org/package=MARSS) fits constrained and
    unconstrained multivariate autoregressive state-space models using
    an EM algorithm. [mbsts](https://cran.r-project.org/package=mbsts) provides
    tools for multivariate Bayesian structural time series models. All
    of these packages assume the observational and state error terms are
    uncorrelated.
-   *Partially-observed Markov processes* are a generalization of the
    usual linear multivariate state space models, allowing non-Gaussian
    and nonlinear models. These are implemented in the
    [pomp](https://cran.r-project.org/package=pomp) package.
-   Multivariate stochastic volatility models (using latent factors) are
    provided by [factorstochvol](https://cran.r-project.org/package=factorstochvol).

**Analysis of large groups of time series**

-   *Time series features* are computed in
    [feasts](https://cran.r-project.org/package=feasts) for time series in `tsibble`
    format. They are computed using
    [tsfeatures](https://cran.r-project.org/package=tsfeatures) for a list or matrix
    of time series in `ts` format. In both packages, many built-in
    feature functions are included, and users can add their own.
    [Rcatch22](https://cran.r-project.org/package=Rcatch22) provides fast
    computation of 22 features identified as particularly useful.
    [fsMTS](https://cran.r-project.org/package=fsMTS) implements feature selection
    routines for multivariate time series.
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
-   [rucrdtw](https://cran.r-project.org/package=rucrdtw) provides R bindings for
    functions from the UCR Suite to enable ultrafast subsequence search
    for a best match under Dynamic Time Warping and Euclidean Distance.
    [IncDTW](https://cran.r-project.org/package=IncDTW) provides incremental
    calculation of dynamic time warping for streaming time series.
-   Methods for plotting and forecasting collections of hierarchical and
    grouped time series are provided by
    [fable](https://cran.r-project.org/package=fable) and
    [hts](https://cran.r-project.org/package=hts).
    [thief](https://cran.r-project.org/package=thief) uses hierarchical methods to
    reconcile forecasts of temporally aggregated time series.
    [FoReco](https://cran.r-project.org/package=FoReco) provides various forecast
    reconciliation methods for cross-sectional, temporal, and
    cross-temporal constrained time series. An alternative approach to
    reconciling forecasts of hierarchical time series is provided by
    [gtop](https://cran.r-project.org/package=gtop).
    [ProbReco](https://cran.r-project.org/package=ProbReco) provides tools to train
    forecast reconciliation weights by optimizing probability scoring
    functions.

**Functional time series**

-   Tools for visualizing, modeling, forecasting and analysing
    functional time series are implemented in
    [ftsa](https://cran.r-project.org/package=ftsa).
    [NTS](https://cran.r-project.org/package=NTS) also implements functional
    autoregressive models. Seasonal functional autoregression models are
    provided by [Rsfar](https://cran.r-project.org/package=Rsfar).
    [fpcb](https://cran.r-project.org/package=fpcb) implements predictive confidence
    bands for functional time series.
-   [fdaACF](https://cran.r-project.org/package=fdaACF) estimates the
    autocorrelation function for functional time series.
-   [freqdom.fda](https://cran.r-project.org/package=freqdom.fda) provides
    implements of dynamical functional principal components for
    functional time series.
-   [STFTS](https://cran.r-project.org/package=STFTS) contains stationarity, trend
    and unit root tests for functional time series.
-   [wwntests](https://cran.r-project.org/package=wwntests) provides an array of
    white noise hypothesis tests for functional data.

**Matrix and tensor-valued time series**

-   [tensorTS](https://cran.r-project.org/package=tensorTS) provides functions for
    estimation, simulation and prediction of factor and autoregressive
    models for matrix and tensor valued time series.

**Continuous time models**

-   [carfima](https://cran.r-project.org/package=carfima) allows for continuous-time
    ARFIMA models.
-   [Sim.DiffProc](https://cran.r-project.org/package=Sim.DiffProc) simulates and
    models stochastic differential equations.
-   Simulation and inference for stochastic differential equations is
    provided by [sde](https://cran.r-project.org/package=sde) and
    [yuima](https://cran.r-project.org/package=yuima).

**Resampling**

-   *Bootstrapping* : The [boot](https://cran.r-project.org/package=boot) package
    provides function `tsboot()` for time series bootstrapping,
    including block bootstrap with several variants.
    [blocklength](https://cran.r-project.org/package=blocklength) allows for
    selecting the optimal block-length for a dependent bootstrap.
    `tsbootstrap()` from [tseries](https://cran.r-project.org/package=tseries)
    provides fast stationary and block bootstrapping. Maximum entropy
    bootstrap for time series is available in
    [meboot](https://cran.r-project.org/package=meboot).
    [timesboot](https://cran.r-project.org/package=timesboot) computes the bootstrap
    CI for the sample ACF and periodogram.
    [BootPR](https://cran.r-project.org/package=BootPR) computes bias-corrected
    forecasting and bootstrap prediction intervals for autoregressive
    time series. [bootUR](https://cran.r-project.org/package=bootUR) implements
    bootstrap unit root tests.

**Time Series Data**

-   Various data sets in [tsibble](https://cran.r-project.org/package=tsibble)
    format are provided by
    [tsibbledata](https://cran.r-project.org/package=tsibbledata).
-   Data from Cryer and Chan (2010, 2nd ed) *Time series analysis with
    applications in R* are in the [TSA](https://cran.r-project.org/package=TSA)
    package.
-   Data from Hyndman and Athanasopoulos (2018, 2nd ed) *Forecasting:
    principles and practice* are in the
    [fpp2](https://cran.r-project.org/package=fpp2) package.
-   Data from Hyndman and Athanasopoulos (2021, 3rd ed) *Forecasting:
    principles and practice* are in the
    [fpp3](https://cran.r-project.org/package=fpp3) package.
-   Data from Hyndman, Koehler, Ord and Snyder (2008) *Forecasting with
    exponential smoothing* are in the
    [expsmooth](https://cran.r-project.org/package=expsmooth) package.
-   Data from Makridakis, Wheelwright and Hyndman (1998, 3rd ed)
    *Forecasting: methods and applications* are in the
    [fma](https://cran.r-project.org/package=fma) package.
-   Data from Shumway and Stoffer (2017, 4th ed) *Time Series Analysis
    and Its Applications: With R Examples* are in the
    [astsa](https://cran.r-project.org/package=astsa) package.
-   Data from Tsay (2005, 2nd ed) *Analysis of Financial Time Series*
    are in the [FinTS](https://cran.r-project.org/package=FinTS) package.
-   Data from Woodward, Gray, and Elliott (2016, 2nd ed) *Applied Time
    Series Analysis with R* are in the
    [tswge](https://cran.r-project.org/package=tswge) package.
-   [AER](https://cran.r-project.org/package=AER) and
    [Ecdat](https://cran.r-project.org/package=Ecdat) both contain many data sets
    (including time series data) from many econometrics text books
-   Data from the M-competition and M3-competition are provided in the
    [Mcomp](https://cran.r-project.org/package=Mcomp) package.
    [Tcomp](https://cran.r-project.org/package=Tcomp) provides data from the 2010
    IJF Tourism Forecasting Competition.
-   *National time series data:*
    [readabs](https://cran.r-project.org/package=readabs) downloads, imports and
    tidies time series data from the [*Australian* Bureau of
    Statistics](https://www.abs.gov.au) .
    [BETS](https://cran.r-project.org/package=BETS) provides access to the most
    important economic time series in *Brazil* .
    [bundesbank](https://cran.r-project.org/package=bundesbank) allows access to the
    time series databases of the *German* central bank. Data from
    *Switzerland* via [dataseries.org](http://dataseries.org) can be
    downloaded and imported using
    [dataseries](https://cran.r-project.org/package=dataseries).
    [ugatsdb](https://cran.r-project.org/package=ugatsdb) provides an API to access
    time series data for *Uganda* .
-   *Time series data bases:* [fame](https://cran.r-project.org/package=fame)
    provides an interface for FAME time series databases. Economic time
    series and other data from FRED (the Federal Reserve Economic Data)
    can be retrieved using [fredr](https://cran.r-project.org/package=fredr).
    [influxdbr](https://cran.r-project.org/package=influxdbr) provides an interface
    to the InfluxDB time series database.
    [pdfetch](https://cran.r-project.org/package=pdfetch) provides facilities for
    downloading economic and financial time series from public sources.
    Data from the [Quandl](http://www.quandl.com) online portal to
    financial, economical and social datasets can be queried
    interactively using the [Quandl](https://cran.r-project.org/package=Quandl)
    package. [TSdbi](https://cran.r-project.org/package=TSdbi) provides a common
    interface to time series databases.
    [tsdb](https://cran.r-project.org/package=tsdb) implements a simple data base
    for numerical time series.
-   *Synthetic data* are produced by `simulate()` in
    [forecast](https://cran.r-project.org/package=forecast) package or `generate()`
    in [fable](https://cran.r-project.org/package=fable), given a specific model.
    [gratis](https://cran.r-project.org/package=gratis) generates new time series
    with diverse and controllable characteristics using mixture
    autoregression models. [synthesis](https://cran.r-project.org/package=synthesis)
    generates synthetic time series from commonly used statistical
    models, including linear, nonlinear and chaotic systems.
    [tssim](https://cran.r-project.org/package=tssim) flexibly simulates daily or
    monthly time series using seasonal, calendar, and outlier
    components.

**Miscellaneous**

-   [dtw](https://cran.r-project.org/package=dtw): Dynamic time warping algorithms
    for computing and plotting pairwise alignments between time series.
-   [EBMAforecast](https://cran.r-project.org/package=EBMAforecast): Ensemble
    Bayesian model averaging forecasts using Gibbs sampling or EM
    algorithms.
-   [ensembleBMA](https://cran.r-project.org/package=ensembleBMA): Bayesian Model
    Averaging to create probabilistic forecasts from ensemble forecasts
    and weather observations.
-   [earlywarnings](https://cran.r-project.org/package=earlywarnings): Early
    warnings signals toolbox for detecting critical transitions in time
    series
-   [FeedbackTS](https://cran.r-project.org/package=FeedbackTS): Analysis of
    fragmented time directionality to investigate feedback in time
    series.
-   [gsignal](https://cran.r-project.org/package=gsignal) is an R implementation of
    the Octave package "signal", containing a variety of signal
    processing tools.
-   [ifultools](https://cran.r-project.org/package=ifultools) is a collection of
    efficient signal processing, image processing, and time series
    modeling routines.
-   [imputePSF](https://cran.r-project.org/package=imputePSF): imputes missing data
    using pattern sequences.
-   [LPStimeSeries](https://cran.r-project.org/package=LPStimeSeries) aims to find
    "learned pattern similarity" for time series.
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
-   [spTimer](https://cran.r-project.org/package=spTimer): Spatio-temporal Bayesian
    modelling.
-   [surveillance](https://cran.r-project.org/package=surveillance): Temporal and
    spatio-temporal modeling and monitoring of epidemic phenomena.
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

</div>

### CRAN packages:

-   [acp](https://cran.r-project.org/package=acp)
-   [AER](https://cran.r-project.org/package=AER)
-   [ARCensReg](https://cran.r-project.org/package=ARCensReg)
-   [ArDec](https://cran.r-project.org/package=ArDec)
-   [arfima](https://cran.r-project.org/package=arfima)
-   [ASSA](https://cran.r-project.org/package=ASSA)
-   [astsa](https://cran.r-project.org/package=astsa)
-   [autostsm](https://cran.r-project.org/package=autostsm)
-   [BayesARIMAX](https://cran.r-project.org/package=BayesARIMAX)
-   [bayesforecast](https://cran.r-project.org/package=bayesforecast)
-   [BAYSTAR](https://cran.r-project.org/package=BAYSTAR)
-   [bentcableAR](https://cran.r-project.org/package=bentcableAR)
-   [BETS](https://cran.r-project.org/package=BETS)
-   [bfast](https://cran.r-project.org/package=bfast)
-   [BGVAR](https://cran.r-project.org/package=BGVAR)
-   [bigtime](https://cran.r-project.org/package=bigtime)
-   [BigVAR](https://cran.r-project.org/package=BigVAR)
-   [biwavelet](https://cran.r-project.org/package=biwavelet)
-   [blocklength](https://cran.r-project.org/package=blocklength)
-   [BMTAR](https://cran.r-project.org/package=BMTAR)
-   [BNPTSclust](https://cran.r-project.org/package=BNPTSclust)
-   [boot](https://cran.r-project.org/package=boot)
-   [BootPR](https://cran.r-project.org/package=BootPR)
-   [bootUR](https://cran.r-project.org/package=bootUR)
-   [breakfast](https://cran.r-project.org/package=breakfast)
-   [brolgar](https://cran.r-project.org/package=brolgar)
-   [bspec](https://cran.r-project.org/package=bspec)
-   [bssm](https://cran.r-project.org/package=bssm)
-   [bsts](https://cran.r-project.org/package=bsts)
-   [bundesbank](https://cran.r-project.org/package=bundesbank)
-   [BVAR](https://cran.r-project.org/package=BVAR)
-   [bvartools](https://cran.r-project.org/package=bvartools)
-   [CADFtest](https://cran.r-project.org/package=CADFtest)
-   [carfima](https://cran.r-project.org/package=carfima)
-   [carx](https://cran.r-project.org/package=carx)
-   [changepoint](https://cran.r-project.org/package=changepoint)
-   [changepoint.geo](https://cran.r-project.org/package=changepoint.geo)
-   [changepoint.mv](https://cran.r-project.org/package=changepoint.mv)
-   [changepoint.np](https://cran.r-project.org/package=changepoint.np)
-   [chron](https://cran.r-project.org/package=chron)
-   [clock](https://cran.r-project.org/package=clock)
-   [cointReg](https://cran.r-project.org/package=cointReg)
-   [collapse](https://cran.r-project.org/package=collapse)
-   [costat](https://cran.r-project.org/package=costat)
-   [data.table](https://cran.r-project.org/package=data.table)
-   [dataseries](https://cran.r-project.org/package=dataseries)
-   [DChaos](https://cran.r-project.org/package=DChaos)
-   [dCovTS](https://cran.r-project.org/package=dCovTS)
-   [depmix](https://cran.r-project.org/package=depmix)
-   [depmixS4](https://cran.r-project.org/package=depmixS4)
-   [deseasonalize](https://cran.r-project.org/package=deseasonalize)
-   [diffusion](https://cran.r-project.org/package=diffusion)
-   [disaggR](https://cran.r-project.org/package=disaggR)
-   [dLagM](https://cran.r-project.org/package=dLagM)
-   [dlm](https://cran.r-project.org/package=dlm)
-   [dlnm](https://cran.r-project.org/package=dlnm)
-   [dsa](https://cran.r-project.org/package=dsa)
-   [dse](https://cran.r-project.org/package=dse)
-   [DTSg](https://cran.r-project.org/package=DTSg)
-   [dtw](https://cran.r-project.org/package=dtw)
-   [dtwclust](https://cran.r-project.org/package=dtwclust)
-   [dygraphs](https://cran.r-project.org/package=dygraphs)
-   [dyn](https://cran.r-project.org/package=dyn)
-   [dynlm](https://cran.r-project.org/package=dynlm)
-   [dynr](https://cran.r-project.org/package=dynr)
-   [earlywarnings](https://cran.r-project.org/package=earlywarnings)
-   [EBMAforecast](https://cran.r-project.org/package=EBMAforecast)
-   [Ecdat](https://cran.r-project.org/package=Ecdat)
-   [ecm](https://cran.r-project.org/package=ecm)
-   [ecp](https://cran.r-project.org/package=ecp)
-   [EMD](https://cran.r-project.org/package=EMD)
-   [ensembleBMA](https://cran.r-project.org/package=ensembleBMA)
-   [esemifar](https://cran.r-project.org/package=esemifar)
-   [expsmooth](https://cran.r-project.org/package=expsmooth)
-   [fable](https://cran.r-project.org/package=fable) (core)
-   [fable.prophet](https://cran.r-project.org/package=fable.prophet)
-   [fabletools](https://cran.r-project.org/package=fabletools)
-   [factorstochvol](https://cran.r-project.org/package=factorstochvol)
-   [fame](https://cran.r-project.org/package=fame)
-   [fanplot](https://cran.r-project.org/package=fanplot)
-   [fdaACF](https://cran.r-project.org/package=fdaACF)
-   [feasts](https://cran.r-project.org/package=feasts) (core)
-   [FeedbackTS](https://cran.r-project.org/package=FeedbackTS)
-   [fGarch](https://cran.r-project.org/package=fGarch)
-   [FinTS](https://cran.r-project.org/package=FinTS)
-   [FitAR](https://cran.r-project.org/package=FitAR)
-   [FitARMA](https://cran.r-project.org/package=FitARMA)
-   [FKF](https://cran.r-project.org/package=FKF)
-   [FKF.SP](https://cran.r-project.org/package=FKF.SP)
-   [fma](https://cran.r-project.org/package=fma)
-   [fNonlinear](https://cran.r-project.org/package=fNonlinear)
-   [ForeCA](https://cran.r-project.org/package=ForeCA)
-   [forecast](https://cran.r-project.org/package=forecast) (core)
-   [ForecastComb](https://cran.r-project.org/package=ForecastComb)
-   [forecastHybrid](https://cran.r-project.org/package=forecastHybrid)
-   [forecastML](https://cran.r-project.org/package=forecastML)
-   [FoReco](https://cran.r-project.org/package=FoReco)
-   [forecTheta](https://cran.r-project.org/package=forecTheta)
-   [fpcb](https://cran.r-project.org/package=fpcb)
-   [fpp2](https://cran.r-project.org/package=fpp2)
-   [fpp3](https://cran.r-project.org/package=fpp3)
-   [fracdiff](https://cran.r-project.org/package=fracdiff)
-   [fredr](https://cran.r-project.org/package=fredr)
-   [freqdom](https://cran.r-project.org/package=freqdom)
-   [freqdom.fda](https://cran.r-project.org/package=freqdom.fda)
-   [fsMTS](https://cran.r-project.org/package=fsMTS)
-   [fts](https://cran.r-project.org/package=fts)
-   [ftsa](https://cran.r-project.org/package=ftsa)
-   [funtimes](https://cran.r-project.org/package=funtimes)
-   [garma](https://cran.r-project.org/package=garma)
-   [GAS](https://cran.r-project.org/package=GAS)
-   [gdpc](https://cran.r-project.org/package=gdpc)
-   [ggdemetra](https://cran.r-project.org/package=ggdemetra)
-   [ggseas](https://cran.r-project.org/package=ggseas)
-   [glarma](https://cran.r-project.org/package=glarma)
-   [GlarmaVarSel](https://cran.r-project.org/package=GlarmaVarSel)
-   [GMDH](https://cran.r-project.org/package=GMDH)
-   [gmvarkit](https://cran.r-project.org/package=gmvarkit)
-   [GNAR](https://cran.r-project.org/package=GNAR)
-   [graphicalVAR](https://cran.r-project.org/package=graphicalVAR)
-   [gratis](https://cran.r-project.org/package=gratis)
-   [gravitas](https://cran.r-project.org/package=gravitas)
-   [greybox](https://cran.r-project.org/package=greybox)
-   [gsarima](https://cran.r-project.org/package=gsarima)
-   [gsignal](https://cran.r-project.org/package=gsignal)
-   [gtop](https://cran.r-project.org/package=gtop)
-   [HarmonicRegression](https://cran.r-project.org/package=HarmonicRegression)
-   [HDTSA](https://cran.r-project.org/package=HDTSA)
-   [hht](https://cran.r-project.org/package=hht)
-   [hts](https://cran.r-project.org/package=hts)
-   [hwwntest](https://cran.r-project.org/package=hwwntest)
-   [ifultools](https://cran.r-project.org/package=ifultools)
-   [imputePSF](https://cran.r-project.org/package=imputePSF)
-   [imputeTestbench](https://cran.r-project.org/package=imputeTestbench)
-   [imputeTS](https://cran.r-project.org/package=imputeTS)
-   [IncDTW](https://cran.r-project.org/package=IncDTW)
-   [influxdbr](https://cran.r-project.org/package=influxdbr)
-   [InspectChangepoint](https://cran.r-project.org/package=InspectChangepoint)
-   [itsmr](https://cran.r-project.org/package=itsmr)
-   [KFAS](https://cran.r-project.org/package=KFAS)
-   [KFKSDS](https://cran.r-project.org/package=KFKSDS)
-   [kza](https://cran.r-project.org/package=kza)
-   [legion](https://cran.r-project.org/package=legion)
-   [locits](https://cran.r-project.org/package=locits)
-   [lomb](https://cran.r-project.org/package=lomb)
-   [LongMemoryTS](https://cran.r-project.org/package=LongMemoryTS)
-   [LPStimeSeries](https://cran.r-project.org/package=LPStimeSeries)
-   [LSTS](https://cran.r-project.org/package=LSTS)
-   [ltsa](https://cran.r-project.org/package=ltsa)
-   [lubridate](https://cran.r-project.org/package=lubridate)
-   [MAPA](https://cran.r-project.org/package=MAPA)
-   [mAr](https://cran.r-project.org/package=mAr)
-   [mar1s](https://cran.r-project.org/package=mar1s)
-   [MARSS](https://cran.r-project.org/package=MARSS)
-   [mbsts](https://cran.r-project.org/package=mbsts)
-   [mclcar](https://cran.r-project.org/package=mclcar)
-   [Mcomp](https://cran.r-project.org/package=Mcomp)
-   [meboot](https://cran.r-project.org/package=meboot)
-   [mfbvar](https://cran.r-project.org/package=mfbvar)
-   [mFilter](https://cran.r-project.org/package=mFilter)
-   [mgm](https://cran.r-project.org/package=mgm)
-   [mixAR](https://cran.r-project.org/package=mixAR)
-   [mlVAR](https://cran.r-project.org/package=mlVAR)
-   [modeltime](https://cran.r-project.org/package=modeltime)
-   [modeltime.ensemble](https://cran.r-project.org/package=modeltime.ensemble)
-   [mondate](https://cran.r-project.org/package=mondate)
-   [mrf](https://cran.r-project.org/package=mrf)
-   [mssm](https://cran.r-project.org/package=mssm)
-   [MSwM](https://cran.r-project.org/package=MSwM)
-   [MTS](https://cran.r-project.org/package=MTS)
-   [mtsdi](https://cran.r-project.org/package=mtsdi)
-   [multDM](https://cran.r-project.org/package=multDM)
-   [MultipleBubbles](https://cran.r-project.org/package=MultipleBubbles)
-   [multitaper](https://cran.r-project.org/package=multitaper)
-   [mvLSW](https://cran.r-project.org/package=mvLSW)
-   [nardl](https://cran.r-project.org/package=nardl)
-   [nets](https://cran.r-project.org/package=nets)
-   [NGSSEML](https://cran.r-project.org/package=NGSSEML)
-   [NlinTS](https://cran.r-project.org/package=NlinTS)
-   [nlts](https://cran.r-project.org/package=nlts)
-   [nnfor](https://cran.r-project.org/package=nnfor)
-   [nonlinearTseries](https://cran.r-project.org/package=nonlinearTseries)
-   [npst](https://cran.r-project.org/package=npst)
-   [nsarfima](https://cran.r-project.org/package=nsarfima)
-   [NTS](https://cran.r-project.org/package=NTS)
-   [odpc](https://cran.r-project.org/package=odpc)
-   [onlineVAR](https://cran.r-project.org/package=onlineVAR)
-   [opera](https://cran.r-project.org/package=opera)
-   [otsad](https://cran.r-project.org/package=otsad)
-   [paleoTS](https://cran.r-project.org/package=paleoTS)
-   [partsm](https://cran.r-project.org/package=partsm)
-   [pastecs](https://cran.r-project.org/package=pastecs)
-   [PCA4TS](https://cran.r-project.org/package=PCA4TS)
-   [pcdpca](https://cran.r-project.org/package=pcdpca)
-   [pcts](https://cran.r-project.org/package=pcts)
-   [pdc](https://cran.r-project.org/package=pdc)
-   [pdfetch](https://cran.r-project.org/package=pdfetch)
-   [perARMA](https://cran.r-project.org/package=perARMA)
-   [pomp](https://cran.r-project.org/package=pomp)
-   [portes](https://cran.r-project.org/package=portes)
-   [ProbReco](https://cran.r-project.org/package=ProbReco)
-   [profoc](https://cran.r-project.org/package=profoc)
-   [prophet](https://cran.r-project.org/package=prophet)
-   [psd](https://cran.r-project.org/package=psd)
-   [PSF](https://cran.r-project.org/package=PSF)
-   [ptw](https://cran.r-project.org/package=ptw)
-   [Quandl](https://cran.r-project.org/package=Quandl)
-   [quantspec](https://cran.r-project.org/package=quantspec)
-   [Rbeast](https://cran.r-project.org/package=Rbeast)
-   [Rcatch22](https://cran.r-project.org/package=Rcatch22)
-   [readabs](https://cran.r-project.org/package=readabs)
-   [regspec](https://cran.r-project.org/package=regspec)
-   [RGENERATE](https://cran.r-project.org/package=RGENERATE)
-   [rhosa](https://cran.r-project.org/package=rhosa)
-   [RJDemetra](https://cran.r-project.org/package=RJDemetra)
-   [Rlgt](https://cran.r-project.org/package=Rlgt)
-   [Rlibeemd](https://cran.r-project.org/package=Rlibeemd)
-   [rmaf](https://cran.r-project.org/package=rmaf)
-   [RMAWGEN](https://cran.r-project.org/package=RMAWGEN)
-   [robets](https://cran.r-project.org/package=robets)
-   [robfilter](https://cran.r-project.org/package=robfilter)
-   [RobKF](https://cran.r-project.org/package=RobKF)
-   [robustarima](https://cran.r-project.org/package=robustarima)
-   [roll](https://cran.r-project.org/package=roll)
-   [RSEIS](https://cran.r-project.org/package=RSEIS)
-   [Rsfar](https://cran.r-project.org/package=Rsfar)
-   [Rssa](https://cran.r-project.org/package=Rssa)
-   [RTransferEntropy](https://cran.r-project.org/package=RTransferEntropy)
-   [rts](https://cran.r-project.org/package=rts)
-   [rucrdtw](https://cran.r-project.org/package=rucrdtw)
-   [rugarch](https://cran.r-project.org/package=rugarch)
-   [runner](https://cran.r-project.org/package=runner)
-   [runstats](https://cran.r-project.org/package=runstats)
-   [rwt](https://cran.r-project.org/package=rwt)
-   [sazedR](https://cran.r-project.org/package=sazedR)
-   [scoringRules](https://cran.r-project.org/package=scoringRules)
-   [scoringutils](https://cran.r-project.org/package=scoringutils)
-   [SDD](https://cran.r-project.org/package=SDD)
-   [sde](https://cran.r-project.org/package=sde)
-   [seas](https://cran.r-project.org/package=seas)
-   [season](https://cran.r-project.org/package=season)
-   [seasonal](https://cran.r-project.org/package=seasonal)
-   [seasonalview](https://cran.r-project.org/package=seasonalview)
-   [seer](https://cran.r-project.org/package=seer)
-   [Sim.DiffProc](https://cran.r-project.org/package=Sim.DiffProc)
-   [sleekts](https://cran.r-project.org/package=sleekts)
-   [slider](https://cran.r-project.org/package=slider)
-   [smooth](https://cran.r-project.org/package=smooth)
-   [smoots](https://cran.r-project.org/package=smoots)
-   [sparsevar](https://cran.r-project.org/package=sparsevar)
-   [spectral](https://cran.r-project.org/package=spectral)
-   [spTimer](https://cran.r-project.org/package=spTimer)
-   [statespacer](https://cran.r-project.org/package=statespacer)
-   [STFTS](https://cran.r-project.org/package=STFTS)
-   [stlplus](https://cran.r-project.org/package=stlplus)
-   [stochvol](https://cran.r-project.org/package=stochvol)
-   [stR](https://cran.r-project.org/package=stR)
-   [strucchange](https://cran.r-project.org/package=strucchange)
-   [strucchangeRcpp](https://cran.r-project.org/package=strucchangeRcpp)
-   [stsm](https://cran.r-project.org/package=stsm)
-   [sugrrants](https://cran.r-project.org/package=sugrrants)
-   [surveillance](https://cran.r-project.org/package=surveillance)
-   [svars](https://cran.r-project.org/package=svars)
-   [sweep](https://cran.r-project.org/package=sweep)
-   [sym.arma](https://cran.r-project.org/package=sym.arma)
-   [synthesis](https://cran.r-project.org/package=synthesis)
-   [tbrf](https://cran.r-project.org/package=tbrf)
-   [Tcomp](https://cran.r-project.org/package=Tcomp)
-   [tempdisagg](https://cran.r-project.org/package=tempdisagg)
-   [tensorTS](https://cran.r-project.org/package=tensorTS)
-   [testcorr](https://cran.r-project.org/package=testcorr)
-   [tfarima](https://cran.r-project.org/package=tfarima)
-   [tframe](https://cran.r-project.org/package=tframe)
-   [thief](https://cran.r-project.org/package=thief)
-   [Tides](https://cran.r-project.org/package=Tides)
-   [tiger](https://cran.r-project.org/package=tiger)
-   [timechange](https://cran.r-project.org/package=timechange)
-   [timeDate](https://cran.r-project.org/package=timeDate)
-   [TimeProjection](https://cran.r-project.org/package=TimeProjection)
-   [timesboot](https://cran.r-project.org/package=timesboot)
-   [timeSeries](https://cran.r-project.org/package=timeSeries)
-   [timeseriesdb](https://cran.r-project.org/package=timeseriesdb)
-   [timetk](https://cran.r-project.org/package=timetk)
-   [timsac](https://cran.r-project.org/package=timsac)
-   [tis](https://cran.r-project.org/package=tis)
-   [tpr](https://cran.r-project.org/package=tpr)
-   [trend](https://cran.r-project.org/package=trend)
-   [TSA](https://cran.r-project.org/package=TSA)
-   [tsbox](https://cran.r-project.org/package=tsbox)
-   [tsBSS](https://cran.r-project.org/package=tsBSS)
-   [TSclust](https://cran.r-project.org/package=TSclust)
-   [tscount](https://cran.r-project.org/package=tscount)
-   [tsdb](https://cran.r-project.org/package=tsdb)
-   [TSdbi](https://cran.r-project.org/package=TSdbi)
-   [tsdecomp](https://cran.r-project.org/package=tsdecomp)
-   [tsdisagg2](https://cran.r-project.org/package=tsdisagg2)
-   [TSdist](https://cran.r-project.org/package=TSdist)
-   [tsDyn](https://cran.r-project.org/package=tsDyn)
-   [TSEntropies](https://cran.r-project.org/package=TSEntropies)
-   [tseries](https://cran.r-project.org/package=tseries) (core)
-   [tseriesChaos](https://cran.r-project.org/package=tseriesChaos)
-   [tseriesEntropy](https://cran.r-project.org/package=tseriesEntropy)
-   [tsfeatures](https://cran.r-project.org/package=tsfeatures)
-   [tsfknn](https://cran.r-project.org/package=tsfknn)
-   [tsibble](https://cran.r-project.org/package=tsibble) (core)
-   [tsibbledata](https://cran.r-project.org/package=tsibbledata)
-   [tsibbletalk](https://cran.r-project.org/package=tsibbletalk)
-   [tsintermittent](https://cran.r-project.org/package=tsintermittent)
-   [TSMining](https://cran.r-project.org/package=TSMining)
-   [tsModel](https://cran.r-project.org/package=tsModel)
-   [tsoutliers](https://cran.r-project.org/package=tsoutliers)
-   [tsPI](https://cran.r-project.org/package=tsPI)
-   [TSrepr](https://cran.r-project.org/package=TSrepr)
-   [tsrobprep](https://cran.r-project.org/package=tsrobprep)
-   [tssim](https://cran.r-project.org/package=tssim)
-   [TSstudio](https://cran.r-project.org/package=TSstudio)
-   [TSTutorial](https://cran.r-project.org/package=TSTutorial)
-   [tsutils](https://cran.r-project.org/package=tsutils)
-   [tswge](https://cran.r-project.org/package=tswge)
-   [UComp](https://cran.r-project.org/package=UComp)
-   [ugatsdb](https://cran.r-project.org/package=ugatsdb)
-   [uGMAR](https://cran.r-project.org/package=uGMAR)
-   [urca](https://cran.r-project.org/package=urca)
-   [uroot](https://cran.r-project.org/package=uroot)
-   [VAR.etp](https://cran.r-project.org/package=VAR.etp)
-   [VARDetect](https://cran.r-project.org/package=VARDetect)
-   [vars](https://cran.r-project.org/package=vars)
-   [VARshrink](https://cran.r-project.org/package=VARshrink)
-   [VARsignR](https://cran.r-project.org/package=VARsignR)
-   [WaveletComp](https://cran.r-project.org/package=WaveletComp)
-   [wavelets](https://cran.r-project.org/package=wavelets)
-   [waveslim](https://cran.r-project.org/package=waveslim)
-   [wavethresh](https://cran.r-project.org/package=wavethresh)
-   [wavScalogram](https://cran.r-project.org/package=wavScalogram)
-   [WeightedPortTest](https://cran.r-project.org/package=WeightedPortTest)
-   [wktmo](https://cran.r-project.org/package=wktmo)
-   [wwntests](https://cran.r-project.org/package=wwntests)
-   [x12](https://cran.r-project.org/package=x12)
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
-   [TISEAN Project](http://www.mpipks-dresden.mpg.de/~tisean/)
-   [GitHub repository for this Task
    View](https://github.com/robjhyndman/ctv-TimeSeries)
