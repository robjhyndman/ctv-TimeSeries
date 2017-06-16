CRAN Task View: High-Performance and Parallel Computing with R
--------------------------------------------------------------

|                 |                                                              
|-----------------|------------------------------------------------------------  
| **Maintainer:** | Dirk Eddelbuettel                                            
| **Contact:**    | Dirk.Eddelbuettel at R-project.org                           
| **Version:**    | 2017-04-25                                                   
| **URL:**        | <https://CRAN.R-project.org/view=HighPerformanceComputing>   

This CRAN task view contains a list of packages, grouped by topic, that are useful for high-performance computing (HPC) with R. In this context, we are defining 'high-performance computing' rather loosely as just about anything related to pushing R a little further: using compiled code, parallel computing (in both explicit and implicit modes), working with large objects as well as profiling.

Unless otherwise mentioned, all packages presented with hyperlinks are available from CRAN, the Comprehensive R Archive Network.

Several of the areas discussed in this Task View are undergoing rapid change. Please send suggestions for additions and extensions for this task view to the [task view maintainer](mailto:Dirk.Eddelbuettel@R-project.org).

Suggestions and corrections by Achim Zeileis, Markus Schmidberger, Martin Morgan, Max Kuhn, Tomas Radivoyevitch, Jochen Knaus, Tobias Verbeke, Hao Yu, David Rosenberg, Marco Enea, Ivo Welch, Jay Emerson, Wei-Chen Chen, Bill Cleveland, Ross Boylan, Ramon Diaz-Uriarte, Mark Zeligman, and Kevin Ushey (as well as others I may have forgotten to add here) are gratefully acknowledged.

Contributions are always welcome, and encouraged. Since the start of this CRAN task view in October 2008, most contributions have arrived as email suggestions. The source file for this particular task view file now also reside in a GitHub repository (see below) so that pull requests are also possible.

The `ctv` package supports these Task Views. Its functions `install.views` and `update.views` allow, respectively, installation or update of packages from a given Task View; the option `coreOnly` can restrict operations to packages labeled as *core* below.

**Direct support in R started with release 2.14.0** which includes a new package **parallel** incorporating (slightly revised) copies of packages multicore and [snow](https://cran.r-project.org/package=snow/index.html). Some types of clusters are not handled directly by the base package 'parallel'. However, and as explained in the package vignette, the parts of parallel which provide [snow](../packages/snow/index.html) -like functions will accept [snow](../packages/snow) clusters including MPI clusters.
The **parallel** package also contains support for multiple RNG streams following L'Ecuyer et al (2002), with support for both mclapply and snow clusters.
The version released for R 2.14.0 contains base functionality: higher-level convenience functions are planned for later R releases.

**Parallel computing: Explicit parallelism**

-   Several packages provide the communications layer required for parallel computing. The first package in this area was rpvm by Li and Rossini which uses the PVM (Parallel Virtual Machine) standard and libraries. rpvm is no longer actively maintained, but available from its CRAN archive directory.
-   In recent years, the alternative MPI (Message Passing Interface) standard has become the de facto standard in parallel computing. It is supported in R via the [Rmpi](https://cran.r-project.org/package=Rmpi/index.html) by Yu. [Rmpi](../packages/Rmpi/index.html) package is mature yet actively maintained and offers access to numerous functions from the MPI API, as well as a number of R-specific extensions. [Rmpi](../packages/Rmpi) can be used with the LAM/MPI, MPICH / MPICH2, Open MPI, and Deino MPI implementations. It should be noted that LAM/MPI is now in maintenance mode, and new development is focussed on Open MPI.
-   The [pbdMPI](https://cran.r-project.org/package=pbdMPI/index.html) package provides S4 classes to directly interface MPI in order to support the Single Program/Multiple Data (SPMD) parallel programming style which is particularly useful for batch parallel execution. The [pbdSLAP](../packages/pbdSLAP/index.html) builds on this and uses scalable linear algebra packages (namely BLACS, PBLAS, and ScaLAPACK) in double precision based on ScaLAPACK version 2.0.2. The [pbdBASE](../packages/pbdBASE/index.html) builds on these and provides the core classes and methods for distributed data types upon which the [pbdDMAT](../packages/pbdDMAT/index.html) builds to provide distributed dense matrices for "Programming with Big Data". The [pbdNCDF4](../packages/pbdNCDF4/index.html) package permits multiple processes to write to the same file (without manual synchronization) and supports terabyte-sized files. The [pbdDEMO](../packages/pbdDEMO/index.html) package provides examples for these packages, and a detailed vignette. The [pbdPROF](../packages/pbdPROF) package profiles MPI communication SPMD code via MPI profiling libraries, such as fpmpi, mpiP, or TAU.
-   An alternative is provided by the [nws](https://cran.r-project.org/package=nws) (NetWorkSpaces) packages from REvolution Computing. It is the successor to the earlier LindaSpaces approach to parallel computing, and is implemented on top of the Twisted networking toolkit for Python.
-   The [snow](https://cran.r-project.org/package=snow/index.html) (Simple Network of Workstations) package by Tierney et al. can use PVM, MPI, NWS as well as direct networking sockets. It provides an abstraction layer by hiding the communications details. The [snowFT](../packages/snowFT/index.html) package provides fault-tolerance extensions to [snow](../packages/snow).
-   The [snowfall](https://cran.r-project.org/package=snowfall/index.html) package by Knaus provides a more recent alternative to [snow](../packages/snow). Functions can be used in sequential or parallel mode.
-   The [foreach](https://cran.r-project.org/package=foreach/index.html) package allows general iteration over elements in a collection without the use of an explicit loop counter. Using foreach without side effects also facilitates executing the loop in parallel which is possible via the [doMC](../packages/doMC/index.html) (using parallel/multicore on single workstations), [doSNOW](../packages/doSNOW/index.html) (using [snow](../packages/snow/index.html), see above), [doMPI](../packages/doMPI/index.html) (using [Rmpi](../packages/Rmpi/index.html)) packages, [doFuture](../packages/doFuture/index.html) (using [future](../packages/future/index.html) or [future.BatchJobs](../packages/future.BatchJobs/index.html)), and [doRedis](../packages/doRedis/index.html) (using [rredis](../packages/rredis)) packages.
-   The [future](https://cran.r-project.org/package=future) package allows for synchroneous (sequential) and asynchronous (parallel) evaluations via abstraction of futures, either via function calls or implicitly via promises. Global variables are automatically identified. Iteration over elements in a collection is supported.
-   The [Rborist](https://cran.r-project.org/package=Rborist) package employs OpenMP pragmas to exploit predictor-level parallelism in the Random Forest algorithm which promotes efficient use of multicore hardware in restaging data and in determining splitting criteria, both of which are performance bottlenecks in the algorithm.
-   The [h2o](https://cran.r-project.org/package=h2o) package connects to the h2o open source machine learning environment which has scalable implementations of random forests, GBM, GLM (with elastic net regularization), and deep learning.
-   The [randomForestSRC](https://cran.r-project.org/package=randomForestSRC) package can use both OpenMP as well as MPI for random forest extensions suitable for survival analysis, competing risks analysis, classification as well as regression

**Parallel computing: Implicit parallelism**

-   The pnmath package by Tierney ( [link](http://www.stat.uiowa.edu/~luke/R/experimental/)) uses the Open MP parallel processing directives of recent compilers (such gcc 4.2 or later) for implicit parallelism by replacing a number of internal R functions with replacements that can make use of multiple cores --- without any explicit requests from the user. The alternate pnmath0 package offers the same functionality using Pthreads for environments in which the newer compilers are not available. Similar functionality is expected to become integrated into R 'eventually'.
-   The romp package by Jamitzky was presented at useR! 2008 ( [slides](http://www.statistik.tu-dortmund.de/useR-2008/slides/Jamitzky.pdf)) and offers another interface to Open MP using Fortran. The code is still pre-alpha and available from the Google Code project [<span class="Gcode">romp</span>](http://code.google.com/p/romp/). An R-Forge project [<span class="Rforge">romp</span>](https://R-Forge.R-project.org/projects/romp/) was initiated but there is no package, yet.
-   The R/parallel package by Vera, Jansen and Suppi offers a C++-based master-slave dispatch mechanism for parallel execution ( [link](http://www.rparallel.org/))
-   The [Rdsm](https://cran.r-project.org/package=Rdsm) package provides a threads-like parallel computing environment, both on multicore machine and across the network by providing facilities inspired from distributed shared memory programming.
-   The [RhpcBLASctl](https://cran.r-project.org/package=RhpcBLASctl) detects the number of available BLAS cores, and permits explicit selection of the number of cores.
-   The [Rhpc](https://cran.r-project.org/package=Rhpc) permits `*apply()` style dispatch via MPI.
- The [drake](https://cran.r-project.org/package=drake) package is an R-focused reproducible build system. Similarly to [Make](https://www.gnu.org/software/make), it arranges the intermediate steps of a workflow and executes them in parallelizable stages. With three alternative backends (powered by `parallel::mclapply()`, `parallel::parLapply()`, and [Makefiles](https://www.gnu.org/software/make), respectivly) [drake](https://cran.r-project.org/package=drake) supports implicit parallelism for both low-overhead single-node applications and true distributed computing workloads.

**Parallel computing: Grid computing**

-   The multiR package by Grose was presented at useR! 2008 but has not been released. It may offer a snow-style framework on a grid computing platform.
-   The [<span class="Rforge">biocep-distrib</span>](https://R-Forge.R-project.org/projects/biocep-distrib/) project by Chine offers a Java-based framework for local, Grid, or Cloud computing. It is under active development.

**Parallel computing: Hadoop**

-   The RHIPE package, started by Saptarshi Guha and now developed by a core team via [GitHub](https://github.com/saptarshiguha/RHIPE/), provides an interface between R and Hadoop for analysis of large complex data wholly from within R using the Divide and Recombine approach to big data.
-   The rmr package by Revolution Analytics also provides an interface between R and Hadoop for a Map/Reduce programming framework. ( [link](https://github.com/RevolutionAnalytics/RHadoop/wiki/rmr))
-   A related package, segue package by Long, permits easy execution of embarassingly parallel task on Elastic Map Reduce (EMR) at Amazon. ( [link](http://code.google.com/p/segue/))
-   The [RProtoBuf](https://cran.r-project.org/package=RProtoBuf) package provides an interface to Google's language-neutral, platform-neutral, extensible mechanism for serializing structured data. This package can be used in R code to read data streams from other systems in a distributed MapReduce setting where data is serialized and passed back and forth between tasks.
-   The [HistogramTools](https://cran.r-project.org/package=HistogramTools) package provides a number of routines useful for the construction, aggregation, manipulation, and plotting of large numbers of Histograms such as those created by Mappers in a MapReduce application.
-   The [toaster](https://cran.r-project.org/package=toaster) package performs in-database computations utilizing the parallel / distributed Teradata Aster analytical platform

**Parallel computing: Random numbers**

-   Random-number generators for parallel computing are available via the [rlecuyer](https://cran.r-project.org/package=rlecuyer) package by Sevcikova and Rossini.
-   The [doRNG](https://cran.r-project.org/package=doRNG) package provides functions to perform reproducible parallel foreach loops, using independent random streams as generated by the package rstream, suitable for the different foreach backends.

**Parallel computing: Resource managers and batch schedulers**

-   Job-scheduling toolkits permit management of parallel computing resources and tasks. The slurm (Simple Linux Utility for Resource Management) set of programs works well with MPI and slurm jobs can be submitted from R using the [rslurm](https://cran.r-project.org/package=rslurm) package. ( [link](http://slurm.schedmd.com/))
-   The Condor toolkit ( [link](http://www.cs.wisc.edu/condor/)) from the University of Wisconsin-Madison has been used with R as described in this [R News article](http://www.r-project.org/doc/Rnews/Rnews_2005-2.pdf).
-   The sfCluster package by Knaus can be used with [snowfall](https://cran.r-project.org/package=snowfall). ( [link](http://www.imbi.uni-freiburg.de/parallel/)) but is currently limited to LAM/MPI.
-   The [batch](https://cran.r-project.org/package=batch) package by Hoffmann can launch parallel computing requests onto a cluster and gather results.
-   The [BatchJobs](https://cran.r-project.org/package=BatchJobs/index.html) package provides Map, Reduce and Filter variants to manage R jobs and their results on batch computing systems like PBS/Torque, LSF and Sun Grid Engine. Multicore and SSH systems are also supported. The [BatchExperiments](../packages/BatchExperiments/index.html) package extends it with an abstraction layer for running statistical experiments. Package [batchtools](../packages/batchtools) is a successor / extension to both.
-   The [flowr](https://cran.r-project.org/package=flowr) package offers a scatter-gather approach to submit jobs lists (including dependencies) to the computing cluster via simple data.frames as inputs. It supports LSF, SGE, Torque and SLURM.

**Parallel computing: Applications**

-   The [caret](https://cran.r-project.org/package=caret) package by Kuhn can use various frameworks (MPI, NWS etc) to parallelized cross-validation and bootstrap characterizations of predictive models.
-   The [<span class="BioC">maanova</span>](http://www.Bioconductor.ohttps://cran.r-project.org/package=release/bioc/html/maanova.html) package on Bioconductor by Wu can use [snow](../packages/snow/index.html) and [Rmpi](../packages/Rmpi) for the analysis of micro-array experiments.
-   The [pvclust](https://cran.r-project.org/package=pvclust/index.html) package by Suzuki and Shimodaira can use [snow](../packages/snow/index.html) and [Rmpi](../packages/Rmpi) for hierarchical clustering via multiscale bootstraps.
-   The [tm](https://cran.r-project.org/package=tm/index.html) package by Feinerer can use [snow](../packages/snow/index.html) and [Rmpi](../packages/Rmpi) for parallelized text mining.
-   The [varSelRF](https://cran.r-project.org/package=varSelRF/index.html) package by Diaz-Uriarte can use [snow](../packages/snow/index.html) and [Rmpi](../packages/Rmpi) for parallelized use of variable selection via random forests.
-   The [bcp](https://cran.r-project.org/package=bcp/index.html) package by Erdman and Emerson for the Bayesian analysis of change points can use [foreach](../packages/foreach) for parallelized operations.
-   The [<span class="BioC">multtest</span>](http://www.Bioconductor.ohttps://cran.r-project.org/package=release/bioc/html/multtest.html) package by Pollard et al. on Bioconductor can use [snow](../packages/snow/index.html), [Rmpi](../packages/Rmpi) or rpvm for resampling-based testing of multiple hypothesis.
-   The [GAMBoost](https://cran.r-project.org/package=GAMBoost/index.html) package by Binder for `glm` and `gam` model fitting via boosting using b-splines, the [Matching](../packages/Matching/index.html) package by Sekhon for multivariate and propensity score matching, the [STAR](../packages/STAR/index.html) package by Pouzat for spike train analysis, the [bnlearn](../packages/bnlearn/index.html) package by Scutari for bayesian network structure learning, the [latentnet](../packages/latentnet/index.html) package by Krivitsky and Handcock for latent position and cluster models, the [lga](../packages/lga/index.html) package by Harrington for linear grouping analysis, the [peperr](../packages/peperr/index.html) package by Porzelius and Binder for parallised estimation of prediction error, the [orloca](../packages/orloca/index.html) package by Fernandez-Palacin and Munoz-Marquez for operations research locational analysis, the [rgenoud](../packages/rgenoud/index.html) package by Mebane and Sekhon for genetic optimization using derivatives the [<span class="BioC">affyPara</span>](http://www.Bioconductor.org/packages/release/bioc/html/affyPara.html) package by Schmidberger, Vicedo and Mansmann for parallel normalization of Affymetrix microarrays, and the [<span class="BioC">puma</span>](http://www.Bioconductor.org/packages/release/bioc/html/puma.html) package by Pearson et al. which propagates uncertainty into standard microarray analyses such as differential expression all can use [snow](../packages/snow/index.html) for parallelized operations using either one of the MPI, PVM, NWS or socket protocols supported by [snow](../packages/snow).
-   The [<span class="Gcode">bugsparallel</span>](http://code.google.com/p/bugsparallel/) package uses [Rmpi](https://cran.r-project.org/package=Rmpi) for distributed computing of multiple MCMC chains using WinBUGS.
-   The [partDSA](https://cran.r-project.org/package=partDSA/index.html) package uses [nws](../packages/nws) for generating a piecewise constant estimation list of increasingly complex predictors based on an intensive and comprehensive search over the entire covariate space.
-   The [dclone](https://cran.r-project.org/package=dclone/index.html) package provides a global optimization approach and a variant of simulated annealing which exploits Bayesian MCMC tools to get MLE point estimates and standard errors using low level functions for implementing maximum likelihood estimating procedures for complex models using data cloning and Bayesian Markov chain Monte Carlo methods with support for JAGS, WinBUGS and OpenBUGS; parallel computing is supported via the [snow](../packages/snow) package.
-   The [pmclust](https://cran.r-project.org/package=pmclust/index.html) package utilizes unsupervised model-based clustering for high dimensional (ultra) large data. The package uses [pbdMPI](../packages/pbdMPI) to perform a parallel version of the EM algorithm for finite mixture Gaussian models.
-   The [harvestr](https://cran.r-project.org/package=harvestr) package provides helper functions for (reproducible) simulations.
-   Nowadays, many packages can use the facilities offered by the **parallel** package. One example is [pls](https://cran.r-project.org/package=pls/index.html), another is [PGICA](../packages/PGICA) which can run ICA analysis in parallel on SGE or multicore platforms.
-   The [sprint](https://cran.r-project.org/package=sprint) (an acronym for "Simple Parallel R INTerface") package provides a parallel computing framework for R making High Performance Computing (HPC) accessible to users who are not familiar with parallel programming and the use of HPC architectures. It contains a library of parallelised R functions for correlation, partitioning around medoids, apply, permutation testing, bootstrapping, random forest, rank product and hamming distance.
-   The [pbapply](https://cran.r-project.org/package=pbapply) package offers a progress bar for vectorized R functions in the \`\*apply\` family, and supports several backends.

**Parallel computing: GPUs**

-   The [gputools](https://cran.r-project.org/package=gputools) package by Buckner and Seligman provides several common data-mining algorithms which are implemented using a mixture of nVidia's CUDA langauge and cublas library. Given a computer with an nVidia GPU these functions may be substantially more efficient than native R routines.
-   The [cudaBayesreg](https://cran.r-project.org/package=cudaBayesreg/index.html) package by da Silva implements the `rhierLinearModel` from the [bayesm](../packages/bayesm) package using nVidia's CUDA langauge and tools to provide high-performance statistical analysis of fMRI voxels.
-   The rgpu package (see below for link) aims to speed up bioinformatics analysis by using the GPU.
-   The [gcbd](https://cran.r-project.org/package=gcbd/index.html) package implements a benchmarking framework for BLAS and GPUs (using [gputools](../packages/gputools)).
-   The [OpenCL](https://cran.r-project.org/package=OpenCL) package provides an interface from R to OpenCL permitting hardware- and vendor neutral interfaces to GPU programming.
-   The [HiPLARM](https://cran.r-project.org/package=HiPLARM) package provide High-Performance Linear Algebra for R using multi-core and/or GPU support using the PLASMA / MAGMA libraries from UTK, CUDA, and accelerated BLAS.
-   The [permGPU](https://cran.r-project.org/package=permGPU) package computes permutation resampling inference in the context of RNA microarray studies on the GPU, it uses CUDA (&gt;= 4.5)
-   The [gmatrix](https://cran.r-project.org/package=gmatrix) package enables the evaluation of matrix and vector operations using GPU coprocessors such that intermediate computations may be kept on the coprocessor and reused, with potentially significant performance enhancements by minimizing data movement.
-   The [gpuR](https://cran.r-project.org/package=gpuR) package offers GPU-enabled functions: New gpu\* and vcl\* classes are provided to wrap typical R objects (e.g. vector, matrix) mirroring typical R syntax without the need to know OpenCL.

**Large memory and out-of-memory data**

-   The [biglm](https://cran.r-project.org/package=biglm) package by Lumley uses incremental computations to offer `lm()` and `glm()` functionality to data sets stored outside of R's main memory.
-   The [ff](https://cran.r-project.org/package=ff) package by Adler et al. offers file-based access to data sets that are too large to be loaded into memory, along with a number of higher-level functions.
-   The [bigmemory](https://cran.r-project.org/package=bigmemory) package by Kane and Emerson permits storing large objects such as matrices in memory (as well as via files) and uses external pointer objects to refer to them. This permits transparent access from R without bumping against R's internal memory limits. Several R processes on the same computer can also share big memory objects.
-   A large number of database packages, and database-alike packages (such as [sqldf](https://cran.r-project.org/package=sqldf/index.html) by Grothendieck and [data.table](../packages/data.table) by Dowle) are also of potential interest but not reviewed here.
-   The [HadoopStreaming](https://cran.r-project.org/package=HadoopStreaming) package provides a framework for writing map/reduce scripts for use in Hadoop Streaming; it also facilitates operating on data in a streaming fashion which does not require Hadoop.
-   The [speedglm](https://cran.r-project.org/package=speedglm) package permits to fit (generalised) linear models to large data. For in-memory data sets, speedlm() or speedglm() can be used along with update.speedlm() which can update fitted models with new data. For out-of-memory data sets, shglm() is available; it works in the presence of factors and can check for singular matrices.
-   The [biglars](https://cran.r-project.org/package=biglars/index.html) package by Seligman et al can use the [ff](../packages/ff) to support large-than-memory datasets for least-angle regression, lasso and stepwise regression.
-   The [MonetDB.R](https://cran.r-project.org/package=MonetDB.R) package allows R to access the MonetDB column-oriented, open source database system as a backend.
-   The [ffbase](https://cran.r-project.org/package=ffbase/index.html) package by de Jonge et al adds basic statistical functionality to the [ff](../packages/ff) package.
-   The [LaF](https://cran.r-project.org/package=LaF) package provides methods for fast access to large ASCII files in csv or fixed-width format.

**Easier interfaces for Compiled code**

-   The [inline](https://cran.r-project.org/package=inline) package by Sklyar et al eases adding code in C, C++ or Fortran to R. It takes care of the compilation, linking and loading of embeded code segments that are stored as R strings.
-   The [Rcpp](https://cran.r-project.org/package=Rcpp/index.html) package by Eddelbuettel and Francois offers a number of C++ clases that makes transferring R objects to C++ functions (and back) easier, and the [RInside](../packages/RInside) package by the same authors allows easy embedding of R itself into C++ applications for faster and more direct data transfer.
-   The [RcppParallel](https://cran.r-project.org/package=RcppParallel/index.html) package by Allaire et al. bundles the [Intel Threading Building Blocks](https://www.threadingbuildingblocks.org) and [TinyThread](http://tinythreadpp.bitsnbites.eu) libraries. Together with [Rcpp](../packages/Rcpp), RcppParallel makes it easy to write safe, performant, concurrently-executing C++ code, and use that code within R and R packages.
-   The [rJava](https://cran.r-project.org/package=rJava) package by Urbanek provides a low-level interface to Java similar to the `.Call()` interface for C and C++.

**Profiling tools**

-   The [profr](https://cran.r-project.org/package=profr) package by Wickham can visualize output from the `Rprof` interface for profiling.
-   The [proftools](https://cran.r-project.org/package=proftools/index.html) package by Tierney, and the [aprof](../packages/aprof) package by Visser, can also be used to analyse profiling output.
-   The [GUIProfiler](https://cran.r-project.org/package=GUIProfiler) package visualizes the results of profiling R programs.

### CRAN packages:

-   [aprof](https://cran.r-project.org/package=aprof)
-   [batch](https://cran.r-project.org/package=batch)
-   [BatchExperiments](https://cran.r-project.org/package=BatchExperiments)
-   [BatchJobs](https://cran.r-project.org/package=BatchJobs)
-   [batchtools](https://cran.r-project.org/package=batchtools)
-   [bayesm](https://cran.r-project.org/package=bayesm)
-   [bcp](https://cran.r-project.org/package=bcp)
-   [biglars](https://cran.r-project.org/package=biglars)
-   [biglm](https://cran.r-project.org/package=biglm)
-   [bigmemory](https://cran.r-project.org/package=bigmemory)
-   [bnlearn](https://cran.r-project.org/package=bnlearn)
-   [caret](https://cran.r-project.org/package=caret)
-   [cudaBayesreg](https://cran.r-project.org/package=cudaBayesreg)
-   [data.table](https://cran.r-project.org/package=data.table)
-   [dclone](https://cran.r-project.org/package=dclone)
-   [doFuture](https://cran.r-project.org/package=doFuture)
-   [doMC](https://cran.r-project.org/package=doMC)
-   [doMPI](https://cran.r-project.org/package=doMPI)
-   [doRedis](https://cran.r-project.org/package=doRedis)
-   [doRNG](https://cran.r-project.org/package=doRNG)
-   [doSNOW](https://cran.r-project.org/package=doSNOW)
-   [drake](https://cran.r-project.org/package=drake)
-   [ff](https://cran.r-project.org/package=ff)
-   [ffbase](https://cran.r-project.org/package=ffbase)
-   [flowr](https://cran.r-project.org/package=flowr)
-   [foreach](https://cran.r-project.org/package=foreach)
-   [future](https://cran.r-project.org/package=future)
-   [future.BatchJobs](https://cran.r-project.org/package=future.BatchJobs)
-   [GAMBoost](https://cran.r-project.org/package=GAMBoost)
-   [gcbd](https://cran.r-project.org/package=gcbd)
-   [gmatrix](https://cran.r-project.org/package=gmatrix)
-   [gpuR](https://cran.r-project.org/package=gpuR)
-   [gputools](https://cran.r-project.org/package=gputools)
-   [GUIProfiler](https://cran.r-project.org/package=GUIProfiler)
-   [h2o](https://cran.r-project.org/package=h2o)
-   [HadoopStreaming](https://cran.r-project.org/package=HadoopStreaming)
-   [harvestr](https://cran.r-project.org/package=harvestr)
-   [HiPLARM](https://cran.r-project.org/package=HiPLARM)
-   [HistogramTools](https://cran.r-project.org/package=HistogramTools)
-   [inline](https://cran.r-project.org/package=inline)
-   [LaF](https://cran.r-project.org/package=LaF)
-   [latentnet](https://cran.r-project.org/package=latentnet)
-   [lga](https://cran.r-project.org/package=lga)
-   [Matching](https://cran.r-project.org/package=Matching)
-   [MonetDB.R](https://cran.r-project.org/package=MonetDB.R)
-   [nws](https://cran.r-project.org/package=nws)
-   [OpenCL](https://cran.r-project.org/package=OpenCL)
-   [orloca](https://cran.r-project.org/package=orloca)
-   [partDSA](https://cran.r-project.org/package=partDSA)
-   [pbapply](https://cran.r-project.org/package=pbapply)
-   [pbdBASE](https://cran.r-project.org/package=pbdBASE)
-   [pbdDEMO](https://cran.r-project.org/package=pbdDEMO)
-   [pbdDMAT](https://cran.r-project.org/package=pbdDMAT)
-   [pbdMPI](https://cran.r-project.org/package=pbdMPI)
-   [pbdNCDF4](https://cran.r-project.org/package=pbdNCDF4)
-   [pbdPROF](https://cran.r-project.org/package=pbdPROF)
-   [pbdSLAP](https://cran.r-project.org/package=pbdSLAP)
-   [peperr](https://cran.r-project.org/package=peperr)
-   [permGPU](https://cran.r-project.org/package=permGPU)
-   [PGICA](https://cran.r-project.org/package=PGICA)
-   [pls](https://cran.r-project.org/package=pls)
-   [pmclust](https://cran.r-project.org/package=pmclust)
-   [profr](https://cran.r-project.org/package=profr)
-   [proftools](https://cran.r-project.org/package=proftools)
-   [pvclust](https://cran.r-project.org/package=pvclust)
-   [randomForestSRC](https://cran.r-project.org/package=randomForestSRC)
-   [Rborist](https://cran.r-project.org/package=Rborist)
-   [Rcpp](https://cran.r-project.org/package=Rcpp)
-   [RcppParallel](https://cran.r-project.org/package=RcppParallel)
-   [Rdsm](https://cran.r-project.org/package=Rdsm)
-   [rgenoud](https://cran.r-project.org/package=rgenoud)
-   [Rhpc](https://cran.r-project.org/package=Rhpc)
-   [RhpcBLASctl](https://cran.r-project.org/package=RhpcBLASctl)
-   [RInside](https://cran.r-project.org/package=RInside)
-   [rJava](https://cran.r-project.org/package=rJava)
-   [rlecuyer](https://cran.r-project.org/package=rlecuyer)
-   [Rmpi](https://cran.r-project.org/package=Rmpi) (core)
-   [RProtoBuf](https://cran.r-project.org/package=RProtoBuf)
-   [rredis](https://cran.r-project.org/package=rredis)
-   [rslurm](https://cran.r-project.org/package=rslurm)
-   [snow](https://cran.r-project.org/package=snow) (core)
-   [snowfall](https://cran.r-project.org/package=snowfall)
-   [snowFT](https://cran.r-project.org/package=snowFT)
-   [speedglm](https://cran.r-project.org/package=speedglm)
-   [sprint](https://cran.r-project.org/package=sprint)
-   [sqldf](https://cran.r-project.org/package=sqldf)
-   [STAR](https://cran.r-project.org/package=STAR)
-   [tm](https://cran.r-project.org/package=tm)
-   [toaster](https://cran.r-project.org/package=toaster)
-   [varSelRF](https://cran.r-project.org/package=varSelRF)

### Related links:

-   [HPC computing notes by Luke Tierney for HPC class at University of Iowa](http://www.stat.uiowa.edu/~luke/classes/295-hpc/)
-   [Mailing List: R Special Interest Group High Performance Computing](https://stat.ethz.ch/mailman/listinfo/r-sig-hpc/)
-   [Schmidberger, Morgan, Eddelbuettel, Yu, Tierney and Mansmann (2009) paper on 'State of the Art in Parallel Computing with R'](http://www.jstatsoft.org/v31/i01/)
-   [Luke Tierney's code directory for pnmath and pnmath0](http://www.stat.uiowa.edu/~luke/R/experimental/)
-   R-Forge Project: [<span class="Rforge">biocep-distrib</span>](https://R-Forge.R-project.org/projects/biocep-distrib/)
-   Bioconductor Package: [<span class="BioC">affyPara</span>](http://www.Bioconductor.org/packages/release/bioc/html/affyPara.html)
-   Bioconductor Package: [<span class="BioC">maanova</span>](http://www.Bioconductor.org/packages/release/bioc/html/maanova.html)
-   Bioconductor Package: [<span class="BioC">multtest</span>](http://www.Bioconductor.org/packages/release/bioc/html/multtest.html)
-   Bioconductor Package: [<span class="BioC">puma</span>](http://www.Bioconductor.org/packages/release/bioc/html/puma.html)
-   Google Code Project: [<span class="Gcode">romp</span>](http://code.google.com/p/romp/)
-   Google Code Project: [<span class="Gcode">bugsparallel</span>](http://code.google.com/p/bugsparallel/)
-   [Slurm open-source workload manager](http://slurm.schedmd.com/)
-   [Condor project at University of Wisconsin-Madison](http://www.cs.wisc.edu/condor/)
-   [Parallel Computing in R with sfCluster/snowfall](http://www.imbi.uni-freiburg.de/parallel/)
-   [Wikipedia: Message Passing Interface (MPI)](http://en.wikipedia.org/wiki/Message_Passing_Interface)
-   [Wikipedia: Parallel Virtual Machine (PVM)](http://en.wikipedia.org/wiki/Parallel_Virtual_Machine)
-   [Slides from Introduction to High-Performance Computing with R tutorial help in Nov 2009 at the Institute for Statistical Mathematics, Tokyo, Japan](http://dirk.eddelbuettel.com/papers/ismNov2009introHPCwithR.pdf)
-   [rgpu project at nbic.nl](https://gforge.nbic.nl/projects/rgpu/)
-   [Magma: Matrix Algebra on GPU and Multicore architectures](http://icl.cs.utk.edu/magma/)
-   [Parallel R: Data Analysis in the Distributed World"](http://shop.oreilly.com/product/0636920021421.do)
-   [High Performance Statistical Computing for Data Intensive Research](http://thirteen-01.stat.iastate.edu/snoweye/hpsc/)
-   [Rth: Parallel R through Thrust](http://heather.cs.ucdavis.edu/~matloff/rth.html)
-   [Programming with Big Data in R](http://r-pbd.org)
-   [RHIPE](https://github.com/saptarshiguha/RHIPE/)
-   [Beyond Single Core: Parallel Analysis in R](https://github.com/ljdursi/beyond-single-core-R)
-   [GitHub repository for this Task View](https://github.com/eddelbuettel/ctv-hpc)
