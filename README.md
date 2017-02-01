CRAN Task View: High-Performance and Parallel Computing with R
--------------------------------------------------------------

|                 |                                                              
|-----------------|------------------------------------------------------------  
| **Maintainer:** | Dirk Eddelbuettel                                            
| **Contact:**    | Dirk.Eddelbuettel at R-project.org                           
| **Version:**    | 2017-02-01                                                   
| **URL:**        | <https://CRAN.R-project.org/view=HighPerformanceComputing>   

This CRAN task view contains a list of packages, grouped by topic, that are useful for high-performance computing (HPC) with R. In this context, we are defining 'high-performance computing' rather loosely as just about anything related to pushing R a little further: using compiled code, parallel computing (in both explicit and implicit modes), working with large objects as well as profiling.

Unless otherwise mentioned, all packages presented with hyperlinks are available from CRAN, the Comprehensive R Archive Network.

Several of the areas discussed in this Task View are undergoing rapid change. Please send suggestions for additions and extensions for this task view to the [task view maintainer](mailto:Dirk.Eddelbuettel@R-project.org).

Suggestions and corrections by Achim Zeileis, Markus Schmidberger, Martin Morgan, Max Kuhn, Tomas Radivoyevitch, Jochen Knaus, Tobias Verbeke, Hao Yu, David Rosenberg, Marco Enea, Ivo Welch, Jay Emerson, Wei-Chen Chen, Bill Cleveland, Ross Boylan, Ramon Diaz-Uriarte, Mark Zeligman, and Kevin Ushey (as well as others I may have forgotten to add here) are gratefully acknowledged.

Contributions are always welcome, and encouraged. Since the start of this CRAN task view in October 2008, most contributions have arrived as email suggestions. The source file for this particular task view file now also reside in a GitHub repository (see below) so that pull requests are also possible.

**Direct support in R started with release 2.14.0** which includes a new package **parallel** incorporating (slightly revised) copies of packages multicore and [snow](https://cloud.r-project.org/web/packages/snow/index.html). Some types of clusters are not handled directly by the base package 'parallel'. However, and as explained in the package vignette, the parts of parallel which provide [snow](https://cloud.r-project.org/web/packages/snow/index.html) -like functions will accept [snow](https://cloud.r-project.org/web/packages/snow/index.html) clusters including MPI clusters.
The **parallel** package also contains support for multiple RNG streams following L'Ecuyer et al (2002), with support for both mclapply and snow clusters.
The version released for R 2.14.0 contains base functionality: higher-level convenience functions are planned for later R releases.

**Parallel computing: Explicit parallelism**

-   Several packages provide the communications layer required for parallel computing. The first package in this area was rpvm by Li and Rossini which uses the PVM (Parallel Virtual Machine) standard and libraries. rpvm is no longer actively maintained, but available from its CRAN archive directory.
-   In recent years, the alternative MPI (Message Passing Interface) standard has become the de facto standard in parallel computing. It is supported in R via the [Rmpi](https://cloud.r-project.org/web/packages/Rmpi/index.html) by Yu. [Rmpi](https://cloud.r-project.org/web/packages/Rmpi/index.html) package is mature yet actively maintained and offers access to numerous functions from the MPI API, as well as a number of R-specific extensions. [Rmpi](https://cloud.r-project.org/web/packages/Rmpi/index.html) can be used with the LAM/MPI, MPICH / MPICH2, Open MPI, and Deino MPI implementations. It should be noted that LAM/MPI is now in maintenance mode, and new development is focussed on Open MPI.
-   The [pbdMPI](https://cloud.r-project.org/web/packages/pbdMPI/index.html) package provides S4 classes to directly interface MPI in order to support the Single Program/Multiple Data (SPMD) parallel programming style which is particularly useful for batch parallel execution. The [pbdSLAP](https://cloud.r-project.org/web/packages/pbdSLAP/index.html) builds on this and uses scalable linear algebra packages (namely BLACS, PBLAS, and ScaLAPACK) in double precision based on ScaLAPACK version 2.0.2. The [pbdBASE](https://cloud.r-project.org/web/packages/pbdBASE/index.html) builds on these and provides the core classes and methods for distributed data types upon which the [pbdDMAT](https://cloud.r-project.org/web/packages/pbdDMAT/index.html) builds to provide distributed dense matrices for "Programming with Big Data". The [pbdNCDF4](https://cloud.r-project.org/web/packages/pbdNCDF4/index.html) package permits multiple processes to write to the same file (without manual synchronization) and supports terabyte-sized files. The [pbdDEMO](https://cloud.r-project.org/web/packages/pbdDEMO/index.html) package provides examples for these packages, and a detailed vignette. The [pbdPROF](https://cloud.r-project.org/web/packages/pbdPROF/index.html) package profiles MPI communication SPMD code via MPI profiling libraries, such as fpmpi, mpiP, or TAU.
-   An alternative is provided by the [nws](https://cloud.r-project.org/web/packages/nws/index.html) (NetWorkSpaces) packages from REvolution Computing. It is the successor to the earlier LindaSpaces approach to parallel computing, and is implemented on top of the Twisted networking toolkit for Python.
-   The [snow](https://cloud.r-project.org/web/packages/snow/index.html) (Simple Network of Workstations) package by Tierney et al. can use PVM, MPI, NWS as well as direct networking sockets. It provides an abstraction layer by hiding the communications details. The [snowFT](https://cloud.r-project.org/web/packages/snowFT/index.html) package provides fault-tolerance extensions to [snow](https://cloud.r-project.org/web/packages/snow/index.html).
-   The [snowfall](https://cloud.r-project.org/web/packages/snowfall/index.html) package by Knaus provides a more recent alternative to [snow](https://cloud.r-project.org/web/packages/snow/index.html). Functions can be used in sequential or parallel mode.
-   The [foreach](https://cloud.r-project.org/web/packages/foreach/index.html) package allows general iteration over elements in a collection without the use of an explicit loop counter. Using foreach without side effects also facilitates executing the loop in parallel which is possible via the [doMC](https://cloud.r-project.org/web/packages/doMC/index.html) (using parallel/multicore on single workstations), [doSNOW](https://cloud.r-project.org/web/packages/doSNOW/index.html) (using [snow](https://cloud.r-project.org/web/packages/snow/index.html), see above), [doMPI](https://cloud.r-project.org/web/packages/doMPI/index.html) (using [Rmpi](https://cloud.r-project.org/web/packages/Rmpi/index.html)) packages, [doFuture](https://cloud.r-project.org/web/packages/doFuture/index.html) (using [future](https://cloud.r-project.org/web/packages/future/index.html) or [future.BatchJobs](https://cloud.r-project.org/web/packages/future.BatchJobs/index.html)), and [doRedis](https://cloud.r-project.org/web/packages/doRedis/index.html) (using [rredis](https://cloud.r-project.org/web/packages/rredis/index.html)) packages.
-   The [future](https://cloud.r-project.org/web/packages/future/index.html) package allows for synchroneous (sequential) and asynchronous (parallel) evaluations via abstraction of futures, either via function calls or implicitly via promises. Global variables are automatically identified. Iteration over elements in a collection is supported.
-   The [Rborist](https://cloud.r-project.org/web/packages/Rborist/index.html) package employs OpenMP pragmas to exploit predictor-level parallelism in the Random Forest algorithm which promotes efficient use of multicore hardware in restaging data and in determining splitting criteria, both of which are performance bottlenecks in the algorithm.
-   The [h2o](https://cloud.r-project.org/web/packages/h2o/index.html) package connects to the h2o open source machine learning environment which has scalable implementations of random forests, GBM, GLM (with elastic net regularization), and deep learning.
-   The [randomForestSRC](https://cloud.r-project.org/web/packages/randomForestSRC/index.html) package can use both OpenMP as well as MPI for random forest extensions suitable for survival analysis, competing risks analysis, classification as well as regression

**Parallel computing: Implicit parallelism**

-   The pnmath package by Tierney ( [link](http://www.stat.uiowa.edu/~luke/R/experimental/)) uses the Open MP parallel processing directives of recent compilers (such gcc 4.2 or later) for implicit parallelism by replacing a number of internal R functions with replacements that can make use of multiple cores --- without any explicit requests from the user. The alternate pnmath0 package offers the same functionality using Pthreads for environments in which the newer compilers are not available. Similar functionality is expected to become integrated into R 'eventually'.
-   The romp package by Jamitzky was presented at useR! 2008 ( [slides](http://www.statistik.tu-dortmund.de/useR-2008/slides/Jamitzky.pdf)) and offers another interface to Open MP using Fortran. The code is still pre-alpha and available from the Google Code project [<span class="Gcode">romp</span>](http://code.google.com/p/romp/). An R-Forge project [<span class="Rforge">romp</span>](https://R-Forge.R-project.org/projects/romp/) was initiated but there is no package, yet.
-   The R/parallel package by Vera, Jansen and Suppi offers a C++-based master-slave dispatch mechanism for parallel execution ( [link](http://www.rparallel.org/))
-   The [Rdsm](https://cloud.r-project.org/web/packages/Rdsm/index.html) package provides a threads-like parallel computing environment, both on multicore machine and across the network by providing facilities inspired from distributed shared memory programming.
-   The [RhpcBLASctl](https://cloud.r-project.org/web/packages/RhpcBLASctl/index.html) detects the number of available BLAS cores, and permits explicit selection of the number of cores.
-   The [Rhpc](https://cloud.r-project.org/web/packages/Rhpc/index.html) permits `*apply()` style dispatch via MPI.

**Parallel computing: Grid computing**

-   The multiR package by Grose was presented at useR! 2008 but has not been released. It may offer a snow-style framework on a grid computing platform.
-   The [<span class="Rforge">biocep-distrib</span>](https://R-Forge.R-project.org/projects/biocep-distrib/) project by Chine offers a Java-based framework for local, Grid, or Cloud computing. It is under active development.

**Parallel computing: Hadoop**

-   The RHIPE package, started by Saptarshi Guha and now developed by a core team via [GitHub](https://github.com/saptarshiguha/RHIPE/), provides an interface between R and Hadoop for analysis of large complex data wholly from within R using the Divide and Recombine approach to big data.
-   The rmr package by Revolution Analytics also provides an interface between R and Hadoop for a Map/Reduce programming framework. ( [link](https://github.com/RevolutionAnalytics/RHadoop/wiki/rmr))
-   A related package, segue package by Long, permits easy execution of embarassingly parallel task on Elastic Map Reduce (EMR) at Amazon. ( [link](http://code.google.com/p/segue/))
-   The [RProtoBuf](https://cloud.r-project.org/web/packages/RProtoBuf/index.html) package provides an interface to Google's language-neutral, platform-neutral, extensible mechanism for serializing structured data. This package can be used in R code to read data streams from other systems in a distributed MapReduce setting where data is serialized and passed back and forth between tasks.
-   The [HistogramTools](https://cloud.r-project.org/web/packages/HistogramTools/index.html) package provides a number of routines useful for the construction, aggregation, manipulation, and plotting of large numbers of Histograms such as those created by Mappers in a MapReduce application.
-   The [toaster](https://cloud.r-project.org/web/packages/toaster/index.html) package performs in-database computations utilizing the parallel / distributed Teradata Aster analytical platform

**Parallel computing: Random numbers**

-   Random-number generators for parallel computing are available via the [rlecuyer](https://cloud.r-project.org/web/packages/rlecuyer/index.html) package by Sevcikova and Rossini.
-   The [doRNG](https://cloud.r-project.org/web/packages/doRNG/index.html) package provides functions to perform reproducible parallel foreach loops, using independent random streams as generated by the package rstream, suitable for the different foreach backends.

**Parallel computing: Resource managers and batch schedulers**

-   Job-scheduling toolkits permit management of parallel computing resources and tasks. The slurm (Simple Linux Utility for Resource Management) set of programs works well with MPI and slurm jobs can be submitted from R using the [rslurm](https://cloud.r-project.org/web/packages/rslurm/index.html) package. ( [link](http://slurm.schedmd.com/))
-   The Condor toolkit ( [link](http://www.cs.wisc.edu/condor/)) from the University of Wisconsin-Madison has been used with R as described in this [R News article](http://www.r-project.org/doc/Rnews/Rnews_2005-2.pdf).
-   The sfCluster package by Knaus can be used with [snowfall](https://cloud.r-project.org/web/packages/snowfall/index.html). ( [link](http://www.imbi.uni-freiburg.de/parallel/)) but is currently limited to LAM/MPI.
-   The [batch](https://cloud.r-project.org/web/packages/batch/index.html) package by Hoffmann can launch parallel computing requests onto a cluster and gather results.
-   The [BatchJobs](https://cloud.r-project.org/web/packages/BatchJobs/index.html) package provides Map, Reduce and Filter variants to manage R jobs and their results on batch computing systems like PBS/Torque, LSF and Sun Grid Engine. Multicore and SSH systems are also supported. The [BatchExperiments](https://cloud.r-project.org/web/packages/BatchExperiments/index.html) package extends it with an abstraction layer for running statistical experiments. Package [batchtools](https://cloud.r-project.org/web/packages/batchtools/index.html) is a successor / extension to both.
-   The [flowr](https://cloud.r-project.org/web/packages/flowr/index.html) package offers a scatter-gather approach to submit jobs lists (including dependencies) to the computing cluster via simple data.frames as inputs. It supports LSF, SGE, Torque and SLURM.

**Parallel computing: Applications**

-   The [caret](https://cloud.r-project.org/web/packages/caret/index.html) package by Kuhn can use various frameworks (MPI, NWS etc) to parallelized cross-validation and bootstrap characterizations of predictive models.
-   The [<span class="BioC">maanova</span>](http://www.Bioconductor.ohttps://cloud.r-project.org/web/packages/release/bioc/html/maanova.html) package on Bioconductor by Wu can use [snow](https://cloud.r-project.org/web/packages/snow/index.html) and [Rmpi](https://cloud.r-project.org/web/packages/Rmpi/index.html) for the analysis of micro-array experiments.
-   The [pvclust](https://cloud.r-project.org/web/packages/pvclust/index.html) package by Suzuki and Shimodaira can use [snow](https://cloud.r-project.org/web/packages/snow/index.html) and [Rmpi](https://cloud.r-project.org/web/packages/Rmpi/index.html) for hierarchical clustering via multiscale bootstraps.
-   The [tm](https://cloud.r-project.org/web/packages/tm/index.html) package by Feinerer can use [snow](https://cloud.r-project.org/web/packages/snow/index.html) and [Rmpi](https://cloud.r-project.org/web/packages/Rmpi/index.html) for parallelized text mining.
-   The [varSelRF](https://cloud.r-project.org/web/packages/varSelRF/index.html) package by Diaz-Uriarte can use [snow](https://cloud.r-project.org/web/packages/snow/index.html) and [Rmpi](https://cloud.r-project.org/web/packages/Rmpi/index.html) for parallelized use of variable selection via random forests.
-   The [bcp](https://cloud.r-project.org/web/packages/bcp/index.html) package by Erdman and Emerson for the Bayesian analysis of change points can use [foreach](https://cloud.r-project.org/web/packages/foreach/index.html) for parallelized operations.
-   The [<span class="BioC">multtest</span>](http://www.Bioconductor.ohttps://cloud.r-project.org/web/packages/release/bioc/html/multtest.html) package by Pollard et al. on Bioconductor can use [snow](https://cloud.r-project.org/web/packages/snow/index.html), [Rmpi](https://cloud.r-project.org/web/packages/Rmpi/index.html) or rpvm for resampling-based testing of multiple hypothesis.
-   The [GAMBoost](https://cloud.r-project.org/web/packages/GAMBoost/index.html) package by Binder for `glm` and `gam` model fitting via boosting using b-splines, the [Geneland](https://cloud.r-project.org/web/packages/Geneland/index.html) package by Estoup, Guillot and Santos for structure detection from multilocus genetic data, the [Matching](https://cloud.r-project.org/web/packages/Matching/index.html) package by Sekhon for multivariate and propensity score matching, the [STAR](https://cloud.r-project.org/web/packages/STAR/index.html) package by Pouzat for spike train analysis, the [bnlearn](https://cloud.r-project.org/web/packages/bnlearn/index.html) package by Scutari for bayesian network structure learning, the [latentnet](https://cloud.r-project.org/web/packages/latentnet/index.html) package by Krivitsky and Handcock for latent position and cluster models, the [lga](https://cloud.r-project.org/web/packages/lga/index.html) package by Harrington for linear grouping analysis, the [peperr](https://cloud.r-project.org/web/packages/peperr/index.html) package by Porzelius and Binder for parallised estimation of prediction error, the [orloca](https://cloud.r-project.org/web/packages/orloca/index.html) package by Fernandez-Palacin and Munoz-Marquez for operations research locational analysis, the [rgenoud](https://cloud.r-project.org/web/packages/rgenoud/index.html) package by Mebane and Sekhon for genetic optimization using derivatives the [<span class="BioC">affyPara</span>](http://www.Bioconductor.ohttps://cloud.r-project.org/web/packages/release/bioc/html/affyPara.html) package by Schmidberger, Vicedo and Mansmann for parallel normalization of Affymetrix microarrays, and the [<span class="BioC">puma</span>](http://www.Bioconductor.ohttps://cloud.r-project.org/web/packages/release/bioc/html/puma.html) package by Pearson et al. which propagates uncertainty into standard microarray analyses such as differential expression all can use [snow](https://cloud.r-project.org/web/packages/snow/index.html) for parallelized operations using either one of the MPI, PVM, NWS or socket protocols supported by [snow](https://cloud.r-project.org/web/packages/snow/index.html).
-   The [<span class="Gcode">bugsparallel</span>](http://code.google.com/p/bugsparallel/) package uses [Rmpi](https://cloud.r-project.org/web/packages/Rmpi/index.html) for distributed computing of multiple MCMC chains using WinBUGS.
-   The [partDSA](https://cloud.r-project.org/web/packages/partDSA/index.html) package uses [nws](https://cloud.r-project.org/web/packages/nws/index.html) for generating a piecewise constant estimation list of increasingly complex predictors based on an intensive and comprehensive search over the entire covariate space.
-   The [dclone](https://cloud.r-project.org/web/packages/dclone/index.html) package provides a global optimization approach and a variant of simulated annealing which exploits Bayesian MCMC tools to get MLE point estimates and standard errors using low level functions for implementing maximum likelihood estimating procedures for complex models using data cloning and Bayesian Markov chain Monte Carlo methods with support for JAGS, WinBUGS and OpenBUGS; parallel computing is supported via the [snow](https://cloud.r-project.org/web/packages/snow/index.html) package.
-   The [pmclust](https://cloud.r-project.org/web/packages/pmclust/index.html) package utilizes unsupervised model-based clustering for high dimensional (ultra) large data. The package uses [pbdMPI](https://cloud.r-project.org/web/packages/pbdMPI/index.html) to perform a parallel version of the EM algorithm for finite mixture Gaussian models.
-   The [harvestr](https://cloud.r-project.org/web/packages/harvestr/index.html) package provides helper functions for (reproducible) simulations.
-   Nowadays, many packages can use the facilities offered by the **parallel** package. One example is [pls](https://cloud.r-project.org/web/packages/pls/index.html), another is [PGICA](https://cloud.r-project.org/web/packages/PGICA/index.html) which can run ICA analysis in parallel on SGE or multicore platforms.
-   The [sprint](https://cloud.r-project.org/web/packages/sprint/index.html) (an acronym for "Simple Parallel R INTerface") package provides a parallel computing framework for R making High Performance Computing (HPC) accessible to users who are not familiar with parallel programming and the use of HPC architectures. It contains a library of parallelised R functions for correlation, partitioning around medoids, apply, permutation testing, bootstrapping, random forest, rank product and hamming distance.
-   The [pbapply](https://cloud.r-project.org/web/packages/pbapply/index.html) package offers a progress bar for vectorized R functions in the \`\*apply\` family, and supports several backends.

**Parallel computing: GPUs**

-   The [gputools](https://cloud.r-project.org/web/packages/gputools/index.html) package by Buckner and Seligman provides several common data-mining algorithms which are implemented using a mixture of nVidia's CUDA langauge and cublas library. Given a computer with an nVidia GPU these functions may be substantially more efficient than native R routines.
-   The [cudaBayesreg](https://cloud.r-project.org/web/packages/cudaBayesreg/index.html) package by da Silva implements the `rhierLinearModel` from the [bayesm](https://cloud.r-project.org/web/packages/bayesm/index.html) package using nVidia's CUDA langauge and tools to provide high-performance statistical analysis of fMRI voxels.
-   The rgpu package (see below for link) aims to speed up bioinformatics analysis by using the GPU.
-   The [gcbd](https://cloud.r-project.org/web/packages/gcbd/index.html) package implements a benchmarking framework for BLAS and GPUs (using [gputools](https://cloud.r-project.org/web/packages/gputools/index.html)).
-   The [OpenCL](https://cloud.r-project.org/web/packages/OpenCL/index.html) package provides an interface from R to OpenCL permitting hardware- and vendor neutral interfaces to GPU programming.
-   The [HiPLARM](https://cloud.r-project.org/web/packages/HiPLARM/index.html) package provide High-Performance Linear Algebra for R using multi-core and/or GPU support using the PLASMA / MAGMA libraries from UTK, CUDA, and accelerated BLAS.
-   The [permGPU](https://cloud.r-project.org/web/packages/permGPU/index.html) package computes permutation resampling inference in the context of RNA microarray studies on the GPU, it uses CUDA (&gt;= 4.5)
-   The [gmatrix](https://cloud.r-project.org/web/packages/gmatrix/index.html) package enables the evaluation of matrix and vector operations using GPU coprocessors such that intermediate computations may be kept on the coprocessor and reused, with potentially significant performance enhancements by minimizing data movement.
-   The [gpuR](https://cloud.r-project.org/web/packages/gpuR/index.html) package offers GPU-enabled functions: New gpu\* and vcl\* classes are provided to wrap typical R objects (e.g. vector, matrix) mirroring typical R syntax without the need to know OpenCL.

**Large memory and out-of-memory data**

-   The [biglm](https://cloud.r-project.org/web/packages/biglm/index.html) package by Lumley uses incremental computations to offer `lm()` and `glm()` functionality to data sets stored outside of R's main memory.
-   The [ff](https://cloud.r-project.org/web/packages/ff/index.html) package by Adler et al. offers file-based access to data sets that are too large to be loaded into memory, along with a number of higher-level functions.
-   The [bigmemory](https://cloud.r-project.org/web/packages/bigmemory/index.html) package by Kane and Emerson permits storing large objects such as matrices in memory (as well as via files) and uses external pointer objects to refer to them. This permits transparent access from R without bumping against R's internal memory limits. Several R processes on the same computer can also share big memory objects.
-   A large number of database packages, and database-alike packages (such as [sqldf](https://cloud.r-project.org/web/packages/sqldf/index.html) by Grothendieck and [data.table](https://cloud.r-project.org/web/packages/data.table/index.html) by Dowle) are also of potential interest but not reviewed here.
-   The [HadoopStreaming](https://cloud.r-project.org/web/packages/HadoopStreaming/index.html) package provides a framework for writing map/reduce scripts for use in Hadoop Streaming; it also facilitates operating on data in a streaming fashion which does not require Hadoop.
-   The [speedglm](https://cloud.r-project.org/web/packages/speedglm/index.html) package permits to fit (generalised) linear models to large data. For in-memory data sets, speedlm() or speedglm() can be used along with update.speedlm() which can update fitted models with new data. For out-of-memory data sets, shglm() is available; it works in the presence of factors and can check for singular matrices.
-   The [biglars](https://cloud.r-project.org/web/packages/biglars/index.html) package by Seligman et al can use the [ff](https://cloud.r-project.org/web/packages/ff/index.html) to support large-than-memory datasets for least-angle regression, lasso and stepwise regression.
-   The [MonetDB.R](https://cloud.r-project.org/web/packages/MonetDB.R/index.html) package allows R to access the MonetDB column-oriented, open source database system as a backend.
-   The [ffbase](https://cloud.r-project.org/web/packages/ffbase/index.html) package by de Jonge et al adds basic statistical functionality to the [ff](https://cloud.r-project.org/web/packages/ff/index.html) package.
-   The [LaF](https://cloud.r-project.org/web/packages/LaF/index.html) package provides methods for fast access to large ASCII files in csv or fixed-width format.

**Easier interfaces for Compiled code**

-   The [inline](https://cloud.r-project.org/web/packages/inline/index.html) package by Sklyar et al eases adding code in C, C++ or Fortran to R. It takes care of the compilation, linking and loading of embeded code segments that are stored as R strings.
-   The [Rcpp](https://cloud.r-project.org/web/packages/Rcpp/index.html) package by Eddelbuettel and Francois offers a number of C++ clases that makes transferring R objects to C++ functions (and back) easier, and the [RInside](https://cloud.r-project.org/web/packages/RInside/index.html) package by the same authors allows easy embedding of R itself into C++ applications for faster and more direct data transfer.
-   The [RcppParallel](https://cloud.r-project.org/web/packages/RcppParallel/index.html) package by Allaire et al. bundles the [Intel Threading Building Blocks](https://www.threadingbuildingblocks.org) and [TinyThread](http://tinythreadpp.bitsnbites.eu) libraries. Together with [Rcpp](https://cloud.r-project.org/web/packages/Rcpp/index.html), RcppParallel makes it easy to write safe, performant, concurrently-executing C++ code, and use that code within R and R packages.
-   The [rJava](https://cloud.r-project.org/web/packages/rJava/index.html) package by Urbanek provides a low-level interface to Java similar to the `.Call()` interface for C and C++.

**Profiling tools**

-   The [profr](https://cloud.r-project.org/web/packages/profr/index.html) package by Wickham can visualize output from the `Rprof` interface for profiling.
-   The [proftools](https://cloud.r-project.org/web/packages/proftools/index.html) package by Tierney, and the [aprof](https://cloud.r-project.org/web/packages/aprof/index.html) package by Visser, can also be used to analyse profiling output.
-   The [GUIProfiler](https://cloud.r-project.org/web/packages/GUIProfiler/index.html) package visualizes the results of profiling R programs.

### CRAN packages:

-   [aprof](https://cloud.r-project.org/web/packages/aprof/index.html)
-   [batch](https://cloud.r-project.org/web/packages/batch/index.html)
-   [BatchExperiments](https://cloud.r-project.org/web/packages/BatchExperiments/index.html)
-   [BatchJobs](https://cloud.r-project.org/web/packages/BatchJobs/index.html)
-   [batchtools](https://cloud.r-project.org/web/packages/batchtools/index.html)
-   [bayesm](https://cloud.r-project.org/web/packages/bayesm/index.html)
-   [bcp](https://cloud.r-project.org/web/packages/bcp/index.html)
-   [biglars](https://cloud.r-project.org/web/packages/biglars/index.html)
-   [biglm](https://cloud.r-project.org/web/packages/biglm/index.html)
-   [bigmemory](https://cloud.r-project.org/web/packages/bigmemory/index.html)
-   [bnlearn](https://cloud.r-project.org/web/packages/bnlearn/index.html)
-   [caret](https://cloud.r-project.org/web/packages/caret/index.html)
-   [cudaBayesreg](https://cloud.r-project.org/web/packages/cudaBayesreg/index.html)
-   [data.table](https://cloud.r-project.org/web/packages/data.table/index.html)
-   [dclone](https://cloud.r-project.org/web/packages/dclone/index.html)
-   [doFuture](https://cloud.r-project.org/web/packages/doFuture/index.html)
-   [doMC](https://cloud.r-project.org/web/packages/doMC/index.html)
-   [doMPI](https://cloud.r-project.org/web/packages/doMPI/index.html)
-   [doRedis](https://cloud.r-project.org/web/packages/doRedis/index.html)
-   [doRNG](https://cloud.r-project.org/web/packages/doRNG/index.html)
-   [doSNOW](https://cloud.r-project.org/web/packages/doSNOW/index.html)
-   [ff](https://cloud.r-project.org/web/packages/ff/index.html)
-   [ffbase](https://cloud.r-project.org/web/packages/ffbase/index.html)
-   [flowr](https://cloud.r-project.org/web/packages/flowr/index.html)
-   [foreach](https://cloud.r-project.org/web/packages/foreach/index.html)
-   [future](https://cloud.r-project.org/web/packages/future/index.html)
-   [future.BatchJobs](https://cloud.r-project.org/web/packages/future.BatchJobs/index.html)
-   [GAMBoost](https://cloud.r-project.org/web/packages/GAMBoost/index.html)
-   [gcbd](https://cloud.r-project.org/web/packages/gcbd/index.html)
-   [Geneland](https://cloud.r-project.org/web/packages/Geneland/index.html)
-   [gmatrix](https://cloud.r-project.org/web/packages/gmatrix/index.html)
-   [gpuR](https://cloud.r-project.org/web/packages/gpuR/index.html)
-   [gputools](https://cloud.r-project.org/web/packages/gputools/index.html)
-   [GUIProfiler](https://cloud.r-project.org/web/packages/GUIProfiler/index.html)
-   [h2o](https://cloud.r-project.org/web/packages/h2o/index.html)
-   [HadoopStreaming](https://cloud.r-project.org/web/packages/HadoopStreaming/index.html)
-   [harvestr](https://cloud.r-project.org/web/packages/harvestr/index.html)
-   [HiPLARM](https://cloud.r-project.org/web/packages/HiPLARM/index.html)
-   [HistogramTools](https://cloud.r-project.org/web/packages/HistogramTools/index.html)
-   [inline](https://cloud.r-project.org/web/packages/inline/index.html)
-   [LaF](https://cloud.r-project.org/web/packages/LaF/index.html)
-   [latentnet](https://cloud.r-project.org/web/packages/latentnet/index.html)
-   [lga](https://cloud.r-project.org/web/packages/lga/index.html)
-   [Matching](https://cloud.r-project.org/web/packages/Matching/index.html)
-   [MonetDB.R](https://cloud.r-project.org/web/packages/MonetDB.R/index.html)
-   [nws](https://cloud.r-project.org/web/packages/nws/index.html)
-   [OpenCL](https://cloud.r-project.org/web/packages/OpenCL/index.html)
-   [orloca](https://cloud.r-project.org/web/packages/orloca/index.html)
-   [partDSA](https://cloud.r-project.org/web/packages/partDSA/index.html)
-   [pbapply](https://cloud.r-project.org/web/packages/pbapply/index.html)
-   [pbdBASE](https://cloud.r-project.org/web/packages/pbdBASE/index.html)
-   [pbdDEMO](https://cloud.r-project.org/web/packages/pbdDEMO/index.html)
-   [pbdDMAT](https://cloud.r-project.org/web/packages/pbdDMAT/index.html)
-   [pbdMPI](https://cloud.r-project.org/web/packages/pbdMPI/index.html)
-   [pbdNCDF4](https://cloud.r-project.org/web/packages/pbdNCDF4/index.html)
-   [pbdPROF](https://cloud.r-project.org/web/packages/pbdPROF/index.html)
-   [pbdSLAP](https://cloud.r-project.org/web/packages/pbdSLAP/index.html)
-   [peperr](https://cloud.r-project.org/web/packages/peperr/index.html)
-   [permGPU](https://cloud.r-project.org/web/packages/permGPU/index.html)
-   [PGICA](https://cloud.r-project.org/web/packages/PGICA/index.html)
-   [pls](https://cloud.r-project.org/web/packages/pls/index.html)
-   [pmclust](https://cloud.r-project.org/web/packages/pmclust/index.html)
-   [profr](https://cloud.r-project.org/web/packages/profr/index.html)
-   [proftools](https://cloud.r-project.org/web/packages/proftools/index.html)
-   [pvclust](https://cloud.r-project.org/web/packages/pvclust/index.html)
-   [randomForestSRC](https://cloud.r-project.org/web/packages/randomForestSRC/index.html)
-   [Rborist](https://cloud.r-project.org/web/packages/Rborist/index.html)
-   [Rcpp](https://cloud.r-project.org/web/packages/Rcpp/index.html)
-   [RcppParallel](https://cloud.r-project.org/web/packages/RcppParallel/index.html)
-   [Rdsm](https://cloud.r-project.org/web/packages/Rdsm/index.html)
-   [rgenoud](https://cloud.r-project.org/web/packages/rgenoud/index.html)
-   [Rhpc](https://cloud.r-project.org/web/packages/Rhpc/index.html)
-   [RhpcBLASctl](https://cloud.r-project.org/web/packages/RhpcBLASctl/index.html)
-   [RInside](https://cloud.r-project.org/web/packages/RInside/index.html)
-   [rJava](https://cloud.r-project.org/web/packages/rJava/index.html)
-   [rlecuyer](https://cloud.r-project.org/web/packages/rlecuyer/index.html)
-   [Rmpi](https://cloud.r-project.org/web/packages/Rmpi/index.html) (core)
-   [RProtoBuf](https://cloud.r-project.org/web/packages/RProtoBuf/index.html)
-   [rredis](https://cloud.r-project.org/web/packages/rredis/index.html)
-   [rslurm](https://cloud.r-project.org/web/packages/rslurm/index.html)
-   [snow](https://cloud.r-project.org/web/packages/snow/index.html) (core)
-   [snowfall](https://cloud.r-project.org/web/packages/snowfall/index.html)
-   [snowFT](https://cloud.r-project.org/web/packages/snowFT/index.html)
-   [speedglm](https://cloud.r-project.org/web/packages/speedglm/index.html)
-   [sprint](https://cloud.r-project.org/web/packages/sprint/index.html)
-   [sqldf](https://cloud.r-project.org/web/packages/sqldf/index.html)
-   [STAR](https://cloud.r-project.org/web/packages/STAR/index.html)
-   [tm](https://cloud.r-project.org/web/packages/tm/index.html)
-   [toaster](https://cloud.r-project.org/web/packages/toaster/index.html)
-   [varSelRF](https://cloud.r-project.org/web/packages/varSelRF/index.html)

### Related links:

-   [HPC computing notes by Luke Tierney for HPC class at University of Iowa](http://www.stat.uiowa.edu/~luke/classes/295-hpc/)
-   [Mailing List: R Special Interest Group High Performance Computing](https://stat.ethz.ch/mailman/listinfo/r-sig-hpc/)
-   [Schmidberger, Morgan, Eddelbuettel, Yu, Tierney and Mansmann (2009) paper on 'State of the Art in Parallel Computing with R'](http://www.jstatsoft.org/v31/i01/)
-   [Luke Tierney's code directory for pnmath and pnmath0](http://www.stat.uiowa.edu/~luke/R/experimental/)
-   R-Forge Project: [<span class="Rforge">biocep-distrib</span>](https://R-Forge.R-project.org/projects/biocep-distrib/)
-   Bioconductor Package: [<span class="BioC">affyPara</span>](http://www.Bioconductor.ohttps://cloud.r-project.org/web/packages/release/bioc/html/affyPara.html)
-   Bioconductor Package: [<span class="BioC">maanova</span>](http://www.Bioconductor.ohttps://cloud.r-project.org/web/packages/release/bioc/html/maanova.html)
-   Bioconductor Package: [<span class="BioC">multtest</span>](http://www.Bioconductor.ohttps://cloud.r-project.org/web/packages/release/bioc/html/multtest.html)
-   Bioconductor Package: [<span class="BioC">puma</span>](http://www.Bioconductor.ohttps://cloud.r-project.org/web/packages/release/bioc/html/puma.html)
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
-   [GitHub repository for this Task View](https://github.com/eddelbuettel/ctv-hpc)
