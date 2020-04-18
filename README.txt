Matlab code for extracting latent variables from multi-population data using
DLAG (Delayed Latents Across Groups).

Version 1.02        18 Apr 2020

The DLAG model, along with learning and inference procedures, are described in
detail in the following reference:

“Dissecting feedforward and feedback interactions between populations of
neurons”
by E. A. Gokcen, J. D. Semedo, A. Zandvakili, C. K. Machens*, A. Kohn*,
B. M. Yu*. Cosyne Abstracts, 2020. Talk.

Please read it carefully before using the code, as it describes all of the
terminology and usage modes. Please cite the above reference if using any
portion of this code for your own purposes.

This DLAG implementation started with base code (which has been significantly
changed) from the following source:
https://github.com/karts25/NeuralTraj
Code in this repo implements GPFA and TD-GPFA, which are described in
the following references:

"Gaussian-process factor analysis for low-dimensional single-trial analysis of
neural population activity"
by B. M. Yu, J. P. Cunningham, G. Santhanam, S. I. Ryu, K. V. Shenoy,
and M. Sahani. J Neurophysiol, vol. 102, 2009, pp. 614-635.

"Extracting Low-Dimensional Latent Structure from Time Series in the Presence
of Delays"
by K. C. Lakshmanan, P. T. Sadtler, E. C. Tyler-Kabara, A. P. Batista, B. M. Yu.
Neural Computation 2015

===============
VERSION HISTORY
===============
In notes below, + denotes new feature, - denotes bug fix or removal.

Version 1.01     -- 13 Apr 2020
Version 1.02     -- 18 Apr 2020
    + Added 0-within-group dimension functionality.

====================
HOW TO GET STARTED
====================

Look at example_dlag.m and example_pcca.m
These scripts include comprehensive demonstrations of exploratory data analysis
and model selection / cross-validation.

Note: Matlab should be started in the current directory (DLAG) so
that it executes startup.m automatically. Otherwise execute startup.m before
attempting to run this code.

==================================================
ARE THERE KNOBS I NEED TO TURN WHEN FITTING DLAG?
==================================================

For those who have used the Yu Lab's implementation of GPFA and/or TD-GPFA,
many of the knobs described below should be familiar. These common knobs are
originally described in the repo noted above, and we reproduce their
descriptions here, for completeness. There are a few additional knobs unique to
DLAG, also described below.

Knobs in the preprocessing of data:

1) Data format (continuous-valued or spiking data) ('datFormat')

   The user has the option of inputting continuous-valued data:

   result = fit_dlag(runIdx, dat, 'datFormat', 'seq');

   or spikes in 1ms time bins:

   result = fit_dlag(runIdx, dat, 'datFormat', 'spikes');

   Continuous-valued data can be data of any modality that the user has already
   preprocessed. Spiking data will be further preprocessed and converted to
   binned spike counts.

2) Segment length ('segLength')

   Like GPFA, the time required for fitting a DLAG model scales roughly with
   T^3, where T is the number of timesteps in the longest trial.
   When different trials have the same number of timesteps, expensive
   computations can be reused. This provides motivation for having a
   small number of unique trial lengths in the dataset, as well as
   having those unique trial lengths be small.

   To speed up DLAG model fitting (by roughly 1.5 orders of
   magnitude), we first segment the trials into many segments of the
   same length ('segLength'). The DLAG model parameters are fit
   using those segments.  We then take those fitted parameters and
   the original (non-segmented) trials to extract full-length neural
   trajectories. All of this is done automatically in the code pack,
   but one knob the user can adjust is 'segLength', which can be
   passed as an optional argument into fit_dlag:

   result = fit_dlag(runIdx, dat, 'segLength', 50);

   Note that the optional argument (in this case, 50) is the number of
   timesteps, not units of time (e.g., the number of ms). If the optional
   argument is not passed in, fit_dlag.m will use a default of 20
   timesteps, which we found to be a reasonable value for various
   datasets. Segmenting can be turned off by setting 'segLength' to Inf,
   in which case the original (non-segmented) trials will be used for model
   fitting.

   The choice of 'segLength' should be guided by the following two
   considerations:

   a) Capturing longer timescales and/or delays

   In general, choosing a shorter 'segLength' gives faster DLAG model
   fitting.  However, 'segLength' should not be made so small that
   DLAG cannot find the longer timescales and/or delays, whose adverse
   effect can be quantified by evaluating cross-validated performance.
   Thus, one should typically choose as large of a 'segLength' as
   possible, as long as the run time is acceptable.  From our
   experience, the performance penalty due to segmenting is small and
   the run time gains are huge.

   b) Minimizing overlap

   The segmentation is performed with overlap, where the overlap is
   randomly distributed within a trial. There is always one segment
   that starts at the first timestep and another segment that ends at
   the last timestep.  Take, for example, a trial with 55 timesteps
   and 'segLength'=50.  There will be two segments, each 50 timesteps
   long, but 45 timesteps of overlap.  The overlap creates extra work
   for the fitting algorithm (i.e., it slows down the algorithm) and
   over-represents the central portion of trials (which may mean
   higher cross-validated prediction error). Thus, one should
   typically choose a 'segLength' that minimizes the amount of
   overlap in view of the distribution of original trial lengths.

3) Sample period or spike bin width ('binWidth')

   result = fit_dlag(runIdx, dat, 'binWidth', 20);

   For continuous-valued data, the 'binWidth' argument specifies the sample
   period at which the input data samples were acquired, or the bin width used
   to count spikes. 'binWidth' can be in any unit of time, but be consistent for
   all temporal arguments. The precision at which timescales and delays can be
   resolved will depend on the inherent temporal resolution of the recording
   technique, as well as the sample period at which data was collected, or bin
   width used to count spikes (see below for more on spike count bin width).

   For spiking data, spikes are counted in non-overlapping bins. The bin width
   (in ms) is thus also the time between adjacent timesteps.

   There are two advantages to using a larger bin width:

   - A larger bin width means larger spike counts. The larger the
     counts are, the better the square-root transform stabilizes the
     variance of Poisson-distributed counts (see Kihlberg et al.,
     1972). Recall that all methods here assume stationary observation noise.
     Using a larger bin width is particularly important when firing rates are
     low across the neural population.

   - A larger bin width also means fewer timesteps in each trial. This
     makes the code run faster. Furthermore, a smaller 'segLength'
     (measured in number of timesteps) can be used to cover the same
     amount of absolute time. Running time scales roughly with
     segLength^3, subject to issues of overlap discussed above.

     The drawback of a larger bin width is less temporal resolution.  We
     recommend using the largest bin width possible that still gives
     reasonable temporal resolution. We typically use 'binWidth'
     values of 20, 50, or 100 msec. The default 'binWidth' is 20ms.

Knobs in initialization of the fitting procedure:

1) Initialization method (‘init_method’)
   The EM algorithm for fitting DLAG can be initialized in two ways:
   ‘pCCA’    - 	parameters are initialized using probabilistic canonical
                correlation analysis (pCCA), described in (Bach & Jordan, 2005).
                This is the default option.
   ‘params’  - 	parameters are initialized to a particular realization of DLAG
                model parameters. For example, one may wish to continue training
                a DLAG model for which training was cut short, or performed
                in a different context. In this case, specify the additional
                optional argument 'initParams':
                    fit_dlag(runIdx, dat, 'init_method', 'params', ...
                             'initParams', params)

2) Initial Gaussian process parameters
   If 'init_method' is set to 'pCCA', then the following arguments are relevant:
   'startTau'   - All within- and across-group Gaussian process timescales will
                  be initialized to 'startTau', given in the same units of time
                  as 'binWidth'. The default is 2*binWidth.
   'startEps'   - All within- and across-group Gaussian process noise variances
                  will be initialized to 'startEps'. The default is 1e-3.
   'startDelay' - 'startDelay' is an array with length corresponding to the
                  number of groups (areas) in the data. All delays to a
                  group (given in the same units of time as 'binWidth'), will be
                  initialized to the corresponding entry of 'startDelay'. By
                  default, all delays are initialized at 0.

Knobs that constrain the outputs of the fitting procedure:

1) Learn or fix delays ('learnDelays')
   DLAG automatically learns delays from data. However, if desired, delays
   can be fixed at their initial values throughout the fitting procedure.
   For example, fixing delays at 0 can provide run time benefits during
   cross-validation, when many models will inevitably find spurious across-group
   signals. In these cases, it can take many EM iterations for delays to
   converge, with little benefit. By default, 'learnDelays' is true.

2) Maximum delay (‘maxDelayFrac’)
   This allows the user to bound the largest values the estimated delays can
   take,as described in the Methods section of the reference above.
   ‘maxDelayFrac’ is set by default to 0.5, which means the magnitude of the
   largest allowed delay is 0.5 times the length of the shortest trial. This
   option can be turned off by setting 'maxDelayFrac' to Inf, in which case
   delays will be optimized in an unconstrained manner.

Knobs that determine what/where data gets saved to disk or overwritten

1) Base directory where the directory 'mat_results' will be created ('baseDir')
   By default, results are saved in the user's current directory
   ('./mat_results'). The user can change this location by setting the optional
   'baseDir' argument. Then, results will be saved to 'baseDir/mat_results'.

2) Save train and test data, inferred trajectories ('saveData')
   By default, train and test data, along with trajectories inferred during
   fitting and performance evaluation, will not be saved to disk. By setting
   'saveData' to true, the user can dictate that these data are saved to disk.
   Doing so could be useful for reproducibility. But beware -- one can quickly
   run out of disk space, depending on the dataset.

3) Overwrite or skip existing results files ('overwriteExisting')
   By default, existing results files (that conflict with models to be fitted)
   will be overwritten. To skip over existing results files, set
   'overwriteExisting' to false.

================================================================
WHAT DOES IT MEAN FOR A WITHIN-GROUP DIMENSIONALITY TO BE ZERO?
================================================================
fit_dlag.m gives the option of fitting DLAG models for which one or more groups
has a within-group dimensionality of 0. For example,

    fit_dlag(runIdx, dat, 'xDims_within', {2, 0});

will fit a DLAG model in which group 1 has 2 within-group latents, and group 2
has 0 within-group latents. This model assumes that there exists variance
shared among observations in group 1 that is not shared among observations in
group 2. In contrast, this model assumes that all variance shared among
observations in group 2 arises through its interaction with group 1. All
parameters associated with group 2, such as gamma_within{2} or eps_within{2},
will be empty. The loading matrix, C, will not have any columns associated with
within-group latents for group 2 (since they don't exist).

===========================
NOTES ON CROSS-VALIDATION
===========================

************
THE BASICS
************
DLAG has 1+numGroups hyperparameters, where numGroups is the number of groups
(areas) in the data: xDim_across, the number of across-group latent states; and
xDim_within = [xDim_within(1) ... xDim_within(numGroups)], the number of
within-group latent states for each group.

To fit multiple models, each with different across-group dimensionalities,
set the 'xDims_across' optional argument in fit_dlag.m to a list of desired
values. For example, to fit models with across-group dimensionalities 1, 2,
and 3, run

    fit_dlag(runIdx, dat, 'xDims_across', 1:3);

Or to fit models with across-group dimensionalities 2, 5, and 6, run

    fit_dlag(runIdx, dat, 'xDims_across', [2 5 6]);

To fit multiple models, each with different within-group dimensionalities,
set the 'xDims_within' optional argument in fit_dlag.m. 'xDims_within' can
be specified in two formats:

  - First, suppose we wish to fit a DLAG model with
    two groups. For group 1, we want models with within-group dimensionalities
    1, 2, and 3. And for group 2, we want models with within-group
    dimensionalities 3, 4, and 5. Then we could run

        fit_dlag(runIdx, dat, 'xDims_within', {1:3, 3:5});

    In this case, 'xDims_within' is a (1 x numGroups) cell array, where each
    element is a list of desired within-group dimensionalities for the
    corresponding group. 3x3=9 sets of models will be fitted, with every
    possible combination of within-group dimensionalities.

  - Second, suppose we want to sweep over within-group dimensionalities, but we
    want every group to have the same number of within-group dimensionalities.
    Then we could run

        fit_dlag(runIdx, dat, 'xDims_within', {1:3});

    In this case, 'xDims_within' contains only one element. Then, we will fit
    three sets of models -- i.e., models where all within-group dimensionalities
    are 1, where all within-group dimensionalities are 2, and where all
    within-group dimensionalities are 3.

Putting these examples on 'xDims_across' and 'xDims_within' together, running

    fit_dlag(runIdx, dat, 'xDims_across', 1:3, 'xDims_within', {1:3, 1:3});

will fit 3x3x3=27 models, with all possible combinations of within- and
across-group dimensionalities. On the other hand, running

    fit_dlag(runIdx, dat, 'xDims_across', 1:3, 'xDims_within', {1:3});

will fit 3x3=9 models, where all models are constrained to have the same number
of within-group dimensions for all groups.

************************
UNDERDETERMINED MODELS
************************
fit_dlag.m will NOT train models for which the combined within- and across-group
state dimensionalities is greater than the number of observed features
(neurons). For example, suppose a dataset consists of two groups, each with
10 observed features. If the user runs

    fit_dlag(runIdx, dat, 'xDims_across', 1:6, 'xDims_within', {1:5});

then models for which xDim_across = 6 and xDim_within(i) = 5 will be ignored.
The combined state dimensionality would be 11, which is greater than 10. Such
models are underdetermined and ill-defined. All other models will be fit
normally.

*******************************************
STREAMLINING THE CROSS-VALIDATION PROCESS
*******************************************
Brute force grid search over optimal within- and across-group state
dimensionalities can quickly get out of hand, in terms of scale. For example,

    fit_dlag(runIdx, dat, 'xDims_across', 1:5, 'xDims_within', {1:5, 1:5});

will attempt to fit 5x5x5=125 DLAG models, for every possible combination of
within- and across-group dimensionalities. Thankfully, the search over within-
and across-group dimensionalities can be largely decoupled (see the referenced
paper for an explanation why). And one can narrow the range of potential
optimal models by several approaches. Here are some suggestions (not
necessarily in order):

    1) Apply pCCA first, and cross-validate to find the optimal pCCA latent
       dimensionality (see example_pcca.m and/or example_dlag.m). pCCA is
       relatively fast, and in practice, the estimates for across-group
       dimensionality are similar for pCCA and DLAG. So, for example, if the
       optimal pCCA model has 5 latent dimensions, one could limit the DLAG
       search range to models with, say, 3-8 across-group dimensions.

       Applying pCCA first (or any related static method) is also important
       to establish that there is any significant interaction between groups
       in the first place. If pCCA cross-validation returns 0 latent dimensions,
       then there is no significant interaction, and there is no need to spend
       time fitting insignificant DLAG models.

    2) Apply factor analysis (FA) or GPFA to each group separately. These
       single-area methods will give a rough upper bound on the combined
       within- and across-group dimensionality for each group.

    3) Fit a DLAG model with more within- and across-group dimensions than you
       think might be needed (based on the above suggestions). Then,
       orthonormalize trajectories according to different objectives. This
       codepack provides a number of functions to help with this process (and
       see example_dlag.m):

       - predictiveProjection_dlag.m orthonormalizes and orders across-group
         latent trajectories according to across-group predictive power. These
         "predictive" latents correspond to the inferences behind R2orth
         or MSEorth in pairwise_regress_dlag.m described above.
         Visualization can provide clues as to the number of significant
         across-group latents.

       - orthonormalizeWithinGroups.m orthonormalizes and orders latents
         according to shared variance explained within each area.
         Visualization can provide clues as to the number of significant
         within- and/or across-group latents.

       - findSharedDimCutoff.m finds (for a given DLAG model):
           - the minimum number of across-group dimensions needed to explain
            a certain percentage (e.g., 95%) of across-group shared variance.
           - the minimum number of within-group dimensions needed to explain
             a certain percentage (e.g., 95%) of within-group shared variance,
             for each group individually.
         These values are analogous to 'd_shared', as described in
         “Scaling properties of dimensionality reduction for neural populations
         and network models”
         by R. C. Williamson, et al., PLoS Comput Biol, 2016.

       - pairwise_regress_dlag.m predicts the activity in one group given the
         activity in another group. The metrics R2orth or MSEorth (these metrics
         are also computed during cross-validation) can give an estimate for the
         number across-group dimensions, regardless of the number of
         within-group dimensions included in the model.

    4) Fix DLAG within-group dimensionalities to values larger than might be
       needed (based on the above suggestions). Then sweep over the range of
       across-group dimensionalities based on suggestion 1. So long as
       *at least* as many within-group dimensions are used as needed to capture
       within-group shared variance over time, across-group latents will be
       correctly identified. Hence the optimal across-group dimensionality can
       be identified, largely decoupled from the search for the optimal
       within-group dimensionalities. If within-group states are not interesting
       for your application, then you can even stop at this point.

    5) Once the optimal across-group dimensionality is determined, one can
       fix the across-group dimensionality at this optimal value, and then
       sweep over within-group dimensionalities independently (around
       smaller ranges narrowed by suggestions 2 and/or 3). Hold the other
       within-group dimensionalities fixed (at either larger-than-needed values
       or their optimal values).

Throughout cross-validation, many candidate models will inevitably find
spurious across-group signals. In these cases, it can take many EM iterations
for delays to converge, with little benefit. Thus, fixing delays at 0 (see
the explanation of 'learnDelays' above) can provide run time benefits during
cross-validation. After an optimal model is identified in this manner, sweep
over a small range of across-group dimensionalities around the optimal point,
but allow delays to be learned. To provide further speed up, you can initialize
these models with the parameters of their zero-delay counterparts (see
explanation of 'init_method' above).

(Aside: Assuming delays are truly present, zero-delay models might require
more across-group latent variables to explain the data than their delayed
counterparts. Hence the need to sweep over a small range around the optimal
zero-delay model.)

=====================================
HOW CAN I MAKE THE CODE RUN FASTER?
=====================================

DLAG uses many of the same runtime optimizations as the original GPFA code, and
hence runs about as fast as it can. Runtime is comparable to GPFA, but
inevitably slower due to additional structure and parameters in the DLAG model.
On datasets of the size encountered in the reference noted above, we have found
the time for DLAG model fitting to be on the order of minutes on a standard CPU.
In the absence of parallelization, cross-validation is the ultimate speed
bottleneck. If further speed-up is desired, one can consider the following
options:

1) Adjust 'segLength' and 'binWidth', as described above.

2) Use a smaller state dimensionality, either across or within areas.

   Tradeoff: Might not be able to capture all the structure in data.

   How to do it: Set 'xDim_across' or 'xDim_within', optional arguments to
   fit_dlag.m.

   Run time scales roughly with xDim^3 (either within or across).

3) Use fewer EM iterations (dangerous)

   Tradeoff: EM might not have converged yet.

   How to do it: Set 'maxIters', an optional argument to fit_dlag.m.

   Run time scales linearly with 'emMaxIters'.

4) Use a larger tolerance for delay/timescale convergence.

   Tradeoff: EM might not have converged yet.

   How to do it: Set 'tolLL' or 'tolParam', optional arguments to fit_dlag.m.

   DLAG can determine EM convergence based on two criteria: (1) LL improves by
   less than a set tolerance after a certain number of iterations; and
   (2) Across-area delays and timescales change by less than a set tolerance
   after a certain number of iterations.

5) Adjust how frequently data log-likelihood is evaluated.

   Tradeoff: The EM algorithm might go a few more iterations than it needed to,
   given 'tolLL'. But this is a very minor concern.

   How to do it: By setting 'freqLL', one can specify that data log-likelihood
   is only computed every 'freqLL' EM iterations. By default, 'freqLL' is 10.


6) For cross-validation, use the parallelization option.

   How to do it: Set ‘parallelize’, an optional argument to fit_dlag.m.
   The number of workers can be specified with the 'numWorkers' option.
   Tailor this value to your machine's specifications. Don't use more workers
   than the number of physical cores in your machine.

   Each model and each cross-validation fold will be trained simultaneously
   in an independent thread. If more models need to be fitted than there are
   workers, then models will be trained in batches of 'numWorkers'.

  =============================
  NOTE ON THE USE OF C/MEX CODE
  =============================

   The function invToeplitz(), which inverts Toeplitz matrices such
   as the GP RBF kernels, attains major runtime speedups by exploiting
   special properties of Toeplitz matrices in a C environment better
   suited to such computations (for loops).

   The function makePrecomp(), which calculates posterior covariance
   matrices within learnGPparams.m, attains major runtime speedups by
   offloading parallelizable/pipelineable computations to a C environment
   better suited to such computations (for loops).

   These improvements require the C/MEX MATLAB interface, whereby C
   code can be compiled and used from within MATLAB.  We have compiled
   this code on a number of platforms in the hopes that this will be
   seamless to the user, but some compilation may be required on the
   part of the user if he/she is using an unusual architecture (such as
   Solaris).

   If the proper mex files are not available, this compilation should
   be automatically done by the startup.m file, so please see the notes
   there that call the "mex" command.

   To understand more, please read "help invToeplitz" in
   util/invToeplitz.
   Please also read "help makePrecomp" in util/precomp/makePrecomp.
   Other useful files to read include invToeplitzFast.m,
   invToeplitzFastZohar.c, and makePautoSumFast.c.

   If you have not used MEX previously, try
   "help mex" and some of the online tutorials such as:
   http://www.mathworks.com/support/tech-notes/1600/1605.html .

   Importantly, all of this code is written to default to a native MATLAB
   code implementation if the mex file can not be properly executed.
   Though the user may lose significant runtime performance, the code
   should still work the same and produce identical answers.

================================================
WHAT DOES THE WARNING ABOUT PRIVATE NOISE MEAN?
================================================

The private noise variance (or uniqueness, in Factor Analysis [FA] speak) for
one or more units may be driven to zero.

There are three possible causes:

1) Crosstalk between electrodes.  This can lead to different units
   having the same (or nearly the same) activity on each trial.

   Solution: Remove units affected by crosstalk from 'dat'.

2) The state dimensionality (xDim_across or xDim_within) is too large.
   The extra dimensions in the latent space may be dedicated to explaining
   particular units perfectly, thus giving zero private noise for
   those units.

   Solution: Reduce 'xDim_across' or 'xDim_within', optional arguments to
   fit_dlag.m.

3) You have encountered a Heywood case.  It's an issue with
   maximum-likelihood parameter learning, whereby more likelihood
   can be gained by setting a private noise to 0 than by finding a
   compromise.  This is a corner case of FA that has been known
   since the 1930's.  Various Bayesian FA models have been proposed,
   but here we simply set a minimum private noise variance for each
   unit as a percentage of its raw data variance.  This problem can
   also arise in GPFA and DLAG, but is less common because of the GP
   smoothing.

   Two possible solutions:
   1) Do nothing.  The private variance is automatically capped at
      some minimum non-zero value.
   2) Remove the offending unit from 'dat'.

   For more about Heywood cases, see:

   "Bayesian Estimation in Unrestricted Factor Analysis: A Treatment
   for Heywood Cases"
   J. K. Martin and R. P. McDonald.
   Psychometrika, 40, 4, 505-17, Dec 1975.
