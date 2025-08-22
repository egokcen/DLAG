## DLAG (Delayed Latents Across Groups)

[![][license-img]][license-url]

[license-img]: https://img.shields.io/github/license/mashape/apistatus.svg
[license-url]: https://github.com/egokcen/DLAG/blob/master/LICENSE.md

This project contains a Matlab implementation of DLAG, a dimensionality reduction framework for disentangling the flow of signals between populations of neurons.

### Table of Contents

- [Citing this work](#citing-this-work)
- [System requirements](#system-requirements)
- [Installation guide](#installation-guide)
- [Getting started](#getting-started)
- [Project wiki](#project-wiki)
- [Contact](#contact)
- [License](#license)
- [Acknowledgments](#acknowledgments)

## Citing this work

The DLAG model, along with learning and inference procedures, are described in
detail in the following [paper](https://doi.org/10.1038/s43588-022-00282-5)
(access the
[main text](https://users.ece.cmu.edu/~byronyu/papers/GokcenNatCompSci2022.pdf) and
[supplement](https://users.ece.cmu.edu/~byronyu/papers/GokcenNatCompSci2022_supp.pdf)):

- Gokcen, E., Jasper, A. I., Semedo, J. D., Zandvakili, A., Kohn, A., Machens, C. K.
  & Yu, B. M. Disentangling the flow of signals between populations of neurons. _Nature Computational Science_ 2, 512–525 (2022).

For proper attribution, please cite this reference and [this codepack](CITATION.cff) if
you use any portion of this code in your own work.

This project now also includes a frequency domain approach to dramatically accelerate
the fitting of DLAG models, based on the following [paper](https://doi.org/10.1162/neco.a.22):

- Gokcen, E., Jasper, A. I., Kohn, A., Machens, C. K., & Yu, B. M.
  Fast multigroup Gaussian process factor models. _Neural Computation_, 37, 1709-1782 (2025).

Please cite this reference if you use the frequency domain functionality of this project
in your own work.

## System requirements

This codepack was written in Matlab (The MathWorks, Inc.), and must be run in 
Matlab.

Some functions rely on the C/MEX Matlab interface for speedup, whereby C code 
can be compiled and used from within Matlab. A native C compiler is necessary
to take advantage of this functionality. Windows users may require extra 
installation of, for example, Microsoft Visual C++ or MinGW. The code is written
to default to a native Matlab code implementation if mex files cannot be 
properly executed, and will still work correctly.

## Installation guide

Simply download and extract the [latest release](https://github.com/egokcen/DLAG/releases) of this project to your desired local directory.

Install [Matlab](https://www.mathworks.com/products/matlab.html).

You may need to specifically install the 
[Matlab Bioinformatics Toolbox](https://www.mathworks.com/help/bioinfo/index.html)
before getting started, if it's not already installed in your Matlab build.

For C/MEX Compilation, simply run `startup.m`.

Assuming Matlab is already installed on your machine, setup should not take 
more than a few minutes.

## Getting started

Execute `startup.m` to add all necessary dependencies to the Matlab path.
Then to get familiar with the methods in this codepack and their usage, check out the
[`demo`](demo) directory. There, demo scripts provide comprehensive demonstrations
of exploratory data analysis, model selection / cross validation, and post-selection
inference.

## Project wiki

For additional information on getting started, as well as subtler usage details, see
this project's [wiki](https://github.com/egokcen/DLAG/wiki).

## Contact
For questions, please contact Evren Gokcen at egokcen@cmu.edu.

## License
[MIT](LICENSE.md)

## Acknowledgments

This DLAG implementation started with base code located
[here](https://github.com/karts25/NeuralTraj).

That project contains implementations of GPFA and TD-GPFA, which are described in
the following references:

- Yu, B. M., Cunningham, J. P., Santhanam, G., Ryu, S. I., Shenoy, K. V. & Sahani, M. 
  Gaussian-Process Factor Analysis for Low-Dimensional Single-Trial Analysis of Neural
  Population Activity. _Journal of Neurophysiology_ 102, 614–635 (2009).

- Lakshmanan, K. C., Sadtler, P. T., Tyler-Kabara, E. C., Batista, A. P. & Yu, B. M.
  Extracting Low-Dimensional Latent Structure from Time Series in the Presence of
  Delays. _Neural Computation_ 27, 1825–1856 (2015).
