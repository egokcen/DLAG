## DLAG (Delayed Latents Across Groups)

[![][license-img]][license-url]

[license-img]: https://img.shields.io/github/license/mashape/apistatus.svg
[license-url]: https://github.com/egokcen/DLAG/blob/master/LICENSE.md

This project contains a Matlab implementation of DLAG, a dimensionality reduction framework for dissecting interactions between populations of neurons.

The DLAG model, along with learning and inference procedures, are described in
detail in the following reference:

- “Dissecting feedforward and feedback interactions between populations of
neurons”
by E. A. Gokcen, J. D. Semedo, A. Zandvakili, C. K. Machens*, A. Kohn*,
B. M. Yu*. Cosyne Abstracts, 2020. Talk.

Please read it carefully before using the code, as it describes all of the
terminology and usage modes. Please cite the above reference if using any
portion of this code for your own purposes.

## Installation

Simply download and extract the [current release](https://github.com/egokcen/DLAG/releases/tag/v1.1.0) of this project to your desired local directory.

## Project wiki

For additional information on getting started, as well as subtler usage details, see this project's wiki [here](https://github.com/egokcen/DLAG/wiki).

## Contact
For questions, please contact Evren Gokcen at egokcen@cmu.edu. 

## Acknowledgments

This DLAG implementation started with base code located [here](https://github.com/karts25/NeuralTraj).

That project contains implementations of GPFA and TD-GPFA, which are described in
the following references:

- "Gaussian-process factor analysis for low-dimensional single-trial analysis of
neural population activity"
by B. M. Yu, J. P. Cunningham, G. Santhanam, S. I. Ryu, K. V. Shenoy,
and M. Sahani. J Neurophysiol, vol. 102, 2009, pp. 614-635.

- "Extracting Low-Dimensional Latent Structure from Time Series in the Presence
of Delays"
by K. C. Lakshmanan, P. T. Sadtler, E. C. Tyler-Kabara, A. P. Batista, B. M. Yu.
Neural Computation 2015
