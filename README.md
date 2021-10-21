## DLAG (Delayed Latents Across Groups)

[![][license-img]][license-url]

[license-img]: https://img.shields.io/github/license/mashape/apistatus.svg
[license-url]: https://github.com/egokcen/DLAG/blob/master/LICENSE.md

This project contains a Matlab implementation of DLAG, a dimensionality reduction framework for dissecting interactions between populations of neurons.

The DLAG model, along with learning and inference procedures, are described in
detail in the following reference:

- Gokcen, E., Jasper, A. I., Semedo, J. D., Zandvakili, A., Kohn, A., Machens, C. K. & Yu, B. M. 
Disentangling the flow of signals between populations of neurons. Preprint at https://doi.org/10.1101/2021.08.30.458230 (2021).

Please read it carefully before using the code, as it describes all of the
terminology and usage modes. Please cite the above reference if using any
portion of this code for your own purposes.

## Installation

Simply download and extract the [latest release](https://github.com/egokcen/DLAG/releases) of this project to your desired local directory.

## Project wiki

For additional information on getting started, as well as subtler usage details, see this project's wiki [here](https://github.com/egokcen/DLAG/wiki).

## Contact
For questions, please contact Evren Gokcen at egokcen@cmu.edu. 

## Acknowledgments

This DLAG implementation started with base code located [here](https://github.com/karts25/NeuralTraj).

That project contains implementations of GPFA and TD-GPFA, which are described in
the following references:

- Yu, B. M., Cunningham, J. P., Santhanam, G., Ryu, S. I., Shenoy, K. V. & Sahani, M. 
Gaussian-Process Factor Analysis for Low-Dimensional Single-Trial Analysis of Neural Population Activity. 
Journal of Neurophysiology 102, 614–635 (2009).

- Lakshmanan, K. C., Sadtler, P. T., Tyler-Kabara, E. C., Batista, A. P. & Yu, B. M.
Extracting Low-Dimensional Latent Structure from Time Series in the Presence of Delays.
Neural Computation 27, 1825–1856 (2015).
