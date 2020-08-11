# Comparison of adaptive psychophysical procedures

## License

See the [LICENSE](LICENSE.md) file for license rights and limitations. 
By downloading and/or installing this software and associated files on your computing system you agree to use the software under the terms and condition as specified in the License agreement.

## About

These scripts are used to compare different two-down one-up adaptive procedures in terms of power and duration, by simulating different procedures for different imposed psychometric functions. 

When performing an experiment with an adaptive psychophysical procedure, many choices have to be made, which should be well-considered, but, in the end, are also somewhat arbitrary (or at least a trade-off between precision and duration). These choices are the step size, when and how to change the step size (after a number of trials or reversals), the start value of each run, when a run stops (after a number of trials or reversals), how the just-noticeable difference is calculated (e.g., the average of the last number of trials or reversals), the number of runs per condition and also the number of participants. Here, we perform simulations to verify how precisely (in terms of bias and variance) different procedures are able to predict just-noticeable differences, and how much time (expressed as a number of trials) they take to do so.

NOTE: These scripts and functions were developed and tested in MATLAB R2016b.

## Quick start guide

Run main.m to get started.

All functions are documented in their respective m-files. 
 
## References
 
[1] Dieudonné, B., Van Wilderode, M., and Francart, T. (2020). Temporal quantization deteriorates the discrimination of interaural time differences. _(submitted)_



