# mwf_estimation_external_tools

Estimate Myelin Water Fraction (MWF) maps with: 
1. DECAES toolbox 
https://juliahub.com/ui/Packages/DECAES/LsYv9/0.4.1
2. Anima toolbox 
https://anima.readthedocs.io/en/latest/

See instructions in `EXAMPLE.sh` file 


## About DECAES: 
This is a fast implementation of the original MWF estimation method proposed by Prasloski et al in 2012. 
"Applications of stimulated echo correction to multicomponent T2 analysis"
https://onlinelibrary.wiley.com/doi/full/10.1002/mrm.23157

DECAES is written with Julia language and is ~10x faster (or more) than the original MATLAB implementation.

NB: 
Original MATLAB code was written by Thomas Prasloski (email: tprasloski@gmail.com). Modifications to the MATLAB code were made by Vanessa Wiggermann to enable processing on various MATLAB versions in February 2019. The Julia port is based on this modified version
Prasloski's 2012 algorithm - fast implementation with Julia 

## About Anima: 

This is a C++ implementation of the original MWF estimation method proposed by Chatterjee et al 2018. 
"Multi-compartment model of brain tissues from T2 relaxometry MRI using gamma distribution"
https://hal.archives-ouvertes.fr/hal-01744852/document

Anima software binaries are available from its website for linux Ubuntu distribution. 

We had wrapped these binaries in a docker image, that can be run on any OS. 
https://github.com/sergeicu/anima-docker