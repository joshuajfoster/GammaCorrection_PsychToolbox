# README

The scripts in this repo enable to create and check a gamma-corrected color lookup table (CLUT)

1. Run `runGammaCalibration_v2.m` to create a gamma-corrected CLUT (see top of scripts for instructions for use)

2. Run `checkGammaCalibration_v2.m` to check that your new CLUT results in a linear gamma function (i.e. luminance is a linear function of RGB intensity)

3. See the "Color lookup tables in PsychToolbox" README for instructions on how to work with CLUTs in PsychToolbox. Also see`JJF_17_6_v2.m` - an example task script that uses a gamma correction. 

Note: this gamma correction only applies to greys (i.e. colors when R = G = B). It does not ensure a linear function for each color (e.g. blue) separately. 

`loadIDclut.m` and `readCurrentClut.m` are handy functions that you can run from the Matlab command line to load the identity CLUT or see the current CLUT. 







