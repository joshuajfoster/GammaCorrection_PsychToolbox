# Color lookup tables in PsychToolbox

A color lookup table (CLUT) or gamma table is a 256 x 3 matrix that is loaded onto the graphics card, which specifies the rules for how RGB values (specified in Matlab) are translated before being communicated to the monitor.

Why is the CLUT a 256 x 3 matrix? 

* The 3 columns are for the three colors (red, green, and blue).
* The 256 rows are for the 256 intensity levels (0-255) that are available to us with 8-bit color depth.

PsychToolbox (PTB) makes it really easy to swap in and out different CLUTs. Below are some useful functions.

### Accessing the current CLUT

 Once you’ve opened a PTB window called ‘win’, you can access the current CLUT with this function:

` [gammatable] = Screen('ReadNormalizedGammaTable',win);`

 The default CLUT on our run machines is the *identity CLUT*. The rgb columns of the identity CLUT are all the same, and each column is a linear function starting at 0 (no intensity) to 0.9961 (just short of full intensity).

### Loading your own CLUT

You can load your own CLUT with this line of code:

`Screen('LoadNormalizedGammaTable',win,gammaTable)`

where `win` in the current PTB window and `gammaTable` is your custom CLUT that is already loaded into the matlab workspace. 

It is be best practice to always include a line that specifies the CLUT you want at the start of your experiments. If someone has forgotten to change the CLUT back after their experiment, you can be confident that you have the CLUT you want.

Of course, you can (and should!) restore the default identity CLUT at the end of your experiment (see below).

### Loading the Identity CLUT

In PTB you can load the identity CLUT using the following line of code:

`LoadIdentityClut(win)`

