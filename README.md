WCSTools2
=========


<img style="float: right;" src="https://github.com/markusroellig/WCSTools2/blob/master/wcstools1.jpg" width="400">
GitHub version of WCSTools for Mathematica

WCSTools

WCSTools is a Mathematica package to display astronomical FITS images together with WCS 
coordinate grids. The functions in WCSTools are mostly ported from the respective routines 
in the IDL astrolib library. All credits therefore go to the original authors W. Landsman 
and R. Balsano.

The rountines in astrolib support 26 possible map projections. So far I did not managed to 
implement them all. Especially the cubic map projections have been left out for the moment. 
The details of the various projections have been presented in a series of papers from Greisen, 
E.W. & Calabretta, M.R.. You can download the papers from the following web page:
http://www.atnf.csiro.au/people/mcalabre/WCS/ . If you are 
in desperate need of any of the not yet supported map projections, contact me and I will 
implement them.

The basic functionality of the package is to calculate the celestial coordinates that belong to each {x,y} pixel of an astronomical FITS image. To do so the routines search the FITS headers ("Metadata" in Mathematica-lingo) for certain keywords that specify at which position of the sky the image had been taken. Unfortunately, FITS "standards" are somewhat loosely defined. This means, that a FITS file not necessarily possesses the required keywords for WCSTools to be able to display the coordinate grid. If this is the case send me the FITS header and I will try to update the package.



