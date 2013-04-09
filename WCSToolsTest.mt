(* Mathematica Test File *)

Test[
  WCSXYtoAST[{511, 512}, {"NAXIS" -> 3, "NAXIS1" -> 512, "NAXIS2" -> 512, "NAXIS3" -> 1, 
 "CTYPE1" -> {"RA---TAN", "Decreases in value as sample index"}, 
 "CRVAL1" -> {45.83, "RA at origin (deg)"}, "CRPIX1" -> 256.`, 
 "CDELT1" -> {-0.003, "Pixel Z-width (deg)"}, 
 "CROTA1" -> {0.`, "Twist angle undefined for Z-axis"}, 
 "CTYPE2" -> {"DEC--TAN", "Increases in value as line index"}, 
 "CRVAL2" -> {63.57, "Dec at origin (deg)"}, "CRPIX2" -> 257.`, 
 "CDELT2" -> {0.003, "Pixel Y-width (deg)"}, 
 "CROTA2" -> {0.`, "Rotates +NAXIS2 into Lat axis (angle"}, 
 "CTYPE3" -> "LAMBDA", 
 "CRVAL3" -> {0.00009999999747`, "Wavelength in meters"}, 
 "CRPIX3" -> 1.`, "CROTA3" -> 0.`, "CDELT3" -> 0.`, 
 "BUNIT" -> "MJy/ster", "EQUINOX" -> {1950.`, "EME50"}, 
 "RADECSYS" -> "FK4"}]
	,
	{44.06441861768384, 64.32433165231973}
	,
	TestID->"WCSToolsTest-20100218-P0F4Z9"
]
Test[
	WCSXYtoAST[{511, 512}, {"NAXIS" -> 3, "NAXIS1" -> 512, "NAXIS2" -> 512, "NAXIS3" -> 1, 
 "CTYPE1" -> {"RA---TAN", "Decreases in value as sample index"}, 
 "CRVAL1" -> {45.83, "RA at origin (deg)"}, "CRPIX1" -> 256.`, 
 "CDELT1" -> {-0.003, "Pixel Z-width (deg)"}, 
 "CROTA1" -> {0.`, "Twist angle undefined for Z-axis"}, 
 "CTYPE2" -> {"DEC--TAN", "Increases in value as line index"}, 
 "CRVAL2" -> {63.57, "Dec at origin (deg)"}, "CRPIX2" -> 257.`, 
 "CDELT2" -> {0.003, "Pixel Y-width (deg)"}, 
 "CROTA2" -> {0.`, "Rotates +NAXIS2 into Lat axis (angle"}, 
 "CTYPE3" -> "LAMBDA", 
 "CRVAL3" -> {0.00009999999747`, "Wavelength in meters"}, 
 "CRPIX3" -> 1.`, "CROTA3" -> 0.`, "CDELT3" -> 0.`, 
 "BUNIT" -> "MJy/ster", "EQUINOX" -> {1950.`, "EME50"}, 
 "RADECSYS" -> "FK4"}]//WCSXYtoAST[%, {"NAXIS" -> 3, "NAXIS1" -> 512, "NAXIS2" -> 512, "NAXIS3" -> 1, 
 "CTYPE1" -> {"RA---TAN", "Decreases in value as sample index"}, 
 "CRVAL1" -> {45.83, "RA at origin (deg)"}, "CRPIX1" -> 256.`, 
 "CDELT1" -> {-0.003, "Pixel Z-width (deg)"}, 
 "CROTA1" -> {0.`, "Twist angle undefined for Z-axis"}, 
 "CTYPE2" -> {"DEC--TAN", "Increases in value as line index"}, 
 "CRVAL2" -> {63.57, "Dec at origin (deg)"}, "CRPIX2" -> 257.`, 
 "CDELT2" -> {0.003, "Pixel Y-width (deg)"}, 
 "CROTA2" -> {0.`, "Rotates +NAXIS2 into Lat axis (angle"}, 
 "CTYPE3" -> "LAMBDA", 
 "CRVAL3" -> {0.00009999999747`, "Wavelength in meters"}, 
 "CRPIX3" -> 1.`, "CROTA3" -> 0.`, "CDELT3" -> 0.`, 
 "BUNIT" -> "MJy/ster", "EQUINOX" -> {1950.`, "EME50"}, 
 "RADECSYS" -> "FK4"}]
	,
	result
	,
	TestID->"WCSToolsTest-20100218-W1B5B1"
]est