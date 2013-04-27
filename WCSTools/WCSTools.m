(*
Created on 2010/02/18

Imported from
C:\Users\Markus Roellig\Documents\Projekte\Tools\WCS Plot\wcs2.nb
*)

(* Mathematica Package *)
BeginPackage["WCSTools`"]
(* Exported symbols added here with SymbolName::usage *)  
WCSExtractAstroPars::usage = "WCSExtractAstroPars[FITSHeader] extracts astrometry parameters from a FITS image header. The output has the following form\n {CD,CDELT,CRPIX,CRVAL,CTYPE,LONGPOLE,LATPOLE,PV2,DISTORT}";
WCSXYtoAST::usage = "WCSXYtoAST[{x,y}, FITSHeader] computes R.A. and Dec from X and Y and a FITS astrometry structure. The astrometry structure must first be extracted by EXTAST from a FITS header. The offset from the reference pixel is computed and the CD matrix is applied. If distortion is present then this is corrected. If a WCS projection (Calabretta & Greisen 2002, A&A, 395, 1077) is present, then the procedure wcsxy2sph is used to compute astronomical coordinates. Angles are returned in degrees.";
WCSASTtoXY::usage = "WCSASTtoXY[{RightAscension, Declination}, FITSHeader] computes X and Y from native coordinates and a FITS  astrometry structure. The astrometry structure is extracted by extractAstroPars from a FITS header. The offset from the reference pixel is computed and the CD matrix is applied. If distortion is present then this is corrected. If a WCS projection (Calabretta & Greisen 2002, A&A, 395, 1077) is present, then the procedure wcsxy2sph is used to compute astronomical coordinates. Angles are returned in degrees.";
WCSSPHtoXY::usage = "WCSSPHtoXY[{Longitude, Latitude}, CTYPE, PV2, CRVAL, LONGPOLE, LATPOLE] converts spherical coordinates (longitude and latitude -- sky) to x and y (map) angular coordinates. To convert and y (map) coordinates to spherical (longitude and latitude or sky) coordinates. This procedure is the inverse of wcssph2xy. This is a lower level procedure -- given a FITS header, the user will usually use XYtoAST which will then call wcsxy2sph with the appropriate  parameters.";
WCSXYtoSPH::usage = "WCSXYtoSPH[{x, y},CTYPE, PV2, CRVAL, LONGPOLE, LATPOLE] converts x and y (map) coordinates to spherical coordinates. To convert and y (map) coordinates to spherical (longitude and latitude or sky) coordinates. This procedure is the inverse of wcssph2xy. This is a lower level procedure -- given a FITS header, the user will usually use XYtoAST which will then call wcsxy2sph with the appropriate  parameters.";
WCSGetPole::usage = "WCSGetPole[CRVAL, LONGPOLE, LATPOLE, \[Theta]0] computes the coordinates of the native pole for a non-polar projection. For non-polar (cylindrical or conic) projections, the native pole is not at the reference point, and WCS_GETPOLE is used to determine the position of the native pole. See section 2.4 of the paper \"Representation of Celestial Coordinates in FITS\" by Calabretta Greisen (2002, A&A, 395, 1077, also available at \n htp://fits.gsfc.nasa.gov/fits_wcs.html  \nCalled by WCSRotate";
WCSRotate::usage = "WCSRotate[{\[Phi] or Longitude,\[Theta] or Latitude},CRVAL,LONGPOLE,LATPOLE] rotates between standard (e.g. celestial) and native coordinates. wcsrotate computes a spherical coordinate rotation between native coordinates and  standard celestial coordinate system (celestial, Galactic, or ecliptic).   Applies the equations in Appendix B of the paper \"Representation of Celestial Coordinates in FITS\" by Calabretta  Greisen (2002, A&A, 395, 1077). Also see\n http://fits.gsfc.nasa.gov/fits_wcs.html";
WCSConstRA::usage = "WCSConstRA[RightAscension, y, FITSHeader] obtains the X and Y coordinates of a line of constant right ascension. Returns a set of X pixel values given an image with astrometry, and either\n (1) a set of Y pixel values, and a scalar right ascension (or longitude),  or\n (2) set of right ascension values, and a scalar Y value.\n\n In usage (1), constRA can be used to determine the (X,Y) values of a line of constant right ascension.  In usage (2), constRA can used to determine the X positions of specified RA values, along a line of constant Y.\n You can provide the option ComputeDeclination to give out the corresponding set of declination values. Default is ComputeDeclination->False";
WCSConstDec::usage = "WCSConstDec[Declination, x, FITSHeader] obtains the X and Y coordinates of a line of constant declination. Returns a set of X pixel values given an image with astrometry, and either\n (1) a set of Y pixel values, and a scalar declination (or longitude),  or\n (2) set of declination values, and a scalar Y value.\n\n In usage (1), constDec can be used to determine the (X,Y) values of a line of constant right ascension.  In usage (2), constDec can used to determine the X positions of specified RA values, along a line of constant Y.\n You can provide the option ComputeDeclination to give out the corresponding set of right ascension values. Default is ComputeRightAscension->False";
WCSBPrecess::usage = "WCSBPrecess[{RightAscension, Declination}] calculates the mean place of a star at B1950.0 on the FK4 system from the mean place at J2000.0 on the FK5 system.";
WCSJPrecess::usage = "WCSJPrecess[{RightAscension, Declination},] Calculates the mean place of a star at B1950.0 on the FK4 system from the mean place at J2000.0 on the FK5 system.";
WCSList::usage = "WCSList[{RightAscension, Declination}] Return R.A. (Right Ascension) and Dec (Declination) as list(s) in sexigesimal format. RA and Dec may be entered as either a 2 element vector or as two separate vectors (or scalars).\n WCSList takes three options: RightAscension (True), Declination (False), Hours (False). RightAscension and Declination toggle the wanted output and define what kind of input has been provided if only one argument is given. The default is to interpret a single number as declination.";
WCSDisplay::usage = "WCSDisplay[FITSData, FITSHeader] displays the content of an (astronomical) FITS image. The data is displayed by default via ArrayPlot with an overlayed WCS coordinate grid according to its astronomical position, as extracted from the FITS Header information."
WCSSupportedProjections[]:=Print[{"AZP", "TAN", "SIN", "STG", "ARC", "ZPN","ZEA", "CYP", "CAR", "MER", "CEA", "COP", "COD", "COE", "COO","BON", "SFL", "PAR", "AIT", "MOL", "SZP"}];
WCSNotSupportedProjections[]:=Print[{"AIR", "PCO", "CSC", "QSC", "TSC"}];
WCSSupportedProjections::usage="WCSSupportedProjections[] gives a list of currently supported map projections.";
WCSNotSupportedProjections::usage="WCSNotSupportedProjections[] gives a list of currently not yet supported map projections."
Begin["`Private`"] 
(* Begin Private Context *)

Options[WCSSPHtoXY] = {MapType->"DEF",PoleOffset->{10^(-7),10^(-7)}};
Options[WCSXYtoSPH] = {MapType->"DEF"};
Options[WCSGetPole] = {ORIGIN->False};
Options[WCSRotate] = {REVERSE->True,ORIGIN->False};
Options[WCSConstRA] = {ComputeDeclination->False};
Options[WCSConstDec] = {ComputeRightAscension->False};
Options[WCSBPrecess] = {ProperMotion->{0,0},StellarParallax->{},RadialVelocity->{},Epoch->2000.};
Options[WCSJPrecess] = {ProperMotion->{0,0},StellarParallax->{},RadialVelocity->{},Epoch->1950.};

Options[WCSList] = {RightAscension->True,Declination->True,Hours->False};



SyntaxInformation[WCSExtractAstroPars] = {"ArgumentsPattern"->{_}};
WCSExtractAstroPars::distrx = "Unrecognized distortion acronym: '1'";
WCSExtractAstroPars::keyx = "Keyword `1` missing!";
WCSExtractAstroPars[imgHead_] :=
    Module[ {maptypes,maptypenumber,maptype,dataDim,hasPC,hasCD,hasCDELT,hasCROTA,
    CTYPE,CD,CDELT,CRPIX,CROTA,CRVAL,PROJ,LONGPOLE,LATPOLE,DISTORT,
    theta0,cellat,PV2,CELCOOR,filter,aOrder,apOrder,bOrder,bpOrder,
    a,b,ap,bp,distortFlag},
        filter[keyword_?StringQ,maxNum_?IntegerQ,head_?ListQ] :=
            If[ maxNum>0,
                Map[If[ ListQ[#],
                        First@#,
                        #
                    ]&,
                    Select[Table[keyword<>ToString[i],{i,0,maxNum}]/.imgHead,(ListQ[#]||NumberQ[#]||(StringQ[#]&&!StringMatchQ[#,keyword~~Repeated[DigitCharacter,{0,IntegerLength[maxNum]}]]))&]
                    ],
                Map[If[ ListQ[#],
                        First@#,
                        #
                    ]&,
                    Select[{keyword}/.head,(ListQ[#]||NumberQ[#]||(StringQ[#]&&!StringMatchQ[#,keyword~~Repeated[DigitCharacter,{0,IntegerLength[maxNum]}]]))&]
                    ]
            ];
        maptypes = {"DEF","AZP","TAN","SIN","STG","ARC","ZPN","ZEA","AIR","CYP","CAR","MER",
            "CEA","COP","COD","COE","COO","BON","PCO","SFL","PAR","AIT","MOL","CSC","QSC",
            "TSC","SZP"};
        maptypenumber = Position[maptypes,maptype];
        dataDim = If[ ListQ[#],
                      First@#,
                      #
                  ]&@("NAXIS"/.imgHead);
        hasPC = MemberQ[StringTake[#,3]&/@imgHead[[All,1]],"PC1"];
        hasCD = MemberQ[StringTake[#,3]&/@imgHead[[All,1]],"CD1"];
        hasCDELT = MemberQ[imgHead[[All,1]],"CDELT1"];
        hasCROTA = MemberQ[imgHead[[All,1]],"CROTA1"];
        CTYPE = filter["CTYPE",dataDim,imgHead](* from FITS file*);
        CD = If[ hasCD,
                 Map[First,{{"CD1_1","CD1_2"},{"CD2_1","CD2_2"}}/.imgHead,{2}]
             ];
        CDELT = {1.,1.};
        CRPIX = {1,1};
        If[ hasCDELT,
            CDELT = filter["CDELT",dataDim,imgHead],
            Message[WCSExtractAstroPars::keyx,"CDELT"];
            CDELT = {1.,1.}
        ];

        (* First check CROTA2, if present use this value, if not, take value CROTA1. If none is present set CROTA=0*)
        CROTA = If[ MemberQ[imgHead[[All,1]],"CROTA2"],
                    CROTA = (If[ ArrayQ[#],
                                 First@#,
                                 #
                             ])&@("CROTA2"/.imgHead ),
                    If[ MemberQ[imgHead[[All,1]],"CROTA1"],
                        CROTA = If[ ArrayQ[#],
                                    First@#,
                                    #
                                ]&@("CROTA1"/.imgHead),
                        CROTA = 0.
                    ]
                ];
        CROTA = CROTA*Degree;
        (* calculate rotation matrix CD from CROTA*)
        If[ !hasCD,
            CD = {{Cos[CROTA],-Sin[CROTA]},{Sin[CROTA],Cos[CROTA]}}
        ];
        If[ MemberQ[imgHead[[All,1]],"CRPIX2"],
            CRPIX[[2]] = If[ ArrayQ[#],
                             First@#,
                             #
                         ]&@("CRPIX2"/.imgHead)
        ];
        If[ MemberQ[imgHead[[All,1]],"CRPIX1"],
            CRPIX[[1]] = If[ ArrayQ[#],
                             First@#,
                             #
                         ]&@("CRPIX1"/.imgHead)
        ];
        CRVAL = If[ ArrayQ[#],
                    First@#,
                    #
                ]&/@({"CRVAL1","CRVAL2"}/.imgHead )(* from FITS file*);
        (* find direction of longitude pole, longpole, or calculate it*)
        PROJ = StringTake[CTYPE[[1]],{6,8}];
        PV2 = filter["PV2_",30,imgHead];

        (*
        ;If LONPOLE (or PV1_3) is not defined in the header,then we must determine
        ;its default value.This depends on the value of theta0 (the native;longitude of the fiducial point) of the particular projection)
        *)
        Which[
            MemberQ[imgHead[[All,1]],"PV1_3"],
                LONGPOLE = If[ ArrayQ[#],
                               First@#,
                               #
                           ]&@("PV1_3"/.imgHead),
            MemberQ[imgHead[[All,1]],"LONGPOLE"],
                LONGPOLE = If[ ArrayQ[#],
                               First@#,
                               #
                           ]&@("LONGPOLE"/.imgHead),
            True,
                (
                Which[
                    MemberQ[{"ZAP","SZP","TAN","STG","SIN","ARC","ZPN","ZEA","AIR"},PROJ],
                        theta0 = 90.;,
            MemberQ[{"COP","COE","COD","COO"},PROJ],
                theta0 = PV2[[1]],
            True,
                theta0 = 0.
            ];
                CELCOOR = StringTake[CTYPE[[2]],{1,4}];
                If[ MemberQ[{"RA--","GLON","ELON"},CELCOOR],
                    cellat = CRVAL[[1]],
                    cellat = CRVAL[[2]]
                ];
                If[ cellat>= theta0,
                    LONGPOLE = 0.,
                    LONGPOLE = 180.
                ]
        )
    ];
        If[ MemberQ[imgHead[[All,1]],"LATPOLE"],
            LATPOLE = If[ ArrayQ[#],
                          First@#,
                          #
                      ]&@("LATPOLE"/.imgHead),
            LATPOLE = 180.
        ];
        If[ PROJ== "NCP",
            (StringReplace[CTYPE,"NCP"->"SIN"];
             PV2 = {0.,1/Tan[CRVAL[[2]]/Degree]};
             LONGPOLE = 180.)
        ];
            
        (* Check for DISTORTION Flag , used by SPitzer *)
        If[ StringLength[CTYPE[[1]]]>= 12,
            distortFlag = StringTake[CTYPE[[1]],{10,12}];
            Which[
                ToUpperCase@distortFlag== "SIP",
                    aOrder = First@filter["A_ORDER",0,imgHead];
                    bOrder = First@filter["A_ORDER",0,imgHead];
                    apOrder = First@filter["A_ORDER",0,imgHead];
                    bpOrder = First@filter["A_ORDER",0,imgHead];
                    a = Table[("A_"<>ToString[i]<>"_"<>ToString[j]/.imgHead)/._String->0,{i,1,aOrder+1},{j,1,aOrder+1}];
                    ap = Table[("AP_"<>ToString[i]<>"_"<>ToString[j]/.imgHead)/._String->0,{i,1,apOrder+1},{j,1,apOrder+1}];
                    b = Table[("B_"<>ToString[i]<>"_"<>ToString[j]/.imgHead)/._String->0,{i,1,bOrder+1},{j,1,bOrder+1}];
                    bp = Table[("BP_"<>ToString[i]<>"_"<>ToString[j]/.imgHead)/._String->0,{i,1,bpOrder+1},{j,1,bpOrder+1}];
                    DISTORT = {distortFlag,a,b,ap,bp};,
                True,
                    Message[WCSExtractAstroPars::distrx,distortFlag];
                    DISTORT = {};
                ],
            DISTORT = {}
        ];
        {CD,CDELT,CRPIX,CRVAL,CTYPE,LONGPOLE,LATPOLE,PV2,DISTORT}
    ]



SyntaxInformation[WCSXYtoAST] = {"ArgumentsPattern"->{{_,_},_,_}};
WCSXYtoAST[(*xy:{{_,_}..}*){xin_?NumberQ,y_?NumberQ},imgHead_,astroStructure_:{}] :=
    Module[ {CTYPE,CD,CDELT,CRPIX,CRVAL,LONGPOLE,LATPOLE,PV2,xdiff,ydiff,xsi,eta,
    coord,temp,a,d,x,DISTORT,xdif1,ydif1,b},
(*
;INPUTS:
;X-row position in pixels,scalar or vector
;Y-column position in pixels,scalar or vector
;X and Y should be in the standard IDL convention (first pixel is;0),and not the FITS convention (first pixel is 1).
;ASTR-astrometry structure,output from EXTAST procedure containing:;.CD-2 x 2 array containing the astrometry parameters CD1_1 CD1_2
;in DEGREES/PIXEL CD2_1 CD2_2
;.CDELT-2 element vector giving physical increment at reference pixel
;.CRPIX-2 element vector giving X and Y coordinates of reference pixel
;(def=NAXIS/2)
;.CRVAL-2 element vector giving R.A.and DEC of reference pixel
;in DEGREES
;.CTYPE-2 element vector giving projection types
;.LONGPOLE-scalar longitude of north pole
;.LATPOLE-scalar giving native latitude of the celestial pole
;.PV2-Vector of projection parameter associated with latitude axis
;PV2 will have up to 21 elements for the ZPN projection,up to 3
;for the SIN projection and no more than 2 for any other
;projection
;.DISTORT-Optional substructure specifying distortion parameters
;;
;OUTPUT:
;A-R.A.in DEGREES,same number of elements as X and Y
;D-Dec.in DEGREES,same number of elements as X and Y
;;RESTRICTIONS:
;Note that all angles are in degrees,including CD and CRVAL
;Also note that the CRPIX keyword assumes an FORTRAN type
;array beginning at (1,1),while X and Y give the IDL position
;beginning at (0,0).No parameter checking is performed.;
;NOTES:
;AD2XY tests for presence of WCS coordinates by the presence of a dash
;in the 5th character position in the value of CTYPE (e.g'DEC--SIN').
*)
(* Attention:
FITS data storage and image display is different from Mathematica conventions! 
1st Difference is that ArrayPlot plots the (1,1) position to the top left while astronomical image display conventionally starts to the bottom left. 
2nd difference is the display of images from the southern hemisphere(? only southern?). This leads to an inversion of the right ascension pixels in display, i.e. reversal of image array columns.
*)
        If[ ListQ@astroStructure==={}||Length@astroStructure==9,
            {CD,CDELT,CRPIX,CRVAL,CTYPE,LONGPOLE,LATPOLE,PV2,DISTORT} = astroStructure,
            {CD,CDELT,CRPIX,CRVAL,CTYPE,LONGPOLE,LATPOLE,PV2,DISTORT} = WCSExtractAstroPars[imgHead]
        ];
        (* RA-Dec system and southern hemisphere *)
        (*If[
        ToUpperCase@StringTake[CTYPE[[2]],3]== "DEC"&&CRVAL[[2]]<0,
        x=xin+2(CRPIX[[1]]-xin)(* mirroring at the CRPIX *),
        x=xin];*)
        (* RA-Dec system 
        If[
        ToUpperCase@StringTake[CTYPE[[2]],3]== "DEC",
        x=xin+2(CRPIX[[1]]-xin)(* mirroring at the CRPIX *),
        x=xin]; *)
        x = xin;
        If[ CDELT[[1]]!= 1,
            CD = {
                {CD[[1,1]]CDELT[[1]],CD[[1,2]]CDELT[[2]]},
                {CD[[2,1]]CDELT[[1]],CD[[2,2]]CDELT[[2]]}
                }
        ];
        {xdiff,ydiff} = {x-(CRPIX[[1]]-0),y-(CRPIX[[2]]-0)};

        (* Distortion check implemented, but not yet checked *)
        If[ DISTORT!= {},
            If[ DISTORT[[1]]== "SIP",
                {a,b} = DISTORT[[2;;3]];
                {xdif1,ydif1} = {xdiff,ydiff};
                Do[xdif1+=xdiff^i ydiff^j a[[i,j]],{j,Length@a},{i,Length@a}];
                Do[ydif1+=xdiff^i ydiff^j b[[i,j]],{j,Length@b},{i,Length@b}];
                {xdiff,ydiff} = {xdif1,ydif1}
            ]
        ];

        (*{xdiff,ydiff}=Table[(xy[[i]]-{CRPIX[[1]],CRPIX[[2]]}),{i,1,Length@xy}];*)
        (* subscript (index) reihenfolge unklar. In IDL laufen die indizes anders herum als in mathematica *) 
        (* xsi, und eta = x,y coordinates in Degree!! in the image plane *)
        xsi = CD[[1,1]]*xdiff+CD[[2,1]]*ydiff;
        eta = CD[[1,2]]*xdiff+CD[[2,2]]*ydiff;
        coord = StringTake[CTYPE,4];
        If[ (coord[[1]]== "DEC-"&&coord[[2]]== "RA--")||
            (coord[[1]]== "GLAT"&&coord[[2]]== "GLON")||
            (coord[[1]]== "ELAT"&&coord[[2]]== "ELON"),
            (CRVAL = Reverse[CRVAL];
             temp = xsi;
             xsi = eta;
             eta = temp)
        ];
        If[ StringTake[CTYPE[[1]],{4}]!= "-",
            {a,d} = {xsi,eta},
            {a,d} = If[ Depth@{xsi,eta}== 3,
                        Map[WCSXYtoSPH[#,CTYPE,PV2,CRVAL,LONGPOLE,LATPOLE,{0,0}]&,{xsi,eta}],
                        WCSXYtoSPH[{xsi,eta},CTYPE,PV2,CRVAL,LONGPOLE,LATPOLE,{0,0}]
                    ]
        ];
        Return[{a,d}]
    ]



WCSASTtoXY::projn = "No CTYPE specified, assuming TANgent projection.";
SyntaxInformation[WCSASTtoXY] = {"ArgumentsPattern"->{{_,_},_}};
WCSASTtoXY[(*xy:{{_,_}..}*){RA_,Dec_},imgHead_,astroStructure_:{}] :=
    Module[ {CTYPE,CD,CDELT,CRPIX,CRVAL,LONGPOLE,LATPOLE,PV2,xsi,eta,
    coord,temp,a,d,x,y,reverse,cdinv,xdif,ydif,DISTORT,ap,
    bp, xdif1,ydif1},
(*
;INPUTS:
;A-R.A.or longitude in DEGREES,scalar or vector
;D-Dec.or longitude in DEGREES,scalar or vector
;ASTR-astrometry structure,output from EXTAST procedure containing:;.CD-2 x 2 array containing the astrometry parameters CD1_1 CD1_2
;in DEGREES/PIXEL CD2_1 CD2_2
;.CDELT-2 element vector giving increment at reference point in
;DEGREES/PIXEL
;.CRPIX-2 element vector giving X and Y coordinates of reference pixel
;(def=NAXIS/2) in FITS convention (first pixel is 1,1)
;.CRVAL-2 element vector giving coordinates of the reference pixel
;in DEGREES
;.CTYPE-2 element vector giving projection types
;.LONGPOLE-scalar longitude of north pole (default=180)
;.PV2-Vector of additional parameter (e.g.PV2_1,PV2_2) needed in
;some projections
;.DISTORT-Optional substructure specifying distortion parameters
;;OUTPUTS:
;X-row position in pixels,scalar or vector
;Y-column position in pixels,scalar or vector
;;X,Y will be in the standard IDL convention (first pixel is 0),and
;*not*the FITS convention (first pixel is 1)
;NOTES:
;AD2XY tests for presence of WCS coordinates by the presence of a dash
;in the 5th character position in the value of CTYPE (e.g'DEC--SIN').
*)
        If[ ListQ@astroStructure==={}||Length@astroStructure==9,
            {CD,CDELT,CRPIX,CRVAL,CTYPE,LONGPOLE,LATPOLE,PV2,DISTORT} = astroStructure,
            {CD,CDELT,CRPIX,CRVAL,CTYPE,LONGPOLE,LATPOLE,PV2,DISTORT} = WCSExtractAstroPars[imgHead]
        ];
        coord = StringTake[CTYPE,{1,4}];
        reverse = (
                    coord[[1]]== "DEC-"&&coord[[2]]== "RA--")||
                    (coord[[1]]== "GLAT"&&coord[[2]]== "GLON")||
                    (coord[[1]]== "ELAT"&&coord[[2]]== "ELON"
                );
        If[ reverse,
            CRVAL = Reverse[CRVAL]
        ];
        If[ CTYPE[[1]]== "",
            (CTYPE = {"RA---TAN","DEC--TAN"};
             Message[WCSASTtoXY::projn])
        ];
        If[ StringTake[CTYPE[[1]],{4}]== "-",
            (* spherical proj *)
            {xsi,eta} = WCSSPHtoXY[{RA,Dec},CTYPE,PV2,CRVAL,LONGPOLE,LATPOLE,{0,0}];,
            {xsi,eta} = {a-CRVAL[[1]],d-CRVAL[[2]]}
        ];
        If[ CDELT[[1]]!= 1,
            CD = {
                {CD[[1,1]]CDELT[[1]],CD[[1,2]]CDELT[[2]]},
                {CD[[2,1]]CDELT[[1]],CD[[2,2]]CDELT[[2]]}
                }
        ];
        If[ reverse,
            (temp = xsi;
             xsi = eta;
             eta = temp)
        ];
        cdinv = Inverse@CD;

        (*CRPIX--; I think reason for this is the index convention in IDL *)
        xdif = cdinv[[1,1]]xsi+cdinv[[2,1]]eta;
        ydif = cdinv[[1,2]]xsi+cdinv[[2,2]]eta;

        (* Distortion check implemented, but not yet checked *)
        If[ DISTORT!= {},
            If[ DISTORT[[1]]== "SIP",
                {ap,bp} = DISTORT[[4;;5]];
                {xdif1,ydif1} = {xdif,ydif};
                Do[xdif1+=xdif^i ydif^j ap[[i,j]],{j,Length@ap},{i,Length@ap}];
                Do[ydif1+=xdif^i ydif^j bp[[i,j]],{j,Length@bp},{i,Length@bp}];
                {xdif,ydif} = {xdif1,ydif1}
            ]
        ];
        (* added +1 to xdif and ydif to account for position indices starting with 1 and not with 0 like in IDL *)
        {x,y} = {xdif+CRPIX[[1]],ydif+CRPIX[[2]]};
        Return[{x,y}]
    ]




WCSSPHtoXY::solx = "No valid solution";
WCSSPHtoXY::projx = "Projection type `1` not yet implemented.";
WCSSPHtoXY::defkey = "`1` not set, using default of `1` = `2` for `3` map projection.";
WCSSPHtoXY::keyx = "`1` Projection requires the keyword `2` to be `3` `4`.";

SyntaxInformation[WCSSPHtoXY] = {"ArgumentsPattern"->{{_,_},_,_,_,_,_,_,OptionsPattern[]}};
WCSSPHtoXY[{longitude_?NumberQ,latitude_?NumberQ},CTYPE_,PV2_,CRVAL_,LONGPOLE_,LATPOLE_,CRXY_:{0,0},OptionsPattern[]] :=
    Module[ {maptypes,maptype,theta0,ctype1,ctype2,projtype,g,
            theta,phi,northoffset,southoffset,pv21,pv22,
            zenithal,conic,x,y,rtheta,long,lat,mu,gamma,
            phic,thetac,xp,yp,zp,z,par,poly,sol,xtmp,rlim,
            pos,a,alpha,thetaa,aphi,y0,theta1,theta2,s1,s2,
            sthetaa,const
        },
    (*
    ;INPUT PARAMETERS:
    ;longitude-longitude of data,scalar or vector,in degrees
    ;latitude-latitude of data,same number of elements as longitude,
    ;in degrees
    ;map_type-optional positional parameter,numeric scalar (0-26)
    ;corresponding to a particular map projection.This is not a
    ;FITS standard,it is simply put in to allow function similar
    ;to that of less general map projection procedures (eg AITOFF).
    ;The following list gives the map projection types and their
    ;respective numbers.;
    ;FITS Number Name Comments
    ;code code
    ;--------------------------------------------------------------------
    ;DEF 0 Default=Cartesian
    ;AZP 1 Zenithal perspective PV2_1 required
    ;TAN 2 Gnomic AZP w/mu=0
    ;SIN 3 Orthographic PV2_1,PV2_2 optional
    ;STG 4 Stereographic AZP w/mu=1
    ;ARC 5 Zenithal Equidistant
    ;ZPN 6 Zenithal polynomial PV2_0,PV2_1....PV2_20 possible
    ;ZEA 7 Zenithal equal area
    ;AIR 8 Airy PV2_1 required
    ;CYP 9 Cylindrical perspective PV2_1 and PV2_2 required
    ;CAR 10 Cartesian
    ;MER 11 Mercator
    ;CEA 12 Cylindrical equal area PV2_1 required
    ;COP 13 Conical perspective PV2_1 and PV2_2 required
    ;COD 14 Conical equidistant PV2_1 and PV2_2 required
    ;COE 15 Conical equal area PV2_1 and PV2_2 required
    ;COO 16 Conical orthomorphic PV2_1 and PV2_2 required
    ;BON 17 Bonne's equal area PV2_1 required
    ;PCO 18 Polyconic
    ;SFL 19 Sanson-Flamsteed
    ;PAR 20 Parabolic
    ;AIT 21 Hammer-Aitoff
    ;MOL 22 Mollweide
    ;CSC 23 Cobe Quadrilateralized convergence of inverse is poor
    ;Spherical Cube
    ;QSC 24 Quadrilateralized
    ;Spherical Cube
    ;TSC 25 Tangential Spherical Cube
    ;SZP 26 Slant Zenithal Projection PV2_1,PV2_2,PV2_3 optional
    ;;OPTIONAL INPUT KEYWORD PARAMETERS:;
    ;CTYPE-One,two,or three element vector containing 8 character
    ;strings corresponding to the CTYPE1,CTYPE2,and CTYPE3
    ;FITS keywords:;
    ;CTYPE[0]-first four characters specify standard system
    ;('RA--','GLON' or'ELON' for right ascension,Galactic;longitude or ecliptic longitude respectively),second four
    ;letters specify the type of map projection (eg'-AIT' for;Aitoff projection)
    ;CTYPE[1]-first four characters specify standard system
    ;('DEC-','GLAT' or'ELAT' for declination,galactic latitude;or ecliptic latitude respectively;these must match;the appropriate system of ctype1),second four letters of
    ;ctype2 must match second four letters of ctype1.
    ;CTYPE[2]-if present must be the 8 character string,'CUBEFACE',
    ;only used for spherical cube projections to identify an axis
    ;as containing the face on which each x and y pair of
    ;coordinates lie.
    ;PV2-Vector of projection parameter associated with latitude axis
    ;PV2 will have up to 21 elements for the ZPN projection,up to 3
    ;for the SIN projection and no more than 2 for any other
    ;projection.The first element corresponds to PV2_1,the
    ;second to PV2_2,etc.
    ;CRVAL-2 element vector containing standard system coordinates (the;longitude and latitude) of the reference point
    ;CRXY-2 element vector giving the x and y coordinates of the
    ;reference point,if this is not set the offset is[0,0]
    ;This is not a FITS standard-- it is similar to CRPIX but in
    ;angular X,Y coordinates (degrees) rather than pixel coordinates
    ;LATPOLE-native latitude of the standard system's North Pole
    ;LONGPOLE-native longitude of standard system's North Pole,default
    ;is 180 degrees for Zenithal systems
    ;NORTH_OFFSET-offset (radians) added to input points near north pole.
    ;SOUTH_OFFSET-offset (radians) added to input points near south pole.
    ;BADINDEX-vector,list of transformed points too close to poles.;
    ;;OUTPUT PARAMETERS:;
    ;x-x coordinate of data,same number of elements as longitude,in
    ;degrees;if CRXY is set,then x will be returned offset by
    ;crxy(0).NOTE:x in all map projections increases to the
    ;left,not the right.
    ;y-y coordinate of data,same number of elements as longitude,in
    ;degrees;if CRXY is set,y will be returned offset by crxy[1]
    ;bad-vector returning index to transformed points close to pole.;
    ;OPTIONAL OUTPUT KEYWORD PARAMETERS:
    ;FACE-a output variable used for spherical cube projections to
    ;designate the face of the cube on which the x and y
    ;coordinates lie.Will contain the same number of elements as
    ;X and Y.Must contain at least 1 arbitrary element on input
    ;If FACE is NOT defined on input,it is assumed that the
    ;spherical cube projection is laid out over the whole sky
    ;in the "sideways T" configuration.
    ;NOTES:
    ;The conventions followed here are described in more detail in
    ;"Representations of Celestial Coordinates in FITS" by Calabretta
    ;and Greisen (2002,A&A,395,1077;also see;http://fits.gsfc.nasa.gov/fits_wcs.html).The general
    ;scheme outlined in that article is to first use WCS_ROTATE to convert
    ;coordinates in one of three standard systems (celestial,galactic,
    ;or ecliptic) into a "native system" of latitude and longitude.The
    ;latitude and longitude are then converted into x and y coordinates
    ;which depend on the map projection which is performed.The rotation
    ;from standard to native coordinates can be skipped if one so desires.
    ;This procedure necessitates two basic sections.The first converts
    ;"standard" coordinates to "native" coordinates while the second converts
    ;"native" coordinates to x and y coordinates.The first section is
    ;simply a call to WCS_ROTATE,while the second contains the guts of
    ;the code in which all of the map projection is done.This procedure
    ;can be called in a form similar to AITOFF,EQPOLE,or QDCB by calling
    ;wcssph2xy with a fifth parameter specifying the map projection by
    ;number and by not using any of the keywords related to the map
    ;projection type (e.g.CTYPE).;
    *)
        maptypes = {
            "AZP","TAN","SIN","STG","ARC","ZPN","ZEA","AIR","CYP","CAR","MER",
            "CEA","COP","COD","COE","COO","BON","PCO","SFL","PAR","AIT","MOL",
            "CSC","QSC","TSC","SZP"};
        ctype1 = StringTrim@CTYPE[[1]];
        projtype = ToUpperCase@StringTake[ctype1,{6,8}];

        (* check to see that CTYPE is set correctly *)
        If[ Length@CTYPE>= 2,
            (ctype2 = CTYPE[[2]];
             If[ ToUpperCase@StringTake[ctype2,{6,8}]!= projtype,
                 (Print["CTYPE1 and CTYPE2 must have the same projection type."];
                  Abort[])
             ];
             If[ ((ToUpperCase@StringTake[ctype1,2]== "RA"&&ToUpperCase@StringTake[ctype2,3]!= "DEC")||
                 (ToUpperCase@StringTake[ctype1,4]== "GLON"&&ToUpperCase@StringTake[ctype2,4]!= "GLAT")||
                 (ToUpperCase@StringTake[ctype1,4]== "ELON"&&ToUpperCase@StringTake[ctype2,4]!= "ELAT")),
                 (Print["CTYPE1 and CTYPE2 must have the same standard system."];
                  Abort[])
             ];
            ),
            projtype = "DEF"
        ];
        (* set projection type to default=cartesian if it not yet set or if it is DEF *)
        (* set projection type to default=cartesian if it not yet set or if it is DEF *)
        If[ !MemberQ[maptypes,projtype]||projtype== "DEF",
            projtype = "CAR"
        ];

        (*
        Convert all longitude values into the range-180 to 180 so that equationswork properly.
        *)
        long = longitude;
        If[ !ListQ[longitude],
            long = If[ longitude>= 180,
                       longitude-360,
                       longitude
                   ],(* changed -180 ->-360 typo! *)
            long = Map[If[ #>= 180,
                           #-360,
                           #
                       ]&,longitude]
        ];
        (*
        Make small offsets at poles to allow the transformations to be
        completely invertible.These generally introduce a small fractional error but only at the poles.They are necessary since all maps lose information at the poles when a rotation is applied,because all points within NORTH_ or SOUTH_OFFSET of the poles are mapped to the same points.
        *)
        lat = latitude;
        {northoffset,southoffset} = OptionValue[PoleOffset];
        If[ !ListQ[latitude],
            If[ Abs[latitude-90]<northoffset/Degree,
                lat = 90-northoffset/Degree
            ],
            lat = Map[If[ Abs[#-90]<northoffset/Degree,
                          90-northoffset/Degree,
                          #
                      ]&,latitude]
        ];
        If[ !ListQ[latitude],
            If[ Abs[latitude+90]<southoffset/Degree,
                lat = -90+southoffset/Degree
            ],
            lat = Map[If[ Abs[#+90]<southoffset/Degree,
                          -90+southoffset/Degree,
                          #
                      ]&,latitude]
        ];
            
        (*
        Convert from standard coordinate system to "native" coordinate system if the CRVAL keyword is set.Otherwise assume the latitude and longitude given are in "native" coordinates already (this is essentially what is done;in the procedure AITOFF).
        *)
        If[ ListQ@PV2&&Length@PV2>0,
            pv21 = PV2[[1]],
            pv21 = 0
        ];
        If[ ListQ@PV2&&Length@PV2>1,
            pv22 = PV2[[2]],
            pv22 = 0
        ];
        If[ VectorQ@CRVAL,
            (If[ !NumberQ@maptype,
                 maptype = Position[maptypes,projtype][[1,1]]
             ];
             conic = 16>= maptype>= 13;
             zenithal = 8>= maptype>= 1||maptype== 26;
             (*Rotate from standard celestial coordinates into the native system.*)
             Which[
                 conic,
                     theta0 = pv21,
                 zenithal,
                     theta0 = 90,
                 True,
                     theta0 = 0
                     ];
             {phi,theta} = WCSRotate[{long,lat},CRVAL,LONGPOLE,LATPOLE,theta0,REVERSE->False];
             phi*=Degree;
             theta*=Degree;),
            (phi = long*Degree;
             theta = long*Degree;)
        ];
        (*BRANCH BY MAP PROJECTION TYPE*)
        Which[
            ToUpperCase@projtype== "AZP",
                If[ pv21<0,
                    Message[WCSSPHtoXY::keyx,"AZP","PV2_1","\[GreaterEqual]",0]
                ];
                gamma = pv22 Degree;
                mu = pv21;
                rtheta = 1/Degree Cos[theta](mu+1)/((mu+Sin[theta])+Cos[theta]Cos[phi]Tan[gamma]);
                {x,y} = {rtheta Sin[phi],-rtheta (Cos[phi])/(Cos[gamma])};
                ,
            ToUpperCase@projtype== "SZP",
                mu = If[ ListQ@PV2&&Length@PV2>0,
                         PV2[[1]],
                         0
                     ];
                phic = If[ ListQ@PV2&&Length@PV2>1,
                           PV2[[2]],
                           0
                       ];
                thetac = If[ ListQ@PV2&&Length@PV2>1,
                             PV2[[3]],
                             90
                         ];
                {phic,thetac} = {phic,thetac}Degree;
                xp = -mu Cos[thetac]Sin[phic];
                yp = mu Cos[thetac]Cos[phic];
                zp = mu Sin[thetac]+1;
                x = 1/Degree(zp Cos[theta]Sin[phi]-xp (1-Sin[theta]))/(zp-(1-Sin[theta]));
                y = (-1)/Degree(zp Cos[theta]Cos[phi]+yp (1-Sin[theta]))/(zp-(1-Sin[theta]));
                ,
            ToUpperCase@projtype== "TAN",
                (
                rtheta = 1/(Tan[theta]Degree);
                x = rtheta Sin[phi];
                y = -rtheta Cos[phi];
                (* FUDGE FACTOR! somewhere there is a switched sign, but I cannot find it :-(, hence the post-correction
                {x,y}=-{x,y}; *)
                (* seems to be correct, error appears for xy-> spherical conversion*)
                ),
            ToUpperCase@projtype== "SIN",
                (
                If[ !NumberQ[pv21],
                    pv21 = 0
                ];
                If[ !NumberQ[pv22],
                    pv22 = 0
                ];
                If[ pv21== 0&&pv22== 0,
                    (
                    rtheta = Cos[theta]/Degree;
                    x = rtheta Sin[phi];
                    y = -rtheta Cos[phi];
                    ),
                    (* else NCP projection *)
                    (
                    x = (Cos[theta]Sin[phi]+pv21(1-Sin[theta]))/Degree;
                    y = (Cos[theta]Cos[phi]-pv22(1-Sin[theta]))/(-Degree)
                    );
                ]
                ),
            ToUpperCase@projtype== "STG",
                rtheta = 2/Degree Tan[(\[Pi]/2-theta)/2];
                {x,y} = {rtheta Sin[phi],-rtheta Cos[phi]};
                ,
            ToUpperCase@projtype== "ARC",
                rtheta = 1/Degree(\[Pi]/2-theta);
                {x,y} = {rtheta Sin[phi],-rtheta Cos[phi]};
                ,
            ToUpperCase@projtype== "ZPN",
                z = \[Pi]/2-theta;
                If[ Count[PV2,a_/;a== 0]>0,
                    par = FixedPoint[If[ Last[#]== 0,
                                         Drop[#,-1],
                                         #
                                     ]&,PV2](* keep PV2 up to the largest non-zero position *),
                    par = PV2[[1]]
                ];
                poly = FromDigits[Reverse[(Range[Length@par-1])*par[[2;;Length@par]]],xtmp];
                sol = Flatten[NSolve[poly== 0,xtmp]/.Rule[_,q_]:> q];
                If[ Position[sol,a_/;a \[Element] Reals]== {},
                    Message[WCSSPHtoXY::realrootx,"ZPN"]
                ];
                g = Extract[sol,Position[sol,a_/;(a \[Element] Reals&&a>0)]];
                rlim = Min[g];
                rtheta = FromDigits[Reverse[par],z]/Degree;
                {x,y} = {rtheta Sin[phi],-rtheta Cos[phi]};
                If[ ListQ[x]&&Length@g<Length@sol,
                    pos = Position[sol,#]&/@Complement[sol,g];
                    x[[Flatten[pos]]]== NaN;
                    y[[Flatten[pos]]]== NaN;
                ];
                ,
            ToUpperCase@projtype== "ZEA",
                rtheta = 2/Degree Sin[(\[Pi]/2-theta)/2];
                {x,y} = {rtheta Sin[phi],-rtheta Cos[phi]};
                ,
            (*ToUpperCase@projtype== "AIR",
            If[PV2== {},Message[wcssph2xy::defkey,"PV2_1",90,"AIR"];pv21=90];
            thetab=90;
            xi=(\[Pi]/2-theta)/2;
            (* When theta_b (aka PV2_1 in radians) is equal to pi/2 the normal equations for the AIR projection produce infinities.To avoid the problem values of theta_b equal to pi/2 cause a different set of equations to be used.*)
            If[thetab== \[Pi]/2,
            (* AIR produces the same radii for different latitudes,causing some overlap.To avoid this problem,if latitudes which are far enough south to be a problem are included in the data,the routine will stop.*)
            If[Min[theta]<-36 Degree,
            Print["AIR produces overlap of native latitudes south of -36 with the PV2_1 = 90."];Abort[]];
            (* points with xi too small are labelled as bad to prevent poor behavior of the equation for r_theta*)
            rtheta=long*0;
            If[ListQ[xi],
            pos=Flatten[Position[Abs[xi],a_/;a>= 10^(-10)]];
            If[pos!= {},
            rtheta[[pos]]=(-2)/Degree((Log[Cos[xi[[pos]]]])/(Tan[xi[[pos]]])-0.5 Tan[xi[[pos]]]);,
            xib=(\[Pi]/2-thetab)/2;
            a=Log[Cos[xib]]/Tan[xib]/Tan[xib];];,
            If[Abs[xi]>= 10^(-10),rtheta=(-2)/Degree((Log[Cos[xi]])/(Tan[xi])-0.5 Tan[xi]),
            xib=(\[Pi]/2-thetab)/2;
            a=Log[Cos[xib]]/Tan[xib]/Tan[xib];]];
            xitemp=Range[1,90]Degree;
            radius=(-1)/Degree((Log[Cos[xitemp]])/(Tan[xitemp])+Log[(Cos[xib])/(Tan[xib]Tan[xitemp])]);
            i=0;
            While[(radius[[i+1]]>= radius[[i]])||(i== Length[radius]-2),
            i++];
            If[i<Length[radius]-2,minlat=90-2/Degreexitemp[[i]],minlat=-90];
            If[Min[theta]<If[ListQ[#],First@#,#]&@minlat,
            Print["AIR produces overlap of native latitudes south of"<>ToString[If[ListQ[#],First@#,#]&@minlat]<>" with the PV2_1 = "<>ToString[pv21]];Abort[]];
            pos=Flatten[Position[Abs[xi],a_/;a>= 10^(-10)]];],*)
                
            ToUpperCase@projtype== "CYP",
                If[ PV2== {},
                    Message[WCSSPHtoXY::defkey,"PV2_1",0,"CYP"];
                    pv21 = 0
                ];
                If[ Length@PV2<2,
                    Message[WCSSPHtoXY::defkey,"PV2_2",1,"CYP"],
                    pv21 = 1
                ];
                If[ pv21== -pv22,
                    Print["PV2_1=-PV2_2 is not allowed fo the CYP map projection."]
                ];
                {x,y} = {pv22/Degree phi,1/Degree(pv21+pv22)(Sin[theta])/(pv21+Cos[theta])};,
                ToUpperCase@projtype== "CAR",
                (
                {x,y} = {phi,theta}/Degree;
        );
        ,
    ToUpperCase@projtype== "MER",
        {x,y} = {phi/Degree,Log[Tan[(\[Pi]/2+theta)/2]]/Degree};
        ,
    ToUpperCase@projtype== "CEA",
        If[ Length[PV2]== 0,
            Message[WCSSPHtoXY::keyx,"CEA","PV2_1","set",""]
        ];
        If[ (pv21<= 0)||(pv21>1),
            Message[WCSSPHtoXY::keyx,"CEA","PV2_1","0<PV2_1","\[LessEqual]1"]
        ];
        {x,y} = {phi/Degree,(Sin[theta])/(pv21 Degree)};
        ,
    ToUpperCase@projtype== "COP",
        If[ Length@PV2<1,
            Message[WCSSPHtoXY::keyx,"COP","PV2_1" ,"set",""],
            pv21 = PV2[[1]]
        ];
        If[ Length@PV2<2,
            Message[WCSSPHtoXY::defkey,"PV2_1" ,0,"COP"];
            pv22 = 0,
            pv22 = PV2[[2]]
        ];
        If[ pv21<-90||pv22>90||pv21>90,
            Print["pv2_1 and pv2_2 must satisfy -90\[LessEqual]PV2_1\[LessEqual]90, PV2_2\[LessEqual]90 for COP projection."]
        ];
        If[ pv21== -pv22,
            Print["COP projection with PV2_1=-PV2_2 is better done as a cylindrical projection."]
        ];
        thetaa = pv21 Degree;
        alpha = pv22 Degree;
        If[ (ListQ[theta]&&Position[theta,a_/;(a>thetaa+\[Pi]/2||a<= thetaa-\[Pi]/2)]!= {})||(theta>thetaa+\[Pi]/2||theta<= thetaa-\[Pi]/2),
            Print["COP map projection diverges for native latitude = PV2_1 \[PlusMinus] 90. Remove these points and try again."]
        ];
        rtheta = (Cos[alpha])/Degree(1/(Tan[thetaa])-Tan[theta-thetaa]);
        aphi = phi Sin[thetaa];
        y0 = (Cos[alpha])/(Degree Tan[thetaa]);
        {x,y} = {rtheta Sin[aphi],y0-rtheta Cos[aphi]};
        ,
    ToUpperCase@projtype== "COD",
        If[ Length@PV2<1,
            Message[WCSSPHtoXY::keyx,"COD","PV2_1" ,"set",""],
            pv21 = PV2[[1]]
        ];
        If[ Length@PV2<2,
            Message[WCSSPHtoXY::defkey,"PV2_1" ,0,"COD"];
            pv22 = 0,
            pv22 = PV2[[2]]
        ];
        If[ pv21<-90||pv22>90||pv21>90,
            Print["pv2_1 and pv2_2 must satisfy -90\[LessEqual]PV2_1\[LessEqual]90, PV2_2\[LessEqual]90 for COD projection."]
        ];
        If[ pv21== -pv22,
            Print["COD projection with PV2_1=-PV2_2 is better done as a cylindrical projection."]
        ];
        thetaa = pv21 Degree;
        If[ pv21!= 0,
            alpha = pv22 Degree;
            rtheta = thetaa-theta+alpha/(Tan[alpha]Tan[thetaa]);
            aphi = Sin[thetaa]Sin[alpha]phi/alpha;
            y0 = alpha/Degree 1/(Tan[alpha]Tan[thetaa]);,
            rtheta = thetaa-theta+1/(Tan[thetaa]);
            aphi = phi Sin[thetaa];
            y0 = 1/(Degree Tan[thetaa]);
        ];
        {x,y} = {rtheta/Degree Sin[aphi],y0-rtheta/Degree Cos[aphi]};
        ,
    ToUpperCase@projtype== "COE",
        If[ Length@PV2<1,
            Message[WCSSPHtoXY::keyx,"COE","PV2_1" ,"set",""],
            pv21 = PV2[[1]]
        ];
        If[ Length@PV2<2,
            Message[WCSSPHtoXY::defkey,"PV2_1" ,0,"COE"];
            pv22 = 0,
            pv22 = PV2[[2]]
        ];
        If[ pv21<-90||pv22>90||pv21>90,
            Print["pv2_1 and pv2_2 must satisfy -90\[LessEqual]PV2_1\[LessEqual]90, PV2_2\[LessEqual]90 for COE projection."]
        ];
        If[ pv21== -pv22,
            Print["COE gives divergent equations for PV2_1=-PV2_2."]
        ];
        theta1 = (pv21-pv22) Degree;
        theta2 = (pv21+pv22) Degree;
        s1 = Sin[theta1];
        s2 = Sin[theta2];
        sthetaa = Sin[pv21 Degree];
        gamma = s1+s2;
        rtheta = 2/Degree(Sqrt[1+s1 s2-gamma Sin[theta]])/gamma;
        aphi = phi gamma/2;
        y0 = 2/Degree(Sqrt[1+s1 s2-gamma sthetaa])/gamma;
        {x,y} = {rtheta Sin[aphi],y0-rtheta Cos[aphi]};
        ,
    ToUpperCase@projtype== "COO",
        If[ Length@PV2<1,
            Message[WCSSPHtoXY::keyx,"COO","PV2_1" ,"set",""],
            pv21 = PV2[[1]]
        ];
        If[ Length@PV2<2,
            Message[WCSSPHtoXY::defkey,"PV2_1" ,0,"COO"];
            pv22 = 0,
            pv22 = PV2[[2]]
        ];
        If[ pv21<-90||pv22>90||pv21>90,
            Print["pv2_1 and pv2_2 must satisfy -90\[LessEqual]PV2_1\[LessEqual]90, PV2_2\[LessEqual]90 for COO projection."]
        ];
        If[ pv21== -pv22,
            Print["COE gives divergent equations for PV2_1=-PV2_2."]
        ];
        theta1 = (pv21-pv22) Degree;
        theta2 = (pv21+pv22) Degree;
        thetaa = pv21 Degree;
        If[ pv22== 0,
            const = Sin[theta1],
            const = (Log[(Cos[theta2])/(Cos[theta1])])/(Log[(Tan[(\[Pi]/2-theta2)/2])/(Tan[(\[Pi]/2-theta1)/2])])
        ];
        alpha = 1/Degree(Cos[theta1])/(const(Tan[(\[Pi]/2-theta1)/2])^const);
        rtheta = alpha(Tan[(\[Pi]/2-theta)/2])^const;
        y0 = alpha (Tan[(\[Pi]/2-thetaa)/2])^const;
        aphi = const phi;
        {x,y} = {rtheta Sin[aphi],y0-rtheta Cos[aphi]};
        ,
    ToUpperCase@projtype== "BON",
        If[ Length@PV2<1,
            Message[WCSSPHtoXY::keyx,"BON","PV2_1" ,"set",""],
            pv21 = PV2[[1]]
        ];
        If[ pv21<-90||pv21>90,
            Print["PV2_1 must satisfy -90\[LessEqual]PV2_1\[LessEqual]90, for BON projection."]
        ];
        If[ pv21== 0,
            Print["PV2_1 = 0 for BON map projection is better done with SFL map projection."]
        ];
        theta1 = pv21 Degree;
        s1 = theta1/(Abs[theta1]);
        y0 = 1/(Tan[theta1])+theta1;
        a = phi (Cos[theta])/(y0-theta);
        {x,y} = {(y0-theta)Sin[a],(y0-(y0-theta)Cos[a])}/Degree
        ,
    ToUpperCase@projtype== "PCO",
        If[ ListQ[theta],
            pos = Flatten[Position[theta,a_/;a== 0]]
        ];
        x = long*0;
        y = x;
        If[ ListQ[theta]&&pos!= {},
            x[[pos]] = (phi[[pos]])/Degree;
            y[[pos]] = 0;
        ];
        If[ NumberQ[theta]&&theta== 0,
            x = phi/Degree;
            y = 0;
        ];
        If[ ListQ[theta],
            pos = Flatten[Position[theta,a_/;a!= 0]]
        ];
        If[ ListQ[theta]&&pos!= {},
            x[[pos]] = (Sin[phi[[pos]]Sin[theta[[pos]]]])/(Degree Tan[theta[[pos]]]);
            y[[pos]] = 1/Degree(theta[[pos]]+(1-Cos[phi[[pos]]Sin[theta[[pos]]]])/(Tan[theta[[pos]]]))
        ];
        If[ NumberQ[theta]&&theta!= 0,
            x = (Sin[phi Sin[theta]])/(Degree Tan[theta]);
            y = 1/Degree(theta+(1-Cos[phi Sin[theta]])/(Tan[theta]))
        ];
        ,
    ToUpperCase@projtype== "SFL",
        {x,y} = {phi Cos[theta]/Degree,theta/Degree};
        ,
    ToUpperCase@projtype== "PAR",
        {x,y} = {1/Degree phi(2 Cos[2 theta/3]-1),180 Sin[theta/3]};
        ,
    ToUpperCase@projtype== "AIT",
        alpha = (1/Degree)Sqrt[2/(1+Cos[theta]Cos[0.5 phi])];
        {x,y} = {2 alpha Cos[theta]Sin[0.5 phi],alpha Sin[theta]};
        ,
    ToUpperCase@projtype== "MOL",
        (* find a numerical solution to the equation: alpha+1/2*sin(2*alpha)-1/2*pi*sin(theta)=0 *)
        (* only scalar form yet*)
        alpha = a/.FindRoot[a+1/2 Sin[2 a]-1/2 \[Pi] Sin[theta],{a,0}];
        {x,y} = {2^1.5phi/Degree Cos[alpha]/\[Pi],(Sqrt[2])/Degree Sin[alpha]};
        ,
    True,
        (Message[WCSSPHtoXY::projx,ToUpperCase@projtype];
         Abort[])
    ];
        x = x-CRXY[[1]];
        y = y-CRXY[[2]];
        Return[{x,y}]
    ]


WCSXYtoSPH::solx = "No valid solution";
WCSXYtoSPH::realrootx = "No real root in `1` computation found.";
WCSXYtoSPH::defkey = "`1` not set, using default of `1` = `2` for `3` map projection.";
WCSXYtoSPH::keyx = "`1` Projection requires the keyword `2` to be `3` `4`.";
SyntaxInformation[WCSXYtoSPH] = {"ArgumentsPattern"->{{_,_},_,_,_,_,_,_,OptionsPattern[]}};
WCSXYtoSPH[{x_?NumberQ,y_?NumberQ},CTYPE_,PV2_,CRVAL_,LONGPOLE_,LATPOLE_,CRXY_:{0,0}] :=
    Module[ {
        maptypes,maptype,theta0,xsi,eta,ctype1,ctype2,projtype,xx,yy,r,
        theta,phi,isconic,iszenithal,longitude,latitude,pv21,pv22,a,b,
        c,rad,rad1,rad2,arad1,arad2,gamma,mu,rho,omega,theta1,theta2,
        diff1,diff2,thetac,phic,xp,yp,zp,xb,yb,rtheta,pv,pv2,sol,poly,
        xtmp,thetaa,alpha,y0,dgmin,gamma1,z2,xi,thetab,zetab,test,testPos,
        xiold,const,s1,s2,sthetaa},
    (*
    NOTES: The conventions followed here are described in more detail in the paper "Representations of Celestial Coordinates in FITS" by Calabretta & Greisen 2002, A&A, 395, 1077 also see http://fits.gsfc.nasa.gov/fits_wcs.html).The general scheme outlined in that article is to convert x and y coordinates into a "native" longitude and latitude and then rotate the system into one of three generally recognized systems (celestial,galactic or ecliptic). 
    This procedure necessitates two basic sections.The first converts x and y coordinates to "native" coordinates while the second converts "native" to "standard" coordinates.The first section contains the guts of the code in which all of the map projection is done.The second step is performed by WCS_ROTATE and only involves rotation of coordinate systems.WCSXY2SPH can be called in a form similar to 
    AITOFF,EQPOLE or QDCB by calling wcsxy2sph with a fifth parameter specifying the map projection by number and by not using any of the keywords related to the map projection type (eg ctype1 and ctyp2).
    PROCEDURE:
    The first task of the procedure is to do general error-checking to make sure the procedure was called correctly and none of the parameters or keywords conflict. This is particularly important because the procedure can be called in two ways (either using FITS-type keywords or using a number corresponding 
    a map projection type).All variables are converted into double precision values.

    The second task of the procedure is to take x and y coordinates and convert them into "native" latitude and longitude coordinates.
    Map-specific error-checking is done at this time.All of the equations were obtained from "Representations of Celestial Coordinates in FITS" and cases needing special attention are handled appropriately (see the comments with individual map projections for more information on special cases).WCS_ROTATE is then called to convert the "native" coordinates to "standard" coordinates by rotating the coordinate system.This rotation is governed by the keywords CRVAL,and LONGPOLE.The transformation is a straightforward application of euler angles. Finally, longitude values are converted into the range from 0 to 360 degrees.
    *)
        maptypes = {
            "AZP","TAN","SIN","STG","ARC","ZPN","ZEA","AIR","CYP","CAR","MER",
            "CEA","COP","COD","COE","COO","BON","PCO","SFL","PAR","AIT","MOL",
            "CSC","QSC","TSC","SZP"};
        ctype1 = StringTrim@CTYPE[[1]];
        projtype = ToUpperCase@StringTake[ctype1,{6,8}];

        (* check to see that CTYPE is set correctly *)
        If[ Length@CTYPE>= 2,
            (ctype2 = CTYPE[[2]];
             If[ ToUpperCase@StringTake[ctype2,{6,8}]!= projtype,
                 (Print["CTYPE1 and CTYPE2 must have the same projection type."];
                  Abort[])
             ];
             If[ ((ToUpperCase@StringTake[ctype1,2]== "RA"&&ToUpperCase@StringTake[ctype2,3]!= "DEC")||
                 (ToUpperCase@StringTake[ctype1,4]== "GLON"&&ToUpperCase@StringTake[ctype2,4]!= "GLAT")||
                 (ToUpperCase@StringTake[ctype1,4]== "ELON"&&ToUpperCase@StringTake[ctype2,4]!= "ELAT")),
                 (Print["CTYPE1 and CTYPE2 must have the same standard system."];
                  Abort[])
             ];
            ),
            projtype = "DEF"
        ];
        (* set projection type to default=cartesian if it not yet set or if it is DEF *)
        If[ !MemberQ[maptypes,projtype]||projtype== "DEF",
            projtype = "CAR"
        ];
        {xx,yy} = {x,y}-CRXY;
        Which[
            ToUpperCase@projtype== "AZP",
                pv21 = If[ NumberQ[PV2],
                           PV2,
                           0
                       ];
                pv21 = If[ (ListQ[PV2]&&Length@PV2>= 1),
                           PV2[[1]],
                           0
                       ]; (* PV2_1=mu, spherical radii *)
                pv22 = If[ (ListQ[PV2]&&Length@PV2>= 2),
                           PV2[[2]],
                           0
                       ]; (* pv2_2=gamma, degree *)
                If[ pv21<0,
                    Message[WCSXYtoSPH::keyx,"AZP","PV2_1",">=",0]
                ];
                gamma = pv22 Degree;
                mu = pv21;
                r = Sqrt[xx^2+yy^2(Cos[gamma])^2];
                rho = r/((mu+1)/Degree+yy Sin[gamma]);
                omega = ArcSin[(rho mu)/(Sqrt[rho^2+1])];
                xsi = ArcTan[rho,1];
                phi = ArcTan[-yy Cos[gamma],xx];
                theta1 = xsi-omega;
                theta2 = xsi+omega+\[Pi];
                theta = theta1*0;
                If[ Abs[mu]<1,
                    If[ ListQ[theta1],
                        theta = Map[If[ Abs[#]<\[Pi]/2,
                                        #,
                                        0
                                    ]&,theta1];
                        theta = MapIndexed[If[ Abs[#1]<\[Pi]/2,
                                               #1,
                                               theta[[First@#2]]
                                           ]&,theta2];,
                        If[ Abs[theta1]<\[Pi]/2,
                            theta = theta1
                        ];
                        If[ Abs[theta2]<\[Pi]/2,
                            theta = theta2
                        ];
                    ],
                    {diff1,diff2} = {Abs[\[Pi]/2-theta1],Abs[\[Pi]/2-theta2]};
                    If[ ListQ[theta1],
                        theta = MapIndexed[If[ diff1[[First@#2]]<= diff2[[First@#2]],
                                               #1,
                                               0
                                           ]&,theta1];
                        theta = MapIndexed[If[ diff2[[First@#2]]<diff1[[First@#2]],
                                               #1,
                                               theta[[First@#2]]
                                           ]&,theta2];,
                        If[ diff1<= diff2,
                            theta = theta1
                        ];
                        If[ diff2<diff1,
                            theta = theta2
                        ];
                    ]
                ];
                ,
            ToUpperCase@projtype== "SZP",
                mu = If[ PV2!= {},
                         PV2[[1]],
                         0
                     ];
                phic = If[ Length@PV2>1,
                           PV2[[2]],
                           0
                       ];
                thetac = If[ Length@PV2>2,
                             PV2[[3]],
                             90
                         ];
                {phic,thetac} = {phic,thetac}*Degree;
                xp = -mu Cos[thetac]Sin[phic];
                yp = mu Cos[thetac]Cos[phic];
                zp = mu Sin[thetac]+1;
                {xx,yy} = {xx,yy}Degree;
                {xb,yb} = {(xx-xp)/zp,(yy-yp)/zp};
                a = xb^2+yb^2+1;
                b = xb(xx-xb)+yb(yy-yb);
                c = (xx-xb)^2+(yy-yb)^2-1;
                rad = Sqrt[b^2-a c];
                {rad1,rad2} = {(-b+rad)/a,(-b-rad)/a};
                {arad1,arad2} = Abs@{rad1,rad2};
                rad = rad 0;
                If[ ListQ[rad1],
                    rad = MapIndexed[If[ arad1[[First@#2]]<= \[Pi]/2&&arad2[[First@#2]]>\[Pi]/2,
                                         #,
                                         0
                                     ]&,rad1];
                    rad = MapIndexed[If[ arad2[[First@#2]]<= \[Pi]/2&&arad1[[First@#2]]>\[Pi]/2,
                                         #,
                                         rad[[First@#2]]
                                     ]&,rad2];
                    rad = MapIndexed[If[ arad2[[First@#2]]<= \[Pi]/2&&arad1[[First@#2]]<= \[Pi]/2,
                                         Max[rad2[[First@#2]],rad1[[First@#2]]],
                                         rad[[First@#2]]
                                     ]&,rad2];,
                    If[ arad1<= \[Pi]/2&&arad2>\[Pi]/2,
                        rad = rad1
                    ];
                    If[ arad2<= \[Pi]/2&&arad1>\[Pi]/2,
                        rad = rad2
                    ];
                    If[ arad2<= \[Pi]/2&&arad1<= \[Pi]/2,
                        rad = Max[rad1,rad2]
                    ]
                ];
                theta = ArcSin[rad];
                phi = ArcTan[-(yy-yb(1-Sin[theta])),xx-xb(1-Sin[theta])];
                ,
            ToUpperCase@projtype== "TAN",
                r = Sqrt[xx^2+yy^2];
                theta = If[ NumberQ[xx],
                            If[ r>0,
                                ArcTan[1/r/Degree],
                                \[Pi]/2
                            ],
                            Map[If[ Evaluate[#>0],
                                    ArcTan[1/#/Degree],
                                    \[Pi]/2
                                ]&,r]
                        ];
                phi = If[ yy== xx== 0,
                          0,
                          ArcTan[-yy,xx]
                      (*Changed -yy to yy to account for the Mathematica convention of having {1,1} at the top left corner of the picture nad not at the bottom left. -> check on northern hemisphere*)
         (* Mathematica ArcTan ist umgekehrt zu IDL *)];
                 ,
            ToUpperCase@projtype== "CAR",
                {phi,theta} = {xx,yy}*Degree
                ,
            ToUpperCase@projtype== "SIN",
                pv21 = If[ NumberQ[PV2],
                           PV2,
                           0
                       ];
                pv21 = If[ (ListQ[PV2]&&Length@PV2>= 1),
                           PV2[[1]],
                           0
                       ];
                pv22 = If[ (ListQ[PV2]&&Length@PV2>= 2),
                           PV2[[2]],
                           0
                       ];
                If[ pv21== pv22== 0,
                    theta = ArcCos[Sqrt[xx^2+yy^2]Degree];
                    phi = ArcTan[-yy,xx],
                    (*Changed xx to -xx to account for the Mathematica convention of having {1,1} at the top left corner of the picture nad not at the bottom left. -> check on northern hemisphere*)
                    {x,y} = {xx,yy} Degree;
                    a = pv21^2+pv22^2+1;
                    b = pv21(x-pv21)+pv22(y-pv22);
                    c = (x-pv21)^2+(y-pv22)^2-1;
                    rad = Sqrt[b^2-a c];
                    rad1 = (-b+rad)/a;
                    rad2 = (-b-rad)/a;
                    arad1 = Abs[rad1];
                    arad2 = Abs[rad2];
                    rad = rad*0;
                    If[ NumberQ@rad,
                        Which[
                        arad1<= \[Pi]/2&&arad2>\[Pi]/2,rad = rad1,
                        arad2<= \[Pi]/2&&arad1>\[Pi]/2,rad = rad2,
                        arad2<= \[Pi]/2&&arad1<= \[Pi]/2,rad = Max[rad2,rad2]]
                    ];
                    If[ ListQ[rad],
                        MapIndexed[
                            Which[
                                arad1[[First@#2]]<= \[Pi]/2&&arad2[[First@#2]]>\[Pi]/2,
                                    rad[[First@#2]] = rad1[[First@#2]],
                                arad2[[First@#2]]<= \[Pi]/2&&arad1[[First@#2]]>\[Pi]/2,
                                    rad[[First@#2]] = rad2[[First@#2]],
                                arad2[[First@#2]]<= \[Pi]/2&&arad1[[First@#2]]<= \[Pi]/2,
                                    rad[[First@#2]] = Max[rad2[[First@#2]],rad2[[First@#2]]]
                                    ]&
                            ,rad]
                    ];
                    theta = ArcSin[rad];
                    phi = ArcTan[-(y-pv22(1-Sin[theta])),(x-pv21(1-Sin[theta]))];
                (*Changed xx to -xx to account for the Mathematica convention of having {1,1} at the top left corner of the picture nad not at the bottom left. -> check on northern hemisphere*)
                    ];
                ,
            ToUpperCase@projtype== "STG",
                {theta,phi} = {\[Pi]/2-2ArcTan[(Sqrt[xx^2+yy^2])/(2/Degree)],ArcTan[-yy,xx]};
        ,
    ToUpperCase@projtype== "ARC",
        {theta,phi} = {\[Pi]/2-Sqrt[xx^2+yy^2]Degree,ArcTan[-yy,xx]};
        ,
    ToUpperCase@projtype== "ZPN",
        {rtheta,phi} = {Sqrt[xx^2+yy^2]Degree,ArcTan[-yy,xx]};
        If[ Count[PV2,a_/;a== 0]>0,
            pv2 = FixedPoint[If[ Last[#]== 0,
                                 Drop[#,-1],
                                 #
                             ]&,PV2](* keep PV2 up to the largest non-zero position *),
            pv2 = PV2[[1]]
        ];
        If[ ListQ[xx],
            Do[
                pv = pv2;
                pv[[1]] = pv[[1]]-rtheta[[i]];
                poly = FromDigits[Reverse[pv],xtmp];
                sol = Flatten[NSolve[poly== 0,xtmp]/.Rule[_,q_]:> q];
                If[ Position[sol,a_/;a \[Element] Reals]== {},
                    Message[WCSXYtoSPH::realrootx,"ZPN"]
                ];
                gamma = Extract[sol,Position[sol,a_/;a \[Element] Reals]];
                If[ Length@sol>1,
                    If[ Position[sol,a_/;(\[Pi]/2>= a>= -\[Pi]/2)]== {},
                        gamma = gamma[[0]],
                        gamma = Extract[sol,Position[sol,a_/;(\[Pi]/2>= a>= -\[Pi]/2)]]
                    ]
                ];
                theta[[i]] = \[Pi]/2-gamma;
                ,{i,Length@xx}
                ];,
            pv = pv2;
            pv[[1]] = pv[[1]]-rtheta;
            poly = FromDigits[Reverse[pv],xtmp];
            sol = Flatten[NSolve[poly== 0,xtmp]/.Rule[_,q_]:> q];
            If[ Position[sol,a_/;a \[Element] Reals]== {},
                Message[WCSXYtoSPH::realrootx,"ZPN"]
            ];
            gamma = Extract[sol,Position[sol,a_/;a \[Element] Reals]];
            (* If multiple real roots are found,then we seek the value closest to the approximate linear solution*)
            If[ Length@sol>1,
                gamma1 = -pv[[1]]/pv[[2]];
                dgmin = Position[Abs[gamma-gamma1],Min[Abs[gamma-gamma1]]];
                gamma = Extract[gamma,dgmin];
                If[ Position[gamma,a_/;(\[Pi]/2>= a>= -\[Pi]/2)]== {},
                    gamma = gamma[[1]],
                    gamma = Extract[gamma,Position[gamma,a_/;(\[Pi]/2>= a>= -\[Pi]/2)]]
                ];
            ];
            If[ Length@gamma== 1,
                gamma = First@gamma
            ];
            theta = \[Pi]/2-gamma;
        ]
        ,
    ToUpperCase@projtype== "ZEA",
        {theta,phi} = {\[Pi]/2-2ArcSin[(Sqrt[xx^2+yy^2])/(2/Degree)],ArcTan[-yy,xx]};
        ,
    ToUpperCase@projtype== "AIR",
        If[ PV2== {},
            Message[WCSXYtoSPH::defkey,"PV2_1",90,"AIR"];
            pv21 = 90,
            pv21 = PV2[[1]]
        ];
        xi = thetab = pv21 Degree;
        zetab = (\[Pi]/2-thetab)/2;
        If[ thetab!= \[Pi]/2,
            a = Log[Cos[zetab]]/(Tan[zetab])^2,
            a = 0
        ];
        rtheta = (Sqrt[xx^2+yy^2])/2Degree;
        test = 0;
        While[test== 0,
            If[ ListQ[xi],
                test = Count[Abs[Exp[(-rtheta-a Tan[xi])Tan[xi]]],a_/;a>1];
                testPos = Flatten[Position[Abs[Exp[(-rtheta-a Tan[xi])Tan[xi]]],a_/;a>1]];
                If[ test== 0,
                    xi[[testPos]] = xi[[testPos]]/2
                ];,
                (* else *)
                If[ Abs[Exp[(-rtheta-a Tan[xi])Tan[xi]]]>1,
                    test = 0,
                    test = 1
                ];
                If[ test== 0,
                    xi = xi/2
                ];
            ]
            ];
        xiold = xi+1;
        While[
            Max[Abs[xiold-xi]]>10^(-12),
            xiold = xi;
            xi = ArcCos[Exp[(-rtheta-a Tan[xi])Tan[xi]]];
            ];
        {theta,phi} = {\[Pi]/2-2 xi,ArcTan[-yy,xx]};
        ,
    ToUpperCase@projtype== "CYP",
        If[ Length@PV2== 0,
            Message[WCSXYtoSPH::defkey,"PV2_1",0,"CYP"];
            pv21 = 0,
            pv21 = PV2[[1]]
        ];
        If[ Length@PV2<2,
            Message[WCSXYtoSPH::defkey,"PV2_2",1,"CYP"];
            pv22 = 1,
            pv22 = PV2[[2]]
        ];
        If[ pv21== -pv22,
            Print["PV2_1 = -PV2_2 is not allowed for CYP map projection."]
        ];
        eta = yy/((pv21+pv22)/Degree);
        theta = ArcTan[1,eta]+ArcSin[eta pv21/(Sqrt[eta^2+1])];
        phi = xx/pv22 Degree;
        ,
    ToUpperCase@projtype== "MER",
        {phi,theta} = {xx Degree,2 ArcTan[Exp[yy Degree]]-\[Pi]/2};
        ,
    ToUpperCase@projtype== "CEA",
        If[ Length@PV2<1,
            Message[WCSXYtoSPH::keyx,"CEA","PV2_1" ,"set",""],
            pv21 = PV2[[1]]
        ];
        If[ pv21<= 0||pv21>1,
            Message[WCSXYtoSPH::keyx,"CEA","PV2_1" ,"0<PV2_1\[LessEqual]","1"]
        ];
        {phi,theta} = {xx Degree,ArcSin[yy pv21 Degree]};
        ,
    ToUpperCase@projtype== "COP",
        If[ Length@PV2<1,
            Message[WCSXYtoSPH::keyx,"COP","PV2_1" ,"set",""],
            pv21 = PV2[[1]]
        ];
        If[ Length@PV2<2,
            Message[WCSXYtoSPH::defkey,"PV2_1" ,0,"COP"];
            pv22 = 0,
            pv22 = PV2[[2]]
        ];
        If[ pv21<-90||pv22>90||pv21>90,
            Print["pv2_1 and pv2_2 must satisfy -90\[LessEqual]PV2_1\[LessEqual]90, PV2_2\[LessEqual]90 for COP projection."]
        ];
        If[ pv21== -pv22,
            Print["COP projection with PV2_1=-PV2_2 is better done as a cylindrical projection."]
        ];
        thetaa = pv21 Degree;
        alpha = pv22 Degree;
        y0 = 1/Degree(Cos[alpha])/(Tan[thetaa]);
        rtheta = Sqrt[xx^2+(y0-yy)^2];
        If[ pv21<0,
            rtheta = -rtheta
        ];
        theta = thetaa+ArcTan[1/(Tan[thetaa])-(rtheta Degree)/(Cos[alpha])];
        phi = ArcTan[(y0-yy)/rtheta,xx/rtheta]/Sin[thetaa];,
    ToUpperCase@projtype== "COD",
        If[ Length@PV2<1,
            Message[WCSXYtoSPH::keyx,"COD","PV2_1" ,"set",""],
            pv21 = PV2[[1]]
        ];
        If[ Length@PV2<2,
            Message[WCSXYtoSPH::defkey,"PV2_1" ,0,"COD"];
            pv22 = 0,
            pv22 = PV2[[2]]
        ];
        If[ pv21<-90||pv22>90||pv21>90,
            Print["pv2_1 and pv2_2 must satisfy -90\[LessEqual]PV2_1\[LessEqual]90, PV2_2\[LessEqual]90 for COD projection."]
        ];
        (* use general set of equations for  the case pv21!= pv22 *)
        thetaa = pv21 Degree;
        If[ pv22!= 0,
            alpha = pv22 Degree;
            const = Sin[thetaa]Sin[alpha]/alpha;
            y0 = 1/Degree alpha/Tan[alpha]/Tan[thetaa];
            rtheta = Sqrt[xx^2+(y0-yy)^2];
            If[ pv21<0,
                rtheta = -rtheta
            ];
            theta = thetaa+alpha/(Tan[alpha]Tan[thetaa])-rtheta Degree;,
            (* special case pv21=pv22*)
            const = Sin[thetaa];
            y0 = 1/(Tan[thetaa]Degree);
            rtheta = Sqrt[xx^2+(y0-yy)^2];
            If[ pv21<0,
                rtheta = -rtheta
            ];
            theta = thetaa+1/Tan[thetaa]-rtheta Degree
        ];
        phi = ArcTan[(y0-yy)/rtheta,xx/rtheta]/const;,
    ToUpperCase@projtype== "COE",
        If[ Length@PV2<1,
            Message[WCSXYtoSPH::keyx,"COE","PV2_1" ,"set",""],
            pv21 = PV2[[1]]
        ];
        If[ Length@PV2<2,
            Message[WCSXYtoSPH::defkey,"PV2_1" ,0,"COE"];
            pv22 = 0,
            pv22 = PV2[[2]]
        ];
        If[ pv21<-90||pv22>90||pv21>90,
            Print["pv2_1 and pv2_2 must satisfy -90\[LessEqual]PV2_1\[LessEqual]90, PV2_2\[LessEqual]90 for COE projection."]
        ];
        theta2 = (pv21+pv22)Degree;
        s1 = Sin[(pv21-pv22)Degree];
        s2 = Sin[theta2];
        sthetaa = Sin[pv21 Degree];
        gamma = s1+s2;
        const = gamma/2;
        y0 = 2/Degree(Sqrt[1+s1 s2-gamma sthetaa])/gamma;
        rtheta = (xx^2+(y0-yy)^2);
        If[ pv21<0,
            rtheta = -rtheta
        ];
        phi = 2 ArcTan[(y0-yy)/rtheta,xx/rtheta]/gamma;
        theta = ArcSin[(1+s1 s2-(xx^2+(y0-yy)^2)(gamma/2Degree)^2)/gamma];,
    ToUpperCase@projtype== "COO",
        If[ Length@PV2<1,
            Message[WCSXYtoSPH::keyx,"COO","PV2_1" ,"set",""],
            pv21 = PV2[[1]]
        ];
        If[ Length@PV2<2,
            Message[WCSXYtoSPH::defkey,"PV2_1" ,0,"COO"];
            pv22 = 0,
            pv22 = PV2[[2]]
        ];
        If[ pv21<-90||pv22>90||pv21>90,
            Print["pv2_1 and pv2_2 must satisfy -90\[LessEqual]PV2_1\[LessEqual]90, PV2_2\[LessEqual]90 for COO projection."]
        ];
        theta1 = (pv21-pv22)Degree;
        theta2 = (pv21+pv22)Degree;
        thetaa = pv21 Degree;
        If[ theta1== theta2,
            const = Sin[theta1],
            const = Log[(Cos[theta2])/(Cos[theta1])]/Log[(Tan[(\[Pi]/2-theta2)/2])/(Tan[(\[Pi]/2-theta1)/2])]
        ];
        alpha = 1/Degree(Cos[theta1])/(const(Tan[(\[Pi]/2-theta1)/2])^const);
        y0 = alpha(Tan[(\[Pi]/2-thetaa)/2])^const;
        rtheta = Sqrt[xx^2+(y0-yy)^2];
        If[ pv21<0,
            rtheta = -rtheta
        ];
        phi = ArcTan[(y0-yy)/rtheta,xx/rtheta]/const;
        theta = \[Pi]/2-2ArcTan[(rtheta/alpha)^(1/const)];,
    ToUpperCase@projtype== "BON",
        If[ Length@PV2<1,
            Message[WCSXYtoSPH::keyx,"BON","PV2_1" ,"set",""],
            pv21 = PV2[[1]]
        ];
        If[ pv21<-90||pv21>90,
            Print["pv2_1 must satisfy -90\[LessEqual]PV2_1\[LessEqual]90 for BON projection."]
        ];
        If[ pv21== 0,
            Print["PV2_1 = 0 for BON projection is better done with SFL map projection."]
        ];
        theta1 = pv21 Degree;
        y0 = 1/(Tan[theta1])+theta1;
        s1 = theta1/(Abs[theta1]);
        theta = y0-s1 Sqrt[xx^2+(y0/Degree-yy)^2]Degree;
        phi = s1(y0-theta)ArcTan[(y0/Degree-yy)/(y0/Degree-theta),s1 xx/(y0/Degree-theta)]/Cos[theta];,
    ToUpperCase@projtype== "SFL",
        {phi,theta} = {(xx Degree)/(Cos[yy Degree]),yy Degree};
        ,
    ToUpperCase@projtype== "PAR",
        {phi,theta} = {(xx )/(1-4(yy/\[Pi] Degree)^2)Degree,3ArcSin[yy/\[Pi]Degree]};
        ,
    ToUpperCase@projtype== "AIT",
        z2 = 1-(xx/4Degree)^2-(yy/2Degree)^2;
        phi = 2ArcTan[2 z2-1,(Sqrt[z2]xx)/2Degree];
        theta = ArcSin[yy Sqrt[z2] Degree];
        If[ ListQ[z2],
            phi = ReplacePart[phi,Position[z2,a_/;a<0.5]->NaN];
            theta = ReplacePart[theta,Position[z2,a_/;a<0.5]->NaN];,
            If[ z2<0.5,
                phi = theta = NaN
            ]
        ];
        ,
    ToUpperCase@projtype== "MOL",
        phi = \[Pi] (xx )/(2/Degree Sqrt[2-(yy Degree)^2]);
        theta = ArcSin[(2 ArcSin[(yy Degree)/(Sqrt[2])])/\[Pi]+yy(Sqrt[2-(yy Degree)^2])/180];,
    True,
        (Print["Projection ",projtype," not yet implemented"];
         Abort[])
    ];
(* 
Convert from "native" coordinate system to "standard" coordinate system if the CRVAL keyword is set. Otherwise assume the map projection is complete. 
*)
        {phi,theta}/=Degree;
        If[ Length@CRVAL>= 2,
            If[ !NumberQ[LONGPOLE],
                LONGPOLE = 180
            ];
            maptype = Position[maptypes,projtype][[1,1]];
            isconic = (16>= maptype>= 13);
            iszenithal = ((1<= maptype<= 8)||maptype== 26);
            Which[
                isconic,
                    theta0 = PV2[[1]],
                iszenithal,
                    theta0 = 90,
                True,
                    theta0 = 0
                ];
            (*Print["wcsxy2sph:",{{phi,theta},CRVAL,LONGPOLE,LATPOLE,theta0}];*)
            {longitude,latitude} = WCSRotate[{phi,theta},CRVAL,LONGPOLE,LATPOLE,theta0,REVERSE->True]
        ];
        If[ !ListQ[longitude],
            If[ longitude<0,
                longitude+=360
            ],
            longitude = Map[(#+360)&,longitude]
        ];
        If[ !ListQ[longitude],
            If[ longitude>= (360-0.01),
                longitude-=360
            ],
            longitude = Map[(#-360)&,longitude]
        ];
        Return[{longitude,latitude}]
    ]







WCSGetPole::solx = "No valid solution";
SyntaxInformation[WCSGetPole] = {"ArgumentsPattern"->{_,_,_,_,OptionsPattern[]}};
WCSGetPole[crval_,lonpole_,latpole_,theta0_,opts:OptionsPattern[]] :=
    Module[ {
    sp,cp,sd,cd,tand,alpha0,delta0,phip,alphap,deltap,ctheta,stheta,
    term1,term2,deltap1,deltap2,sdelt},
(*
INPUT PARAMETERS:
crval-2 element vector containing standard system coordinates (the;longitude and latitude) of the reference point in degrees
lonpole-native longitude of the celestial North Pole (degrees)
LATPOLE-native latitude of the celestial North Pole (degrees)
theta0-native latitude of the fiducial point
OUTPUT PARAMETERS:
;alpha_p,delta_p-celestial longitude and latitude of the native pole
;(degrees)
*)
        alpha0 = crval[[1]]*Degree;
        delta0 = crval[[2]]*Degree;
        If[ theta0== 90,
            Return[{alpha0,delta0}]
        ];
        (*
        ;Longpole is the longitude in the native system of the North Pole in the
        ;standard system (default=180 degrees).
        *)
        phip = lonpole*Degree;
        sp = Sin[phip];
        cp = Cos[phip];
        sd = Sin[delta0];
        cd = Cos[delta0];
        tand = Tan[delta0];
        If[ theta0== 0,
            (If[ delta0== 0&&lonpole== 90,
                 deltap = latpole,
                 deltap = ArcCos[sd/cp]
             ];
             If[ latpole!= 90,
                 If[ Abs[latpole+deltap]<Abs[latpole-deltap],
                     deltap = -deltap
                 ]
             ];
             If[ lonpole== 180||cd== 0,
                 alphap = alpha0,
                 alphap = alpha0-ArcTan[-Tan[deltap]*tand,sp/cd]
             ]),
            (*  ;General case for arbitary theta0 *)
            (
            ctheta = Cos[theta0 Degree];
            stheta = Sin[theta0 Degree];
            term1 = ArcTan[ctheta cp,stheta];
            term2 = ArcCos[sd/Sqrt[1-ctheta^2 sp^2]];
            If[ term2== 0,
                deltap = term1,
                deltap1 = Abs[(term1+term2)/Degree];
                deltap2 = Abs[(term1-term2)/Degree];
                Which[
                    (deltap1>90&&deltap2>90),
                        Message[WCSGetPole::solx];
                        Abort[],
                    (deltap1<= 90&&deltap2>90),
                        deltap = term1+term2,
                    (deltap1>90&&deltap2<= 90),
                        deltap = term1-term2,
                    True,
                        deltap1 = (term1+term2)/Degree;
                        deltap2 = (term1-term2)/Degree;
                        If[ Abs[latpole-deltap1]<Abs[latpole-deltap2],
                            deltap = term1+term2,
                            deltap = term1-term2
                        ];
                ];
                If[ cd== 0,
                    alphap = alpha0,
                    (sdelt = Sin[deltap];
                     Which[
                         sdelt== 1,
                             alphap = alpha0-phip-\[Pi],
                         sdelt== -1,
                             alphap = alpha0-phip,
                         True,
                             alphap = alpha0-ArcTan[(sp ctheta)/cd,(stheta-Sin[deltap]sd)/(Cos[deltap]cd)]
                         ]
                    )
                ];
            ])
        ];
        Return[{alphap,deltap}];
    (*
        ;alpha_p,delta_p-celestial longitude and latitude of the native pole
        ;(degrees)
        *)
        ]



SyntaxInformation[WCSRotate] = {"ArgumentsPattern"->{{_,_},_,_,_,_,OptionsPattern[]}};
WCSRotate[{in1_,in2_},CRVAL_,LONGPOLE_,LATPOLE_,theta0_,opts:OptionsPattern[]] :=
    Module[ {
    longitude,latitude,alphap,deltap,phip,cp,sp,sa,ca,sd,cd,r,
    isreverse,phi1,theta1,l,m,n,b0,b1,b2,origin,phi,theta},
        isreverse = OptionValue[REVERSE];
        If[ !isreverse,
            longitude = in1;
            latitude = in2,
            phi = in1;
            theta = in2
        ];
        origin = OptionValue[ORIGIN];
        (* 
        Longpole is the longitude in the native system of the North Pole in the standard system (default=180 degrees).
        *)
        If[ !NumberQ[LONGPOLE],
            LONGPOLE = 180
        ];
        phip = LONGPOLE*Degree;
        sp = Sin[phip];
        cp = Cos[phip];
        (*
        If Theta0=90 then CRVAL gives the coordinates of the origin in the native system. This must be converted using Eq.7 in Greisen & Calabretta with theta0=0) to give the coordinates of the North pole (alpha_p,delta_p)
        *)
        If[ theta0== 90,
            {alphap,deltap} = {CRVAL[[1]],CRVAL[[2]]}*Degree,
            {alphap,deltap} = WCSGetPole[CRVAL,LONGPOLE,LATPOLE,theta0,ORIGIN->origin]
        ];
        {sa,ca,sd,cd} = {Sin[alphap],Cos[alphap],Sin[deltap],Cos[deltap]};
        (* rotation matrix r *)
        r =
            {
                {-sa sp-ca cp sd,ca sp-sa cp sd,cp cd},
                {sa cp-ca sp sd,-ca cp-sa sp sd,sp cd},
                {ca cd,sa cd,sd}
            };
        (* 
        solve the set of equations for each datum point
         *)
        If[ isreverse,
            (
            {latitude,longitude} = {phi,theta};
            (* should check for NaNs, later *)
            {phi1,theta1} = {phi,theta}*Degree;
            r = Transpose@r
            ),
            (
            {phi,phi1,theta1} = {longitude,longitude*Degree,latitude* Degree}
            )
        ];

        (* define the right-hand side of the equations *)
        {l,m,n} = {Cos[theta1]Cos[phi1],Cos[theta1]Sin[phi1],Sin[theta1]};
        (* 
        ;find solution to the system of equations and put it in b
        ;Can't use matrix notation in case l,m,n are vectors
        *)
        {b0,b1,b2} =
        {
            r[[1,1]]l+r[[1,2]]m+r[[1,3]]n,
            r[[2,1]]l+r[[2,2]]m+r[[2,3]]n,
            r[[3,1]]l+r[[3,2]]m+r[[3,3]]n
        };
        b2 = Min[Max[b2,-1],1];
        (* 
        use b0,b1,b2 to compute "native" latitude and longitude
        *)
        If[ isreverse,
            (
            {latitude,longitude} = {ArcSin[b2]/Degree,ArcTan[b0,b1]/Degree};
            Return[{longitude,latitude}]
            ),
            (
            {theta,phi} = {ArcSin[b2]/Degree,ArcTan[b0,b1]/Degree};
            Return[{phi,theta}]
            )
        ]
    ]



WCSConstRA::projx = "Program only works for TAN and SIN projections.";
SyntaxInformation[WCSConstRA] = {"ArgumentsPattern"->{_,_,_,_,OptionsPattern[]}};
WCSConstRA[RA_,y_,FITSHeader_,astroStructure_:{},opts:OptionsPattern[]] :=
    Module[ {
    CD,CTYPE,CDELT,CRPIX,PV2,CRVAL,LONGPOLE,LATPOLE,cdelta,cdel0
    ,sdel0,delra,cdelra,ctype,sdelra,cdi,yy,delta,a,b,c,d,aa,bb,
    denom,x,dectrue,DISTORT},
(*
;INPUTS:
;RA-Right Ascension value in DEGREES (0<RA<360.).If Y is a
;vector,then RA must be a scalar
;Y-Specified Y pixel value(s) for line of constant right ascension
;If RA is a vector,then Y must be a scalar
;ASTR-Astrometry structure as extracted from a FITS header by the
;procedure EXTAST
;OUTPUTS
;X-Computed set of X pixel values.The number of elements of X
;is the maximum of the number of elements of RA and Y.
;OPTIONAL OUTPUT:
;DEC-Computed set of declinations (in DEGREES) for X,Y,coordinates
;NOTES:
;The algorithm (and notation) is based on AIPS Memo 27 by Eric Greisen,
;with modifications for a coordinate description (CD) matrix as
;described in Paper II of Calabretta& Greisen (2002,A&A,395,1077).
;These documents are available from
;http://www.cv.nrao.edu/fits/documents/wcs/wcs.html
;;RESTRICTIONS:
;Implemented only for the TANgent,SIN and CARtesian projections
;
*)
(* Ported from IDL function FUNCTION CONS_RA written by Written,Wayne Landsman STX Co.April,1988, Algorithm adapted from AIPS memo No.27 by Eric Greisen
*)
        If[ ListQ@astroStructure==={}||Length@astroStructure==9,
            {CD,CDELT,CRPIX,CRVAL,CTYPE,LONGPOLE,LATPOLE,PV2,DISTORT} = astroStructure,
            {CD,CDELT,CRPIX,CRVAL,CTYPE,LONGPOLE,LATPOLE,PV2,DISTORT} = WCSExtractAstroPars[FITSHeader]
        ];

        (*    {CD,CDELT,CRPIX,CRVAL,CTYPE,LONGPOLE,LATPOLE,PV2,DISTORT}=WCSExtractAstroPars[FITSHeader];*)
        CRVAL = CRVAL*Degree;
        cdelta = {{CDELT[[1]],0.},{0.,CDELT[[2]]}};
        CD = CD*Degree;
        {cdel0,sdel0} = {Cos[CRVAL[[2]]],Sin[CRVAL[[2]]]};
        delra = RA*Degree-CRVAL[[1]];
        {cdelra,sdelra} = {Cos[delra],Sin[delra]};
        ctype = ToUpperCase@StringTake[CTYPE[[1]],{6,8}];
        Which[
            ctype== "TAN",
                ((*cdi=Inverse[cdelta.CD];*)
                (* IDL # is the invers matrix multiplication, i.e. A#B multiplies the columns of A with the rows of B. Regular dot product act oppositely. Hence, A#B \[LongRightArrow]B.A in Mathematica*)
                cdi = Inverse[CD.cdelta];
                yy = y-CRPIX[[2]];
                delta = ArcTan[(sdel0 cdelra cdi[[2,2]]-Sin[delra]cdi[[1,2]]+yy cdelra cdel0)/(cdel0 cdi[[2,2]]-yy sdel0)];
                ),
            ctype=="SIN",
                (
                {a,b} = -CD[[All,1]]CDELT[[1]];
                {c,d} = CD[[All,2]]CDELT[[2]];
                yy = (y-CRPIX[[2]])(b c-a d);(* new coordinate origin *)
                aa = cdel0 d;
                bb = sdel0 cdelra d+sdelra b;
                denom = Sqrt[aa^2+bb^2];
                delta = ArcTan[bb/aa]+ArcSin[yy/denom];
                );,
    ctype== "CAR",
        (
        {a,b} = -CD[[All,1]]CDELT[[1]];
        {c,d} = CD[[All,2]]CDELT[[2]];
        delta = (y-CRPIX[[2]])(b c-a d)+CRVAL[[2]];(* new coordinate origin *)
        If[ Length@delta== 0&&Length@RA>0,
            delta = ConstantArray[delta,Length@RA]
        ]
        );,
    True,
        (Message[WCSConstRA::projx];
         Abort[]);
    ];
        delta = delta/Degree;
        x = WCSASTtoXY[{RA,delta},FITSHeader];
        (* x contains now the coordinate pair x,y *)
        dectrue = OptionValue[ComputeDeclination];
        If[ dectrue,
            Return[x,delta],
            Return[x]
        ]
    ]



WCSConstDec::projx = "Program only works for TAN and SIN projections.";
SyntaxInformation[WCSConstDec] = {"ArgumentsPattern"->{_,_,_,OptionsPattern[]}};
WCSConstDec[Dec_,x_,FITSHeader_,astroStructure_:{},opts:OptionsPattern[]] :=
    Module[ {
    CD,CDELT,CRPIX,CTYPE,PV2,CRVAL,LONGPOLE,LATPOLE,cdelta,cdel0,
    sdel0,ctype,delta,a,b,c,d,aa,bb,denom,dectrue,xx,sign,alpha,
    y,DISTORT},
(*
;INPUTS:
;DEC-Declination value(s) in DEGREES (-!PI/2<DEC<!PI/2).
;If X is a vector,then DEC must be a scalar.
;X-Specified X pixel value(s) for line of constant declination
;If DEC is a vector,then X must be a scalar.
;ASTR-Astrometry structure,as extracted from a FITS header by the
;procedure EXTAST
;OUTPUT:
;Y-Computed set of Y pixel values.The number of Y values is the
;same as either DEC or X,whichever is greater.;
;OPTIONAL OUTPUT:
;ALPHA-the right ascensions (DEGREES) associated with the (X,Y) points
;;RESTRICTIONS:
;Implemented only for the TANgent,SIN and CAR projections
;;NOTES:
;The algorithm (and notation) is based on AIPS Memo 27 by Eric Greisen,
;with modifications for a coordinate description (CD) matrix as
;described in Paper II of Greisen& Calabretta (2002,A&A,395,1077).
;These documents are available from
;http://www.cv.nrao.edu/fits/documents/wcs/wcs.html

*)
(* Ported from IDL function FUNCTION CONS_DEC written by Written,Wayne Landsman STX Co.April,1988, Algorithm adapted from AIPS memo No.27 by Eric Greisen
*)
        If[ ListQ@astroStructure==={}||Length@astroStructure==9,
            {CD,CDELT,CRPIX,CRVAL,CTYPE,LONGPOLE,LATPOLE,PV2,DISTORT} = astroStructure,
            {CD,CDELT,CRPIX,CRVAL,CTYPE,LONGPOLE,LATPOLE,PV2,DISTORT} = WCSExtractAstroPars[FITSHeader]
        ];

        (*    {CD,CDELT,CRPIX,CRVAL,CTYPE,LONGPOLE,LATPOLE,PV2,DISTORT}=WCSExtractAstroPars[FITSHeader];*)
        CRVAL = CRVAL*Degree;
        cdelta = {{CDELT[[1]],0.},{0.,CDELT[[2]]}};
        CD = CD*Degree;
        a = -CD[[1,1]]CDELT[[1]];
        b = -CD[[2,1]]CDELT[[1]];
        c = CD[[1,2]]CDELT[[2]];
        d = CD[[2,2]]CDELT[[2]];

        (* a=-CD[[1,1]]CDELT[[1]];*)
        xx = x-CRPIX[[1]];
        sdel0 = Sin[CRVAL[[2]]];
        cdel0 = Cos[CRVAL[[2]]];
        ctype = ToUpperCase@StringTake[CTYPE[[1]],{6,8}];
        Which[
            ctype== "TAN",
                (
                aa = d;
                bb = (b c-d a) xx cdel0+sdel0 b;
                sign = If[ aa>0,
                           1,
                           -1
                       ];
                alpha = CRVAL[[1]]+ArcTan[bb/aa]+sign ArcSin[Tan[Dec Degree]((b c-d a)xx sdel0-b cdel0)/(Sqrt[aa^2+bb^2])];
                ),
            ctype=="SIN",
                (
                aa = d;
                bb = b sdel0;
                sign = If[ aa>0,
                           1,
                           -1
                       ];
                denom = Cos[Dec Degree]Sqrt[aa^2+bb^2];
                alpha = CRVAL[[1]]+ArcTan[bb/aa]+sign ArcSin[((b c-a d)xx-Sin[Dec Degree]cdel0 b)/denom];
                );,
    ctype== "CAR",
        (
        alpha = CRVAL[[1]]+(b c-a d)xx;
        If[ AtomQ[alpha]&&ListQ[Dec],
            alpha = ConstantArray[alpha,Length@Dec]
        ];
        );,
    True,
        (Message[WCSConstDec::projx];
         Abort[]);
    ];
        alpha = alpha/Degree;
        y = WCSASTtoXY[{alpha,Dec},FITSHeader];
        (* y contains now the coordinate pair x,y *)
        dectrue = OptionValue[ComputeRightAscension];
        If[ dectrue,
            Return[y,delta],
            Return[y]
        ]
    ]



WCSBPrecess::projn = "No CTYPE specified, assuming TANgent projection.";

SyntaxInformation[WCSBPrecess] = {"ArgumentsPattern"->{{_,_},OptionsPattern[]}};
WCSBPrecess[(*xy:{{_,_}..}*){RAJ2000_,DecJ2000_},opts:OptionsPattern[]] :=
    Module[ {
    MURADEC,PARALLAX,RADVEL,m,sectoradian,rarad,decrad,sinra,cosra,sindec,
    cosdec,adot,dec1950,ra1950,mua,mud,a,r0,r0dot,R0,R1,x1,rmag,s1,s,x,r2,
    r1,rdot,xdot,EPOCH,r1dot,s1dot,y1,z1,r,y,z,muradecset,ydot,zdot,neg},
        MURADEC = OptionValue[ProperMotion]; (* proper motion in Ra Dec in Arcsec per Century *)
        PARALLAX = OptionValue[StellarParallax]; (*parallax in arc sec *)
        RADVEL = OptionValue[RadialVelocity]; (*radial velocity in km/s *)
        EPOCH = OptionValue[Epoch];(* default is 2000. *)
        If[ RADVEL== {},
            If[ NumberQ[RAJ2000],
                RADVEL = {0}
            ];
            If[ ListQ[RAJ2000],
                RADVEL = ConstantArray[0,Length@RAJ2000]
            ],
            If[ Head[RAJ2000]!= Head[RADVEL],
                Print["ERROR - RAD_VEL keyword must be of same type as RA/Dec"];
                Abort[]
            ];
            If[ ListQ[RAJ2000]&&ListQ[RADVEL]&&Length[RAJ2000]!= Length@RADVEL,
                Print["ERROR - RAD_VEL  must be of same dimension as RA/Dec"];
                Abort[]
            ];
        ];
        If[ NumberQ[RADVEL],
            RADVEL = {RADVEL}
        ];
        muradecset = False;
        If[ MURADEC!= {},
            muradecset = True,
            If[ ListQ[RAJ2000],
                If[ Dimensions@MURADEC!= {Length@RAJ2000,Length@RAJ2000},
                    Print["ERROR - MU_RADEC (proper motion) must be dimensioned "<>ToString[Length@RAJ2000]<>" x "<>ToString[Length@RAJ2000]"."];
                ],
                If[ Dimensions@MURADEC!= {2},
                    Print["ERROR - MU_RADEC (proper motion) must be a 2-dim vector."]
                ];
            ],
            MURADEC = {{0,0},{0,0}};
        ];
        If[ VectorQ[MURADEC],
            MURADEC = {MURADEC}
        ];
        If[ PARALLAX== {},
            If[ NumberQ[RAJ2000],
                PARALLAX = {0}
            ];
            If[ ListQ[RAJ2000],
                PARALLAX = ConstantArray[0,Length@RAJ2000]
            ]
        ];
        If[ NumberQ[PARALLAX],
            PARALLAX = {PARALLAX}
        ];
        sectoradian = Degree/3600;
        m = {
            {0.9999256795,-0.0111814828,-0.0048590040,-0.000551,-0.238560,0.435730},
            {0.0111814828,0.9999374849,-0.0000271557,0.238509,-0.002667,-0.008541},
            {0.0048590039,-0.0000271771,0.9999881946,-0.435614,0.012254,0.002117},
            {-0.00000242389840,0.00000002710544,+0.00000001177742,0.99990432,-0.01118145,-0.00485852},
            {-0.00000002710544,-0.00000242392702,+0.00000000006585,0.01118145,0.99991613,-0.00002716},
            {-0.00000001177742,+0.00000000006585,-0.00000242404995,0.00485852,-0.00002717,+0.99996684}
            };
        adot = 10^(-3){1.244,-1.579,-0.660};(*in arc seconds per century*)
        (* switched to vector based computation! *)
        rarad = If[ !ListQ[RAJ2000],
                    {RAJ2000 },
                    RAJ2000
                ]Degree;
        decrad = If[ !ListQ[DecJ2000],
                     {DecJ2000},
                     DecJ2000
                 ] Degree;
        cosra = Cos[rarad];
        sinra = Sin[rarad];
        cosdec = Cos[decrad];
        sindec = Sin[decrad];
        dec1950 = If[ !ListQ[RAJ2000],
                      {RAJ2000},
                      RAJ2000
                  ]*0;
        ra1950 = If[ !ListQ[DecJ2000],
                     {DecJ2000},
                     DecJ2000
                 ]*0;
        Do[
            a = 10^(-6){-1.62557,-0.31919,-0.13843}; (*in radians*)
            r0 = {cosra[[i]]cosdec[[i]],sinra[[i]]cosdec[[i]],sindec[[i]]};
            If[ MURADEC!= {},
                mua = MURADEC[[i,1]];
                mud = MURADEC[[i,2]];
                r0dot = {
                    -mua sinra[[i]]cosdec[[i]]-mud cosra[[i]]sindec[[i]],
                    mua cosra[[i]]cosdec[[i]]-mud sinra[[i]]sindec[[i]],
                    mud cosdec[[i]]
                    }+21.095 RADVEL[[i]] PARALLAX[[i]] r0;(* velocity vector*),
                r0dot = {0,0,0}
            ];
            R0 = Join[r0,r0dot];
            R1 = R0.m;
            (*Include the effects of the E-terms of aberration to form r and r_dot.*)
            r1 = R1[[1;;3]];
            r1dot = R1[[4;;6]];
            If[ MURADEC== {{0,0},{0,0}},
                r1 = r1+sectoradian r1dot (EPOCH-1950)/100;
                a = a+sectoradian adot (EPOCH-1950)/100;
            ];
            x1 = R1[[1]];
            y1 = R1[[2]];
            z1 = R1[[3]];
            rmag = Sqrt[x1^2+y1^2+z1^2];
            s1 = r1/rmag;
            s1dot = r1dot/rmag;
            s = s1;
            Do[
                r = s1+a-Total[s a]s;
                s = r/rmag;
                ,{3}];
            {x,y,z} = r[[1;;3]];
            r2 = x^2+y^2+z^2;
            rmag = Sqrt[r2];
            If[ muradecset,
                rdot = s1dot+adot-Total[s adot]s;
                {xdot,ydot,zdot} = rdot[[1;;3]];
                MURADEC[[i,1]] = (x ydot-y xdot)/(x^2+y^2);
                MURADEC[[i,2]] = (zdot (x^2+y^2)-z(x xdot +y ydot))/(r2 Sqrt[x^2+y^2])
            ];
            dec1950[[i]] = ArcSin[z/rmag];
            ra1950[[i]] = ArcTan[x,y];
            If[ PARALLAX[[i]]>0,
                RADVEL[[i]] = (x xdot+y ydot +z zdot)/(21.095 PARALLAX[[i]] rmag);
                PARALLAX[[i]] = PARALLAX[[i]]/rmag
            ];
            ,{i,Length@rarad}
            ];
        neg = Flatten[Position[ra1950,a_/;a<0]];
        If[ Length@neg>0,
            ra1950[[neg]] = ra1950[[neg]]+2\[Pi]
        ];
        (* Make output scalar if input was scalar *)
        If[ !ListQ[RAJ2000],
            ra1950 = First@ra1950;
            dec1950 = First@dec1950
        ];
        Return[{ra1950,dec1950}/Degree]
    ]





WCSJPrecess::projn = "No CTYPE specified, assuming TANgent projection.";
SyntaxInformation[WCSJPrecess] = {"ArgumentsPattern"->{{_,_},OptionsPattern[]}};
WCSJPrecess[(*xy:{{_,_}..}*){RAB1950_,DecB1950_},opts:OptionsPattern[]] :=
    Module[ {
    MURADEC,PARALLAX,RADVEL,m,sectoradian,rarad,decrad,sinra,cosra,sindec,cosdec,adot,
    dec2000,ra2000,mua,mud,a,rr,r0,r0dot,R1,rr1,R,v,t,rmag,x,r2,r1,xdot,EPOCH,r1dot,y,z,muradecset,
    ydot,zdot,neg},
        MURADEC = OptionValue[ProperMotion]; (* proper motion in Ra Dec in Arcsec per Century *)
        PARALLAX = OptionValue[StellarParallax]; (*parallax in arc sec *)
        RADVEL = OptionValue[RadialVelocity]; (*radial velocity in km/s *)
        EPOCH = OptionValue[Epoch];(* default is 2000. *)
        If[ RADVEL== {},
            If[ NumberQ[RAB1950],
                RADVEL = {0}
            ];
            If[ ListQ[RAB1950],
                RADVEL = ConstantArray[0,Length@RAB1950]
            ],
            If[ Head[RAB1950]!= Head[RADVEL],
                Print["ERROR - RAD_VEL keyword must be of same type as RA/Dec"];
                Abort[]
            ];
            If[ ListQ[RAB1950]&&ListQ[RADVEL]&&Length[RAB1950]!= Length@RADVEL,
                Print["ERROR - RAD_VEL  must be of same dimension as RA/Dec"];
                Abort[]
            ];
        ];
        If[ NumberQ[RADVEL],
            RADVEL = {RADVEL}
        ];
        muradecset = False;
        If[ MURADEC!= {},
            muradecset = True,
            If[ ListQ[RAB1950],
                If[ Dimensions@MURADEC!= {Length@RAB1950,Length@RAB1950},
                    Print["ERROR - MU_RADEC (proper motion) must be dimensioned "<>ToString[Length@RAB1950]<>" x "<>ToString[Length@RAB1950]"."];
                ],
                If[ Dimensions@MURADEC!= {2},
                    Print["ERROR - MU_RADEC (proper motion) must be a 2-dim vector."]
                ];
            ],
            MURADEC = {{0,0},{0,0}};
        ];
        If[ VectorQ[MURADEC],
            MURADEC = {MURADEC}
        ];
        If[ PARALLAX== {},
            If[ NumberQ[RAB1950],
                PARALLAX = {0}
            ];
            If[ ListQ[RAB1950],
                PARALLAX = ConstantArray[0,Length@RAB1950]
            ]
        ];
        If[ NumberQ[PARALLAX],
            PARALLAX = {PARALLAX}
        ];
        sectoradian = Degree/3600;
        m = {
            {0.9999256782,+0.0111820610,+0.0048579479,-0.000551,+0.238514,-0.435623},
            {-0.0111820611,+0.9999374784,-0.0000271474,-0.238565,-0.002667,+0.012254},
            {-0.0048579477,-0.0000271765,+0.9999881997,+0.435739,-0.008541,+0.002117},
            {0.00000242395018,+0.00000002710663,+0.00000001177656,+0.99994704,+0.01118251,+0.00485767},
            {-0.00000002710663,+0.00000242397878,-0.00000000006582,-0.01118251,+0.99995883,-0.00002714},
            {-0.00000001177656,-0.00000000006587,0.00000242410173,-0.00485767,-0.00002718,1.00000956}
            };
        adot = 10^(-3){1.244,-1.579,-0.660};(*in arc seconds per century*)
        a = 10^(-6){-1.62557,-0.31919,-0.13843}; (*in radians*)
        If[ EPOCH!= 1950,
            a = a+sectoradian adot(EPOCH-1950)/100
        ];
        (* switched to vector based computation! *)
        rarad = If[ !ListQ[RAB1950],
                    {RAB1950 },
                    RAB1950
                ]Degree;
        decrad = If[ !ListQ[DecB1950],
                     {DecB1950},
                     DecB1950
                 ] Degree;
        cosra = Cos[rarad];
        sinra = Sin[rarad];
        cosdec = Cos[decrad];
        sindec = Sin[decrad];
        dec2000 = If[ !ListQ[RAB1950],
                      {RAB1950},
                      RAB1950
                  ]*0;
        ra2000 = If[ !ListQ[DecB1950],
                     {DecB1950},
                     DecB1950
                 ]*0;
        Do[    
            r0 = {cosra[[i]]cosdec[[i]],sinra[[i]]cosdec[[i]],sindec[[i]]};
            If[ MURADEC!= {},
                mua = MURADEC[[i,1]];
                mud = MURADEC[[i,2]];
                r0dot = {
                    -mua sinra[[i]]cosdec[[i]]-mud cosra[[i]]sindec[[i]],
                    mua cosra[[i]]cosdec[[i]]-mud sinra[[i]]sindec[[i]],
                    mud cosdec[[i]]
                    }+21.095 RADVEL[[i]] PARALLAX[[i]] r0;(* velocity vector*),
                r0dot = {0,0,0}
            ];
            (* Remove the effects of the E-terms of aberration to form r1 and r1_dot. *)
            r1 = r0-a+Total[r0 a]r0;
            r1dot = r0dot-adot+Total[r0 adot]r0;
            R1 = Join[r1, r1dot];
            R = R1.m;
            If[ !muradecset,
                rr = R[[1;;3]];
                v = R[[4;;6]];
                t = ((EPOCH-1950)-50.00021)/100;
                rr1 = rr+sectoradian v t;
                {x,y,z} = rr1[[1;;3]];,
                {x,y,z} = R[[1;;3]];
                {xdot,ydot,zdot} = R[[4;;6]];
            ];
            r2 = x^2+y^2+z^2;
            rmag = Sqrt[r2];
            dec2000[[i]] = ArcSin[z/rmag];
            ra2000[[i]] = ArcTan[x,y];
            If[ muradecset,
                MURADEC[[i,1]] = (x ydot-y xdot)/(x^2+y^2);
                MURADEC[[i,2]] = (zdot (x^2+y^2)-z(x xdot +y ydot))/(r2 Sqrt[x^2+y^2])
            ];
            If[ PARALLAX[[i]]>0,
                RADVEL[[i]] = (x xdot+y ydot +z zdot)/(21.095 PARALLAX[[i]] rmag);
                PARALLAX[[i]] = PARALLAX[[i]]/rmag
            ];
            ,{i,Length@rarad}
            ];
        neg = Flatten[Position[ra2000,a_/;a<0]];
        If[ Length@neg>0,
            ra2000[[neg]] = ra2000[[neg]]+2\[Pi]
        ];
        (* Make output scalar if input was scalar *)
        If[ !ListQ[RAB1950],
            ra2000 = First@ra2000;
            dec2000 = First@dec2000
        ];
        Return[{ra2000,dec2000}/Degree]
    ]



SyntaxInformation[WCSList] = {"ArgumentsPattern"->{{_,_},OptionsPattern[]}};
WCSList[RA_?NumberQ,Dec_?NumberQ,opts:OptionsPattern[]] :=
    WCSList[{RA,Dec},opts];
WCSList[Dec_?NumberQ,opts:OptionsPattern[]] :=
    Module[ {},
        If[ OptionValue[Declination],
            WCSList[{0,Dec},RightAscension->False],
            WCSList[{Dec,0},Declination->False]
        ]
    ]
WCSList[(*xy:{{_,_}..}*){RA_,Dec_},OptionsPattern[]] :=
    Module[ {ra,dec,imin,xsec,ideg,imn,xsc,isra,isdec,ishours,ihr,xmin,outra,outdec},
        isra = OptionValue[RightAscension];
        isdec = OptionValue[Declination];
        ishours = OptionValue[Hours];
        ra = Mod[RA,360];
        dec = Dec;
        (* Compute RA *)
        If[ ishours,
            ra = Mod[ra,24];
            ra = ra+24Min[ra,0];
            ihr = IntegerPart[ra];
            xmin = Abs[ra 60-ihr 60];,
            ra = Mod[ra,360];
            ra = ra+360 Min[ra,0];
            ihr = IntegerPart[ra/15];
            xmin = Abs[ra 4-ihr 60];
        ];
        imin = IntegerPart[xmin];
        xsec = (xmin-imin)60;
        (* Compute Dec *)
        ideg = IntegerPart[dec];
        xmin = Abs[dec-ideg]60;
        imn = IntegerPart[xmin];
        xsc = (xmin-imn)60;
        If[ ideg== 0&&dec<0,
            imn = imn-2 imn;
            If[ imn== 0,
                xsc = xsc-2 xsc
            ]
        ];
        outra = {ihr,imin,xsec};
        outdec = {ideg,imn,xsc};
        If[ isra&&isdec,
            Return[{outra,outdec}]
        ];
        If[ isra&&!isdec,
            Return[outra]
        ];
        If[ !isra&&isdec,
            Return[outdec]
        ]
    ];

(*
XYtoWCS::keyx="FITS header keyword `1` is missing. Using default setting of `2`.";
XYtoWCS[(*xy:{{_,_}..}*){xin_,yin_},imgHead_]:=Module[{CTYPE,CD,CDELT,CRPIX,CROTA,CRVAL,PROJ,LONGPOLE,LATPOLE,PV2,xdiff,ydiff,xsi,eta,coord,temp,a,d,x,y,xFITS,yFITS,DISTORT,xdif1,ydif1,b,naxis1,naxis2},
{CD,CDELT,CRPIX,CRVAL,CTYPE,LONGPOLE,LATPOLE,PV2,DISTORT}=extractAstroPars[imgHead];
(* RA-Dec system and southern hemisphere *)
(* mirroring at the CRPIX *)
{naxis1,naxis2}={"NAXIS1","NAXIS2"}/.imgHead;
naxis1=If[ListQ@naxis1,First@naxis1,naxis1];
naxis2=If[ListQ@naxis2,First@naxis2,naxis2];
If[CRPIX== {},Message[XYtoWCS::keyx,"CRPIX1",(naxis1-1)/2];Message[XYtoWCS::keyx,"CRPIX2",(naxis2-1)/2];CRPIX={(naxis1-1)/2,(naxis2-1)/2}];
If[Length@CRPIX== 1,Message[XYtoWCS::keyx,"CRPIX2",(naxis2-1)/2];CRPIX={CRPIX[[1]],(naxis2-1)/2}];
x=xin+2(CRPIX[[1]]-xin);
y=yin+2(CRPIX[[2]]-yin);
{a,d}=XYtoAST[{x,y},imgHead];
Return[{a,d}]
]


Options[WCSconstRA]={ComputeDeclination->False};
WCSconstRA[RA_,yin_,FITSHeader_,opts:OptionsPattern[]]:=Module[{longitude,latitude,alphap,deltap,phip,cp,sp,sa,ca,sd,cd,r,isreverse,phi1,theta1,l,m,n,b0,b1,b2,CD,CDELT,CRPIX,CTYPE,origin,PV2,CRVAL,LONGPOLE,LATPOLE,cdelta,cdel0,sdel0,delra,cdelra,ctype,sdelra,cdi,yy,delta,a,b,c,d,aa,bb,denom,x,dectrue,DISTORT,naxis1,naxis2,y},

{CD,CDELT,CRPIX,CRVAL,CTYPE,LONGPOLE,LATPOLE,PV2,DISTORT}=extractAstroPars[FITSHeader];
(* mirroring at the CRPIX *)
{naxis1,naxis2}={"NAXIS1","NAXIS2"}/.FITSHeader;
naxis1=If[ListQ@naxis1,First@naxis1,naxis1];
naxis2=If[ListQ@naxis2,First@naxis2,naxis2];
If[CRPIX== {},Message[XYtoWCS::keyx,"CRPIX1",(naxis1-1)/2];Message[XYtoWCS::keyx,"CRPIX2",(naxis2-1)/2];CRPIX={(naxis1-1)/2,(naxis2-1)/2}];
If[Length@CRPIX== 1,Message[XYtoWCS::keyx,"CRPIX2",(naxis2-1)/2];CRPIX={CRPIX[[1]],(naxis2-1)/2}];

y=yin+2(CRPIX[[2]]-yin);
(* ------------------------*)
CRVAL=CRVAL*Degree;
cdelta={{CDELT[[1]],0.},{0.,CDELT[[2]]}};
CD=CD*Degree;
{cdel0,sdel0}={Cos[CRVAL[[2]]],Sin[CRVAL[[2]]]};
delra=RA*Degree-CRVAL[[1]];
{cdelra,sdelra}={Cos[delra],Sin[delra]};
ctype=ToUpperCase@StringTake[CTYPE[[1]],{6,8}];
Which[
ctype== "TAN",
((*cdi=Inverse[cdelta.CD];*)
(* IDL # is the invers matrix multiplication, i.e. A#B multiplies the columns of A with the rows of B. Regular dot product act oppositely. Hence, A#B \[LongRightArrow]B.A in Mathematica*)
cdi=Inverse[CD.cdelta];
yy=y-CRPIX[[2]];
delta=ArcTan[(sdel0 cdelra cdi[[2,2]]-Sin[delra]cdi[[1,2]]+yy cdelra cdel0)/(cdel0 cdi[[2,2]]-yy sdel0)];
),
ctype=="SIN",
(
{a,b}=-CD[[All,1]]CDELT[[1]];
{c,d}=CD[[All,2]]CDELT[[2]];
yy=(y-CRPIX[[2]])(b c-a d);(* new coordinate origin *)
aa=cdel0 d;
bb=sdel0 cdelra d+sdelra b;
denom=Sqrt[aa^2+bb^2];
delta=ArcTan[bb/aa]+ArcSin[yy/denom];
);,
ctype== "CAR",
(
{a,b}=-CD[[All,1]]CDELT[[1]];
{c,d}=CD[[All,2]]CDELT[[2]];
delta=(y-CRPIX[[2]])(b c-a d)+CRVAL[[2]];(* new coordinate origin *)
If[Length@delta== 0&&Length@RA>0,delta=ConstantArray[delta,Length@RA]]
);,
True,
(Message[constRA::projx];Abort[]);
];
delta=delta/Degree;
x=ASTtoXY[{RA,delta},FITSHeader];
(* x contains now the coordinate pair x,y *)
dectrue=OptionValue[ComputeDeclination];
If[dectrue, Return[x,delta],Return[x]]
]

Clear[WCSconstDec];
WCSconstDec::usage="WCSconstDec[Dec,x,FITSHeader_] obtains the X and Y coordinates of a line of constant declination. Returns a set of X pixel values given an image with astrometry, and either\n (1) a set of Y pixel values, and a scalar declination (or longitude),  or\n (2) set of declination values, and a scalar Y value.\n\n In usage (1), constDec can be used to determine the (X,Y) values of a line of constant right ascension.  In usage (2), constDec can used to determine the X positions of specified RA values, along a line of constant Y.\n You can provide the option ComputeDeclination to give out the corresponding set of right ascension values. Default is ComputeRightAscension->False";
WCSconstDec::projx="Program only works for TAN and SIN projections.";
Options[WCSconstDec]={ComputeRightAscension->False};
SyntaxInformation[constDec]={"ArgumentsPattern"->{_,_,_,OptionsPattern[]}};
WCSconstDec[Dec_,xin_,FITSHeader_,opts:OptionsPattern[]]:=Module[{longitude,latitude,alphap,deltap,phip,cp,sp,sa,ca,sd,cd,r,isreverse,phi1,theta1,l,m,n,b0,b1,b2,CD,CDELT,CRPIX,CTYPE,origin,PV2,CRVAL,LONGPOLE,LATPOLE,cdelta,cdel0,sdel0,delra,cdelra,ctype,sdelra,cdi,yy,delta,a,b,c,d,aa,bb,denom,dectrue,xx,sign,alpha,y,DISTORT,x,naxis1,naxis2},
{CD,CDELT,CRPIX,CRVAL,CTYPE,LONGPOLE,LATPOLE,PV2,DISTORT}=extractAstroPars[FITSHeader];
(* mirroring at the CRPIX *)
{naxis1,naxis2}={"NAXIS1","NAXIS2"}/.FITSHeader;
naxis1=If[ListQ@naxis1,First@naxis1,naxis1];
naxis2=If[ListQ@naxis2,First@naxis2,naxis2];
If[CRPIX== {},Message[XYtoWCS::keyx,"CRPIX1",(naxis1-1)/2];Message[XYtoWCS::keyx,"CRPIX2",(naxis2-1)/2];CRPIX={(naxis1-1)/2,(naxis2-1)/2}];
If[Length@CRPIX== 1,Message[XYtoWCS::keyx,"CRPIX2",(naxis2-1)/2];CRPIX={CRPIX[[1]],(naxis2-1)/2}];

x=xin+2(CRPIX[[2]]-xin);
(* ------------------------*)
CRVAL=CRVAL*Degree;
cdelta={{CDELT[[1]],0.},{0.,CDELT[[2]]}};
CD=CD*Degree;
a=-CD[[1,1]]CDELT[[1]];
b=-CD[[2,1]]CDELT[[1]];
c=CD[[1,2]]CDELT[[2]];
d=CD[[2,2]]CDELT[[2]];
(* a=-CD[[1,1]]CDELT[[1]];*)
xx=x-CRPIX[[1]];
sdel0=Sin[CRVAL[[2]]];
cdel0=Cos[CRVAL[[2]]];
ctype=ToUpperCase@StringTake[CTYPE[[1]],{6,8}];
Which[
ctype== "TAN",
(
aa=d;
bb=(b c-d a) xx cdel0+sdel0 b;
sign=If[aa>0,1,-1];
alpha=CRVAL[[1]]+ArcTan[bb/aa]+sign ArcSin[Tan[Dec Degree]((b c-d a)xx sdel0-b cdel0)/(Sqrt[aa^2+bb^2])];
),
ctype=="SIN",
(
aa=d;
bb=b sdel0;
sign=If[aa>0,1,-1];
denom=Cos[Dec Degree]Sqrt[aa^2+bb^2];
alpha=CRVAL[[1]]+ArcTan[bb/aa]+sign ArcSin[((b c-a d)xx-Sin[Dec Degree]cdel0 b)/denom];
);,
ctype== "CAR",
(
alpha=CRVAL[[1]]+(b c-a d)xx;
If[AtomQ[alpha]&&ListQ[Dec],alpha=ConstantArray[alpha,Length@Dec]];
);,
True,
(Message[constDec::projx];Abort[]);
];
alpha=alpha/Degree;
y=ASTtoXY[{alpha,Dec},FITSHeader];
(* y contains now the coordinate pair x,y *)
dectrue=OptionValue[ComputeRightAscension];
If[dectrue, Return[y,delta],Return[y]]
]
*)



FITSListContourOptions = {
    FITSAlignmentPoint->Center,FITSAspectRatio->1,FITSAxes->False,FITSAxesLabel->None,
    FITSAxesOrigin->Automatic,FITSAxesStyle->{},FITSBackground->None,FITSBaselinePosition->Automatic,
    FITSBaseStyle->{},FITSBoundaryStyle->None,FITSBoxRatios->Automatic,FITSClippingStyle->None,
    FITSColorFunction->Automatic,FITSColorFunctionScaling->True,FITSColorOutput->Automatic,
    FITSContentSelectable->Automatic,FITSContourLabels->Automatic,FITSContourLines->True,
    FITSContours->Automatic,FITSContourShading->Automatic,FITSContourStyle->Automatic,
    FITSCoordinatesToolOptions->Automatic,FITSDataRange->Automatic,
    FITSDisplayFunction:> $DisplayFunction,FITSEpilog->{},FITSFormatType:> TraditionalForm,
    FITSFrame->True,FITSFrameLabel->None,FITSFrameStyle->{},FITSFrameTicks->Automatic,
    FITSFrameTicksStyle->{},FITSGridLines->None,FITSGridLinesStyle->{},FITSImageMargins->0.`,
    FITSImagePadding->All,FITSImageSize->Automatic,FITSImageSizeRaw->Automatic,
    FITSInterpolationOrder->None,FITSLabelStyle->{},FITSLightingAngle->None,
    FITSMaxPlotPoints->Automatic,FITSMesh->None,FITSMeshFunctions->{},
    FITSMeshStyle->Automatic,FITSMethod->Automatic,FITSPerformanceGoal:> $PerformanceGoal,
    FITSPlotLabel->None,FITSPlotRange->{Full,Full,Automatic},FITSPlotRangeClipping->True,
    FITSPlotRangePadding->Automatic,FITSPlotRegion->Automatic,FITSPreserveImageOptions->Automatic,
    FITSProlog->{},FITSRegionFunction->(True&),FITSRotateLabel->True,FITSTicks->Automatic,
    FITSTicksStyle->{}};
WCSOptions = {
    WCSAlignmentPoint->Center,WCSAspectRatio->1,WCSAxes->False,WCSAxesLabel->None,
    WCSAxesOrigin->Automatic,WCSAxesStyle->{},WCSBackground->None,WCSBaselinePosition->Automatic,
    WCSBaseStyle->{},WCSBoundaryStyle->None,WCSBoxRatios->Automatic,WCSClippingStyle->None,
    WCSColorFunction->Automatic,WCSColorFunctionScaling->True,WCSColorOutput->Automatic,
    WCSContentSelectable->Automatic,WCSContourLabels->Automatic,WCSContourLabelsRA->Automatic,
    WCSContourLabelsDEC->Automatic,WCSContourLines->True,WCSContoursRA->Automatic,
    WCSContoursDEC->Automatic,WCSContourShading->Automatic,WCSContourStyle->Automatic,
    WCSContourStyleRA->Automatic,WCSContourStyleDEC->Automatic,WCSCoordinatesToolOptions->Automatic,
    WCSDisplayFunction:> $DisplayFunction,WCSEpilog->{},WCSEvaluated->Automatic,
    WCSEvaluationMonitor->None,WCSExclusions->Automatic,WCSExclusionsStyle->None,
    WCSFormatType:> TraditionalForm,WCSFrame->True,WCSFrameLabel->None,WCSFrameStyle->{},
    WCSFrameTicks->Automatic,WCSFrameTicksStyle->{},WCSGridLines->None,
    WCSGridLinesStyle->{},WCSImageMargins->0.`,WCSImagePadding->All,WCSImageSize->Automatic,
    WCSImageSizeRaw->Automatic,WCSLabelStyle->{},WCSLightingAngle->None,
    WCSMaxRecursion->Automatic,WCSMesh->None,WCSMeshFunctions->{},WCSMeshStyle->Automatic,
    WCSMethod->Automatic,WCSPerformanceGoal:> $PerformanceGoal,WCSPlotLabel->None,
    WCSPlotPoints->Automatic,WCSPlotRange->{Full,Full,Automatic},WCSPlotRangeClipping->True,
    WCSPlotRangePadding->Automatic,WCSPlotRegion->Automatic,WCSPreserveImageOptions->Automatic,
    WCSProlog->{},WCSRegionFunction->(True&),WCSRotateLabel->True,WCSTicks->Automatic,
    WCSTicksStyle->{},WCSWorkingPrecision->MachinePrecision
    };
Options[WCSDisplay] = Join[
    {DataReversed->True},
    Options[ArrayPlot],
    Options[ContourPlot],
    WCSOptions,
    FITSListContourOptions,
        {
            WCSCoordinateGrid->True,
            FITSContourPlot->False,
            FITSArrayPlot->True
            }];


WCSDisplay[FITSdata_,FITShead_,opts:OptionsPattern[]] :=
    Quiet[
    Module[ {
    arrayPlotOptions,xrange,yrange,gridtab,arrayPlot,gridplot,rulePlotRange,
    zrange,subarray,xdim1,xdim2,ydim1,ydim2,clippedData,showFITS,showFITSContours,
    showWCS,FITScontours,wcsoptions,fitscontouroptions,RAoptions,DECoptions,astroStructure},
        astroStructure = WCSExtractAstroPars[FITShead];
        arrayPlotOptions = Cases[FilterRules[{opts}, Options[ArrayPlot]],Except[Rule[DataReversed,_]]];
        rulePlotRange = FilterRules[{opts},PlotRange];
        {xdim1,xdim2} = {1,Dimensions[FITSdata][[2]]};
        {ydim1,ydim2} = {1,Dimensions[FITSdata][[1]]};
        If[ rulePlotRange!= {},
            Which[
                MatrixQ[rulePlotRange[[1,2]]]&&Length@rulePlotRange[[1,2]]== 3,
                    {xrange,yrange,zrange} = rulePlotRange[[1,2]],
                ListQ[rulePlotRange[[1,2]]]&&Length@rulePlotRange[[1,2]]== 3,
                    {xrange,yrange,zrange} = rulePlotRange[[1,2]],
                MatrixQ[rulePlotRange[[1,2]]]&&Length@rulePlotRange[[1,2]]== 2,
                    {xrange,yrange} = rulePlotRange[[1,2]];
                    zrange = Automatic;,
                VectorQ[rulePlotRange[[1,2]]]&&Length@rulePlotRange[[1,2]]== 2,
                    {xrange,yrange} = {{xdim1,xdim2},{ydim1,ydim2}};
                    zrange = rulePlotRange[[1,2]]
                ],
            {xrange,yrange} = {{xdim1,xdim2},{ydim1,ydim2}}
        ];
        clippedData = FITSdata;
        If[ !MemberQ[{Automatic,Full},zrange]&&zrange!= {},
            clippedData = Clip[FITSdata,zrange,{Missing[],Missing[]}]
        ];
        (* BUG in ArrayPlot: cannot specify PlotRange and DataRange together.
        If both options are given, DataRange is taken and PlotRange is reversed to full array display.
        Workaround: clip the img array*)
        If[ xrange[[1]]>1,
            xdim1 = xrange[[1]],
            xdim1 = 1
        ];
        If[ xrange[[2]]<Dimensions[FITSdata][[2]],
            xdim2 = xrange[[2]],
            xdim2 = Dimensions[FITSdata][[2]]
        ];
        If[ yrange[[1]]>1,
            ydim1 = yrange[[1]]
        ];
        If[ yrange[[2]]<Dimensions[FITSdata][[1]],
            ydim2 = yrange[[2]]
        ];
        subarray = clippedData[[ydim1;;ydim2,xdim1;;xdim2]];(* ydim is the 1st array dimensions!*)
        wcsoptions = FilterRules[Cases[{opts},Rule[a_/;(StringTake[ToString[a],3]== "WCS"),b_]:>Rule[ToExpression@StringDrop[ToString[a],3],b]],Options[ContourPlot]];
        RAoptions = FilterRules[Cases[{opts},Rule[a_/;(StringTake[ToString[a],3]== "WCS")&&(StringTake[ToString[a],-2]== "RA"),b_]:>Rule[ToExpression@StringDrop[StringDrop[ToString[a],3],-2],b]],Options[ContourPlot]];
        DECoptions = FilterRules[Cases[{opts},Rule[a_/;(StringTake[ToString[a],3]== "WCS")&&(StringTake[ToString[a],-3]== "DEC"),b_]:>Rule[ToExpression@StringDrop[StringDrop[ToString[a],3],-3],b]],Options[ContourPlot]];
        fitscontouroptions = FilterRules[Cases[{opts},Rule[a_/;(StringTake[ToString[a],4]== "FITS"),b_]:>Rule[ToExpression@StringDrop[ToString[a],4],b]],Options[ListContourPlot]];
        showFITSContours = OptionValue[FITSContourPlot];
        showWCS = OptionValue[WCSCoordinateGrid];
        showFITS = OptionValue[FITSArrayPlot];
        If[ showFITS,
            arrayPlot = ArrayPlot[subarray,
                DataRange->{{xdim1,xdim2},{ydim1,ydim2}},
                Cases[Evaluate[arrayPlotOptions],Except[Rule[PlotRange,_]]],DataReversed->True]
        ];
        If[ showFITSContours,
            FITScontours = ListContourPlot[subarray,
                PerformanceGoal->"Speed",
                Method -> {"DelaunayDomainScaling" -> True},
                Evaluate[fitscontouroptions],
                ContourShading->None,
                Method -> {"DelaunayDomainScaling" -> True},
                MaxPlotPoints->(Max[xdim2,ydim2]-Min[xdim1,ydim1]),
                DataRange->{{xdim1,xdim2},{ydim1,ydim2}}]
        ];
        If[ showWCS,
            gridtab = Table[{i,j,WCSXYtoAST[{i,j},FITShead,astroStructure]},{j,yrange[[1]],yrange[[2]],(yrange[[2]]-yrange[[1]])/10},{i,xrange[[1]],xrange[[2]],(xrange[[2]]-xrange[[1]])/10}];
            gridplot = Show[{
                ListContourPlot[gridtab[[All,All,3,2]],
                    Evaluate[Join[wcsoptions,DECoptions]],
                    ContourShading->None,
                    Method -> {"DelaunayDomainScaling" -> True},
                    DataRange->{xrange,yrange}],
                ListContourPlot[gridtab[[All,All,3,1]],
                    Evaluate[Join[wcsoptions,RAoptions]],
                    ContourShading->None,
                    Method -> {"DelaunayDomainScaling" -> True},
                    DataRange->{xrange,yrange}]}]
        ];
        
        Show[
            Which[
                showFITS&&showWCS&&showFITSContours,
                    {arrayPlot,FITScontours,gridplot},
                showFITS&&showWCS&&!showFITSContours,
                    {arrayPlot,gridplot},
                !showFITS&&showWCS&&showFITSContours,
                    {FITScontours,gridplot},
                showFITS&&!showWCS&&showFITSContours,
                    {arrayPlot,FITScontours},
                showFITS&&!showWCS&&!showFITSContours,
                    {arrayPlot},
                !showFITS&&!showWCS&&showFITSContours,
                    {FITScontours},
                !showFITS&&showWCS&&!showFITSContours,
                    {gridplot},
                True,
                    {arrayPlot}],
                Frame->True,
                FrameTicks->Automatic]
    (*gridtab=Table[{i,j,XYtoAST[{i,j},FITShead]},{j,-400,881,40},{i,-400,881,40}];*)],
    OptionValue::"nodef"];



End[] 
(* End Private Context *)

EndPackage[]
