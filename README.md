## Fortran codes for approach of residual load and synthesis of residual load effects using SRBFs
https://www.zcyphygeodesy.com/en/h-nd-136.html
## [Algorithm purpose]
    From the regional residual equivalent water height (EWH) variation grid (cm), approach the regional residual surface loads using spherical radial basis functions (SRBFs) and then calculate the residual EWH estimation (cm) and residual load effects on the geoid or height anomaly (mm), ground gravity (μGal), gravity disturbance (μGal), ground tilt (SW, to the south and to the west, mas), vertical deflection (SW, to the south and to the west, mas), horizontal displacement (EN, to the east and to the north, mm), ground radial displacement (mm), ground normal or orthometric height (mm), radial gravity gradient (mE) or horizontal gravity gradient (NW, to the north and to the west, E) using SRBF synthesis.
    When computing the load effects of sea level variations, the height of the calculation point is the normal or orthometric height. When computing the load effects of land water variations, the height of the calculation point is the height relative to the Earth’s surface.
## [Geophysical models]
    The Earth’s Load Love number file love_load_cm.dat from a Regional EIAstic Rebound calculator (REAR1.0, 2015).
## [Main program for test entrance]
    SurfaceEwhSRBFloadeffect.f90
    The record format of the input calculation point file calcpntfl: ID (point no / point name), longitude (decimal degrees), latitude (decimal degrees), height (m) relative to the Earth’s surface …
    Input the regional residual equivalent water height variation grid file ewhgridfl and Earth’s Load Love number file loadlovfl.
    Input parameters: para(1) - the cumulative SRBF approach times. para(1) - method of the solution of normal equation (=1 LU triangular decomposition method, =2 Cholesky decomposition, =3 least square QR decomposition, =4 Minimum norm singular value decomposition).
    Input main SRBF parameters: para(3:9) - the spherical radial basis functions (=0 radial multipole kernel function, =1 Poisson wavelet kernel function), the order number m, minimum and maximum degree of SRBF Legendre expansion, Bjerhammar sphere burial depth (km), action distance (km) of SRBF center and Reuter network level K.
    Input cumulative approach SRBF parameters: para(10:16).
## (1) Module for approach of residual load and synthesis of load effects using SRBFs
    EwhSRBFLoadeffectpnt(calcpntfl,ewhgridfl,loadlovfl,para)
    Output the residual EWH estimation and residual load effect file reslt.txt.
    The output file header is the same as the input file. Behind the input file record, adds the residual EWH estimation (cm) and residual load effects on the geoid or height anomaly (mm), ground gravity (μGal), gravity disturbance (μGal), ground tilt (SW, to the south and to the west, mas), vertical deflection (SW, to the south and to the west, mas), horizontal displacement (EN, to the east and to the north, mm), ground radial displacement (mm), ground normal or orthometric height (mm), radial gravity gradient (mE) or horizontal gravity gradient (NW, to the north and to the west, E).
    Input parameters: calcpntfl - The space calculation point file name.
    Input parameters: ewhgridfl - The regional residual equivalent water height variation grid file name.
    Input parameters: loadlovfl - The Earth’s Load Love number file name.
    The module also outputs the SRBF spatial curve file SRBFspc.rbf and spectral curve file SRBFdgr.rbf of 11 kinds of geodetic variations into the current directory.
    SRBFspc.rbf file header format: SRBF type (0-radial multipole kernel function, 1-Poisson wavelet kernel function), order of SRBF, Minimum and maximum degree of SRBF Legendre expansion, buried depth (km). The record format: spherical distance (km), the normalized SRBF values from the load EWH, height anomaly, ground gravity, gravity disturbance, ground tilt, vertical deflection, horizontal displacement, radial displacement, orthometric height, radial gravity gradient and horizontal gradient variations.
    The file header of SRBFdgr.rbf is the same as SRBFspc.rbf. The record format: degree n of SRBF Legendre expansion, the degree-n normalized SRBF values from the load EWH, height anomaly, ground gravity, gravity disturbance, ground tilt, vertical deflection, horizontal displacement, radial displacement, orthometric height, radial gravity gradient and horizontal gradient variations.
## (2) Computation module for the Reuter network parameters
    ReuterGrid(rhd,lvl,Kt,blat,nn,mm,nln,sr,dl,nrd,lon)
    Input parameters: rhd(4) - minimum and maximum longitude, minimum and maximum geocentric latitude of the Reuter network.
    Input parameters: lvl, nn, mm - the Reuter network level, maximum number of Reuter centers in the meridian direction and that in the parallel direction.
    Return parameter: Blat - the geocentric latitude (degree decimal) of Reuter centers in the first parallel direction.
    Return parameters: Kt - the number of Reuter centers, equal to the number of unknowns to be estimated.
    Return parameters: nln(nn) - the number of the Reuter centers in the parallel direction.
    Return parameters: sr(nn) - the percentage of the difference between the area of the Reuter cell-grid in the parallel direction and the area of the equatorial cell-grid.
    Return parameters: dl(nn) - the longitude interval (degree decimal) of the Reuter centers in the parallel direction.
    Return parameters: nrd(nn,mm) - ordinal number value of the Reuter centers.
    Return parameters: lon(nn,mm) - longitude value of the Reuter centers.
## (3) Module for the best match between the Reuter centers and observation points
    Edgnode(enode,rlatlon,lvl,edgn,lon,blat,nln,gpnt,nn,mm)
    The module calculates the number of observation points in the Reuter cell-grid and the number gpnt of Reuter centers corrected.
    Return parameter: edgn - the number of Reuter centers around the edge of the Reuter grid.
    Return parameters: enode(edgn) - the ordinal number value of Reuter centers in the edge of the Reuter grid.
    Return parameter: rlatlon(edgn,2) - the geocentric latitude and longitude (degree decimal) of Reuter centers in the edge of the Reuter grid.
## (4) Computation module for the SRBF curves of all 11 elements
    SRBF11all(RBF,flv,order,krbf,mpn,mdp,mp2,minN,maxN,NF,nta)
    Return parameters: RBF(NF+1,11) - The SRBF curves of all 11 elements, which is calculated by the action distance of SRBF center and Reuter grid level.
    Where RBF(NF+1,knd): knd=1 EWH, =2 height anomaly, =3 ground gravity, =4 gravity disturbance, =5 ground tilt, =6 vertical deflection, =7 horizontal displacement, =8 ground radial displacement, =9 ground normal or orthometric height, =10 radial gravity gradient, =11 horizontal gravity gradient.
    Input parameters: flv(:,3) - Load love numbers.
    Input parameters: nta - the bandwidth parameter and nta = (r0-dpth)/r0, here dpth is the Bjerhammar sphere burial depth and r0 is the average geocentric distance of observation points.
    Input parameters: krbf,order - krbf=0 radial multipole kernel function, =1 Poisson wavelet kernel function and order is the order number of SRBF.
    Input parameters: mpn(maxN-minN+1, NF+1), mdp(maxN-minN+1, NF+1), mp2(maxN-minN+1, NF+1) - all minN to maxN degree Legendre functions and their first and second derivatives.
## (5) Calculation module for the position of the calculation point in the Reuter grid
    RtGridij(rln,ki,kj,blat,lvl,nn,mm,nln,dl,lon)
    Input parameters: rln(3) - the spherical coordinates of the calculation point。
    Return parameters: ki,kj - the position of the calculation point rln(3) in the Reuter grid, It is represented by the element of the 2-D ordinal number array of the Reuter grid and ki>0, kj>0.
## (6) Solution module of large positive definite symmetric equations
    Equsolve(BB,xx,nn,BL,knd,bf)
    The mkl_lapack95_ilp64.lib library is called to solve the large equations BB.xx = BL. bf(8) is the property of the solution.
    Input parameter: knd - method of the solution of normal equation and  knd =1 LU triangular decomposition method, =2 Cholesky decomposition, =3 least square QR decomposition, =4 Minimum norm singular value decomposition.
## (7) Calculation module for the normal gravity field
    normdjn(GRS,djn); GNormalfd(BLH,NFD,GRS)
    Return parameters: NFD(5) - the normal geopotential (m2/s2), normal gravity (mGal), normal gravity gradient (E), normal gravity line direction (', expressed by its north declination relative to the center of the Earth center of mass) or normal gravity gradient direction (', expressed by its north declination relative to the Earth center of mass).
## (8) Calculation module for Legendre functions and their derivatives to ψ
    LegPn_dt2(pn,dp1,dp2,n,t) ! t=cos ψ
## (9) Algorithm library for transforming of geodetic coordinates
    BLH_RLAT(GRS, BLH, RLAT); BLH_XYZ(GRS, BLH, XYZ)
    RLAT_BLH(GRS, RLAT, BLH)
## (10) Other auxiliary modules
    LegPn02(mpn,mdp,mp2,minN,maxN,NF,dr); PickRecord(str0, kln, rec, nn)
    RBFvalue(RBF(:,1),NF,dr,dln(2),tmp); drln(rln,rlnk,dln); Stat1d(dt,nn,rst)
## [For compile and link]
    Fortran90, 132 Columns fixed format. Fortran compiler. mkl_lapack95_ilp64.lib link library required.
## [Algorithmic formula] PAGravf4.5 User Reference https://www.zcyphygeodesy.com/en/
    7.2.3 The geodetic numerical grid file
    8.7 Load deformation field approach from heterogeneous variations using SRBFs
The zip compression package includes the test project in visual studio 2017 - intel fortran integrated environment, DOS executable test file and all input and output data.
