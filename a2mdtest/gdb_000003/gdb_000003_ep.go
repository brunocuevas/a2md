 Entering Gaussian System, Link 0=g09
 Initial command:
 /usr/local/g09/l1.exe /mnt/f/scratch/Gau-172.inp -scrdir=/mnt/f/scratch/
 Entering Link 1 = /usr/local/g09/l1.exe PID=       173.
  
 Copyright (c) 1988,1990,1992,1993,1995,1998,2003,2009, Gaussian, Inc.
                  All Rights Reserved.
  
 This is part of the Gaussian(R) 09 program.  It is based on
 the Gaussian(R) 03 system (copyright 2003, Gaussian, Inc.),
 the Gaussian(R) 98 system (copyright 1998, Gaussian, Inc.),
 the Gaussian(R) 94 system (copyright 1995, Gaussian, Inc.),
 the Gaussian 92(TM) system (copyright 1992, Gaussian, Inc.),
 the Gaussian 90(TM) system (copyright 1990, Gaussian, Inc.),
 the Gaussian 88(TM) system (copyright 1988, Gaussian, Inc.),
 the Gaussian 86(TM) system (copyright 1986, Carnegie Mellon
 University), and the Gaussian 82(TM) system (copyright 1983,
 Carnegie Mellon University). Gaussian is a federally registered
 trademark of Gaussian, Inc.
  
 This software contains proprietary and confidential information,
 including trade secrets, belonging to Gaussian, Inc.
  
 This software is provided under written license and may be
 used, copied, transmitted, or stored only in accord with that
 written license.
  
 The following legend is applicable only to US Government
 contracts under FAR:
  
                    RESTRICTED RIGHTS LEGEND
  
 Use, reproduction and disclosure by the US Government is
 subject to restrictions as set forth in subparagraphs (a)
 and (c) of the Commercial Computer Software - Restricted
 Rights clause in FAR 52.227-19.
  
 Gaussian, Inc.
 340 Quinnipiac St., Bldg. 40, Wallingford CT 06492
  
  
 ---------------------------------------------------------------
 Warning -- This program may not be used in any manner that
 competes with the business of Gaussian, Inc. or will provide
 assistance to any competitor of Gaussian, Inc.  The licensee
 of this program is prohibited from giving any competitor of
 Gaussian, Inc. access to this program.  By using this program,
 the user acknowledges that Gaussian, Inc. is engaged in the
 business of creating and licensing software in the field of
 computational chemistry and represents and warrants to the
 licensee that it is not a competitor of Gaussian, Inc. and that
 it will not use this program in any manner prohibited above.
 ---------------------------------------------------------------
  

 Cite this work as:
 Gaussian 09, Revision A.02,
 M. J. Frisch, G. W. Trucks, H. B. Schlegel, G. E. Scuseria, 
 M. A. Robb, J. R. Cheeseman, G. Scalmani, V. Barone, B. Mennucci, 
 G. A. Petersson, H. Nakatsuji, M. Caricato, X. Li, H. P. Hratchian, 
 A. F. Izmaylov, J. Bloino, G. Zheng, J. L. Sonnenberg, M. Hada, 
 M. Ehara, K. Toyota, R. Fukuda, J. Hasegawa, M. Ishida, T. Nakajima, 
 Y. Honda, O. Kitao, H. Nakai, T. Vreven, J. A. Montgomery, Jr., 
 J. E. Peralta, F. Ogliaro, M. Bearpark, J. J. Heyd, E. Brothers, 
 K. N. Kudin, V. N. Staroverov, R. Kobayashi, J. Normand, 
 K. Raghavachari, A. Rendell, J. C. Burant, S. S. Iyengar, J. Tomasi, 
 M. Cossi, N. Rega, J. M. Millam, M. Klene, J. E. Knox, J. B. Cross, 
 V. Bakken, C. Adamo, J. Jaramillo, R. Gomperts, R. E. Stratmann, 
 O. Yazyev, A. J. Austin, R. Cammi, C. Pomelli, J. W. Ochterski, 
 R. L. Martin, K. Morokuma, V. G. Zakrzewski, G. A. Voth, 
 P. Salvador, J. J. Dannenberg, S. Dapprich, A. D. Daniels, 
 O. Farkas, J. B. Foresman, J. V. Ortiz, J. Cioslowski, 
 and D. J. Fox, Gaussian, Inc., Wallingford CT, 2009.
 
 ******************************************
 Gaussian 09:  AM64L-G09RevA.02 11-Jun-2009
                30-May-2019 
 ******************************************
 --------------------------------------------------------
 #mp2/6-311++G(d,p) density=current prop=(read,potential)
 --------------------------------------------------------
 1/30=1,38=1/1;
 2/12=2,17=6,18=5,40=1/2;
 3/5=4,6=6,7=1111,11=9,16=1,25=1,30=1,71=1/1,2,3;
 4//1;
 5/5=2,38=5/2;
 8/6=4,10=2/1;
 9/15=1,16=-3/6;
 10/5=1,13=10/2;
 6/7=2,8=2,9=2,10=2,14=3,15=1,22=-1/1,2;
 99/5=1,9=1/99;
 -------------
 water pruebas
 -------------
 Symbolic Z-matrix:
 Charge =  0 Multiplicity = 1
 O                     0.        0.11086   0. 
 H                     0.7979   -0.44344   0. 
 H                    -0.7979   -0.44344   0. 
 
                          Input orientation:                          
 ---------------------------------------------------------------------
 Center     Atomic      Atomic             Coordinates (Angstroms)
 Number     Number       Type             X           Y           Z
 ---------------------------------------------------------------------
      1          8           0        0.000000    0.110861    0.000000
      2          1           0        0.797899   -0.443444    0.000000
      3          1           0       -0.797899   -0.443443    0.000000
 ---------------------------------------------------------------------
                    Distance matrix (angstroms):
                    1          2          3
     1  O    0.000000
     2  H    0.971544   0.000000
     3  H    0.971543   1.595798   0.000000
 Stoichiometry    H2O
 Framework group  C2V[C2(O),SGV(H2)]
 Deg. of freedom     2
 Full point group                 C2V     NOp   4
 Largest Abelian subgroup         C2V     NOp   4
 Largest concise Abelian subgroup C2      NOp   2
                         Standard orientation:                         
 ---------------------------------------------------------------------
 Center     Atomic      Atomic             Coordinates (Angstroms)
 Number     Number       Type             X           Y           Z
 ---------------------------------------------------------------------
      1          8           0        0.000000    0.000000    0.110861
      2          1           0        0.000000    0.797899   -0.443443
      3          1           0        0.000000   -0.797899   -0.443443
 ---------------------------------------------------------------------
 Rotational constants (GHZ):    918.8638053    393.8272137    275.6730769
 Standard basis: 6-311++G(d,p) (5D, 7F)
 There are    17 symmetry adapted basis functions of A1  symmetry.
 There are     2 symmetry adapted basis functions of A2  symmetry.
 There are     6 symmetry adapted basis functions of B1  symmetry.
 There are    11 symmetry adapted basis functions of B2  symmetry.
 Integral buffers will be    131072 words long.
 Raffenetti 1 integral format.
 Two-electron integral symmetry is turned on.
    36 basis functions,    54 primitive gaussians,    37 cartesian basis functions
     5 alpha electrons        5 beta electrons
       nuclear repulsion energy         9.0464358024 Hartrees.
 NAtoms=    3 NActive=    3 NUniq=    2 SFac= 1.69D+00 NAtFMM=   80 NAOKFM=F Big=F
 One-electron integrals computed using PRISM.
 NBasis=    36 RedAO= T  NBF=    17     2     6    11
 NBsUse=    36 1.00D-06 NBFU=    17     2     6    11
 Harris functional with IExCor=  205 diagonalized for initial guess.
 ExpMin= 3.60D-02 ExpMax= 8.59D+03 ExpMxC= 1.30D+03 IAcc=3 IRadAn=         0 AccDes= 0.00D+00
 HarFok:  IExCor=  205 AccDes= 0.00D+00 IRadAn=         0 IDoV= 1
 ScaDFX=  1.000000  1.000000  1.000000  1.000000
 FoFCou: FMM=F IPFlag=           0 FMFlag=      100000 FMFlg1=           0
         NFxFlg=           0 DoJE=T BraDBF=F KetDBF=T FulRan=T
         Omega=  0.000000  0.000000  1.000000  0.000000  0.000000 ICntrl=     500 IOpCl=  0
         NMat0=    1 NMatS0=    1 NMatT0=    0 NMatD0=    1 NMtDS0=    0 NMtDT0=    0
         I1Cent=           4 NGrid=           0.
 Petite list used in FoFCou.
 Initial guess orbital symmetries:
       Occupied  (A1) (A1) (B2) (A1) (B1)
       Virtual   (A1) (B2) (A1) (B2) (B1) (A1) (B2) (A1) (B2) (A1)
                 (A1) (B1) (B2) (A2) (A1) (A1) (B1) (B2) (B2) (A1)
                 (A1) (B2) (B1) (A2) (A1) (A1) (B2) (B1) (A1) (B2)
                 (A1)
 The electronic state of the initial guess is 1-A1.
 Requested convergence on RMS density matrix=1.00D-08 within 128 cycles.
 Requested convergence on MAX density matrix=1.00D-06.
 Requested convergence on             energy=1.00D-06.
 No special actions if energy rises.
 Keep R1 ints in memory in canonical form, NReq=1090309.
 SCF Done:  E(RHF) =  -76.0508474560     A.U. after   11 cycles
             Convg  =    0.3345D-08             -V/T =  2.0017
 ExpMin= 3.60D-02 ExpMax= 8.59D+03 ExpMxC= 1.30D+03 IAcc=3 IRadAn=         5 AccDes= 0.00D+00
 HarFok:  IExCor=  205 AccDes= 0.00D+00 IRadAn=         5 IDoV=-2
 ScaDFX=  1.000000  1.000000  1.000000  1.000000
 Range of M.O.s used for correlation:     2    36
 NBasis=    36 NAE=     5 NBE=     5 NFC=     1 NFV=     0
 NROrb=     35 NOA=     4 NOB=     4 NVA=    31 NVB=    31
 Fully in-core method, ICMem=     6710399.
 JobTyp=1 Pass  1 fully in-core, NPsUse=  1.
 Spin components of T(2) and E(2):
     alpha-alpha T2 =       0.7101864464D-02 E2=     -0.2810777840D-01
     alpha-beta  T2 =       0.4171906820D-01 E2=     -0.1663294894D+00
     beta-beta   T2 =       0.7101864464D-02 E2=     -0.2810777840D-01
 ANorm=    0.1027581042D+01
 E2 =    -0.2225450462D+00 EUMP2 =    -0.76273392502226D+02
          Differentiating once with respect to electric field.
                with respect to dipole field.
          Electric field/nuclear overlap derivatives assumed to be zero.
          Keep R1 ints in memory in canonical form, NReq=1068080.
          There are     1 degrees of freedom in the 1st order CPHF.  IDoFFX=0.
 LinEq1:  Iter=  0 NonCon=     1 RMS=7.91D-03 Max=6.54D-02
 AX will form     1 AO Fock derivatives at one time.
 LinEq1:  Iter=  1 NonCon=     1 RMS=8.38D-04 Max=4.61D-03
 LinEq1:  Iter=  2 NonCon=     1 RMS=1.17D-04 Max=5.31D-04
 LinEq1:  Iter=  3 NonCon=     1 RMS=4.13D-05 Max=3.07D-04
 LinEq1:  Iter=  4 NonCon=     1 RMS=8.39D-06 Max=3.77D-05
 LinEq1:  Iter=  5 NonCon=     1 RMS=1.46D-06 Max=6.75D-06
 LinEq1:  Iter=  6 NonCon=     1 RMS=2.78D-07 Max=1.27D-06
 LinEq1:  Iter=  7 NonCon=     1 RMS=3.47D-08 Max=1.63D-07
 LinEq1:  Iter=  8 NonCon=     1 RMS=4.16D-09 Max=2.92D-08
 LinEq1:  Iter=  9 NonCon=     1 RMS=4.54D-10 Max=1.90D-09
 LinEq1:  Iter= 10 NonCon=     0 RMS=3.57D-11 Max=1.51D-10
 Linear equations converged to 1.000D-10 1.000D-09 after    10 iterations.
 End of Minotr Frequency-dependent properties file   721 does not exist.
 End of Minotr Frequency-dependent properties file   722 does not exist.

 **********************************************************************

            Population analysis using the MP2 density.

 **********************************************************************

 Orbital symmetries:
       Occupied  (A1) (A1) (B2) (A1) (B1)
       Virtual   (A1) (B2) (A1) (A1) (B1) (B2) (B2) (A1) (B2) (A1)
                 (A1) (B2) (B1) (A2) (A1) (A1) (B1) (B2) (B2) (A1)
                 (A1) (B2) (B1) (A2) (A1) (A1) (B2) (B1) (A1) (B2)
                 (A1)
 The electronic state is 1-A1.
 Alpha  occ. eigenvalues --  -20.56633  -1.34417  -0.72367  -0.57003  -0.50712
 Alpha virt. eigenvalues --    0.04335   0.07162   0.23989   0.24487   0.24596
 Alpha virt. eigenvalues --    0.26036   0.31729   0.32267   0.68864   0.70978
 Alpha virt. eigenvalues --    1.20336   1.24099   1.24513   1.48421   1.52352
 Alpha virt. eigenvalues --    1.64712   1.74117   1.87619   2.31342   2.41698
 Alpha virt. eigenvalues --    2.74307   2.86139   3.44189   3.54301   3.66968
 Alpha virt. eigenvalues --    4.12055   4.28090   5.49006   5.82051   6.19804
 Alpha virt. eigenvalues --   51.55991
          Condensed to atoms (all electrons):
              1          2          3
     1  O    7.913128   0.307659   0.307659
     2  H    0.307659   0.452234  -0.024116
     3  H    0.307659  -0.024116   0.452234
 Mulliken atomic charges:
              1
     1  O   -0.528446
     2  H    0.264223
     3  H    0.264223
 Sum of Mulliken atomic charges =   0.00000
 Mulliken charges with hydrogens summed into heavy atoms:
              1
     1  O    0.000000
 Sum of Mulliken charges with hydrogens summed into heavy atoms =   0.00000
 Electronic spatial extent (au):  <R**2>=             20.2825
 Charge=              0.0000 electrons
 Dipole moment (field-independent basis, Debye):
    X=              0.0000    Y=              0.0000    Z=             -2.0963  Tot=              2.0963
 Quadrupole moment (field-independent basis, Debye-Ang):
   XX=             -7.8903   YY=             -4.1837   ZZ=             -6.7295
   XY=              0.0000   XZ=              0.0000   YZ=              0.0000
 Traceless Quadrupole moment (field-independent basis, Debye-Ang):
   XX=             -1.6225   YY=              2.0842   ZZ=             -0.4617
   XY=              0.0000   XZ=              0.0000   YZ=              0.0000
 Octapole moment (field-independent basis, Debye-Ang**2):
  XXX=              0.0000  YYY=              0.0000  ZZZ=             -1.4230  XYY=              0.0000
  XXY=              0.0000  XXZ=             -0.4286  XZZ=              0.0000  YZZ=              0.0000
  YYZ=             -1.4064  XYZ=              0.0000
 Hexadecapole moment (field-independent basis, Debye-Ang**3):
 XXXX=             -7.7725 YYYY=             -6.4164 ZZZZ=             -8.1553 XXXY=              0.0000
 XXXZ=              0.0000 YYYX=              0.0000 YYYZ=              0.0000 ZZZX=              0.0000
 ZZZY=              0.0000 XXYY=             -2.7295 XXZZ=             -2.6994 YYZZ=             -2.1733
 XXYZ=              0.0000 YYXZ=              0.0000 ZZXY=              0.0000
 N-N= 9.046435802398D+00 E-N=-1.987137020826D+02  KE= 7.618265313738D+01
 Symmetry A1   KE= 6.803258833210D+01
 Symmetry A2   KE= 1.776214917528D-02
 Symmetry B1   KE= 4.480413887993D+00
 Symmetry B2   KE= 3.651888768107D+00

 **********************************************************************

          Electrostatic Properties Using The MP2 Density

 **********************************************************************

       Atomic Center    1 is at   0.000000  0.000000  0.110861
       Atomic Center    2 is at   0.000000  0.797899 -0.443443
       Atomic Center    3 is at   0.000000 -0.797899 -0.443443
      Read-in Center    4 is at   0.000000 -2.645900 -2.645900
      Read-in Center    5 is at   0.000000 -2.593000 -2.645900
      Read-in Center    6 is at   0.000000 -2.540000 -2.645900
      Read-in Center    7 is at   0.000000 -2.487100 -2.645900
      Read-in Center    8 is at   0.000000 -2.434200 -2.645900
      Read-in Center    9 is at   0.000000 -2.381300 -2.645900
      Read-in Center   10 is at   0.000000 -2.328400 -2.645900
      Read-in Center   11 is at   0.000000 -2.275500 -2.645900
      Read-in Center   12 is at   0.000000 -2.222500 -2.645900
      Read-in Center   13 is at   0.000000 -2.169600 -2.645900
      Read-in Center   14 is at   0.000000 -2.116700 -2.645900
      Read-in Center   15 is at   0.000000 -2.063800 -2.645900
      Read-in Center   16 is at   0.000000 -2.010900 -2.645900
      Read-in Center   17 is at   0.000000 -1.958000 -2.645900
      Read-in Center   18 is at   0.000000 -1.905000 -2.645900
      Read-in Center   19 is at   0.000000 -1.852100 -2.645900
      Read-in Center   20 is at   0.000000 -1.799200 -2.645900
      Read-in Center   21 is at   0.000000 -1.746300 -2.645900
      Read-in Center   22 is at   0.000000 -1.693400 -2.645900
      Read-in Center   23 is at   0.000000 -1.640400 -2.645900
      Read-in Center   24 is at   0.000000 -1.587500 -2.645900
      Read-in Center   25 is at   0.000000 -1.534600 -2.645900
      Read-in Center   26 is at   0.000000 -1.481700 -2.645900
      Read-in Center   27 is at   0.000000 -1.428800 -2.645900
      Read-in Center   28 is at   0.000000 -1.375900 -2.645900
      Read-in Center   29 is at   0.000000 -1.322900 -2.645900
      Read-in Center   30 is at   0.000000 -1.270000 -2.645900
      Read-in Center   31 is at   0.000000 -1.217100 -2.645900
      Read-in Center   32 is at   0.000000 -1.164200 -2.645900
      Read-in Center   33 is at   0.000000 -1.111300 -2.645900
      Read-in Center   34 is at   0.000000 -1.058400 -2.645900
      Read-in Center   35 is at   0.000000 -1.005400 -2.645900
      Read-in Center   36 is at   0.000000 -0.952500 -2.645900
      Read-in Center   37 is at   0.000000 -0.899600 -2.645900
      Read-in Center   38 is at   0.000000 -0.846700 -2.645900
      Read-in Center   39 is at   0.000000 -0.793800 -2.645900
      Read-in Center   40 is at   0.000000 -0.740800 -2.645900
      Read-in Center   41 is at   0.000000 -0.687900 -2.645900
      Read-in Center   42 is at   0.000000 -0.635000 -2.645900
      Read-in Center   43 is at   0.000000 -0.582100 -2.645900
      Read-in Center   44 is at   0.000000 -0.529200 -2.645900
      Read-in Center   45 is at   0.000000 -0.476300 -2.645900
      Read-in Center   46 is at   0.000000 -0.423300 -2.645900
      Read-in Center   47 is at   0.000000 -0.370400 -2.645900
      Read-in Center   48 is at   0.000000 -0.317500 -2.645900
      Read-in Center   49 is at   0.000000 -0.264600 -2.645900
      Read-in Center   50 is at   0.000000 -0.211700 -2.645900
      Read-in Center   51 is at   0.000000 -0.158800 -2.645900
      Read-in Center   52 is at   0.000000 -0.105800 -2.645900
      Read-in Center   53 is at   0.000000 -0.052900 -2.645900
      Read-in Center   54 is at   0.000000  0.000000 -2.645900
      Read-in Center   55 is at   0.000000  0.052900 -2.645900
      Read-in Center   56 is at   0.000000  0.105800 -2.645900
      Read-in Center   57 is at   0.000000  0.158800 -2.645900
      Read-in Center   58 is at   0.000000  0.211700 -2.645900
      Read-in Center   59 is at   0.000000  0.264600 -2.645900
      Read-in Center   60 is at   0.000000  0.317500 -2.645900
      Read-in Center   61 is at   0.000000  0.370400 -2.645900
      Read-in Center   62 is at   0.000000  0.423300 -2.645900
      Read-in Center   63 is at   0.000000  0.476300 -2.645900
      Read-in Center   64 is at   0.000000  0.529200 -2.645900
      Read-in Center   65 is at   0.000000  0.582100 -2.645900
      Read-in Center   66 is at   0.000000  0.635000 -2.645900
      Read-in Center   67 is at   0.000000  0.687900 -2.645900
      Read-in Center   68 is at   0.000000  0.740800 -2.645900
      Read-in Center   69 is at   0.000000  0.793800 -2.645900
      Read-in Center   70 is at   0.000000  0.846700 -2.645900
      Read-in Center   71 is at   0.000000  0.899600 -2.645900
      Read-in Center   72 is at   0.000000  0.952500 -2.645900
      Read-in Center   73 is at   0.000000  1.005400 -2.645900
      Read-in Center   74 is at   0.000000  1.058400 -2.645900
      Read-in Center   75 is at   0.000000  1.111300 -2.645900
      Read-in Center   76 is at   0.000000  1.164200 -2.645900
      Read-in Center   77 is at   0.000000  1.217100 -2.645900
      Read-in Center   78 is at   0.000000  1.270000 -2.645900
      Read-in Center   79 is at   0.000000  1.322900 -2.645900
      Read-in Center   80 is at   0.000000  1.375900 -2.645900
      Read-in Center   81 is at   0.000000  1.428800 -2.645900
      Read-in Center   82 is at   0.000000  1.481700 -2.645900
      Read-in Center   83 is at   0.000000  1.534600 -2.645900
      Read-in Center   84 is at   0.000000  1.587500 -2.645900
      Read-in Center   85 is at   0.000000  1.640400 -2.645900
      Read-in Center   86 is at   0.000000  1.693400 -2.645900
      Read-in Center   87 is at   0.000000  1.746300 -2.645900
      Read-in Center   88 is at   0.000000  1.799200 -2.645900
      Read-in Center   89 is at   0.000000  1.852100 -2.645900
      Read-in Center   90 is at   0.000000  1.905000 -2.645900
      Read-in Center   91 is at   0.000000  1.958000 -2.645900
      Read-in Center   92 is at   0.000000  2.010900 -2.645900
      Read-in Center   93 is at   0.000000  2.063800 -2.645900
      Read-in Center   94 is at   0.000000  2.116700 -2.645900
      Read-in Center   95 is at   0.000000  2.169600 -2.645900
      Read-in Center   96 is at   0.000000  2.222500 -2.645900
      Read-in Center   97 is at   0.000000  2.275500 -2.645900
      Read-in Center   98 is at   0.000000  2.328400 -2.645900
      Read-in Center   99 is at   0.000000  2.381300 -2.645900
      Read-in Center  100 is at   0.000000  2.434200 -2.645900
      Read-in Center  101 is at   0.000000  2.487100 -2.645900
      Read-in Center  102 is at   0.000000  2.540000 -2.645900
      Read-in Center  103 is at   0.000000  2.593000 -2.645900
      Read-in Center  104 is at   0.000000 -2.645900 -2.593000
      Read-in Center  105 is at   0.000000 -2.593000 -2.593000
      Read-in Center  106 is at   0.000000 -2.540000 -2.593000
      Read-in Center  107 is at   0.000000 -2.487100 -2.593000
      Read-in Center  108 is at   0.000000 -2.434200 -2.593000
      Read-in Center  109 is at   0.000000 -2.381300 -2.593000
      Read-in Center  110 is at   0.000000 -2.328400 -2.593000
      Read-in Center  111 is at   0.000000 -2.275500 -2.593000
      Read-in Center  112 is at   0.000000 -2.222500 -2.593000
      Read-in Center  113 is at   0.000000 -2.169600 -2.593000
      Read-in Center  114 is at   0.000000 -2.116700 -2.593000
      Read-in Center  115 is at   0.000000 -2.063800 -2.593000
      Read-in Center  116 is at   0.000000 -2.010900 -2.593000
      Read-in Center  117 is at   0.000000 -1.958000 -2.593000
      Read-in Center  118 is at   0.000000 -1.905000 -2.593000
      Read-in Center  119 is at   0.000000 -1.852100 -2.593000
      Read-in Center  120 is at   0.000000 -1.799200 -2.593000
      Read-in Center  121 is at   0.000000 -1.746300 -2.593000
      Read-in Center  122 is at   0.000000 -1.693400 -2.593000
      Read-in Center  123 is at   0.000000 -1.640400 -2.593000
      Read-in Center  124 is at   0.000000 -1.587500 -2.593000
      Read-in Center  125 is at   0.000000 -1.534600 -2.593000
      Read-in Center  126 is at   0.000000 -1.481700 -2.593000
      Read-in Center  127 is at   0.000000 -1.428800 -2.593000
      Read-in Center  128 is at   0.000000 -1.375900 -2.593000
      Read-in Center  129 is at   0.000000 -1.322900 -2.593000
      Read-in Center  130 is at   0.000000 -1.270000 -2.593000
      Read-in Center  131 is at   0.000000 -1.217100 -2.593000
      Read-in Center  132 is at   0.000000 -1.164200 -2.593000
      Read-in Center  133 is at   0.000000 -1.111300 -2.593000
      Read-in Center  134 is at   0.000000 -1.058400 -2.593000
      Read-in Center  135 is at   0.000000 -1.005400 -2.593000
      Read-in Center  136 is at   0.000000 -0.952500 -2.593000
      Read-in Center  137 is at   0.000000 -0.899600 -2.593000
      Read-in Center  138 is at   0.000000 -0.846700 -2.593000
      Read-in Center  139 is at   0.000000 -0.793800 -2.593000
      Read-in Center  140 is at   0.000000 -0.740800 -2.593000
      Read-in Center  141 is at   0.000000 -0.687900 -2.593000
      Read-in Center  142 is at   0.000000 -0.635000 -2.593000
      Read-in Center  143 is at   0.000000 -0.582100 -2.593000
      Read-in Center  144 is at   0.000000 -0.529200 -2.593000
      Read-in Center  145 is at   0.000000 -0.476300 -2.593000
      Read-in Center  146 is at   0.000000 -0.423300 -2.593000
      Read-in Center  147 is at   0.000000 -0.370400 -2.593000
      Read-in Center  148 is at   0.000000 -0.317500 -2.593000
      Read-in Center  149 is at   0.000000 -0.264600 -2.593000
      Read-in Center  150 is at   0.000000 -0.211700 -2.593000
      Read-in Center  151 is at   0.000000 -0.158800 -2.593000
      Read-in Center  152 is at   0.000000 -0.105800 -2.593000
      Read-in Center  153 is at   0.000000 -0.052900 -2.593000
      Read-in Center  154 is at   0.000000  0.000000 -2.593000
      Read-in Center  155 is at   0.000000  0.052900 -2.593000
      Read-in Center  156 is at   0.000000  0.105800 -2.593000
      Read-in Center  157 is at   0.000000  0.158800 -2.593000
      Read-in Center  158 is at   0.000000  0.211700 -2.593000
      Read-in Center  159 is at   0.000000  0.264600 -2.593000
      Read-in Center  160 is at   0.000000  0.317500 -2.593000
      Read-in Center  161 is at   0.000000  0.370400 -2.593000
      Read-in Center  162 is at   0.000000  0.423300 -2.593000
      Read-in Center  163 is at   0.000000  0.476300 -2.593000
      Read-in Center  164 is at   0.000000  0.529200 -2.593000
      Read-in Center  165 is at   0.000000  0.582100 -2.593000
      Read-in Center  166 is at   0.000000  0.635000 -2.593000
      Read-in Center  167 is at   0.000000  0.687900 -2.593000
      Read-in Center  168 is at   0.000000  0.740800 -2.593000
      Read-in Center  169 is at   0.000000  0.793800 -2.593000
      Read-in Center  170 is at   0.000000  0.846700 -2.593000
      Read-in Center  171 is at   0.000000  0.899600 -2.593000
      Read-in Center  172 is at   0.000000  0.952500 -2.593000
      Read-in Center  173 is at   0.000000  1.005400 -2.593000
      Read-in Center  174 is at   0.000000  1.058400 -2.593000
      Read-in Center  175 is at   0.000000  1.111300 -2.593000
      Read-in Center  176 is at   0.000000  1.164200 -2.593000
      Read-in Center  177 is at   0.000000  1.217100 -2.593000
      Read-in Center  178 is at   0.000000  1.270000 -2.593000
      Read-in Center  179 is at   0.000000  1.322900 -2.593000
      Read-in Center  180 is at   0.000000  1.375900 -2.593000
      Read-in Center  181 is at   0.000000  1.428800 -2.593000
      Read-in Center  182 is at   0.000000  1.481700 -2.593000
      Read-in Center  183 is at   0.000000  1.534600 -2.593000
      Read-in Center  184 is at   0.000000  1.587500 -2.593000
      Read-in Center  185 is at   0.000000  1.640400 -2.593000
      Read-in Center  186 is at   0.000000  1.693400 -2.593000
      Read-in Center  187 is at   0.000000  1.746300 -2.593000
      Read-in Center  188 is at   0.000000  1.799200 -2.593000
      Read-in Center  189 is at   0.000000  1.852100 -2.593000
      Read-in Center  190 is at   0.000000  1.905000 -2.593000
      Read-in Center  191 is at   0.000000  1.958000 -2.593000
      Read-in Center  192 is at   0.000000  2.010900 -2.593000
      Read-in Center  193 is at   0.000000  2.063800 -2.593000
      Read-in Center  194 is at   0.000000  2.116700 -2.593000
      Read-in Center  195 is at   0.000000  2.169600 -2.593000
      Read-in Center  196 is at   0.000000  2.222500 -2.593000
      Read-in Center  197 is at   0.000000  2.275500 -2.593000
      Read-in Center  198 is at   0.000000  2.328400 -2.593000
      Read-in Center  199 is at   0.000000  2.381300 -2.593000
      Read-in Center  200 is at   0.000000  2.434200 -2.593000
      Read-in Center  201 is at   0.000000  2.487100 -2.593000
      Read-in Center  202 is at   0.000000  2.540000 -2.593000
      Read-in Center  203 is at   0.000000  2.593000 -2.593000
      Read-in Center  204 is at   0.000000 -2.645900 -2.540000
      Read-in Center  205 is at   0.000000 -2.593000 -2.540000
      Read-in Center  206 is at   0.000000 -2.540000 -2.540000
      Read-in Center  207 is at   0.000000 -2.487100 -2.540000
      Read-in Center  208 is at   0.000000 -2.434200 -2.540000
      Read-in Center  209 is at   0.000000 -2.381300 -2.540000
      Read-in Center  210 is at   0.000000 -2.328400 -2.540000
      Read-in Center  211 is at   0.000000 -2.275500 -2.540000
      Read-in Center  212 is at   0.000000 -2.222500 -2.540000
      Read-in Center  213 is at   0.000000 -2.169600 -2.540000
      Read-in Center  214 is at   0.000000 -2.116700 -2.540000
      Read-in Center  215 is at   0.000000 -2.063800 -2.540000
      Read-in Center  216 is at   0.000000 -2.010900 -2.540000
      Read-in Center  217 is at   0.000000 -1.958000 -2.540000
      Read-in Center  218 is at   0.000000 -1.905000 -2.540000
      Read-in Center  219 is at   0.000000 -1.852100 -2.540000
      Read-in Center  220 is at   0.000000 -1.799200 -2.540000
      Read-in Center  221 is at   0.000000 -1.746300 -2.540000
      Read-in Center  222 is at   0.000000 -1.693400 -2.540000
      Read-in Center  223 is at   0.000000 -1.640400 -2.540000
      Read-in Center  224 is at   0.000000 -1.587500 -2.540000
      Read-in Center  225 is at   0.000000 -1.534600 -2.540000
      Read-in Center  226 is at   0.000000 -1.481700 -2.540000
      Read-in Center  227 is at   0.000000 -1.428800 -2.540000
      Read-in Center  228 is at   0.000000 -1.375900 -2.540000
      Read-in Center  229 is at   0.000000 -1.322900 -2.540000
      Read-in Center  230 is at   0.000000 -1.270000 -2.540000
      Read-in Center  231 is at   0.000000 -1.217100 -2.540000
      Read-in Center  232 is at   0.000000 -1.164200 -2.540000
      Read-in Center  233 is at   0.000000 -1.111300 -2.540000
      Read-in Center  234 is at   0.000000 -1.058400 -2.540000
      Read-in Center  235 is at   0.000000 -1.005400 -2.540000
      Read-in Center  236 is at   0.000000 -0.952500 -2.540000
      Read-in Center  237 is at   0.000000 -0.899600 -2.540000
      Read-in Center  238 is at   0.000000 -0.846700 -2.540000
      Read-in Center  239 is at   0.000000 -0.793800 -2.540000
      Read-in Center  240 is at   0.000000 -0.740800 -2.540000
      Read-in Center  241 is at   0.000000 -0.687900 -2.540000
      Read-in Center  242 is at   0.000000 -0.635000 -2.540000
      Read-in Center  243 is at   0.000000 -0.582100 -2.540000
      Read-in Center  244 is at   0.000000 -0.529200 -2.540000
      Read-in Center  245 is at   0.000000 -0.476300 -2.540000
      Read-in Center  246 is at   0.000000 -0.423300 -2.540000
      Read-in Center  247 is at   0.000000 -0.370400 -2.540000
      Read-in Center  248 is at   0.000000 -0.317500 -2.540000
      Read-in Center  249 is at   0.000000 -0.264600 -2.540000
      Read-in Center  250 is at   0.000000 -0.211700 -2.540000
      Read-in Center  251 is at   0.000000 -0.158800 -2.540000
      Read-in Center  252 is at   0.000000 -0.105800 -2.540000
      Read-in Center  253 is at   0.000000 -0.052900 -2.540000
      Read-in Center  254 is at   0.000000  0.000000 -2.540000
      Read-in Center  255 is at   0.000000  0.052900 -2.540000
      Read-in Center  256 is at   0.000000  0.105800 -2.540000
      Read-in Center  257 is at   0.000000  0.158800 -2.540000
      Read-in Center  258 is at   0.000000  0.211700 -2.540000
      Read-in Center  259 is at   0.000000  0.264600 -2.540000
      Read-in Center  260 is at   0.000000  0.317500 -2.540000
      Read-in Center  261 is at   0.000000  0.370400 -2.540000
      Read-in Center  262 is at   0.000000  0.423300 -2.540000
      Read-in Center  263 is at   0.000000  0.476300 -2.540000
      Read-in Center  264 is at   0.000000  0.529200 -2.540000
      Read-in Center  265 is at   0.000000  0.582100 -2.540000
      Read-in Center  266 is at   0.000000  0.635000 -2.540000
      Read-in Center  267 is at   0.000000  0.687900 -2.540000
      Read-in Center  268 is at   0.000000  0.740800 -2.540000
      Read-in Center  269 is at   0.000000  0.793800 -2.540000
      Read-in Center  270 is at   0.000000  0.846700 -2.540000
      Read-in Center  271 is at   0.000000  0.899600 -2.540000
      Read-in Center  272 is at   0.000000  0.952500 -2.540000
      Read-in Center  273 is at   0.000000  1.005400 -2.540000
      Read-in Center  274 is at   0.000000  1.058400 -2.540000
      Read-in Center  275 is at   0.000000  1.111300 -2.540000
      Read-in Center  276 is at   0.000000  1.164200 -2.540000
      Read-in Center  277 is at   0.000000  1.217100 -2.540000
      Read-in Center  278 is at   0.000000  1.270000 -2.540000
      Read-in Center  279 is at   0.000000  1.322900 -2.540000
      Read-in Center  280 is at   0.000000  1.375900 -2.540000
      Read-in Center  281 is at   0.000000  1.428800 -2.540000
      Read-in Center  282 is at   0.000000  1.481700 -2.540000
      Read-in Center  283 is at   0.000000  1.534600 -2.540000
      Read-in Center  284 is at   0.000000  1.587500 -2.540000
      Read-in Center  285 is at   0.000000  1.640400 -2.540000
      Read-in Center  286 is at   0.000000  1.693400 -2.540000
      Read-in Center  287 is at   0.000000  1.746300 -2.540000
      Read-in Center  288 is at   0.000000  1.799200 -2.540000
      Read-in Center  289 is at   0.000000  1.852100 -2.540000
      Read-in Center  290 is at   0.000000  1.905000 -2.540000
      Read-in Center  291 is at   0.000000  1.958000 -2.540000
      Read-in Center  292 is at   0.000000  2.010900 -2.540000
      Read-in Center  293 is at   0.000000  2.063800 -2.540000
      Read-in Center  294 is at   0.000000  2.116700 -2.540000
      Read-in Center  295 is at   0.000000  2.169600 -2.540000
      Read-in Center  296 is at   0.000000  2.222500 -2.540000
      Read-in Center  297 is at   0.000000  2.275500 -2.540000
      Read-in Center  298 is at   0.000000  2.328400 -2.540000
      Read-in Center  299 is at   0.000000  2.381300 -2.540000
      Read-in Center  300 is at   0.000000  2.434200 -2.540000
      Read-in Center  301 is at   0.000000  2.487100 -2.540000
      Read-in Center  302 is at   0.000000  2.540000 -2.540000
      Read-in Center  303 is at   0.000000  2.593000 -2.540000
      Read-in Center  304 is at   0.000000 -2.645900 -2.487100
      Read-in Center  305 is at   0.000000 -2.593000 -2.487100
      Read-in Center  306 is at   0.000000 -2.540000 -2.487100
      Read-in Center  307 is at   0.000000 -2.487100 -2.487100
      Read-in Center  308 is at   0.000000 -2.434200 -2.487100
      Read-in Center  309 is at   0.000000 -2.381300 -2.487100
      Read-in Center  310 is at   0.000000 -2.328400 -2.487100
      Read-in Center  311 is at   0.000000 -2.275500 -2.487100
      Read-in Center  312 is at   0.000000 -2.222500 -2.487100
      Read-in Center  313 is at   0.000000 -2.169600 -2.487100
      Read-in Center  314 is at   0.000000 -2.116700 -2.487100
      Read-in Center  315 is at   0.000000 -2.063800 -2.487100
      Read-in Center  316 is at   0.000000 -2.010900 -2.487100
      Read-in Center  317 is at   0.000000 -1.958000 -2.487100
      Read-in Center  318 is at   0.000000 -1.905000 -2.487100
      Read-in Center  319 is at   0.000000 -1.852100 -2.487100
      Read-in Center  320 is at   0.000000 -1.799200 -2.487100
      Read-in Center  321 is at   0.000000 -1.746300 -2.487100
      Read-in Center  322 is at   0.000000 -1.693400 -2.487100
      Read-in Center  323 is at   0.000000 -1.640400 -2.487100
      Read-in Center  324 is at   0.000000 -1.587500 -2.487100
      Read-in Center  325 is at   0.000000 -1.534600 -2.487100
      Read-in Center  326 is at   0.000000 -1.481700 -2.487100
      Read-in Center  327 is at   0.000000 -1.428800 -2.487100
      Read-in Center  328 is at   0.000000 -1.375900 -2.487100
      Read-in Center  329 is at   0.000000 -1.322900 -2.487100
      Read-in Center  330 is at   0.000000 -1.270000 -2.487100
      Read-in Center  331 is at   0.000000 -1.217100 -2.487100
      Read-in Center  332 is at   0.000000 -1.164200 -2.487100
      Read-in Center  333 is at   0.000000 -1.111300 -2.487100
      Read-in Center  334 is at   0.000000 -1.058400 -2.487100
      Read-in Center  335 is at   0.000000 -1.005400 -2.487100
      Read-in Center  336 is at   0.000000 -0.952500 -2.487100
      Read-in Center  337 is at   0.000000 -0.899600 -2.487100
      Read-in Center  338 is at   0.000000 -0.846700 -2.487100
      Read-in Center  339 is at   0.000000 -0.793800 -2.487100
      Read-in Center  340 is at   0.000000 -0.740800 -2.487100
      Read-in Center  341 is at   0.000000 -0.687900 -2.487100
      Read-in Center  342 is at   0.000000 -0.635000 -2.487100
      Read-in Center  343 is at   0.000000 -0.582100 -2.487100
      Read-in Center  344 is at   0.000000 -0.529200 -2.487100
      Read-in Center  345 is at   0.000000 -0.476300 -2.487100
      Read-in Center  346 is at   0.000000 -0.423300 -2.487100
      Read-in Center  347 is at   0.000000 -0.370400 -2.487100
      Read-in Center  348 is at   0.000000 -0.317500 -2.487100
      Read-in Center  349 is at   0.000000 -0.264600 -2.487100
      Read-in Center  350 is at   0.000000 -0.211700 -2.487100
      Read-in Center  351 is at   0.000000 -0.158800 -2.487100
      Read-in Center  352 is at   0.000000 -0.105800 -2.487100
      Read-in Center  353 is at   0.000000 -0.052900 -2.487100
      Read-in Center  354 is at   0.000000  0.000000 -2.487100
      Read-in Center  355 is at   0.000000  0.052900 -2.487100
      Read-in Center  356 is at   0.000000  0.105800 -2.487100
      Read-in Center  357 is at   0.000000  0.158800 -2.487100
      Read-in Center  358 is at   0.000000  0.211700 -2.487100
      Read-in Center  359 is at   0.000000  0.264600 -2.487100
      Read-in Center  360 is at   0.000000  0.317500 -2.487100
      Read-in Center  361 is at   0.000000  0.370400 -2.487100
      Read-in Center  362 is at   0.000000  0.423300 -2.487100
      Read-in Center  363 is at   0.000000  0.476300 -2.487100
      Read-in Center  364 is at   0.000000  0.529200 -2.487100
      Read-in Center  365 is at   0.000000  0.582100 -2.487100
      Read-in Center  366 is at   0.000000  0.635000 -2.487100
      Read-in Center  367 is at   0.000000  0.687900 -2.487100
      Read-in Center  368 is at   0.000000  0.740800 -2.487100
      Read-in Center  369 is at   0.000000  0.793800 -2.487100
      Read-in Center  370 is at   0.000000  0.846700 -2.487100
      Read-in Center  371 is at   0.000000  0.899600 -2.487100
      Read-in Center  372 is at   0.000000  0.952500 -2.487100
      Read-in Center  373 is at   0.000000  1.005400 -2.487100
      Read-in Center  374 is at   0.000000  1.058400 -2.487100
      Read-in Center  375 is at   0.000000  1.111300 -2.487100
      Read-in Center  376 is at   0.000000  1.164200 -2.487100
      Read-in Center  377 is at   0.000000  1.217100 -2.487100
      Read-in Center  378 is at   0.000000  1.270000 -2.487100
      Read-in Center  379 is at   0.000000  1.322900 -2.487100
      Read-in Center  380 is at   0.000000  1.375900 -2.487100
      Read-in Center  381 is at   0.000000  1.428800 -2.487100
      Read-in Center  382 is at   0.000000  1.481700 -2.487100
      Read-in Center  383 is at   0.000000  1.534600 -2.487100
      Read-in Center  384 is at   0.000000  1.587500 -2.487100
      Read-in Center  385 is at   0.000000  1.640400 -2.487100
      Read-in Center  386 is at   0.000000  1.693400 -2.487100
      Read-in Center  387 is at   0.000000  1.746300 -2.487100
      Read-in Center  388 is at   0.000000  1.799200 -2.487100
      Read-in Center  389 is at   0.000000  1.852100 -2.487100
      Read-in Center  390 is at   0.000000  1.905000 -2.487100
      Read-in Center  391 is at   0.000000  1.958000 -2.487100
      Read-in Center  392 is at   0.000000  2.010900 -2.487100
      Read-in Center  393 is at   0.000000  2.063800 -2.487100
      Read-in Center  394 is at   0.000000  2.116700 -2.487100
      Read-in Center  395 is at   0.000000  2.169600 -2.487100
      Read-in Center  396 is at   0.000000  2.222500 -2.487100
      Read-in Center  397 is at   0.000000  2.275500 -2.487100
      Read-in Center  398 is at   0.000000  2.328400 -2.487100
      Read-in Center  399 is at   0.000000  2.381300 -2.487100
      Read-in Center  400 is at   0.000000  2.434200 -2.487100
      Read-in Center  401 is at   0.000000  2.487100 -2.487100
      Read-in Center  402 is at   0.000000  2.540000 -2.487100
      Read-in Center  403 is at   0.000000  2.593000 -2.487100
      Read-in Center  404 is at   0.000000 -2.645900 -2.434200
      Read-in Center  405 is at   0.000000 -2.593000 -2.434200
      Read-in Center  406 is at   0.000000 -2.540000 -2.434200
      Read-in Center  407 is at   0.000000 -2.487100 -2.434200
      Read-in Center  408 is at   0.000000 -2.434200 -2.434200
      Read-in Center  409 is at   0.000000 -2.381300 -2.434200
      Read-in Center  410 is at   0.000000 -2.328400 -2.434200
      Read-in Center  411 is at   0.000000 -2.275500 -2.434200
      Read-in Center  412 is at   0.000000 -2.222500 -2.434200
      Read-in Center  413 is at   0.000000 -2.169600 -2.434200
      Read-in Center  414 is at   0.000000 -2.116700 -2.434200
      Read-in Center  415 is at   0.000000 -2.063800 -2.434200
      Read-in Center  416 is at   0.000000 -2.010900 -2.434200
      Read-in Center  417 is at   0.000000 -1.958000 -2.434200
      Read-in Center  418 is at   0.000000 -1.905000 -2.434200
      Read-in Center  419 is at   0.000000 -1.852100 -2.434200
      Read-in Center  420 is at   0.000000 -1.799200 -2.434200
      Read-in Center  421 is at   0.000000 -1.746300 -2.434200
      Read-in Center  422 is at   0.000000 -1.693400 -2.434200
      Read-in Center  423 is at   0.000000 -1.640400 -2.434200
      Read-in Center  424 is at   0.000000 -1.587500 -2.434200
      Read-in Center  425 is at   0.000000 -1.534600 -2.434200
      Read-in Center  426 is at   0.000000 -1.481700 -2.434200
      Read-in Center  427 is at   0.000000 -1.428800 -2.434200
      Read-in Center  428 is at   0.000000 -1.375900 -2.434200
      Read-in Center  429 is at   0.000000 -1.322900 -2.434200
      Read-in Center  430 is at   0.000000 -1.270000 -2.434200
      Read-in Center  431 is at   0.000000 -1.217100 -2.434200
      Read-in Center  432 is at   0.000000 -1.164200 -2.434200
      Read-in Center  433 is at   0.000000 -1.111300 -2.434200
      Read-in Center  434 is at   0.000000 -1.058400 -2.434200
      Read-in Center  435 is at   0.000000 -1.005400 -2.434200
      Read-in Center  436 is at   0.000000 -0.952500 -2.434200
      Read-in Center  437 is at   0.000000 -0.899600 -2.434200
      Read-in Center  438 is at   0.000000 -0.846700 -2.434200
      Read-in Center  439 is at   0.000000 -0.793800 -2.434200
      Read-in Center  440 is at   0.000000 -0.740800 -2.434200
      Read-in Center  441 is at   0.000000 -0.687900 -2.434200
      Read-in Center  442 is at   0.000000 -0.635000 -2.434200
      Read-in Center  443 is at   0.000000 -0.582100 -2.434200
      Read-in Center  444 is at   0.000000 -0.529200 -2.434200
      Read-in Center  445 is at   0.000000 -0.476300 -2.434200
      Read-in Center  446 is at   0.000000 -0.423300 -2.434200
      Read-in Center  447 is at   0.000000 -0.370400 -2.434200
      Read-in Center  448 is at   0.000000 -0.317500 -2.434200
      Read-in Center  449 is at   0.000000 -0.264600 -2.434200
      Read-in Center  450 is at   0.000000 -0.211700 -2.434200
      Read-in Center  451 is at   0.000000 -0.158800 -2.434200
      Read-in Center  452 is at   0.000000 -0.105800 -2.434200
      Read-in Center  453 is at   0.000000 -0.052900 -2.434200
      Read-in Center  454 is at   0.000000  0.000000 -2.434200
      Read-in Center  455 is at   0.000000  0.052900 -2.434200
      Read-in Center  456 is at   0.000000  0.105800 -2.434200
      Read-in Center  457 is at   0.000000  0.158800 -2.434200
      Read-in Center  458 is at   0.000000  0.211700 -2.434200
      Read-in Center  459 is at   0.000000  0.264600 -2.434200
      Read-in Center  460 is at   0.000000  0.317500 -2.434200
      Read-in Center  461 is at   0.000000  0.370400 -2.434200
      Read-in Center  462 is at   0.000000  0.423300 -2.434200
      Read-in Center  463 is at   0.000000  0.476300 -2.434200
      Read-in Center  464 is at   0.000000  0.529200 -2.434200
      Read-in Center  465 is at   0.000000  0.582100 -2.434200
      Read-in Center  466 is at   0.000000  0.635000 -2.434200
      Read-in Center  467 is at   0.000000  0.687900 -2.434200
      Read-in Center  468 is at   0.000000  0.740800 -2.434200
      Read-in Center  469 is at   0.000000  0.793800 -2.434200
      Read-in Center  470 is at   0.000000  0.846700 -2.434200
      Read-in Center  471 is at   0.000000  0.899600 -2.434200
      Read-in Center  472 is at   0.000000  0.952500 -2.434200
      Read-in Center  473 is at   0.000000  1.005400 -2.434200
      Read-in Center  474 is at   0.000000  1.058400 -2.434200
      Read-in Center  475 is at   0.000000  1.111300 -2.434200
      Read-in Center  476 is at   0.000000  1.164200 -2.434200
      Read-in Center  477 is at   0.000000  1.217100 -2.434200
      Read-in Center  478 is at   0.000000  1.270000 -2.434200
      Read-in Center  479 is at   0.000000  1.322900 -2.434200
      Read-in Center  480 is at   0.000000  1.375900 -2.434200
      Read-in Center  481 is at   0.000000  1.428800 -2.434200
      Read-in Center  482 is at   0.000000  1.481700 -2.434200
      Read-in Center  483 is at   0.000000  1.534600 -2.434200
      Read-in Center  484 is at   0.000000  1.587500 -2.434200
      Read-in Center  485 is at   0.000000  1.640400 -2.434200
      Read-in Center  486 is at   0.000000  1.693400 -2.434200
      Read-in Center  487 is at   0.000000  1.746300 -2.434200
      Read-in Center  488 is at   0.000000  1.799200 -2.434200
      Read-in Center  489 is at   0.000000  1.852100 -2.434200
      Read-in Center  490 is at   0.000000  1.905000 -2.434200
      Read-in Center  491 is at   0.000000  1.958000 -2.434200
      Read-in Center  492 is at   0.000000  2.010900 -2.434200
      Read-in Center  493 is at   0.000000  2.063800 -2.434200
      Read-in Center  494 is at   0.000000  2.116700 -2.434200
      Read-in Center  495 is at   0.000000  2.169600 -2.434200
      Read-in Center  496 is at   0.000000  2.222500 -2.434200
      Read-in Center  497 is at   0.000000  2.275500 -2.434200
      Read-in Center  498 is at   0.000000  2.328400 -2.434200
      Read-in Center  499 is at   0.000000  2.381300 -2.434200
      Read-in Center  500 is at   0.000000  2.434200 -2.434200
      Read-in Center  501 is at   0.000000  2.487100 -2.434200
      Read-in Center  502 is at   0.000000  2.540000 -2.434200
      Read-in Center  503 is at   0.000000  2.593000 -2.434200
      Read-in Center  504 is at   0.000000 -2.645900 -2.381300
      Read-in Center  505 is at   0.000000 -2.593000 -2.381300
      Read-in Center  506 is at   0.000000 -2.540000 -2.381300
      Read-in Center  507 is at   0.000000 -2.487100 -2.381300
      Read-in Center  508 is at   0.000000 -2.434200 -2.381300
      Read-in Center  509 is at   0.000000 -2.381300 -2.381300
      Read-in Center  510 is at   0.000000 -2.328400 -2.381300
      Read-in Center  511 is at   0.000000 -2.275500 -2.381300
      Read-in Center  512 is at   0.000000 -2.222500 -2.381300
      Read-in Center  513 is at   0.000000 -2.169600 -2.381300
      Read-in Center  514 is at   0.000000 -2.116700 -2.381300
      Read-in Center  515 is at   0.000000 -2.063800 -2.381300
      Read-in Center  516 is at   0.000000 -2.010900 -2.381300
      Read-in Center  517 is at   0.000000 -1.958000 -2.381300
      Read-in Center  518 is at   0.000000 -1.905000 -2.381300
      Read-in Center  519 is at   0.000000 -1.852100 -2.381300
      Read-in Center  520 is at   0.000000 -1.799200 -2.381300
      Read-in Center  521 is at   0.000000 -1.746300 -2.381300
      Read-in Center  522 is at   0.000000 -1.693400 -2.381300
      Read-in Center  523 is at   0.000000 -1.640400 -2.381300
      Read-in Center  524 is at   0.000000 -1.587500 -2.381300
      Read-in Center  525 is at   0.000000 -1.534600 -2.381300
      Read-in Center  526 is at   0.000000 -1.481700 -2.381300
      Read-in Center  527 is at   0.000000 -1.428800 -2.381300
      Read-in Center  528 is at   0.000000 -1.375900 -2.381300
      Read-in Center  529 is at   0.000000 -1.322900 -2.381300
      Read-in Center  530 is at   0.000000 -1.270000 -2.381300
      Read-in Center  531 is at   0.000000 -1.217100 -2.381300
      Read-in Center  532 is at   0.000000 -1.164200 -2.381300
      Read-in Center  533 is at   0.000000 -1.111300 -2.381300
      Read-in Center  534 is at   0.000000 -1.058400 -2.381300
      Read-in Center  535 is at   0.000000 -1.005400 -2.381300
      Read-in Center  536 is at   0.000000 -0.952500 -2.381300
      Read-in Center  537 is at   0.000000 -0.899600 -2.381300
      Read-in Center  538 is at   0.000000 -0.846700 -2.381300
      Read-in Center  539 is at   0.000000 -0.793800 -2.381300
      Read-in Center  540 is at   0.000000 -0.740800 -2.381300
      Read-in Center  541 is at   0.000000 -0.687900 -2.381300
      Read-in Center  542 is at   0.000000 -0.635000 -2.381300
      Read-in Center  543 is at   0.000000 -0.582100 -2.381300
      Read-in Center  544 is at   0.000000 -0.529200 -2.381300
      Read-in Center  545 is at   0.000000 -0.476300 -2.381300
      Read-in Center  546 is at   0.000000 -0.423300 -2.381300
      Read-in Center  547 is at   0.000000 -0.370400 -2.381300
      Read-in Center  548 is at   0.000000 -0.317500 -2.381300
      Read-in Center  549 is at   0.000000 -0.264600 -2.381300
      Read-in Center  550 is at   0.000000 -0.211700 -2.381300
      Read-in Center  551 is at   0.000000 -0.158800 -2.381300
      Read-in Center  552 is at   0.000000 -0.105800 -2.381300
      Read-in Center  553 is at   0.000000 -0.052900 -2.381300
      Read-in Center  554 is at   0.000000  0.000000 -2.381300
      Read-in Center  555 is at   0.000000  0.052900 -2.381300
      Read-in Center  556 is at   0.000000  0.105800 -2.381300
      Read-in Center  557 is at   0.000000  0.158800 -2.381300
      Read-in Center  558 is at   0.000000  0.211700 -2.381300
      Read-in Center  559 is at   0.000000  0.264600 -2.381300
      Read-in Center  560 is at   0.000000  0.317500 -2.381300
      Read-in Center  561 is at   0.000000  0.370400 -2.381300
      Read-in Center  562 is at   0.000000  0.423300 -2.381300
      Read-in Center  563 is at   0.000000  0.476300 -2.381300
      Read-in Center  564 is at   0.000000  0.529200 -2.381300
      Read-in Center  565 is at   0.000000  0.582100 -2.381300
      Read-in Center  566 is at   0.000000  0.635000 -2.381300
      Read-in Center  567 is at   0.000000  0.687900 -2.381300
      Read-in Center  568 is at   0.000000  0.740800 -2.381300
      Read-in Center  569 is at   0.000000  0.793800 -2.381300
      Read-in Center  570 is at   0.000000  0.846700 -2.381300
      Read-in Center  571 is at   0.000000  0.899600 -2.381300
      Read-in Center  572 is at   0.000000  0.952500 -2.381300
      Read-in Center  573 is at   0.000000  1.005400 -2.381300
      Read-in Center  574 is at   0.000000  1.058400 -2.381300
      Read-in Center  575 is at   0.000000  1.111300 -2.381300
      Read-in Center  576 is at   0.000000  1.164200 -2.381300
      Read-in Center  577 is at   0.000000  1.217100 -2.381300
      Read-in Center  578 is at   0.000000  1.270000 -2.381300
      Read-in Center  579 is at   0.000000  1.322900 -2.381300
      Read-in Center  580 is at   0.000000  1.375900 -2.381300
      Read-in Center  581 is at   0.000000  1.428800 -2.381300
      Read-in Center  582 is at   0.000000  1.481700 -2.381300
      Read-in Center  583 is at   0.000000  1.534600 -2.381300
      Read-in Center  584 is at   0.000000  1.587500 -2.381300
      Read-in Center  585 is at   0.000000  1.640400 -2.381300
      Read-in Center  586 is at   0.000000  1.693400 -2.381300
      Read-in Center  587 is at   0.000000  1.746300 -2.381300
      Read-in Center  588 is at   0.000000  1.799200 -2.381300
      Read-in Center  589 is at   0.000000  1.852100 -2.381300
      Read-in Center  590 is at   0.000000  1.905000 -2.381300
      Read-in Center  591 is at   0.000000  1.958000 -2.381300
      Read-in Center  592 is at   0.000000  2.010900 -2.381300
      Read-in Center  593 is at   0.000000  2.063800 -2.381300
      Read-in Center  594 is at   0.000000  2.116700 -2.381300
      Read-in Center  595 is at   0.000000  2.169600 -2.381300
      Read-in Center  596 is at   0.000000  2.222500 -2.381300
      Read-in Center  597 is at   0.000000  2.275500 -2.381300
      Read-in Center  598 is at   0.000000  2.328400 -2.381300
      Read-in Center  599 is at   0.000000  2.381300 -2.381300
      Read-in Center  600 is at   0.000000  2.434200 -2.381300
      Read-in Center  601 is at   0.000000  2.487100 -2.381300
      Read-in Center  602 is at   0.000000  2.540000 -2.381300
      Read-in Center  603 is at   0.000000  2.593000 -2.381300
      Read-in Center  604 is at   0.000000 -2.645900 -2.328400
      Read-in Center  605 is at   0.000000 -2.593000 -2.328400
      Read-in Center  606 is at   0.000000 -2.540000 -2.328400
      Read-in Center  607 is at   0.000000 -2.487100 -2.328400
      Read-in Center  608 is at   0.000000 -2.434200 -2.328400
      Read-in Center  609 is at   0.000000 -2.381300 -2.328400
      Read-in Center  610 is at   0.000000 -2.328400 -2.328400
      Read-in Center  611 is at   0.000000 -2.275500 -2.328400
      Read-in Center  612 is at   0.000000 -2.222500 -2.328400
      Read-in Center  613 is at   0.000000 -2.169600 -2.328400
      Read-in Center  614 is at   0.000000 -2.116700 -2.328400
      Read-in Center  615 is at   0.000000 -2.063800 -2.328400
      Read-in Center  616 is at   0.000000 -2.010900 -2.328400
      Read-in Center  617 is at   0.000000 -1.958000 -2.328400
      Read-in Center  618 is at   0.000000 -1.905000 -2.328400
      Read-in Center  619 is at   0.000000 -1.852100 -2.328400
      Read-in Center  620 is at   0.000000 -1.799200 -2.328400
      Read-in Center  621 is at   0.000000 -1.746300 -2.328400
      Read-in Center  622 is at   0.000000 -1.693400 -2.328400
      Read-in Center  623 is at   0.000000 -1.640400 -2.328400
      Read-in Center  624 is at   0.000000 -1.587500 -2.328400
      Read-in Center  625 is at   0.000000 -1.534600 -2.328400
      Read-in Center  626 is at   0.000000 -1.481700 -2.328400
      Read-in Center  627 is at   0.000000 -1.428800 -2.328400
      Read-in Center  628 is at   0.000000 -1.375900 -2.328400
      Read-in Center  629 is at   0.000000 -1.322900 -2.328400
      Read-in Center  630 is at   0.000000 -1.270000 -2.328400
      Read-in Center  631 is at   0.000000 -1.217100 -2.328400
      Read-in Center  632 is at   0.000000 -1.164200 -2.328400
      Read-in Center  633 is at   0.000000 -1.111300 -2.328400
      Read-in Center  634 is at   0.000000 -1.058400 -2.328400
      Read-in Center  635 is at   0.000000 -1.005400 -2.328400
      Read-in Center  636 is at   0.000000 -0.952500 -2.328400
      Read-in Center  637 is at   0.000000 -0.899600 -2.328400
      Read-in Center  638 is at   0.000000 -0.846700 -2.328400
      Read-in Center  639 is at   0.000000 -0.793800 -2.328400
      Read-in Center  640 is at   0.000000 -0.740800 -2.328400
      Read-in Center  641 is at   0.000000 -0.687900 -2.328400
      Read-in Center  642 is at   0.000000 -0.635000 -2.328400
      Read-in Center  643 is at   0.000000 -0.582100 -2.328400
      Read-in Center  644 is at   0.000000 -0.529200 -2.328400
      Read-in Center  645 is at   0.000000 -0.476300 -2.328400
      Read-in Center  646 is at   0.000000 -0.423300 -2.328400
      Read-in Center  647 is at   0.000000 -0.370400 -2.328400
      Read-in Center  648 is at   0.000000 -0.317500 -2.328400
      Read-in Center  649 is at   0.000000 -0.264600 -2.328400
      Read-in Center  650 is at   0.000000 -0.211700 -2.328400
      Read-in Center  651 is at   0.000000 -0.158800 -2.328400
      Read-in Center  652 is at   0.000000 -0.105800 -2.328400
      Read-in Center  653 is at   0.000000 -0.052900 -2.328400
      Read-in Center  654 is at   0.000000  0.000000 -2.328400
      Read-in Center  655 is at   0.000000  0.052900 -2.328400
      Read-in Center  656 is at   0.000000  0.105800 -2.328400
      Read-in Center  657 is at   0.000000  0.158800 -2.328400
      Read-in Center  658 is at   0.000000  0.211700 -2.328400
      Read-in Center  659 is at   0.000000  0.264600 -2.328400
      Read-in Center  660 is at   0.000000  0.317500 -2.328400
      Read-in Center  661 is at   0.000000  0.370400 -2.328400
      Read-in Center  662 is at   0.000000  0.423300 -2.328400
      Read-in Center  663 is at   0.000000  0.476300 -2.328400
      Read-in Center  664 is at   0.000000  0.529200 -2.328400
      Read-in Center  665 is at   0.000000  0.582100 -2.328400
      Read-in Center  666 is at   0.000000  0.635000 -2.328400
      Read-in Center  667 is at   0.000000  0.687900 -2.328400
      Read-in Center  668 is at   0.000000  0.740800 -2.328400
      Read-in Center  669 is at   0.000000  0.793800 -2.328400
      Read-in Center  670 is at   0.000000  0.846700 -2.328400
      Read-in Center  671 is at   0.000000  0.899600 -2.328400
      Read-in Center  672 is at   0.000000  0.952500 -2.328400
      Read-in Center  673 is at   0.000000  1.005400 -2.328400
      Read-in Center  674 is at   0.000000  1.058400 -2.328400
      Read-in Center  675 is at   0.000000  1.111300 -2.328400
      Read-in Center  676 is at   0.000000  1.164200 -2.328400
      Read-in Center  677 is at   0.000000  1.217100 -2.328400
      Read-in Center  678 is at   0.000000  1.270000 -2.328400
      Read-in Center  679 is at   0.000000  1.322900 -2.328400
      Read-in Center  680 is at   0.000000  1.375900 -2.328400
      Read-in Center  681 is at   0.000000  1.428800 -2.328400
      Read-in Center  682 is at   0.000000  1.481700 -2.328400
      Read-in Center  683 is at   0.000000  1.534600 -2.328400
      Read-in Center  684 is at   0.000000  1.587500 -2.328400
      Read-in Center  685 is at   0.000000  1.640400 -2.328400
      Read-in Center  686 is at   0.000000  1.693400 -2.328400
      Read-in Center  687 is at   0.000000  1.746300 -2.328400
      Read-in Center  688 is at   0.000000  1.799200 -2.328400
      Read-in Center  689 is at   0.000000  1.852100 -2.328400
      Read-in Center  690 is at   0.000000  1.905000 -2.328400
      Read-in Center  691 is at   0.000000  1.958000 -2.328400
      Read-in Center  692 is at   0.000000  2.010900 -2.328400
      Read-in Center  693 is at   0.000000  2.063800 -2.328400
      Read-in Center  694 is at   0.000000  2.116700 -2.328400
      Read-in Center  695 is at   0.000000  2.169600 -2.328400
      Read-in Center  696 is at   0.000000  2.222500 -2.328400
      Read-in Center  697 is at   0.000000  2.275500 -2.328400
      Read-in Center  698 is at   0.000000  2.328400 -2.328400
      Read-in Center  699 is at   0.000000  2.381300 -2.328400
      Read-in Center  700 is at   0.000000  2.434200 -2.328400
      Read-in Center  701 is at   0.000000  2.487100 -2.328400
      Read-in Center  702 is at   0.000000  2.540000 -2.328400
      Read-in Center  703 is at   0.000000  2.593000 -2.328400
      Read-in Center  704 is at   0.000000 -2.645900 -2.275500
      Read-in Center  705 is at   0.000000 -2.593000 -2.275500
      Read-in Center  706 is at   0.000000 -2.540000 -2.275500
      Read-in Center  707 is at   0.000000 -2.487100 -2.275500
      Read-in Center  708 is at   0.000000 -2.434200 -2.275500
      Read-in Center  709 is at   0.000000 -2.381300 -2.275500
      Read-in Center  710 is at   0.000000 -2.328400 -2.275500
      Read-in Center  711 is at   0.000000 -2.275500 -2.275500
      Read-in Center  712 is at   0.000000 -2.222500 -2.275500
      Read-in Center  713 is at   0.000000 -2.169600 -2.275500
      Read-in Center  714 is at   0.000000 -2.116700 -2.275500
      Read-in Center  715 is at   0.000000 -2.063800 -2.275500
      Read-in Center  716 is at   0.000000 -2.010900 -2.275500
      Read-in Center  717 is at   0.000000 -1.958000 -2.275500
      Read-in Center  718 is at   0.000000 -1.905000 -2.275500
      Read-in Center  719 is at   0.000000 -1.852100 -2.275500
      Read-in Center  720 is at   0.000000 -1.799200 -2.275500
      Read-in Center  721 is at   0.000000 -1.746300 -2.275500
      Read-in Center  722 is at   0.000000 -1.693400 -2.275500
      Read-in Center  723 is at   0.000000 -1.640400 -2.275500
      Read-in Center  724 is at   0.000000 -1.587500 -2.275500
      Read-in Center  725 is at   0.000000 -1.534600 -2.275500
      Read-in Center  726 is at   0.000000 -1.481700 -2.275500
      Read-in Center  727 is at   0.000000 -1.428800 -2.275500
      Read-in Center  728 is at   0.000000 -1.375900 -2.275500
      Read-in Center  729 is at   0.000000 -1.322900 -2.275500
      Read-in Center  730 is at   0.000000 -1.270000 -2.275500
      Read-in Center  731 is at   0.000000 -1.217100 -2.275500
      Read-in Center  732 is at   0.000000 -1.164200 -2.275500
      Read-in Center  733 is at   0.000000 -1.111300 -2.275500
      Read-in Center  734 is at   0.000000 -1.058400 -2.275500
      Read-in Center  735 is at   0.000000 -1.005400 -2.275500
      Read-in Center  736 is at   0.000000 -0.952500 -2.275500
      Read-in Center  737 is at   0.000000 -0.899600 -2.275500
      Read-in Center  738 is at   0.000000 -0.846700 -2.275500
      Read-in Center  739 is at   0.000000 -0.793800 -2.275500
      Read-in Center  740 is at   0.000000 -0.740800 -2.275500
      Read-in Center  741 is at   0.000000 -0.687900 -2.275500
      Read-in Center  742 is at   0.000000 -0.635000 -2.275500
      Read-in Center  743 is at   0.000000 -0.582100 -2.275500
      Read-in Center  744 is at   0.000000 -0.529200 -2.275500
      Read-in Center  745 is at   0.000000 -0.476300 -2.275500
      Read-in Center  746 is at   0.000000 -0.423300 -2.275500
      Read-in Center  747 is at   0.000000 -0.370400 -2.275500
      Read-in Center  748 is at   0.000000 -0.317500 -2.275500
      Read-in Center  749 is at   0.000000 -0.264600 -2.275500
      Read-in Center  750 is at   0.000000 -0.211700 -2.275500
      Read-in Center  751 is at   0.000000 -0.158800 -2.275500
      Read-in Center  752 is at   0.000000 -0.105800 -2.275500
      Read-in Center  753 is at   0.000000 -0.052900 -2.275500
      Read-in Center  754 is at   0.000000  0.000000 -2.275500
      Read-in Center  755 is at   0.000000  0.052900 -2.275500
      Read-in Center  756 is at   0.000000  0.105800 -2.275500
      Read-in Center  757 is at   0.000000  0.158800 -2.275500
      Read-in Center  758 is at   0.000000  0.211700 -2.275500
      Read-in Center  759 is at   0.000000  0.264600 -2.275500
      Read-in Center  760 is at   0.000000  0.317500 -2.275500
      Read-in Center  761 is at   0.000000  0.370400 -2.275500
      Read-in Center  762 is at   0.000000  0.423300 -2.275500
      Read-in Center  763 is at   0.000000  0.476300 -2.275500
      Read-in Center  764 is at   0.000000  0.529200 -2.275500
      Read-in Center  765 is at   0.000000  0.582100 -2.275500
      Read-in Center  766 is at   0.000000  0.635000 -2.275500
      Read-in Center  767 is at   0.000000  0.687900 -2.275500
      Read-in Center  768 is at   0.000000  0.740800 -2.275500
      Read-in Center  769 is at   0.000000  0.793800 -2.275500
      Read-in Center  770 is at   0.000000  0.846700 -2.275500
      Read-in Center  771 is at   0.000000  0.899600 -2.275500
      Read-in Center  772 is at   0.000000  0.952500 -2.275500
      Read-in Center  773 is at   0.000000  1.005400 -2.275500
      Read-in Center  774 is at   0.000000  1.058400 -2.275500
      Read-in Center  775 is at   0.000000  1.111300 -2.275500
      Read-in Center  776 is at   0.000000  1.164200 -2.275500
      Read-in Center  777 is at   0.000000  1.217100 -2.275500
      Read-in Center  778 is at   0.000000  1.270000 -2.275500
      Read-in Center  779 is at   0.000000  1.322900 -2.275500
      Read-in Center  780 is at   0.000000  1.375900 -2.275500
      Read-in Center  781 is at   0.000000  1.428800 -2.275500
      Read-in Center  782 is at   0.000000  1.481700 -2.275500
      Read-in Center  783 is at   0.000000  1.534600 -2.275500
      Read-in Center  784 is at   0.000000  1.587500 -2.275500
      Read-in Center  785 is at   0.000000  1.640400 -2.275500
      Read-in Center  786 is at   0.000000  1.693400 -2.275500
      Read-in Center  787 is at   0.000000  1.746300 -2.275500
      Read-in Center  788 is at   0.000000  1.799200 -2.275500
      Read-in Center  789 is at   0.000000  1.852100 -2.275500
      Read-in Center  790 is at   0.000000  1.905000 -2.275500
      Read-in Center  791 is at   0.000000  1.958000 -2.275500
      Read-in Center  792 is at   0.000000  2.010900 -2.275500
      Read-in Center  793 is at   0.000000  2.063800 -2.275500
      Read-in Center  794 is at   0.000000  2.116700 -2.275500
      Read-in Center  795 is at   0.000000  2.169600 -2.275500
      Read-in Center  796 is at   0.000000  2.222500 -2.275500
      Read-in Center  797 is at   0.000000  2.275500 -2.275500
      Read-in Center  798 is at   0.000000  2.328400 -2.275500
      Read-in Center  799 is at   0.000000  2.381300 -2.275500
      Read-in Center  800 is at   0.000000  2.434200 -2.275500
      Read-in Center  801 is at   0.000000  2.487100 -2.275500
      Read-in Center  802 is at   0.000000  2.540000 -2.275500
      Read-in Center  803 is at   0.000000  2.593000 -2.275500
      Read-in Center  804 is at   0.000000 -2.645900 -2.222500
      Read-in Center  805 is at   0.000000 -2.593000 -2.222500
      Read-in Center  806 is at   0.000000 -2.540000 -2.222500
      Read-in Center  807 is at   0.000000 -2.487100 -2.222500
      Read-in Center  808 is at   0.000000 -2.434200 -2.222500
      Read-in Center  809 is at   0.000000 -2.381300 -2.222500
      Read-in Center  810 is at   0.000000 -2.328400 -2.222500
      Read-in Center  811 is at   0.000000 -2.275500 -2.222500
      Read-in Center  812 is at   0.000000 -2.222500 -2.222500
      Read-in Center  813 is at   0.000000 -2.169600 -2.222500
      Read-in Center  814 is at   0.000000 -2.116700 -2.222500
      Read-in Center  815 is at   0.000000 -2.063800 -2.222500
      Read-in Center  816 is at   0.000000 -2.010900 -2.222500
      Read-in Center  817 is at   0.000000 -1.958000 -2.222500
      Read-in Center  818 is at   0.000000 -1.905000 -2.222500
      Read-in Center  819 is at   0.000000 -1.852100 -2.222500
      Read-in Center  820 is at   0.000000 -1.799200 -2.222500
      Read-in Center  821 is at   0.000000 -1.746300 -2.222500
      Read-in Center  822 is at   0.000000 -1.693400 -2.222500
      Read-in Center  823 is at   0.000000 -1.640400 -2.222500
      Read-in Center  824 is at   0.000000 -1.587500 -2.222500
      Read-in Center  825 is at   0.000000 -1.534600 -2.222500
      Read-in Center  826 is at   0.000000 -1.481700 -2.222500
      Read-in Center  827 is at   0.000000 -1.428800 -2.222500
      Read-in Center  828 is at   0.000000 -1.375900 -2.222500
      Read-in Center  829 is at   0.000000 -1.322900 -2.222500
      Read-in Center  830 is at   0.000000 -1.270000 -2.222500
      Read-in Center  831 is at   0.000000 -1.217100 -2.222500
      Read-in Center  832 is at   0.000000 -1.164200 -2.222500
      Read-in Center  833 is at   0.000000 -1.111300 -2.222500
      Read-in Center  834 is at   0.000000 -1.058400 -2.222500
      Read-in Center  835 is at   0.000000 -1.005400 -2.222500
      Read-in Center  836 is at   0.000000 -0.952500 -2.222500
      Read-in Center  837 is at   0.000000 -0.899600 -2.222500
      Read-in Center  838 is at   0.000000 -0.846700 -2.222500
      Read-in Center  839 is at   0.000000 -0.793800 -2.222500
      Read-in Center  840 is at   0.000000 -0.740800 -2.222500
      Read-in Center  841 is at   0.000000 -0.687900 -2.222500
      Read-in Center  842 is at   0.000000 -0.635000 -2.222500
      Read-in Center  843 is at   0.000000 -0.582100 -2.222500
      Read-in Center  844 is at   0.000000 -0.529200 -2.222500
      Read-in Center  845 is at   0.000000 -0.476300 -2.222500
      Read-in Center  846 is at   0.000000 -0.423300 -2.222500
      Read-in Center  847 is at   0.000000 -0.370400 -2.222500
      Read-in Center  848 is at   0.000000 -0.317500 -2.222500
      Read-in Center  849 is at   0.000000 -0.264600 -2.222500
      Read-in Center  850 is at   0.000000 -0.211700 -2.222500
      Read-in Center  851 is at   0.000000 -0.158800 -2.222500
      Read-in Center  852 is at   0.000000 -0.105800 -2.222500
      Read-in Center  853 is at   0.000000 -0.052900 -2.222500
      Read-in Center  854 is at   0.000000  0.000000 -2.222500
      Read-in Center  855 is at   0.000000  0.052900 -2.222500
      Read-in Center  856 is at   0.000000  0.105800 -2.222500
      Read-in Center  857 is at   0.000000  0.158800 -2.222500
      Read-in Center  858 is at   0.000000  0.211700 -2.222500
      Read-in Center  859 is at   0.000000  0.264600 -2.222500
      Read-in Center  860 is at   0.000000  0.317500 -2.222500
      Read-in Center  861 is at   0.000000  0.370400 -2.222500
      Read-in Center  862 is at   0.000000  0.423300 -2.222500
      Read-in Center  863 is at   0.000000  0.476300 -2.222500
      Read-in Center  864 is at   0.000000  0.529200 -2.222500
      Read-in Center  865 is at   0.000000  0.582100 -2.222500
      Read-in Center  866 is at   0.000000  0.635000 -2.222500
      Read-in Center  867 is at   0.000000  0.687900 -2.222500
      Read-in Center  868 is at   0.000000  0.740800 -2.222500
      Read-in Center  869 is at   0.000000  0.793800 -2.222500
      Read-in Center  870 is at   0.000000  0.846700 -2.222500
      Read-in Center  871 is at   0.000000  0.899600 -2.222500
      Read-in Center  872 is at   0.000000  0.952500 -2.222500
      Read-in Center  873 is at   0.000000  1.005400 -2.222500
      Read-in Center  874 is at   0.000000  1.058400 -2.222500
      Read-in Center  875 is at   0.000000  1.111300 -2.222500
      Read-in Center  876 is at   0.000000  1.164200 -2.222500
      Read-in Center  877 is at   0.000000  1.217100 -2.222500
      Read-in Center  878 is at   0.000000  1.270000 -2.222500
      Read-in Center  879 is at   0.000000  1.322900 -2.222500
      Read-in Center  880 is at   0.000000  1.375900 -2.222500
      Read-in Center  881 is at   0.000000  1.428800 -2.222500
      Read-in Center  882 is at   0.000000  1.481700 -2.222500
      Read-in Center  883 is at   0.000000  1.534600 -2.222500
      Read-in Center  884 is at   0.000000  1.587500 -2.222500
      Read-in Center  885 is at   0.000000  1.640400 -2.222500
      Read-in Center  886 is at   0.000000  1.693400 -2.222500
      Read-in Center  887 is at   0.000000  1.746300 -2.222500
      Read-in Center  888 is at   0.000000  1.799200 -2.222500
      Read-in Center  889 is at   0.000000  1.852100 -2.222500
      Read-in Center  890 is at   0.000000  1.905000 -2.222500
      Read-in Center  891 is at   0.000000  1.958000 -2.222500
      Read-in Center  892 is at   0.000000  2.010900 -2.222500
      Read-in Center  893 is at   0.000000  2.063800 -2.222500
      Read-in Center  894 is at   0.000000  2.116700 -2.222500
      Read-in Center  895 is at   0.000000  2.169600 -2.222500
      Read-in Center  896 is at   0.000000  2.222500 -2.222500
      Read-in Center  897 is at   0.000000  2.275500 -2.222500
      Read-in Center  898 is at   0.000000  2.328400 -2.222500
      Read-in Center  899 is at   0.000000  2.381300 -2.222500
      Read-in Center  900 is at   0.000000  2.434200 -2.222500
      Read-in Center  901 is at   0.000000  2.487100 -2.222500
      Read-in Center  902 is at   0.000000  2.540000 -2.222500
      Read-in Center  903 is at   0.000000  2.593000 -2.222500
      Read-in Center  904 is at   0.000000 -2.645900 -2.169600
      Read-in Center  905 is at   0.000000 -2.593000 -2.169600
      Read-in Center  906 is at   0.000000 -2.540000 -2.169600
      Read-in Center  907 is at   0.000000 -2.487100 -2.169600
      Read-in Center  908 is at   0.000000 -2.434200 -2.169600
      Read-in Center  909 is at   0.000000 -2.381300 -2.169600
      Read-in Center  910 is at   0.000000 -2.328400 -2.169600
      Read-in Center  911 is at   0.000000 -2.275500 -2.169600
      Read-in Center  912 is at   0.000000 -2.222500 -2.169600
      Read-in Center  913 is at   0.000000 -2.169600 -2.169600
      Read-in Center  914 is at   0.000000 -2.116700 -2.169600
      Read-in Center  915 is at   0.000000 -2.063800 -2.169600
      Read-in Center  916 is at   0.000000 -2.010900 -2.169600
      Read-in Center  917 is at   0.000000 -1.958000 -2.169600
      Read-in Center  918 is at   0.000000 -1.905000 -2.169600
      Read-in Center  919 is at   0.000000 -1.852100 -2.169600
      Read-in Center  920 is at   0.000000 -1.799200 -2.169600
      Read-in Center  921 is at   0.000000 -1.746300 -2.169600
      Read-in Center  922 is at   0.000000 -1.693400 -2.169600
      Read-in Center  923 is at   0.000000 -1.640400 -2.169600
      Read-in Center  924 is at   0.000000 -1.587500 -2.169600
      Read-in Center  925 is at   0.000000 -1.534600 -2.169600
      Read-in Center  926 is at   0.000000 -1.481700 -2.169600
      Read-in Center  927 is at   0.000000 -1.428800 -2.169600
      Read-in Center  928 is at   0.000000 -1.375900 -2.169600
      Read-in Center  929 is at   0.000000 -1.322900 -2.169600
      Read-in Center  930 is at   0.000000 -1.270000 -2.169600
      Read-in Center  931 is at   0.000000 -1.217100 -2.169600
      Read-in Center  932 is at   0.000000 -1.164200 -2.169600
      Read-in Center  933 is at   0.000000 -1.111300 -2.169600
      Read-in Center  934 is at   0.000000 -1.058400 -2.169600
      Read-in Center  935 is at   0.000000 -1.005400 -2.169600
      Read-in Center  936 is at   0.000000 -0.952500 -2.169600
      Read-in Center  937 is at   0.000000 -0.899600 -2.169600
      Read-in Center  938 is at   0.000000 -0.846700 -2.169600
      Read-in Center  939 is at   0.000000 -0.793800 -2.169600
      Read-in Center  940 is at   0.000000 -0.740800 -2.169600
      Read-in Center  941 is at   0.000000 -0.687900 -2.169600
      Read-in Center  942 is at   0.000000 -0.635000 -2.169600
      Read-in Center  943 is at   0.000000 -0.582100 -2.169600
      Read-in Center  944 is at   0.000000 -0.529200 -2.169600
      Read-in Center  945 is at   0.000000 -0.476300 -2.169600
      Read-in Center  946 is at   0.000000 -0.423300 -2.169600
      Read-in Center  947 is at   0.000000 -0.370400 -2.169600
      Read-in Center  948 is at   0.000000 -0.317500 -2.169600
      Read-in Center  949 is at   0.000000 -0.264600 -2.169600
      Read-in Center  950 is at   0.000000 -0.211700 -2.169600
      Read-in Center  951 is at   0.000000 -0.158800 -2.169600
      Read-in Center  952 is at   0.000000 -0.105800 -2.169600
      Read-in Center  953 is at   0.000000 -0.052900 -2.169600
      Read-in Center  954 is at   0.000000  0.000000 -2.169600
      Read-in Center  955 is at   0.000000  0.052900 -2.169600
      Read-in Center  956 is at   0.000000  0.105800 -2.169600
      Read-in Center  957 is at   0.000000  0.158800 -2.169600
      Read-in Center  958 is at   0.000000  0.211700 -2.169600
      Read-in Center  959 is at   0.000000  0.264600 -2.169600
      Read-in Center  960 is at   0.000000  0.317500 -2.169600
      Read-in Center  961 is at   0.000000  0.370400 -2.169600
      Read-in Center  962 is at   0.000000  0.423300 -2.169600
      Read-in Center  963 is at   0.000000  0.476300 -2.169600
      Read-in Center  964 is at   0.000000  0.529200 -2.169600
      Read-in Center  965 is at   0.000000  0.582100 -2.169600
      Read-in Center  966 is at   0.000000  0.635000 -2.169600
      Read-in Center  967 is at   0.000000  0.687900 -2.169600
      Read-in Center  968 is at   0.000000  0.740800 -2.169600
      Read-in Center  969 is at   0.000000  0.793800 -2.169600
      Read-in Center  970 is at   0.000000  0.846700 -2.169600
      Read-in Center  971 is at   0.000000  0.899600 -2.169600
      Read-in Center  972 is at   0.000000  0.952500 -2.169600
      Read-in Center  973 is at   0.000000  1.005400 -2.169600
      Read-in Center  974 is at   0.000000  1.058400 -2.169600
      Read-in Center  975 is at   0.000000  1.111300 -2.169600
      Read-in Center  976 is at   0.000000  1.164200 -2.169600
      Read-in Center  977 is at   0.000000  1.217100 -2.169600
      Read-in Center  978 is at   0.000000  1.270000 -2.169600
      Read-in Center  979 is at   0.000000  1.322900 -2.169600
      Read-in Center  980 is at   0.000000  1.375900 -2.169600
      Read-in Center  981 is at   0.000000  1.428800 -2.169600
      Read-in Center  982 is at   0.000000  1.481700 -2.169600
      Read-in Center  983 is at   0.000000  1.534600 -2.169600
      Read-in Center  984 is at   0.000000  1.587500 -2.169600
      Read-in Center  985 is at   0.000000  1.640400 -2.169600
      Read-in Center  986 is at   0.000000  1.693400 -2.169600
      Read-in Center  987 is at   0.000000  1.746300 -2.169600
      Read-in Center  988 is at   0.000000  1.799200 -2.169600
      Read-in Center  989 is at   0.000000  1.852100 -2.169600
      Read-in Center  990 is at   0.000000  1.905000 -2.169600
      Read-in Center  991 is at   0.000000  1.958000 -2.169600
      Read-in Center  992 is at   0.000000  2.010900 -2.169600
      Read-in Center  993 is at   0.000000  2.063800 -2.169600
      Read-in Center  994 is at   0.000000  2.116700 -2.169600
      Read-in Center  995 is at   0.000000  2.169600 -2.169600
      Read-in Center  996 is at   0.000000  2.222500 -2.169600
      Read-in Center  997 is at   0.000000  2.275500 -2.169600
      Read-in Center  998 is at   0.000000  2.328400 -2.169600
      Read-in Center  999 is at   0.000000  2.381300 -2.169600
      Read-in Center 1000 is at   0.000000  2.434200 -2.169600
      Read-in Center 1001 is at   0.000000  2.487100 -2.169600
      Read-in Center 1002 is at   0.000000  2.540000 -2.169600
      Read-in Center 1003 is at   0.000000  2.593000 -2.169600
      Read-in Center 1004 is at   0.000000 -2.645900 -2.116700
      Read-in Center 1005 is at   0.000000 -2.593000 -2.116700
      Read-in Center 1006 is at   0.000000 -2.540000 -2.116700
      Read-in Center 1007 is at   0.000000 -2.487100 -2.116700
      Read-in Center 1008 is at   0.000000 -2.434200 -2.116700
      Read-in Center 1009 is at   0.000000 -2.381300 -2.116700
      Read-in Center 1010 is at   0.000000 -2.328400 -2.116700
      Read-in Center 1011 is at   0.000000 -2.275500 -2.116700
      Read-in Center 1012 is at   0.000000 -2.222500 -2.116700
      Read-in Center 1013 is at   0.000000 -2.169600 -2.116700
      Read-in Center 1014 is at   0.000000 -2.116700 -2.116700
      Read-in Center 1015 is at   0.000000 -2.063800 -2.116700
      Read-in Center 1016 is at   0.000000 -2.010900 -2.116700
      Read-in Center 1017 is at   0.000000 -1.958000 -2.116700
      Read-in Center 1018 is at   0.000000 -1.905000 -2.116700
      Read-in Center 1019 is at   0.000000 -1.852100 -2.116700
      Read-in Center 1020 is at   0.000000 -1.799200 -2.116700
      Read-in Center 1021 is at   0.000000 -1.746300 -2.116700
      Read-in Center 1022 is at   0.000000 -1.693400 -2.116700
      Read-in Center 1023 is at   0.000000 -1.640400 -2.116700
      Read-in Center 1024 is at   0.000000 -1.587500 -2.116700
      Read-in Center 1025 is at   0.000000 -1.534600 -2.116700
      Read-in Center 1026 is at   0.000000 -1.481700 -2.116700
      Read-in Center 1027 is at   0.000000 -1.428800 -2.116700
      Read-in Center 1028 is at   0.000000 -1.375900 -2.116700
      Read-in Center 1029 is at   0.000000 -1.322900 -2.116700
      Read-in Center 1030 is at   0.000000 -1.270000 -2.116700
      Read-in Center 1031 is at   0.000000 -1.217100 -2.116700
      Read-in Center 1032 is at   0.000000 -1.164200 -2.116700
      Read-in Center 1033 is at   0.000000 -1.111300 -2.116700
      Read-in Center 1034 is at   0.000000 -1.058400 -2.116700
      Read-in Center 1035 is at   0.000000 -1.005400 -2.116700
      Read-in Center 1036 is at   0.000000 -0.952500 -2.116700
      Read-in Center 1037 is at   0.000000 -0.899600 -2.116700
      Read-in Center 1038 is at   0.000000 -0.846700 -2.116700
      Read-in Center 1039 is at   0.000000 -0.793800 -2.116700
      Read-in Center 1040 is at   0.000000 -0.740800 -2.116700
      Read-in Center 1041 is at   0.000000 -0.687900 -2.116700
      Read-in Center 1042 is at   0.000000 -0.635000 -2.116700
      Read-in Center 1043 is at   0.000000 -0.582100 -2.116700
      Read-in Center 1044 is at   0.000000 -0.529200 -2.116700
      Read-in Center 1045 is at   0.000000 -0.476300 -2.116700
      Read-in Center 1046 is at   0.000000 -0.423300 -2.116700
      Read-in Center 1047 is at   0.000000 -0.370400 -2.116700
      Read-in Center 1048 is at   0.000000 -0.317500 -2.116700
      Read-in Center 1049 is at   0.000000 -0.264600 -2.116700
      Read-in Center 1050 is at   0.000000 -0.211700 -2.116700
      Read-in Center 1051 is at   0.000000 -0.158800 -2.116700
      Read-in Center 1052 is at   0.000000 -0.105800 -2.116700
      Read-in Center 1053 is at   0.000000 -0.052900 -2.116700
      Read-in Center 1054 is at   0.000000  0.000000 -2.116700
      Read-in Center 1055 is at   0.000000  0.052900 -2.116700
      Read-in Center 1056 is at   0.000000  0.105800 -2.116700
      Read-in Center 1057 is at   0.000000  0.158800 -2.116700
      Read-in Center 1058 is at   0.000000  0.211700 -2.116700
      Read-in Center 1059 is at   0.000000  0.264600 -2.116700
      Read-in Center 1060 is at   0.000000  0.317500 -2.116700
      Read-in Center 1061 is at   0.000000  0.370400 -2.116700
      Read-in Center 1062 is at   0.000000  0.423300 -2.116700
      Read-in Center 1063 is at   0.000000  0.476300 -2.116700
      Read-in Center 1064 is at   0.000000  0.529200 -2.116700
      Read-in Center 1065 is at   0.000000  0.582100 -2.116700
      Read-in Center 1066 is at   0.000000  0.635000 -2.116700
      Read-in Center 1067 is at   0.000000  0.687900 -2.116700
      Read-in Center 1068 is at   0.000000  0.740800 -2.116700
      Read-in Center 1069 is at   0.000000  0.793800 -2.116700
      Read-in Center 1070 is at   0.000000  0.846700 -2.116700
      Read-in Center 1071 is at   0.000000  0.899600 -2.116700
      Read-in Center 1072 is at   0.000000  0.952500 -2.116700
      Read-in Center 1073 is at   0.000000  1.005400 -2.116700
      Read-in Center 1074 is at   0.000000  1.058400 -2.116700
      Read-in Center 1075 is at   0.000000  1.111300 -2.116700
      Read-in Center 1076 is at   0.000000  1.164200 -2.116700
      Read-in Center 1077 is at   0.000000  1.217100 -2.116700
      Read-in Center 1078 is at   0.000000  1.270000 -2.116700
      Read-in Center 1079 is at   0.000000  1.322900 -2.116700
      Read-in Center 1080 is at   0.000000  1.375900 -2.116700
      Read-in Center 1081 is at   0.000000  1.428800 -2.116700
      Read-in Center 1082 is at   0.000000  1.481700 -2.116700
      Read-in Center 1083 is at   0.000000  1.534600 -2.116700
      Read-in Center 1084 is at   0.000000  1.587500 -2.116700
      Read-in Center 1085 is at   0.000000  1.640400 -2.116700
      Read-in Center 1086 is at   0.000000  1.693400 -2.116700
      Read-in Center 1087 is at   0.000000  1.746300 -2.116700
      Read-in Center 1088 is at   0.000000  1.799200 -2.116700
      Read-in Center 1089 is at   0.000000  1.852100 -2.116700
      Read-in Center 1090 is at   0.000000  1.905000 -2.116700
      Read-in Center 1091 is at   0.000000  1.958000 -2.116700
      Read-in Center 1092 is at   0.000000  2.010900 -2.116700
      Read-in Center 1093 is at   0.000000  2.063800 -2.116700
      Read-in Center 1094 is at   0.000000  2.116700 -2.116700
      Read-in Center 1095 is at   0.000000  2.169600 -2.116700
      Read-in Center 1096 is at   0.000000  2.222500 -2.116700
      Read-in Center 1097 is at   0.000000  2.275500 -2.116700
      Read-in Center 1098 is at   0.000000  2.328400 -2.116700
      Read-in Center 1099 is at   0.000000  2.381300 -2.116700
      Read-in Center 1100 is at   0.000000  2.434200 -2.116700
      Read-in Center 1101 is at   0.000000  2.487100 -2.116700
      Read-in Center 1102 is at   0.000000  2.540000 -2.116700
      Read-in Center 1103 is at   0.000000  2.593000 -2.116700
      Read-in Center 1104 is at   0.000000 -2.645900 -2.063800
      Read-in Center 1105 is at   0.000000 -2.593000 -2.063800
      Read-in Center 1106 is at   0.000000 -2.540000 -2.063800
      Read-in Center 1107 is at   0.000000 -2.487100 -2.063800
      Read-in Center 1108 is at   0.000000 -2.434200 -2.063800
      Read-in Center 1109 is at   0.000000 -2.381300 -2.063800
      Read-in Center 1110 is at   0.000000 -2.328400 -2.063800
      Read-in Center 1111 is at   0.000000 -2.275500 -2.063800
      Read-in Center 1112 is at   0.000000 -2.222500 -2.063800
      Read-in Center 1113 is at   0.000000 -2.169600 -2.063800
      Read-in Center 1114 is at   0.000000 -2.116700 -2.063800
      Read-in Center 1115 is at   0.000000 -2.063800 -2.063800
      Read-in Center 1116 is at   0.000000 -2.010900 -2.063800
      Read-in Center 1117 is at   0.000000 -1.958000 -2.063800
      Read-in Center 1118 is at   0.000000 -1.905000 -2.063800
      Read-in Center 1119 is at   0.000000 -1.852100 -2.063800
      Read-in Center 1120 is at   0.000000 -1.799200 -2.063800
      Read-in Center 1121 is at   0.000000 -1.746300 -2.063800
      Read-in Center 1122 is at   0.000000 -1.693400 -2.063800
      Read-in Center 1123 is at   0.000000 -1.640400 -2.063800
      Read-in Center 1124 is at   0.000000 -1.587500 -2.063800
      Read-in Center 1125 is at   0.000000 -1.534600 -2.063800
      Read-in Center 1126 is at   0.000000 -1.481700 -2.063800
      Read-in Center 1127 is at   0.000000 -1.428800 -2.063800
      Read-in Center 1128 is at   0.000000 -1.375900 -2.063800
      Read-in Center 1129 is at   0.000000 -1.322900 -2.063800
      Read-in Center 1130 is at   0.000000 -1.270000 -2.063800
      Read-in Center 1131 is at   0.000000 -1.217100 -2.063800
      Read-in Center 1132 is at   0.000000 -1.164200 -2.063800
      Read-in Center 1133 is at   0.000000 -1.111300 -2.063800
      Read-in Center 1134 is at   0.000000 -1.058400 -2.063800
      Read-in Center 1135 is at   0.000000 -1.005400 -2.063800
      Read-in Center 1136 is at   0.000000 -0.952500 -2.063800
      Read-in Center 1137 is at   0.000000 -0.899600 -2.063800
      Read-in Center 1138 is at   0.000000 -0.846700 -2.063800
      Read-in Center 1139 is at   0.000000 -0.793800 -2.063800
      Read-in Center 1140 is at   0.000000 -0.740800 -2.063800
      Read-in Center 1141 is at   0.000000 -0.687900 -2.063800
      Read-in Center 1142 is at   0.000000 -0.635000 -2.063800
      Read-in Center 1143 is at   0.000000 -0.582100 -2.063800
      Read-in Center 1144 is at   0.000000 -0.529200 -2.063800
      Read-in Center 1145 is at   0.000000 -0.476300 -2.063800
      Read-in Center 1146 is at   0.000000 -0.423300 -2.063800
      Read-in Center 1147 is at   0.000000 -0.370400 -2.063800
      Read-in Center 1148 is at   0.000000 -0.317500 -2.063800
      Read-in Center 1149 is at   0.000000 -0.264600 -2.063800
      Read-in Center 1150 is at   0.000000 -0.211700 -2.063800
      Read-in Center 1151 is at   0.000000 -0.158800 -2.063800
      Read-in Center 1152 is at   0.000000 -0.105800 -2.063800
      Read-in Center 1153 is at   0.000000 -0.052900 -2.063800
      Read-in Center 1154 is at   0.000000  0.000000 -2.063800
      Read-in Center 1155 is at   0.000000  0.052900 -2.063800
      Read-in Center 1156 is at   0.000000  0.105800 -2.063800
      Read-in Center 1157 is at   0.000000  0.158800 -2.063800
      Read-in Center 1158 is at   0.000000  0.211700 -2.063800
      Read-in Center 1159 is at   0.000000  0.264600 -2.063800
      Read-in Center 1160 is at   0.000000  0.317500 -2.063800
      Read-in Center 1161 is at   0.000000  0.370400 -2.063800
      Read-in Center 1162 is at   0.000000  0.423300 -2.063800
      Read-in Center 1163 is at   0.000000  0.476300 -2.063800
      Read-in Center 1164 is at   0.000000  0.529200 -2.063800
      Read-in Center 1165 is at   0.000000  0.582100 -2.063800
      Read-in Center 1166 is at   0.000000  0.635000 -2.063800
      Read-in Center 1167 is at   0.000000  0.687900 -2.063800
      Read-in Center 1168 is at   0.000000  0.740800 -2.063800
      Read-in Center 1169 is at   0.000000  0.793800 -2.063800
      Read-in Center 1170 is at   0.000000  0.846700 -2.063800
      Read-in Center 1171 is at   0.000000  0.899600 -2.063800
      Read-in Center 1172 is at   0.000000  0.952500 -2.063800
      Read-in Center 1173 is at   0.000000  1.005400 -2.063800
      Read-in Center 1174 is at   0.000000  1.058400 -2.063800
      Read-in Center 1175 is at   0.000000  1.111300 -2.063800
      Read-in Center 1176 is at   0.000000  1.164200 -2.063800
      Read-in Center 1177 is at   0.000000  1.217100 -2.063800
      Read-in Center 1178 is at   0.000000  1.270000 -2.063800
      Read-in Center 1179 is at   0.000000  1.322900 -2.063800
      Read-in Center 1180 is at   0.000000  1.375900 -2.063800
      Read-in Center 1181 is at   0.000000  1.428800 -2.063800
      Read-in Center 1182 is at   0.000000  1.481700 -2.063800
      Read-in Center 1183 is at   0.000000  1.534600 -2.063800
      Read-in Center 1184 is at   0.000000  1.587500 -2.063800
      Read-in Center 1185 is at   0.000000  1.640400 -2.063800
      Read-in Center 1186 is at   0.000000  1.693400 -2.063800
      Read-in Center 1187 is at   0.000000  1.746300 -2.063800
      Read-in Center 1188 is at   0.000000  1.799200 -2.063800
      Read-in Center 1189 is at   0.000000  1.852100 -2.063800
      Read-in Center 1190 is at   0.000000  1.905000 -2.063800
      Read-in Center 1191 is at   0.000000  1.958000 -2.063800
      Read-in Center 1192 is at   0.000000  2.010900 -2.063800
      Read-in Center 1193 is at   0.000000  2.063800 -2.063800
      Read-in Center 1194 is at   0.000000  2.116700 -2.063800
      Read-in Center 1195 is at   0.000000  2.169600 -2.063800
      Read-in Center 1196 is at   0.000000  2.222500 -2.063800
      Read-in Center 1197 is at   0.000000  2.275500 -2.063800
      Read-in Center 1198 is at   0.000000  2.328400 -2.063800
      Read-in Center 1199 is at   0.000000  2.381300 -2.063800
      Read-in Center 1200 is at   0.000000  2.434200 -2.063800
      Read-in Center 1201 is at   0.000000  2.487100 -2.063800
      Read-in Center 1202 is at   0.000000  2.540000 -2.063800
      Read-in Center 1203 is at   0.000000  2.593000 -2.063800
      Read-in Center 1204 is at   0.000000 -2.645900 -2.010900
      Read-in Center 1205 is at   0.000000 -2.593000 -2.010900
      Read-in Center 1206 is at   0.000000 -2.540000 -2.010900
      Read-in Center 1207 is at   0.000000 -2.487100 -2.010900
      Read-in Center 1208 is at   0.000000 -2.434200 -2.010900
      Read-in Center 1209 is at   0.000000 -2.381300 -2.010900
      Read-in Center 1210 is at   0.000000 -2.328400 -2.010900
      Read-in Center 1211 is at   0.000000 -2.275500 -2.010900
      Read-in Center 1212 is at   0.000000 -2.222500 -2.010900
      Read-in Center 1213 is at   0.000000 -2.169600 -2.010900
      Read-in Center 1214 is at   0.000000 -2.116700 -2.010900
      Read-in Center 1215 is at   0.000000 -2.063800 -2.010900
      Read-in Center 1216 is at   0.000000 -2.010900 -2.010900
      Read-in Center 1217 is at   0.000000 -1.958000 -2.010900
      Read-in Center 1218 is at   0.000000 -1.905000 -2.010900
      Read-in Center 1219 is at   0.000000 -1.852100 -2.010900
      Read-in Center 1220 is at   0.000000 -1.799200 -2.010900
      Read-in Center 1221 is at   0.000000 -1.746300 -2.010900
      Read-in Center 1222 is at   0.000000 -1.693400 -2.010900
      Read-in Center 1223 is at   0.000000 -1.640400 -2.010900
      Read-in Center 1224 is at   0.000000 -1.587500 -2.010900
      Read-in Center 1225 is at   0.000000 -1.534600 -2.010900
      Read-in Center 1226 is at   0.000000 -1.481700 -2.010900
      Read-in Center 1227 is at   0.000000 -1.428800 -2.010900
      Read-in Center 1228 is at   0.000000 -1.375900 -2.010900
      Read-in Center 1229 is at   0.000000 -1.322900 -2.010900
      Read-in Center 1230 is at   0.000000 -1.270000 -2.010900
      Read-in Center 1231 is at   0.000000 -1.217100 -2.010900
      Read-in Center 1232 is at   0.000000 -1.164200 -2.010900
      Read-in Center 1233 is at   0.000000 -1.111300 -2.010900
      Read-in Center 1234 is at   0.000000 -1.058400 -2.010900
      Read-in Center 1235 is at   0.000000 -1.005400 -2.010900
      Read-in Center 1236 is at   0.000000 -0.952500 -2.010900
      Read-in Center 1237 is at   0.000000 -0.899600 -2.010900
      Read-in Center 1238 is at   0.000000 -0.846700 -2.010900
      Read-in Center 1239 is at   0.000000 -0.793800 -2.010900
      Read-in Center 1240 is at   0.000000 -0.740800 -2.010900
      Read-in Center 1241 is at   0.000000 -0.687900 -2.010900
      Read-in Center 1242 is at   0.000000 -0.635000 -2.010900
      Read-in Center 1243 is at   0.000000 -0.582100 -2.010900
      Read-in Center 1244 is at   0.000000 -0.529200 -2.010900
      Read-in Center 1245 is at   0.000000 -0.476300 -2.010900
      Read-in Center 1246 is at   0.000000 -0.423300 -2.010900
      Read-in Center 1247 is at   0.000000 -0.370400 -2.010900
      Read-in Center 1248 is at   0.000000 -0.317500 -2.010900
      Read-in Center 1249 is at   0.000000 -0.264600 -2.010900
      Read-in Center 1250 is at   0.000000 -0.211700 -2.010900
      Read-in Center 1251 is at   0.000000 -0.158800 -2.010900
      Read-in Center 1252 is at   0.000000 -0.105800 -2.010900
      Read-in Center 1253 is at   0.000000 -0.052900 -2.010900
      Read-in Center 1254 is at   0.000000  0.000000 -2.010900
      Read-in Center 1255 is at   0.000000  0.052900 -2.010900
      Read-in Center 1256 is at   0.000000  0.105800 -2.010900
      Read-in Center 1257 is at   0.000000  0.158800 -2.010900
      Read-in Center 1258 is at   0.000000  0.211700 -2.010900
      Read-in Center 1259 is at   0.000000  0.264600 -2.010900
      Read-in Center 1260 is at   0.000000  0.317500 -2.010900
      Read-in Center 1261 is at   0.000000  0.370400 -2.010900
      Read-in Center 1262 is at   0.000000  0.423300 -2.010900
      Read-in Center 1263 is at   0.000000  0.476300 -2.010900
      Read-in Center 1264 is at   0.000000  0.529200 -2.010900
      Read-in Center 1265 is at   0.000000  0.582100 -2.010900
      Read-in Center 1266 is at   0.000000  0.635000 -2.010900
      Read-in Center 1267 is at   0.000000  0.687900 -2.010900
      Read-in Center 1268 is at   0.000000  0.740800 -2.010900
      Read-in Center 1269 is at   0.000000  0.793800 -2.010900
      Read-in Center 1270 is at   0.000000  0.846700 -2.010900
      Read-in Center 1271 is at   0.000000  0.899600 -2.010900
      Read-in Center 1272 is at   0.000000  0.952500 -2.010900
      Read-in Center 1273 is at   0.000000  1.005400 -2.010900
      Read-in Center 1274 is at   0.000000  1.058400 -2.010900
      Read-in Center 1275 is at   0.000000  1.111300 -2.010900
      Read-in Center 1276 is at   0.000000  1.164200 -2.010900
      Read-in Center 1277 is at   0.000000  1.217100 -2.010900
      Read-in Center 1278 is at   0.000000  1.270000 -2.010900
      Read-in Center 1279 is at   0.000000  1.322900 -2.010900
      Read-in Center 1280 is at   0.000000  1.375900 -2.010900
      Read-in Center 1281 is at   0.000000  1.428800 -2.010900
      Read-in Center 1282 is at   0.000000  1.481700 -2.010900
      Read-in Center 1283 is at   0.000000  1.534600 -2.010900
      Read-in Center 1284 is at   0.000000  1.587500 -2.010900
      Read-in Center 1285 is at   0.000000  1.640400 -2.010900
      Read-in Center 1286 is at   0.000000  1.693400 -2.010900
      Read-in Center 1287 is at   0.000000  1.746300 -2.010900
      Read-in Center 1288 is at   0.000000  1.799200 -2.010900
      Read-in Center 1289 is at   0.000000  1.852100 -2.010900
      Read-in Center 1290 is at   0.000000  1.905000 -2.010900
      Read-in Center 1291 is at   0.000000  1.958000 -2.010900
      Read-in Center 1292 is at   0.000000  2.010900 -2.010900
      Read-in Center 1293 is at   0.000000  2.063800 -2.010900
      Read-in Center 1294 is at   0.000000  2.116700 -2.010900
      Read-in Center 1295 is at   0.000000  2.169600 -2.010900
      Read-in Center 1296 is at   0.000000  2.222500 -2.010900
      Read-in Center 1297 is at   0.000000  2.275500 -2.010900
      Read-in Center 1298 is at   0.000000  2.328400 -2.010900
      Read-in Center 1299 is at   0.000000  2.381300 -2.010900
      Read-in Center 1300 is at   0.000000  2.434200 -2.010900
      Read-in Center 1301 is at   0.000000  2.487100 -2.010900
      Read-in Center 1302 is at   0.000000  2.540000 -2.010900
      Read-in Center 1303 is at   0.000000  2.593000 -2.010900
      Read-in Center 1304 is at   0.000000 -2.645900 -1.958000
      Read-in Center 1305 is at   0.000000 -2.593000 -1.958000
      Read-in Center 1306 is at   0.000000 -2.540000 -1.958000
      Read-in Center 1307 is at   0.000000 -2.487100 -1.958000
      Read-in Center 1308 is at   0.000000 -2.434200 -1.958000
      Read-in Center 1309 is at   0.000000 -2.381300 -1.958000
      Read-in Center 1310 is at   0.000000 -2.328400 -1.958000
      Read-in Center 1311 is at   0.000000 -2.275500 -1.958000
      Read-in Center 1312 is at   0.000000 -2.222500 -1.958000
      Read-in Center 1313 is at   0.000000 -2.169600 -1.958000
      Read-in Center 1314 is at   0.000000 -2.116700 -1.958000
      Read-in Center 1315 is at   0.000000 -2.063800 -1.958000
      Read-in Center 1316 is at   0.000000 -2.010900 -1.958000
      Read-in Center 1317 is at   0.000000 -1.958000 -1.958000
      Read-in Center 1318 is at   0.000000 -1.905000 -1.958000
      Read-in Center 1319 is at   0.000000 -1.852100 -1.958000
      Read-in Center 1320 is at   0.000000 -1.799200 -1.958000
      Read-in Center 1321 is at   0.000000 -1.746300 -1.958000
      Read-in Center 1322 is at   0.000000 -1.693400 -1.958000
      Read-in Center 1323 is at   0.000000 -1.640400 -1.958000
      Read-in Center 1324 is at   0.000000 -1.587500 -1.958000
      Read-in Center 1325 is at   0.000000 -1.534600 -1.958000
      Read-in Center 1326 is at   0.000000 -1.481700 -1.958000
      Read-in Center 1327 is at   0.000000 -1.428800 -1.958000
      Read-in Center 1328 is at   0.000000 -1.375900 -1.958000
      Read-in Center 1329 is at   0.000000 -1.322900 -1.958000
      Read-in Center 1330 is at   0.000000 -1.270000 -1.958000
      Read-in Center 1331 is at   0.000000 -1.217100 -1.958000
      Read-in Center 1332 is at   0.000000 -1.164200 -1.958000
      Read-in Center 1333 is at   0.000000 -1.111300 -1.958000
      Read-in Center 1334 is at   0.000000 -1.058400 -1.958000
      Read-in Center 1335 is at   0.000000 -1.005400 -1.958000
      Read-in Center 1336 is at   0.000000 -0.952500 -1.958000
      Read-in Center 1337 is at   0.000000 -0.899600 -1.958000
      Read-in Center 1338 is at   0.000000 -0.846700 -1.958000
      Read-in Center 1339 is at   0.000000 -0.793800 -1.958000
      Read-in Center 1340 is at   0.000000 -0.740800 -1.958000
      Read-in Center 1341 is at   0.000000 -0.687900 -1.958000
      Read-in Center 1342 is at   0.000000 -0.635000 -1.958000
      Read-in Center 1343 is at   0.000000 -0.582100 -1.958000
      Read-in Center 1344 is at   0.000000 -0.529200 -1.958000
      Read-in Center 1345 is at   0.000000 -0.476300 -1.958000
      Read-in Center 1346 is at   0.000000 -0.423300 -1.958000
      Read-in Center 1347 is at   0.000000 -0.370400 -1.958000
      Read-in Center 1348 is at   0.000000 -0.317500 -1.958000
      Read-in Center 1349 is at   0.000000 -0.264600 -1.958000
      Read-in Center 1350 is at   0.000000 -0.211700 -1.958000
      Read-in Center 1351 is at   0.000000 -0.158800 -1.958000
      Read-in Center 1352 is at   0.000000 -0.105800 -1.958000
      Read-in Center 1353 is at   0.000000 -0.052900 -1.958000
      Read-in Center 1354 is at   0.000000  0.000000 -1.958000
      Read-in Center 1355 is at   0.000000  0.052900 -1.958000
      Read-in Center 1356 is at   0.000000  0.105800 -1.958000
      Read-in Center 1357 is at   0.000000  0.158800 -1.958000
      Read-in Center 1358 is at   0.000000  0.211700 -1.958000
      Read-in Center 1359 is at   0.000000  0.264600 -1.958000
      Read-in Center 1360 is at   0.000000  0.317500 -1.958000
      Read-in Center 1361 is at   0.000000  0.370400 -1.958000
      Read-in Center 1362 is at   0.000000  0.423300 -1.958000
      Read-in Center 1363 is at   0.000000  0.476300 -1.958000
      Read-in Center 1364 is at   0.000000  0.529200 -1.958000
      Read-in Center 1365 is at   0.000000  0.582100 -1.958000
      Read-in Center 1366 is at   0.000000  0.635000 -1.958000
      Read-in Center 1367 is at   0.000000  0.687900 -1.958000
      Read-in Center 1368 is at   0.000000  0.740800 -1.958000
      Read-in Center 1369 is at   0.000000  0.793800 -1.958000
      Read-in Center 1370 is at   0.000000  0.846700 -1.958000
      Read-in Center 1371 is at   0.000000  0.899600 -1.958000
      Read-in Center 1372 is at   0.000000  0.952500 -1.958000
      Read-in Center 1373 is at   0.000000  1.005400 -1.958000
      Read-in Center 1374 is at   0.000000  1.058400 -1.958000
      Read-in Center 1375 is at   0.000000  1.111300 -1.958000
      Read-in Center 1376 is at   0.000000  1.164200 -1.958000
      Read-in Center 1377 is at   0.000000  1.217100 -1.958000
      Read-in Center 1378 is at   0.000000  1.270000 -1.958000
      Read-in Center 1379 is at   0.000000  1.322900 -1.958000
      Read-in Center 1380 is at   0.000000  1.375900 -1.958000
      Read-in Center 1381 is at   0.000000  1.428800 -1.958000
      Read-in Center 1382 is at   0.000000  1.481700 -1.958000
      Read-in Center 1383 is at   0.000000  1.534600 -1.958000
      Read-in Center 1384 is at   0.000000  1.587500 -1.958000
      Read-in Center 1385 is at   0.000000  1.640400 -1.958000
      Read-in Center 1386 is at   0.000000  1.693400 -1.958000
      Read-in Center 1387 is at   0.000000  1.746300 -1.958000
      Read-in Center 1388 is at   0.000000  1.799200 -1.958000
      Read-in Center 1389 is at   0.000000  1.852100 -1.958000
      Read-in Center 1390 is at   0.000000  1.905000 -1.958000
      Read-in Center 1391 is at   0.000000  1.958000 -1.958000
      Read-in Center 1392 is at   0.000000  2.010900 -1.958000
      Read-in Center 1393 is at   0.000000  2.063800 -1.958000
      Read-in Center 1394 is at   0.000000  2.116700 -1.958000
      Read-in Center 1395 is at   0.000000  2.169600 -1.958000
      Read-in Center 1396 is at   0.000000  2.222500 -1.958000
      Read-in Center 1397 is at   0.000000  2.275500 -1.958000
      Read-in Center 1398 is at   0.000000  2.328400 -1.958000
      Read-in Center 1399 is at   0.000000  2.381300 -1.958000
      Read-in Center 1400 is at   0.000000  2.434200 -1.958000
      Read-in Center 1401 is at   0.000000  2.487100 -1.958000
      Read-in Center 1402 is at   0.000000  2.540000 -1.958000
      Read-in Center 1403 is at   0.000000  2.593000 -1.958000
      Read-in Center 1404 is at   0.000000 -2.645900 -1.905000
      Read-in Center 1405 is at   0.000000 -2.593000 -1.905000
      Read-in Center 1406 is at   0.000000 -2.540000 -1.905000
      Read-in Center 1407 is at   0.000000 -2.487100 -1.905000
      Read-in Center 1408 is at   0.000000 -2.434200 -1.905000
      Read-in Center 1409 is at   0.000000 -2.381300 -1.905000
      Read-in Center 1410 is at   0.000000 -2.328400 -1.905000
      Read-in Center 1411 is at   0.000000 -2.275500 -1.905000
      Read-in Center 1412 is at   0.000000 -2.222500 -1.905000
      Read-in Center 1413 is at   0.000000 -2.169600 -1.905000
      Read-in Center 1414 is at   0.000000 -2.116700 -1.905000
      Read-in Center 1415 is at   0.000000 -2.063800 -1.905000
      Read-in Center 1416 is at   0.000000 -2.010900 -1.905000
      Read-in Center 1417 is at   0.000000 -1.958000 -1.905000
      Read-in Center 1418 is at   0.000000 -1.905000 -1.905000
      Read-in Center 1419 is at   0.000000 -1.852100 -1.905000
      Read-in Center 1420 is at   0.000000 -1.799200 -1.905000
      Read-in Center 1421 is at   0.000000 -1.746300 -1.905000
      Read-in Center 1422 is at   0.000000 -1.693400 -1.905000
      Read-in Center 1423 is at   0.000000 -1.640400 -1.905000
      Read-in Center 1424 is at   0.000000 -1.587500 -1.905000
      Read-in Center 1425 is at   0.000000 -1.534600 -1.905000
      Read-in Center 1426 is at   0.000000 -1.481700 -1.905000
      Read-in Center 1427 is at   0.000000 -1.428800 -1.905000
      Read-in Center 1428 is at   0.000000 -1.375900 -1.905000
      Read-in Center 1429 is at   0.000000 -1.322900 -1.905000
      Read-in Center 1430 is at   0.000000 -1.270000 -1.905000
      Read-in Center 1431 is at   0.000000 -1.217100 -1.905000
      Read-in Center 1432 is at   0.000000 -1.164200 -1.905000
      Read-in Center 1433 is at   0.000000 -1.111300 -1.905000
      Read-in Center 1434 is at   0.000000 -1.058400 -1.905000
      Read-in Center 1435 is at   0.000000 -1.005400 -1.905000
      Read-in Center 1436 is at   0.000000 -0.952500 -1.905000
      Read-in Center 1437 is at   0.000000 -0.899600 -1.905000
      Read-in Center 1438 is at   0.000000 -0.846700 -1.905000
      Read-in Center 1439 is at   0.000000 -0.793800 -1.905000
      Read-in Center 1440 is at   0.000000 -0.740800 -1.905000
      Read-in Center 1441 is at   0.000000 -0.687900 -1.905000
      Read-in Center 1442 is at   0.000000 -0.635000 -1.905000
      Read-in Center 1443 is at   0.000000 -0.582100 -1.905000
      Read-in Center 1444 is at   0.000000 -0.529200 -1.905000
      Read-in Center 1445 is at   0.000000 -0.476300 -1.905000
      Read-in Center 1446 is at   0.000000 -0.423300 -1.905000
      Read-in Center 1447 is at   0.000000 -0.370400 -1.905000
      Read-in Center 1448 is at   0.000000 -0.317500 -1.905000
      Read-in Center 1449 is at   0.000000 -0.264600 -1.905000
      Read-in Center 1450 is at   0.000000 -0.211700 -1.905000
      Read-in Center 1451 is at   0.000000 -0.158800 -1.905000
      Read-in Center 1452 is at   0.000000 -0.105800 -1.905000
      Read-in Center 1453 is at   0.000000 -0.052900 -1.905000
      Read-in Center 1454 is at   0.000000  0.000000 -1.905000
      Read-in Center 1455 is at   0.000000  0.052900 -1.905000
      Read-in Center 1456 is at   0.000000  0.105800 -1.905000
      Read-in Center 1457 is at   0.000000  0.158800 -1.905000
      Read-in Center 1458 is at   0.000000  0.211700 -1.905000
      Read-in Center 1459 is at   0.000000  0.264600 -1.905000
      Read-in Center 1460 is at   0.000000  0.317500 -1.905000
      Read-in Center 1461 is at   0.000000  0.370400 -1.905000
      Read-in Center 1462 is at   0.000000  0.423300 -1.905000
      Read-in Center 1463 is at   0.000000  0.476300 -1.905000
      Read-in Center 1464 is at   0.000000  0.529200 -1.905000
      Read-in Center 1465 is at   0.000000  0.582100 -1.905000
      Read-in Center 1466 is at   0.000000  0.635000 -1.905000
      Read-in Center 1467 is at   0.000000  0.687900 -1.905000
      Read-in Center 1468 is at   0.000000  0.740800 -1.905000
      Read-in Center 1469 is at   0.000000  0.793800 -1.905000
      Read-in Center 1470 is at   0.000000  0.846700 -1.905000
      Read-in Center 1471 is at   0.000000  0.899600 -1.905000
      Read-in Center 1472 is at   0.000000  0.952500 -1.905000
      Read-in Center 1473 is at   0.000000  1.005400 -1.905000
      Read-in Center 1474 is at   0.000000  1.058400 -1.905000
      Read-in Center 1475 is at   0.000000  1.111300 -1.905000
      Read-in Center 1476 is at   0.000000  1.164200 -1.905000
      Read-in Center 1477 is at   0.000000  1.217100 -1.905000
      Read-in Center 1478 is at   0.000000  1.270000 -1.905000
      Read-in Center 1479 is at   0.000000  1.322900 -1.905000
      Read-in Center 1480 is at   0.000000  1.375900 -1.905000
      Read-in Center 1481 is at   0.000000  1.428800 -1.905000
      Read-in Center 1482 is at   0.000000  1.481700 -1.905000
      Read-in Center 1483 is at   0.000000  1.534600 -1.905000
      Read-in Center 1484 is at   0.000000  1.587500 -1.905000
      Read-in Center 1485 is at   0.000000  1.640400 -1.905000
      Read-in Center 1486 is at   0.000000  1.693400 -1.905000
      Read-in Center 1487 is at   0.000000  1.746300 -1.905000
      Read-in Center 1488 is at   0.000000  1.799200 -1.905000
      Read-in Center 1489 is at   0.000000  1.852100 -1.905000
      Read-in Center 1490 is at   0.000000  1.905000 -1.905000
      Read-in Center 1491 is at   0.000000  1.958000 -1.905000
      Read-in Center 1492 is at   0.000000  2.010900 -1.905000
      Read-in Center 1493 is at   0.000000  2.063800 -1.905000
      Read-in Center 1494 is at   0.000000  2.116700 -1.905000
      Read-in Center 1495 is at   0.000000  2.169600 -1.905000
      Read-in Center 1496 is at   0.000000  2.222500 -1.905000
      Read-in Center 1497 is at   0.000000  2.275500 -1.905000
      Read-in Center 1498 is at   0.000000  2.328400 -1.905000
      Read-in Center 1499 is at   0.000000  2.381300 -1.905000
      Read-in Center 1500 is at   0.000000  2.434200 -1.905000
      Read-in Center 1501 is at   0.000000  2.487100 -1.905000
      Read-in Center 1502 is at   0.000000  2.540000 -1.905000
      Read-in Center 1503 is at   0.000000  2.593000 -1.905000
      Read-in Center 1504 is at   0.000000 -2.645900 -1.852100
      Read-in Center 1505 is at   0.000000 -2.593000 -1.852100
      Read-in Center 1506 is at   0.000000 -2.540000 -1.852100
      Read-in Center 1507 is at   0.000000 -2.487100 -1.852100
      Read-in Center 1508 is at   0.000000 -2.434200 -1.852100
      Read-in Center 1509 is at   0.000000 -2.381300 -1.852100
      Read-in Center 1510 is at   0.000000 -2.328400 -1.852100
      Read-in Center 1511 is at   0.000000 -2.275500 -1.852100
      Read-in Center 1512 is at   0.000000 -2.222500 -1.852100
      Read-in Center 1513 is at   0.000000 -2.169600 -1.852100
      Read-in Center 1514 is at   0.000000 -2.116700 -1.852100
      Read-in Center 1515 is at   0.000000 -2.063800 -1.852100
      Read-in Center 1516 is at   0.000000 -2.010900 -1.852100
      Read-in Center 1517 is at   0.000000 -1.958000 -1.852100
      Read-in Center 1518 is at   0.000000 -1.905000 -1.852100
      Read-in Center 1519 is at   0.000000 -1.852100 -1.852100
      Read-in Center 1520 is at   0.000000 -1.799200 -1.852100
      Read-in Center 1521 is at   0.000000 -1.746300 -1.852100
      Read-in Center 1522 is at   0.000000 -1.693400 -1.852100
      Read-in Center 1523 is at   0.000000 -1.640400 -1.852100
      Read-in Center 1524 is at   0.000000 -1.587500 -1.852100
      Read-in Center 1525 is at   0.000000 -1.534600 -1.852100
      Read-in Center 1526 is at   0.000000 -1.481700 -1.852100
      Read-in Center 1527 is at   0.000000 -1.428800 -1.852100
      Read-in Center 1528 is at   0.000000 -1.375900 -1.852100
      Read-in Center 1529 is at   0.000000 -1.322900 -1.852100
      Read-in Center 1530 is at   0.000000 -1.270000 -1.852100
      Read-in Center 1531 is at   0.000000 -1.217100 -1.852100
      Read-in Center 1532 is at   0.000000 -1.164200 -1.852100
      Read-in Center 1533 is at   0.000000 -1.111300 -1.852100
      Read-in Center 1534 is at   0.000000 -1.058400 -1.852100
      Read-in Center 1535 is at   0.000000 -1.005400 -1.852100
      Read-in Center 1536 is at   0.000000 -0.952500 -1.852100
      Read-in Center 1537 is at   0.000000 -0.899600 -1.852100
      Read-in Center 1538 is at   0.000000 -0.846700 -1.852100
      Read-in Center 1539 is at   0.000000 -0.793800 -1.852100
      Read-in Center 1540 is at   0.000000 -0.740800 -1.852100
      Read-in Center 1541 is at   0.000000 -0.687900 -1.852100
      Read-in Center 1542 is at   0.000000 -0.635000 -1.852100
      Read-in Center 1543 is at   0.000000 -0.582100 -1.852100
      Read-in Center 1544 is at   0.000000 -0.529200 -1.852100
      Read-in Center 1545 is at   0.000000 -0.476300 -1.852100
      Read-in Center 1546 is at   0.000000 -0.423300 -1.852100
      Read-in Center 1547 is at   0.000000 -0.370400 -1.852100
      Read-in Center 1548 is at   0.000000 -0.317500 -1.852100
      Read-in Center 1549 is at   0.000000 -0.264600 -1.852100
      Read-in Center 1550 is at   0.000000 -0.211700 -1.852100
      Read-in Center 1551 is at   0.000000 -0.158800 -1.852100
      Read-in Center 1552 is at   0.000000 -0.105800 -1.852100
      Read-in Center 1553 is at   0.000000 -0.052900 -1.852100
      Read-in Center 1554 is at   0.000000  0.000000 -1.852100
      Read-in Center 1555 is at   0.000000  0.052900 -1.852100
      Read-in Center 1556 is at   0.000000  0.105800 -1.852100
      Read-in Center 1557 is at   0.000000  0.158800 -1.852100
      Read-in Center 1558 is at   0.000000  0.211700 -1.852100
      Read-in Center 1559 is at   0.000000  0.264600 -1.852100
      Read-in Center 1560 is at   0.000000  0.317500 -1.852100
      Read-in Center 1561 is at   0.000000  0.370400 -1.852100
      Read-in Center 1562 is at   0.000000  0.423300 -1.852100
      Read-in Center 1563 is at   0.000000  0.476300 -1.852100
      Read-in Center 1564 is at   0.000000  0.529200 -1.852100
      Read-in Center 1565 is at   0.000000  0.582100 -1.852100
      Read-in Center 1566 is at   0.000000  0.635000 -1.852100
      Read-in Center 1567 is at   0.000000  0.687900 -1.852100
      Read-in Center 1568 is at   0.000000  0.740800 -1.852100
      Read-in Center 1569 is at   0.000000  0.793800 -1.852100
      Read-in Center 1570 is at   0.000000  0.846700 -1.852100
      Read-in Center 1571 is at   0.000000  0.899600 -1.852100
      Read-in Center 1572 is at   0.000000  0.952500 -1.852100
      Read-in Center 1573 is at   0.000000  1.005400 -1.852100
      Read-in Center 1574 is at   0.000000  1.058400 -1.852100
      Read-in Center 1575 is at   0.000000  1.111300 -1.852100
      Read-in Center 1576 is at   0.000000  1.164200 -1.852100
      Read-in Center 1577 is at   0.000000  1.217100 -1.852100
      Read-in Center 1578 is at   0.000000  1.270000 -1.852100
      Read-in Center 1579 is at   0.000000  1.322900 -1.852100
      Read-in Center 1580 is at   0.000000  1.375900 -1.852100
      Read-in Center 1581 is at   0.000000  1.428800 -1.852100
      Read-in Center 1582 is at   0.000000  1.481700 -1.852100
      Read-in Center 1583 is at   0.000000  1.534600 -1.852100
      Read-in Center 1584 is at   0.000000  1.587500 -1.852100
      Read-in Center 1585 is at   0.000000  1.640400 -1.852100
      Read-in Center 1586 is at   0.000000  1.693400 -1.852100
      Read-in Center 1587 is at   0.000000  1.746300 -1.852100
      Read-in Center 1588 is at   0.000000  1.799200 -1.852100
      Read-in Center 1589 is at   0.000000  1.852100 -1.852100
      Read-in Center 1590 is at   0.000000  1.905000 -1.852100
      Read-in Center 1591 is at   0.000000  1.958000 -1.852100
      Read-in Center 1592 is at   0.000000  2.010900 -1.852100
      Read-in Center 1593 is at   0.000000  2.063800 -1.852100
      Read-in Center 1594 is at   0.000000  2.116700 -1.852100
      Read-in Center 1595 is at   0.000000  2.169600 -1.852100
      Read-in Center 1596 is at   0.000000  2.222500 -1.852100
      Read-in Center 1597 is at   0.000000  2.275500 -1.852100
      Read-in Center 1598 is at   0.000000  2.328400 -1.852100
      Read-in Center 1599 is at   0.000000  2.381300 -1.852100
      Read-in Center 1600 is at   0.000000  2.434200 -1.852100
      Read-in Center 1601 is at   0.000000  2.487100 -1.852100
      Read-in Center 1602 is at   0.000000  2.540000 -1.852100
      Read-in Center 1603 is at   0.000000  2.593000 -1.852100
      Read-in Center 1604 is at   0.000000 -2.645900 -1.799200
      Read-in Center 1605 is at   0.000000 -2.593000 -1.799200
      Read-in Center 1606 is at   0.000000 -2.540000 -1.799200
      Read-in Center 1607 is at   0.000000 -2.487100 -1.799200
      Read-in Center 1608 is at   0.000000 -2.434200 -1.799200
      Read-in Center 1609 is at   0.000000 -2.381300 -1.799200
      Read-in Center 1610 is at   0.000000 -2.328400 -1.799200
      Read-in Center 1611 is at   0.000000 -2.275500 -1.799200
      Read-in Center 1612 is at   0.000000 -2.222500 -1.799200
      Read-in Center 1613 is at   0.000000 -2.169600 -1.799200
      Read-in Center 1614 is at   0.000000 -2.116700 -1.799200
      Read-in Center 1615 is at   0.000000 -2.063800 -1.799200
      Read-in Center 1616 is at   0.000000 -2.010900 -1.799200
      Read-in Center 1617 is at   0.000000 -1.958000 -1.799200
      Read-in Center 1618 is at   0.000000 -1.905000 -1.799200
      Read-in Center 1619 is at   0.000000 -1.852100 -1.799200
      Read-in Center 1620 is at   0.000000 -1.799200 -1.799200
      Read-in Center 1621 is at   0.000000 -1.746300 -1.799200
      Read-in Center 1622 is at   0.000000 -1.693400 -1.799200
      Read-in Center 1623 is at   0.000000 -1.640400 -1.799200
      Read-in Center 1624 is at   0.000000 -1.587500 -1.799200
      Read-in Center 1625 is at   0.000000 -1.534600 -1.799200
      Read-in Center 1626 is at   0.000000 -1.481700 -1.799200
      Read-in Center 1627 is at   0.000000 -1.428800 -1.799200
      Read-in Center 1628 is at   0.000000 -1.375900 -1.799200
      Read-in Center 1629 is at   0.000000 -1.322900 -1.799200
      Read-in Center 1630 is at   0.000000 -1.270000 -1.799200
      Read-in Center 1631 is at   0.000000 -1.217100 -1.799200
      Read-in Center 1632 is at   0.000000 -1.164200 -1.799200
      Read-in Center 1633 is at   0.000000 -1.111300 -1.799200
      Read-in Center 1634 is at   0.000000 -1.058400 -1.799200
      Read-in Center 1635 is at   0.000000 -1.005400 -1.799200
      Read-in Center 1636 is at   0.000000 -0.952500 -1.799200
      Read-in Center 1637 is at   0.000000 -0.899600 -1.799200
      Read-in Center 1638 is at   0.000000 -0.846700 -1.799200
      Read-in Center 1639 is at   0.000000 -0.793800 -1.799200
      Read-in Center 1640 is at   0.000000 -0.740800 -1.799200
      Read-in Center 1641 is at   0.000000 -0.687900 -1.799200
      Read-in Center 1642 is at   0.000000 -0.635000 -1.799200
      Read-in Center 1643 is at   0.000000 -0.582100 -1.799200
      Read-in Center 1644 is at   0.000000 -0.529200 -1.799200
      Read-in Center 1645 is at   0.000000 -0.476300 -1.799200
      Read-in Center 1646 is at   0.000000 -0.423300 -1.799200
      Read-in Center 1647 is at   0.000000 -0.370400 -1.799200
      Read-in Center 1648 is at   0.000000 -0.317500 -1.799200
      Read-in Center 1649 is at   0.000000 -0.264600 -1.799200
      Read-in Center 1650 is at   0.000000 -0.211700 -1.799200
      Read-in Center 1651 is at   0.000000 -0.158800 -1.799200
      Read-in Center 1652 is at   0.000000 -0.105800 -1.799200
      Read-in Center 1653 is at   0.000000 -0.052900 -1.799200
      Read-in Center 1654 is at   0.000000  0.000000 -1.799200
      Read-in Center 1655 is at   0.000000  0.052900 -1.799200
      Read-in Center 1656 is at   0.000000  0.105800 -1.799200
      Read-in Center 1657 is at   0.000000  0.158800 -1.799200
      Read-in Center 1658 is at   0.000000  0.211700 -1.799200
      Read-in Center 1659 is at   0.000000  0.264600 -1.799200
      Read-in Center 1660 is at   0.000000  0.317500 -1.799200
      Read-in Center 1661 is at   0.000000  0.370400 -1.799200
      Read-in Center 1662 is at   0.000000  0.423300 -1.799200
      Read-in Center 1663 is at   0.000000  0.476300 -1.799200
      Read-in Center 1664 is at   0.000000  0.529200 -1.799200
      Read-in Center 1665 is at   0.000000  0.582100 -1.799200
      Read-in Center 1666 is at   0.000000  0.635000 -1.799200
      Read-in Center 1667 is at   0.000000  0.687900 -1.799200
      Read-in Center 1668 is at   0.000000  0.740800 -1.799200
      Read-in Center 1669 is at   0.000000  0.793800 -1.799200
      Read-in Center 1670 is at   0.000000  0.846700 -1.799200
      Read-in Center 1671 is at   0.000000  0.899600 -1.799200
      Read-in Center 1672 is at   0.000000  0.952500 -1.799200
      Read-in Center 1673 is at   0.000000  1.005400 -1.799200
      Read-in Center 1674 is at   0.000000  1.058400 -1.799200
      Read-in Center 1675 is at   0.000000  1.111300 -1.799200
      Read-in Center 1676 is at   0.000000  1.164200 -1.799200
      Read-in Center 1677 is at   0.000000  1.217100 -1.799200
      Read-in Center 1678 is at   0.000000  1.270000 -1.799200
      Read-in Center 1679 is at   0.000000  1.322900 -1.799200
      Read-in Center 1680 is at   0.000000  1.375900 -1.799200
      Read-in Center 1681 is at   0.000000  1.428800 -1.799200
      Read-in Center 1682 is at   0.000000  1.481700 -1.799200
      Read-in Center 1683 is at   0.000000  1.534600 -1.799200
      Read-in Center 1684 is at   0.000000  1.587500 -1.799200
      Read-in Center 1685 is at   0.000000  1.640400 -1.799200
      Read-in Center 1686 is at   0.000000  1.693400 -1.799200
      Read-in Center 1687 is at   0.000000  1.746300 -1.799200
      Read-in Center 1688 is at   0.000000  1.799200 -1.799200
      Read-in Center 1689 is at   0.000000  1.852100 -1.799200
      Read-in Center 1690 is at   0.000000  1.905000 -1.799200
      Read-in Center 1691 is at   0.000000  1.958000 -1.799200
      Read-in Center 1692 is at   0.000000  2.010900 -1.799200
      Read-in Center 1693 is at   0.000000  2.063800 -1.799200
      Read-in Center 1694 is at   0.000000  2.116700 -1.799200
      Read-in Center 1695 is at   0.000000  2.169600 -1.799200
      Read-in Center 1696 is at   0.000000  2.222500 -1.799200
      Read-in Center 1697 is at   0.000000  2.275500 -1.799200
      Read-in Center 1698 is at   0.000000  2.328400 -1.799200
      Read-in Center 1699 is at   0.000000  2.381300 -1.799200
      Read-in Center 1700 is at   0.000000  2.434200 -1.799200
      Read-in Center 1701 is at   0.000000  2.487100 -1.799200
      Read-in Center 1702 is at   0.000000  2.540000 -1.799200
      Read-in Center 1703 is at   0.000000  2.593000 -1.799200
      Read-in Center 1704 is at   0.000000 -2.645900 -1.746300
      Read-in Center 1705 is at   0.000000 -2.593000 -1.746300
      Read-in Center 1706 is at   0.000000 -2.540000 -1.746300
      Read-in Center 1707 is at   0.000000 -2.487100 -1.746300
      Read-in Center 1708 is at   0.000000 -2.434200 -1.746300
      Read-in Center 1709 is at   0.000000 -2.381300 -1.746300
      Read-in Center 1710 is at   0.000000 -2.328400 -1.746300
      Read-in Center 1711 is at   0.000000 -2.275500 -1.746300
      Read-in Center 1712 is at   0.000000 -2.222500 -1.746300
      Read-in Center 1713 is at   0.000000 -2.169600 -1.746300
      Read-in Center 1714 is at   0.000000 -2.116700 -1.746300
      Read-in Center 1715 is at   0.000000 -2.063800 -1.746300
      Read-in Center 1716 is at   0.000000 -2.010900 -1.746300
      Read-in Center 1717 is at   0.000000 -1.958000 -1.746300
      Read-in Center 1718 is at   0.000000 -1.905000 -1.746300
      Read-in Center 1719 is at   0.000000 -1.852100 -1.746300
      Read-in Center 1720 is at   0.000000 -1.799200 -1.746300
      Read-in Center 1721 is at   0.000000 -1.746300 -1.746300
      Read-in Center 1722 is at   0.000000 -1.693400 -1.746300
      Read-in Center 1723 is at   0.000000 -1.640400 -1.746300
      Read-in Center 1724 is at   0.000000 -1.587500 -1.746300
      Read-in Center 1725 is at   0.000000 -1.534600 -1.746300
      Read-in Center 1726 is at   0.000000 -1.481700 -1.746300
      Read-in Center 1727 is at   0.000000 -1.428800 -1.746300
      Read-in Center 1728 is at   0.000000 -1.375900 -1.746300
      Read-in Center 1729 is at   0.000000 -1.322900 -1.746300
      Read-in Center 1730 is at   0.000000 -1.270000 -1.746300
      Read-in Center 1731 is at   0.000000 -1.217100 -1.746300
      Read-in Center 1732 is at   0.000000 -1.164200 -1.746300
      Read-in Center 1733 is at   0.000000 -1.111300 -1.746300
      Read-in Center 1734 is at   0.000000 -1.058400 -1.746300
      Read-in Center 1735 is at   0.000000 -1.005400 -1.746300
      Read-in Center 1736 is at   0.000000 -0.952500 -1.746300
      Read-in Center 1737 is at   0.000000 -0.899600 -1.746300
      Read-in Center 1738 is at   0.000000 -0.846700 -1.746300
      Read-in Center 1739 is at   0.000000 -0.793800 -1.746300
      Read-in Center 1740 is at   0.000000 -0.740800 -1.746300
      Read-in Center 1741 is at   0.000000 -0.687900 -1.746300
      Read-in Center 1742 is at   0.000000 -0.635000 -1.746300
      Read-in Center 1743 is at   0.000000 -0.582100 -1.746300
      Read-in Center 1744 is at   0.000000 -0.529200 -1.746300
      Read-in Center 1745 is at   0.000000 -0.476300 -1.746300
      Read-in Center 1746 is at   0.000000 -0.423300 -1.746300
      Read-in Center 1747 is at   0.000000 -0.370400 -1.746300
      Read-in Center 1748 is at   0.000000 -0.317500 -1.746300
      Read-in Center 1749 is at   0.000000 -0.264600 -1.746300
      Read-in Center 1750 is at   0.000000 -0.211700 -1.746300
      Read-in Center 1751 is at   0.000000 -0.158800 -1.746300
      Read-in Center 1752 is at   0.000000 -0.105800 -1.746300
      Read-in Center 1753 is at   0.000000 -0.052900 -1.746300
      Read-in Center 1754 is at   0.000000  0.000000 -1.746300
      Read-in Center 1755 is at   0.000000  0.052900 -1.746300
      Read-in Center 1756 is at   0.000000  0.105800 -1.746300
      Read-in Center 1757 is at   0.000000  0.158800 -1.746300
      Read-in Center 1758 is at   0.000000  0.211700 -1.746300
      Read-in Center 1759 is at   0.000000  0.264600 -1.746300
      Read-in Center 1760 is at   0.000000  0.317500 -1.746300
      Read-in Center 1761 is at   0.000000  0.370400 -1.746300
      Read-in Center 1762 is at   0.000000  0.423300 -1.746300
      Read-in Center 1763 is at   0.000000  0.476300 -1.746300
      Read-in Center 1764 is at   0.000000  0.529200 -1.746300
      Read-in Center 1765 is at   0.000000  0.582100 -1.746300
      Read-in Center 1766 is at   0.000000  0.635000 -1.746300
      Read-in Center 1767 is at   0.000000  0.687900 -1.746300
      Read-in Center 1768 is at   0.000000  0.740800 -1.746300
      Read-in Center 1769 is at   0.000000  0.793800 -1.746300
      Read-in Center 1770 is at   0.000000  0.846700 -1.746300
      Read-in Center 1771 is at   0.000000  0.899600 -1.746300
      Read-in Center 1772 is at   0.000000  0.952500 -1.746300
      Read-in Center 1773 is at   0.000000  1.005400 -1.746300
      Read-in Center 1774 is at   0.000000  1.058400 -1.746300
      Read-in Center 1775 is at   0.000000  1.111300 -1.746300
      Read-in Center 1776 is at   0.000000  1.164200 -1.746300
      Read-in Center 1777 is at   0.000000  1.217100 -1.746300
      Read-in Center 1778 is at   0.000000  1.270000 -1.746300
      Read-in Center 1779 is at   0.000000  1.322900 -1.746300
      Read-in Center 1780 is at   0.000000  1.375900 -1.746300
      Read-in Center 1781 is at   0.000000  1.428800 -1.746300
      Read-in Center 1782 is at   0.000000  1.481700 -1.746300
      Read-in Center 1783 is at   0.000000  1.534600 -1.746300
      Read-in Center 1784 is at   0.000000  1.587500 -1.746300
      Read-in Center 1785 is at   0.000000  1.640400 -1.746300
      Read-in Center 1786 is at   0.000000  1.693400 -1.746300
      Read-in Center 1787 is at   0.000000  1.746300 -1.746300
      Read-in Center 1788 is at   0.000000  1.799200 -1.746300
      Read-in Center 1789 is at   0.000000  1.852100 -1.746300
      Read-in Center 1790 is at   0.000000  1.905000 -1.746300
      Read-in Center 1791 is at   0.000000  1.958000 -1.746300
      Read-in Center 1792 is at   0.000000  2.010900 -1.746300
      Read-in Center 1793 is at   0.000000  2.063800 -1.746300
      Read-in Center 1794 is at   0.000000  2.116700 -1.746300
      Read-in Center 1795 is at   0.000000  2.169600 -1.746300
      Read-in Center 1796 is at   0.000000  2.222500 -1.746300
      Read-in Center 1797 is at   0.000000  2.275500 -1.746300
      Read-in Center 1798 is at   0.000000  2.328400 -1.746300
      Read-in Center 1799 is at   0.000000  2.381300 -1.746300
      Read-in Center 1800 is at   0.000000  2.434200 -1.746300
      Read-in Center 1801 is at   0.000000  2.487100 -1.746300
      Read-in Center 1802 is at   0.000000  2.540000 -1.746300
      Read-in Center 1803 is at   0.000000  2.593000 -1.746300
      Read-in Center 1804 is at   0.000000 -2.645900 -1.693400
      Read-in Center 1805 is at   0.000000 -2.593000 -1.693400
      Read-in Center 1806 is at   0.000000 -2.540000 -1.693400
      Read-in Center 1807 is at   0.000000 -2.487100 -1.693400
      Read-in Center 1808 is at   0.000000 -2.434200 -1.693400
      Read-in Center 1809 is at   0.000000 -2.381300 -1.693400
      Read-in Center 1810 is at   0.000000 -2.328400 -1.693400
      Read-in Center 1811 is at   0.000000 -2.275500 -1.693400
      Read-in Center 1812 is at   0.000000 -2.222500 -1.693400
      Read-in Center 1813 is at   0.000000 -2.169600 -1.693400
      Read-in Center 1814 is at   0.000000 -2.116700 -1.693400
      Read-in Center 1815 is at   0.000000 -2.063800 -1.693400
      Read-in Center 1816 is at   0.000000 -2.010900 -1.693400
      Read-in Center 1817 is at   0.000000 -1.958000 -1.693400
      Read-in Center 1818 is at   0.000000 -1.905000 -1.693400
      Read-in Center 1819 is at   0.000000 -1.852100 -1.693400
      Read-in Center 1820 is at   0.000000 -1.799200 -1.693400
      Read-in Center 1821 is at   0.000000 -1.746300 -1.693400
      Read-in Center 1822 is at   0.000000 -1.693400 -1.693400
      Read-in Center 1823 is at   0.000000 -1.640400 -1.693400
      Read-in Center 1824 is at   0.000000 -1.587500 -1.693400
      Read-in Center 1825 is at   0.000000 -1.534600 -1.693400
      Read-in Center 1826 is at   0.000000 -1.481700 -1.693400
      Read-in Center 1827 is at   0.000000 -1.428800 -1.693400
      Read-in Center 1828 is at   0.000000 -1.375900 -1.693400
      Read-in Center 1829 is at   0.000000 -1.322900 -1.693400
      Read-in Center 1830 is at   0.000000 -1.270000 -1.693400
      Read-in Center 1831 is at   0.000000 -1.217100 -1.693400
      Read-in Center 1832 is at   0.000000 -1.164200 -1.693400
      Read-in Center 1833 is at   0.000000 -1.111300 -1.693400
      Read-in Center 1834 is at   0.000000 -1.058400 -1.693400
      Read-in Center 1835 is at   0.000000 -1.005400 -1.693400
      Read-in Center 1836 is at   0.000000 -0.952500 -1.693400
      Read-in Center 1837 is at   0.000000 -0.899600 -1.693400
      Read-in Center 1838 is at   0.000000 -0.846700 -1.693400
      Read-in Center 1839 is at   0.000000 -0.793800 -1.693400
      Read-in Center 1840 is at   0.000000 -0.740800 -1.693400
      Read-in Center 1841 is at   0.000000 -0.687900 -1.693400
      Read-in Center 1842 is at   0.000000 -0.635000 -1.693400
      Read-in Center 1843 is at   0.000000 -0.582100 -1.693400
      Read-in Center 1844 is at   0.000000 -0.529200 -1.693400
      Read-in Center 1845 is at   0.000000 -0.476300 -1.693400
      Read-in Center 1846 is at   0.000000 -0.423300 -1.693400
      Read-in Center 1847 is at   0.000000 -0.370400 -1.693400
      Read-in Center 1848 is at   0.000000 -0.317500 -1.693400
      Read-in Center 1849 is at   0.000000 -0.264600 -1.693400
      Read-in Center 1850 is at   0.000000 -0.211700 -1.693400
      Read-in Center 1851 is at   0.000000 -0.158800 -1.693400
      Read-in Center 1852 is at   0.000000 -0.105800 -1.693400
      Read-in Center 1853 is at   0.000000 -0.052900 -1.693400
      Read-in Center 1854 is at   0.000000  0.000000 -1.693400
      Read-in Center 1855 is at   0.000000  0.052900 -1.693400
      Read-in Center 1856 is at   0.000000  0.105800 -1.693400
      Read-in Center 1857 is at   0.000000  0.158800 -1.693400
      Read-in Center 1858 is at   0.000000  0.211700 -1.693400
      Read-in Center 1859 is at   0.000000  0.264600 -1.693400
      Read-in Center 1860 is at   0.000000  0.317500 -1.693400
      Read-in Center 1861 is at   0.000000  0.370400 -1.693400
      Read-in Center 1862 is at   0.000000  0.423300 -1.693400
      Read-in Center 1863 is at   0.000000  0.476300 -1.693400
      Read-in Center 1864 is at   0.000000  0.529200 -1.693400
      Read-in Center 1865 is at   0.000000  0.582100 -1.693400
      Read-in Center 1866 is at   0.000000  0.635000 -1.693400
      Read-in Center 1867 is at   0.000000  0.687900 -1.693400
      Read-in Center 1868 is at   0.000000  0.740800 -1.693400
      Read-in Center 1869 is at   0.000000  0.793800 -1.693400
      Read-in Center 1870 is at   0.000000  0.846700 -1.693400
      Read-in Center 1871 is at   0.000000  0.899600 -1.693400
      Read-in Center 1872 is at   0.000000  0.952500 -1.693400
      Read-in Center 1873 is at   0.000000  1.005400 -1.693400
      Read-in Center 1874 is at   0.000000  1.058400 -1.693400
      Read-in Center 1875 is at   0.000000  1.111300 -1.693400
      Read-in Center 1876 is at   0.000000  1.164200 -1.693400
      Read-in Center 1877 is at   0.000000  1.217100 -1.693400
      Read-in Center 1878 is at   0.000000  1.270000 -1.693400
      Read-in Center 1879 is at   0.000000  1.322900 -1.693400
      Read-in Center 1880 is at   0.000000  1.375900 -1.693400
      Read-in Center 1881 is at   0.000000  1.428800 -1.693400
      Read-in Center 1882 is at   0.000000  1.481700 -1.693400
      Read-in Center 1883 is at   0.000000  1.534600 -1.693400
      Read-in Center 1884 is at   0.000000  1.587500 -1.693400
      Read-in Center 1885 is at   0.000000  1.640400 -1.693400
      Read-in Center 1886 is at   0.000000  1.693400 -1.693400
      Read-in Center 1887 is at   0.000000  1.746300 -1.693400
      Read-in Center 1888 is at   0.000000  1.799200 -1.693400
      Read-in Center 1889 is at   0.000000  1.852100 -1.693400
      Read-in Center 1890 is at   0.000000  1.905000 -1.693400
      Read-in Center 1891 is at   0.000000  1.958000 -1.693400
      Read-in Center 1892 is at   0.000000  2.010900 -1.693400
      Read-in Center 1893 is at   0.000000  2.063800 -1.693400
      Read-in Center 1894 is at   0.000000  2.116700 -1.693400
      Read-in Center 1895 is at   0.000000  2.169600 -1.693400
      Read-in Center 1896 is at   0.000000  2.222500 -1.693400
      Read-in Center 1897 is at   0.000000  2.275500 -1.693400
      Read-in Center 1898 is at   0.000000  2.328400 -1.693400
      Read-in Center 1899 is at   0.000000  2.381300 -1.693400
      Read-in Center 1900 is at   0.000000  2.434200 -1.693400
      Read-in Center 1901 is at   0.000000  2.487100 -1.693400
      Read-in Center 1902 is at   0.000000  2.540000 -1.693400
      Read-in Center 1903 is at   0.000000  2.593000 -1.693400
      Read-in Center 1904 is at   0.000000 -2.645900 -1.640400
      Read-in Center 1905 is at   0.000000 -2.593000 -1.640400
      Read-in Center 1906 is at   0.000000 -2.540000 -1.640400
      Read-in Center 1907 is at   0.000000 -2.487100 -1.640400
      Read-in Center 1908 is at   0.000000 -2.434200 -1.640400
      Read-in Center 1909 is at   0.000000 -2.381300 -1.640400
      Read-in Center 1910 is at   0.000000 -2.328400 -1.640400
      Read-in Center 1911 is at   0.000000 -2.275500 -1.640400
      Read-in Center 1912 is at   0.000000 -2.222500 -1.640400
      Read-in Center 1913 is at   0.000000 -2.169600 -1.640400
      Read-in Center 1914 is at   0.000000 -2.116700 -1.640400
      Read-in Center 1915 is at   0.000000 -2.063800 -1.640400
      Read-in Center 1916 is at   0.000000 -2.010900 -1.640400
      Read-in Center 1917 is at   0.000000 -1.958000 -1.640400
      Read-in Center 1918 is at   0.000000 -1.905000 -1.640400
      Read-in Center 1919 is at   0.000000 -1.852100 -1.640400
      Read-in Center 1920 is at   0.000000 -1.799200 -1.640400
      Read-in Center 1921 is at   0.000000 -1.746300 -1.640400
      Read-in Center 1922 is at   0.000000 -1.693400 -1.640400
      Read-in Center 1923 is at   0.000000 -1.640400 -1.640400
      Read-in Center 1924 is at   0.000000 -1.587500 -1.640400
      Read-in Center 1925 is at   0.000000 -1.534600 -1.640400
      Read-in Center 1926 is at   0.000000 -1.481700 -1.640400
      Read-in Center 1927 is at   0.000000 -1.428800 -1.640400
      Read-in Center 1928 is at   0.000000 -1.375900 -1.640400
      Read-in Center 1929 is at   0.000000 -1.322900 -1.640400
      Read-in Center 1930 is at   0.000000 -1.270000 -1.640400
      Read-in Center 1931 is at   0.000000 -1.217100 -1.640400
      Read-in Center 1932 is at   0.000000 -1.164200 -1.640400
      Read-in Center 1933 is at   0.000000 -1.111300 -1.640400
      Read-in Center 1934 is at   0.000000 -1.058400 -1.640400
      Read-in Center 1935 is at   0.000000 -1.005400 -1.640400
      Read-in Center 1936 is at   0.000000 -0.952500 -1.640400
      Read-in Center 1937 is at   0.000000 -0.899600 -1.640400
      Read-in Center 1938 is at   0.000000 -0.846700 -1.640400
      Read-in Center 1939 is at   0.000000 -0.793800 -1.640400
      Read-in Center 1940 is at   0.000000 -0.740800 -1.640400
      Read-in Center 1941 is at   0.000000 -0.687900 -1.640400
      Read-in Center 1942 is at   0.000000 -0.635000 -1.640400
      Read-in Center 1943 is at   0.000000 -0.582100 -1.640400
      Read-in Center 1944 is at   0.000000 -0.529200 -1.640400
      Read-in Center 1945 is at   0.000000 -0.476300 -1.640400
      Read-in Center 1946 is at   0.000000 -0.423300 -1.640400
      Read-in Center 1947 is at   0.000000 -0.370400 -1.640400
      Read-in Center 1948 is at   0.000000 -0.317500 -1.640400
      Read-in Center 1949 is at   0.000000 -0.264600 -1.640400
      Read-in Center 1950 is at   0.000000 -0.211700 -1.640400
      Read-in Center 1951 is at   0.000000 -0.158800 -1.640400
      Read-in Center 1952 is at   0.000000 -0.105800 -1.640400
      Read-in Center 1953 is at   0.000000 -0.052900 -1.640400
      Read-in Center 1954 is at   0.000000  0.000000 -1.640400
      Read-in Center 1955 is at   0.000000  0.052900 -1.640400
      Read-in Center 1956 is at   0.000000  0.105800 -1.640400
      Read-in Center 1957 is at   0.000000  0.158800 -1.640400
      Read-in Center 1958 is at   0.000000  0.211700 -1.640400
      Read-in Center 1959 is at   0.000000  0.264600 -1.640400
      Read-in Center 1960 is at   0.000000  0.317500 -1.640400
      Read-in Center 1961 is at   0.000000  0.370400 -1.640400
      Read-in Center 1962 is at   0.000000  0.423300 -1.640400
      Read-in Center 1963 is at   0.000000  0.476300 -1.640400
      Read-in Center 1964 is at   0.000000  0.529200 -1.640400
      Read-in Center 1965 is at   0.000000  0.582100 -1.640400
      Read-in Center 1966 is at   0.000000  0.635000 -1.640400
      Read-in Center 1967 is at   0.000000  0.687900 -1.640400
      Read-in Center 1968 is at   0.000000  0.740800 -1.640400
      Read-in Center 1969 is at   0.000000  0.793800 -1.640400
      Read-in Center 1970 is at   0.000000  0.846700 -1.640400
      Read-in Center 1971 is at   0.000000  0.899600 -1.640400
      Read-in Center 1972 is at   0.000000  0.952500 -1.640400
      Read-in Center 1973 is at   0.000000  1.005400 -1.640400
      Read-in Center 1974 is at   0.000000  1.058400 -1.640400
      Read-in Center 1975 is at   0.000000  1.111300 -1.640400
      Read-in Center 1976 is at   0.000000  1.164200 -1.640400
      Read-in Center 1977 is at   0.000000  1.217100 -1.640400
      Read-in Center 1978 is at   0.000000  1.270000 -1.640400
      Read-in Center 1979 is at   0.000000  1.322900 -1.640400
      Read-in Center 1980 is at   0.000000  1.375900 -1.640400
      Read-in Center 1981 is at   0.000000  1.428800 -1.640400
      Read-in Center 1982 is at   0.000000  1.481700 -1.640400
      Read-in Center 1983 is at   0.000000  1.534600 -1.640400
      Read-in Center 1984 is at   0.000000  1.587500 -1.640400
      Read-in Center 1985 is at   0.000000  1.640400 -1.640400
      Read-in Center 1986 is at   0.000000  1.693400 -1.640400
      Read-in Center 1987 is at   0.000000  1.746300 -1.640400
      Read-in Center 1988 is at   0.000000  1.799200 -1.640400
      Read-in Center 1989 is at   0.000000  1.852100 -1.640400
      Read-in Center 1990 is at   0.000000  1.905000 -1.640400
      Read-in Center 1991 is at   0.000000  1.958000 -1.640400
      Read-in Center 1992 is at   0.000000  2.010900 -1.640400
      Read-in Center 1993 is at   0.000000  2.063800 -1.640400
      Read-in Center 1994 is at   0.000000  2.116700 -1.640400
      Read-in Center 1995 is at   0.000000  2.169600 -1.640400
      Read-in Center 1996 is at   0.000000  2.222500 -1.640400
      Read-in Center 1997 is at   0.000000  2.275500 -1.640400
      Read-in Center 1998 is at   0.000000  2.328400 -1.640400
      Read-in Center 1999 is at   0.000000  2.381300 -1.640400
      Read-in Center 2000 is at   0.000000  2.434200 -1.640400
      Read-in Center 2001 is at   0.000000  2.487100 -1.640400
      Read-in Center 2002 is at   0.000000  2.540000 -1.640400
      Read-in Center 2003 is at   0.000000  2.593000 -1.640400
      Read-in Center 2004 is at   0.000000 -2.645900 -1.587500
      Read-in Center 2005 is at   0.000000 -2.593000 -1.587500
      Read-in Center 2006 is at   0.000000 -2.540000 -1.587500
      Read-in Center 2007 is at   0.000000 -2.487100 -1.587500
      Read-in Center 2008 is at   0.000000 -2.434200 -1.587500
      Read-in Center 2009 is at   0.000000 -2.381300 -1.587500
      Read-in Center 2010 is at   0.000000 -2.328400 -1.587500
      Read-in Center 2011 is at   0.000000 -2.275500 -1.587500
      Read-in Center 2012 is at   0.000000 -2.222500 -1.587500
      Read-in Center 2013 is at   0.000000 -2.169600 -1.587500
      Read-in Center 2014 is at   0.000000 -2.116700 -1.587500
      Read-in Center 2015 is at   0.000000 -2.063800 -1.587500
      Read-in Center 2016 is at   0.000000 -2.010900 -1.587500
      Read-in Center 2017 is at   0.000000 -1.958000 -1.587500
      Read-in Center 2018 is at   0.000000 -1.905000 -1.587500
      Read-in Center 2019 is at   0.000000 -1.852100 -1.587500
      Read-in Center 2020 is at   0.000000 -1.799200 -1.587500
      Read-in Center 2021 is at   0.000000 -1.746300 -1.587500
      Read-in Center 2022 is at   0.000000 -1.693400 -1.587500
      Read-in Center 2023 is at   0.000000 -1.640400 -1.587500
      Read-in Center 2024 is at   0.000000 -1.587500 -1.587500
      Read-in Center 2025 is at   0.000000 -1.534600 -1.587500
      Read-in Center 2026 is at   0.000000 -1.481700 -1.587500
      Read-in Center 2027 is at   0.000000 -1.428800 -1.587500
      Read-in Center 2028 is at   0.000000 -1.375900 -1.587500
      Read-in Center 2029 is at   0.000000 -1.322900 -1.587500
      Read-in Center 2030 is at   0.000000 -1.270000 -1.587500
      Read-in Center 2031 is at   0.000000 -1.217100 -1.587500
      Read-in Center 2032 is at   0.000000 -1.164200 -1.587500
      Read-in Center 2033 is at   0.000000 -1.111300 -1.587500
      Read-in Center 2034 is at   0.000000 -1.058400 -1.587500
      Read-in Center 2035 is at   0.000000 -1.005400 -1.587500
      Read-in Center 2036 is at   0.000000 -0.952500 -1.587500
      Read-in Center 2037 is at   0.000000 -0.899600 -1.587500
      Read-in Center 2038 is at   0.000000 -0.846700 -1.587500
      Read-in Center 2039 is at   0.000000 -0.793800 -1.587500
      Read-in Center 2040 is at   0.000000 -0.740800 -1.587500
      Read-in Center 2041 is at   0.000000 -0.687900 -1.587500
      Read-in Center 2042 is at   0.000000 -0.635000 -1.587500
      Read-in Center 2043 is at   0.000000 -0.582100 -1.587500
      Read-in Center 2044 is at   0.000000 -0.529200 -1.587500
      Read-in Center 2045 is at   0.000000 -0.476300 -1.587500
      Read-in Center 2046 is at   0.000000 -0.423300 -1.587500
      Read-in Center 2047 is at   0.000000 -0.370400 -1.587500
      Read-in Center 2048 is at   0.000000 -0.317500 -1.587500
      Read-in Center 2049 is at   0.000000 -0.264600 -1.587500
      Read-in Center 2050 is at   0.000000 -0.211700 -1.587500
      Read-in Center 2051 is at   0.000000 -0.158800 -1.587500
      Read-in Center 2052 is at   0.000000 -0.105800 -1.587500
      Read-in Center 2053 is at   0.000000 -0.052900 -1.587500
      Read-in Center 2054 is at   0.000000  0.000000 -1.587500
      Read-in Center 2055 is at   0.000000  0.052900 -1.587500
      Read-in Center 2056 is at   0.000000  0.105800 -1.587500
      Read-in Center 2057 is at   0.000000  0.158800 -1.587500
      Read-in Center 2058 is at   0.000000  0.211700 -1.587500
      Read-in Center 2059 is at   0.000000  0.264600 -1.587500
      Read-in Center 2060 is at   0.000000  0.317500 -1.587500
      Read-in Center 2061 is at   0.000000  0.370400 -1.587500
      Read-in Center 2062 is at   0.000000  0.423300 -1.587500
      Read-in Center 2063 is at   0.000000  0.476300 -1.587500
      Read-in Center 2064 is at   0.000000  0.529200 -1.587500
      Read-in Center 2065 is at   0.000000  0.582100 -1.587500
      Read-in Center 2066 is at   0.000000  0.635000 -1.587500
      Read-in Center 2067 is at   0.000000  0.687900 -1.587500
      Read-in Center 2068 is at   0.000000  0.740800 -1.587500
      Read-in Center 2069 is at   0.000000  0.793800 -1.587500
      Read-in Center 2070 is at   0.000000  0.846700 -1.587500
      Read-in Center 2071 is at   0.000000  0.899600 -1.587500
      Read-in Center 2072 is at   0.000000  0.952500 -1.587500
      Read-in Center 2073 is at   0.000000  1.005400 -1.587500
      Read-in Center 2074 is at   0.000000  1.058400 -1.587500
      Read-in Center 2075 is at   0.000000  1.111300 -1.587500
      Read-in Center 2076 is at   0.000000  1.164200 -1.587500
      Read-in Center 2077 is at   0.000000  1.217100 -1.587500
      Read-in Center 2078 is at   0.000000  1.270000 -1.587500
      Read-in Center 2079 is at   0.000000  1.322900 -1.587500
      Read-in Center 2080 is at   0.000000  1.375900 -1.587500
      Read-in Center 2081 is at   0.000000  1.428800 -1.587500
      Read-in Center 2082 is at   0.000000  1.481700 -1.587500
      Read-in Center 2083 is at   0.000000  1.534600 -1.587500
      Read-in Center 2084 is at   0.000000  1.587500 -1.587500
      Read-in Center 2085 is at   0.000000  1.640400 -1.587500
      Read-in Center 2086 is at   0.000000  1.693400 -1.587500
      Read-in Center 2087 is at   0.000000  1.746300 -1.587500
      Read-in Center 2088 is at   0.000000  1.799200 -1.587500
      Read-in Center 2089 is at   0.000000  1.852100 -1.587500
      Read-in Center 2090 is at   0.000000  1.905000 -1.587500
      Read-in Center 2091 is at   0.000000  1.958000 -1.587500
      Read-in Center 2092 is at   0.000000  2.010900 -1.587500
      Read-in Center 2093 is at   0.000000  2.063800 -1.587500
      Read-in Center 2094 is at   0.000000  2.116700 -1.587500
      Read-in Center 2095 is at   0.000000  2.169600 -1.587500
      Read-in Center 2096 is at   0.000000  2.222500 -1.587500
      Read-in Center 2097 is at   0.000000  2.275500 -1.587500
      Read-in Center 2098 is at   0.000000  2.328400 -1.587500
      Read-in Center 2099 is at   0.000000  2.381300 -1.587500
      Read-in Center 2100 is at   0.000000  2.434200 -1.587500
      Read-in Center 2101 is at   0.000000  2.487100 -1.587500
      Read-in Center 2102 is at   0.000000  2.540000 -1.587500
      Read-in Center 2103 is at   0.000000  2.593000 -1.587500
      Read-in Center 2104 is at   0.000000 -2.645900 -1.534600
      Read-in Center 2105 is at   0.000000 -2.593000 -1.534600
      Read-in Center 2106 is at   0.000000 -2.540000 -1.534600
      Read-in Center 2107 is at   0.000000 -2.487100 -1.534600
      Read-in Center 2108 is at   0.000000 -2.434200 -1.534600
      Read-in Center 2109 is at   0.000000 -2.381300 -1.534600
      Read-in Center 2110 is at   0.000000 -2.328400 -1.534600
      Read-in Center 2111 is at   0.000000 -2.275500 -1.534600
      Read-in Center 2112 is at   0.000000 -2.222500 -1.534600
      Read-in Center 2113 is at   0.000000 -2.169600 -1.534600
      Read-in Center 2114 is at   0.000000 -2.116700 -1.534600
      Read-in Center 2115 is at   0.000000 -2.063800 -1.534600
      Read-in Center 2116 is at   0.000000 -2.010900 -1.534600
      Read-in Center 2117 is at   0.000000 -1.958000 -1.534600
      Read-in Center 2118 is at   0.000000 -1.905000 -1.534600
      Read-in Center 2119 is at   0.000000 -1.852100 -1.534600
      Read-in Center 2120 is at   0.000000 -1.799200 -1.534600
      Read-in Center 2121 is at   0.000000 -1.746300 -1.534600
      Read-in Center 2122 is at   0.000000 -1.693400 -1.534600
      Read-in Center 2123 is at   0.000000 -1.640400 -1.534600
      Read-in Center 2124 is at   0.000000 -1.587500 -1.534600
      Read-in Center 2125 is at   0.000000 -1.534600 -1.534600
      Read-in Center 2126 is at   0.000000 -1.481700 -1.534600
      Read-in Center 2127 is at   0.000000 -1.428800 -1.534600
      Read-in Center 2128 is at   0.000000 -1.375900 -1.534600
      Read-in Center 2129 is at   0.000000 -1.322900 -1.534600
      Read-in Center 2130 is at   0.000000 -1.270000 -1.534600
      Read-in Center 2131 is at   0.000000 -1.217100 -1.534600
      Read-in Center 2132 is at   0.000000 -1.164200 -1.534600
      Read-in Center 2133 is at   0.000000 -1.111300 -1.534600
      Read-in Center 2134 is at   0.000000 -1.058400 -1.534600
      Read-in Center 2135 is at   0.000000 -1.005400 -1.534600
      Read-in Center 2136 is at   0.000000 -0.952500 -1.534600
      Read-in Center 2137 is at   0.000000 -0.899600 -1.534600
      Read-in Center 2138 is at   0.000000 -0.846700 -1.534600
      Read-in Center 2139 is at   0.000000 -0.793800 -1.534600
      Read-in Center 2140 is at   0.000000 -0.740800 -1.534600
      Read-in Center 2141 is at   0.000000 -0.687900 -1.534600
      Read-in Center 2142 is at   0.000000 -0.635000 -1.534600
      Read-in Center 2143 is at   0.000000 -0.582100 -1.534600
      Read-in Center 2144 is at   0.000000 -0.529200 -1.534600
      Read-in Center 2145 is at   0.000000 -0.476300 -1.534600
      Read-in Center 2146 is at   0.000000 -0.423300 -1.534600
      Read-in Center 2147 is at   0.000000 -0.370400 -1.534600
      Read-in Center 2148 is at   0.000000 -0.317500 -1.534600
      Read-in Center 2149 is at   0.000000 -0.264600 -1.534600
      Read-in Center 2150 is at   0.000000 -0.211700 -1.534600
      Read-in Center 2151 is at   0.000000 -0.158800 -1.534600
      Read-in Center 2152 is at   0.000000 -0.105800 -1.534600
      Read-in Center 2153 is at   0.000000 -0.052900 -1.534600
      Read-in Center 2154 is at   0.000000  0.000000 -1.534600
      Read-in Center 2155 is at   0.000000  0.052900 -1.534600
      Read-in Center 2156 is at   0.000000  0.105800 -1.534600
      Read-in Center 2157 is at   0.000000  0.158800 -1.534600
      Read-in Center 2158 is at   0.000000  0.211700 -1.534600
      Read-in Center 2159 is at   0.000000  0.264600 -1.534600
      Read-in Center 2160 is at   0.000000  0.317500 -1.534600
      Read-in Center 2161 is at   0.000000  0.370400 -1.534600
      Read-in Center 2162 is at   0.000000  0.423300 -1.534600
      Read-in Center 2163 is at   0.000000  0.476300 -1.534600
      Read-in Center 2164 is at   0.000000  0.529200 -1.534600
      Read-in Center 2165 is at   0.000000  0.582100 -1.534600
      Read-in Center 2166 is at   0.000000  0.635000 -1.534600
      Read-in Center 2167 is at   0.000000  0.687900 -1.534600
      Read-in Center 2168 is at   0.000000  0.740800 -1.534600
      Read-in Center 2169 is at   0.000000  0.793800 -1.534600
      Read-in Center 2170 is at   0.000000  0.846700 -1.534600
      Read-in Center 2171 is at   0.000000  0.899600 -1.534600
      Read-in Center 2172 is at   0.000000  0.952500 -1.534600
      Read-in Center 2173 is at   0.000000  1.005400 -1.534600
      Read-in Center 2174 is at   0.000000  1.058400 -1.534600
      Read-in Center 2175 is at   0.000000  1.111300 -1.534600
      Read-in Center 2176 is at   0.000000  1.164200 -1.534600
      Read-in Center 2177 is at   0.000000  1.217100 -1.534600
      Read-in Center 2178 is at   0.000000  1.270000 -1.534600
      Read-in Center 2179 is at   0.000000  1.322900 -1.534600
      Read-in Center 2180 is at   0.000000  1.375900 -1.534600
      Read-in Center 2181 is at   0.000000  1.428800 -1.534600
      Read-in Center 2182 is at   0.000000  1.481700 -1.534600
      Read-in Center 2183 is at   0.000000  1.534600 -1.534600
      Read-in Center 2184 is at   0.000000  1.587500 -1.534600
      Read-in Center 2185 is at   0.000000  1.640400 -1.534600
      Read-in Center 2186 is at   0.000000  1.693400 -1.534600
      Read-in Center 2187 is at   0.000000  1.746300 -1.534600
      Read-in Center 2188 is at   0.000000  1.799200 -1.534600
      Read-in Center 2189 is at   0.000000  1.852100 -1.534600
      Read-in Center 2190 is at   0.000000  1.905000 -1.534600
      Read-in Center 2191 is at   0.000000  1.958000 -1.534600
      Read-in Center 2192 is at   0.000000  2.010900 -1.534600
      Read-in Center 2193 is at   0.000000  2.063800 -1.534600
      Read-in Center 2194 is at   0.000000  2.116700 -1.534600
      Read-in Center 2195 is at   0.000000  2.169600 -1.534600
      Read-in Center 2196 is at   0.000000  2.222500 -1.534600
      Read-in Center 2197 is at   0.000000  2.275500 -1.534600
      Read-in Center 2198 is at   0.000000  2.328400 -1.534600
      Read-in Center 2199 is at   0.000000  2.381300 -1.534600
      Read-in Center 2200 is at   0.000000  2.434200 -1.534600
      Read-in Center 2201 is at   0.000000  2.487100 -1.534600
      Read-in Center 2202 is at   0.000000  2.540000 -1.534600
      Read-in Center 2203 is at   0.000000  2.593000 -1.534600
      Read-in Center 2204 is at   0.000000 -2.645900 -1.481700
      Read-in Center 2205 is at   0.000000 -2.593000 -1.481700
      Read-in Center 2206 is at   0.000000 -2.540000 -1.481700
      Read-in Center 2207 is at   0.000000 -2.487100 -1.481700
      Read-in Center 2208 is at   0.000000 -2.434200 -1.481700
      Read-in Center 2209 is at   0.000000 -2.381300 -1.481700
      Read-in Center 2210 is at   0.000000 -2.328400 -1.481700
      Read-in Center 2211 is at   0.000000 -2.275500 -1.481700
      Read-in Center 2212 is at   0.000000 -2.222500 -1.481700
      Read-in Center 2213 is at   0.000000 -2.169600 -1.481700
      Read-in Center 2214 is at   0.000000 -2.116700 -1.481700
      Read-in Center 2215 is at   0.000000 -2.063800 -1.481700
      Read-in Center 2216 is at   0.000000 -2.010900 -1.481700
      Read-in Center 2217 is at   0.000000 -1.958000 -1.481700
      Read-in Center 2218 is at   0.000000 -1.905000 -1.481700
      Read-in Center 2219 is at   0.000000 -1.852100 -1.481700
      Read-in Center 2220 is at   0.000000 -1.799200 -1.481700
      Read-in Center 2221 is at   0.000000 -1.746300 -1.481700
      Read-in Center 2222 is at   0.000000 -1.693400 -1.481700
      Read-in Center 2223 is at   0.000000 -1.640400 -1.481700
      Read-in Center 2224 is at   0.000000 -1.587500 -1.481700
      Read-in Center 2225 is at   0.000000 -1.534600 -1.481700
      Read-in Center 2226 is at   0.000000 -1.481700 -1.481700
      Read-in Center 2227 is at   0.000000 -1.428800 -1.481700
      Read-in Center 2228 is at   0.000000 -1.375900 -1.481700
      Read-in Center 2229 is at   0.000000 -1.322900 -1.481700
      Read-in Center 2230 is at   0.000000 -1.270000 -1.481700
      Read-in Center 2231 is at   0.000000 -1.217100 -1.481700
      Read-in Center 2232 is at   0.000000 -1.164200 -1.481700
      Read-in Center 2233 is at   0.000000 -1.111300 -1.481700
      Read-in Center 2234 is at   0.000000 -1.058400 -1.481700
      Read-in Center 2235 is at   0.000000 -1.005400 -1.481700
      Read-in Center 2236 is at   0.000000 -0.952500 -1.481700
      Read-in Center 2237 is at   0.000000 -0.899600 -1.481700
      Read-in Center 2238 is at   0.000000 -0.846700 -1.481700
      Read-in Center 2239 is at   0.000000 -0.793800 -1.481700
      Read-in Center 2240 is at   0.000000 -0.740800 -1.481700
      Read-in Center 2241 is at   0.000000 -0.687900 -1.481700
      Read-in Center 2242 is at   0.000000 -0.635000 -1.481700
      Read-in Center 2243 is at   0.000000 -0.582100 -1.481700
      Read-in Center 2244 is at   0.000000 -0.529200 -1.481700
      Read-in Center 2245 is at   0.000000 -0.476300 -1.481700
      Read-in Center 2246 is at   0.000000 -0.423300 -1.481700
      Read-in Center 2247 is at   0.000000 -0.370400 -1.481700
      Read-in Center 2248 is at   0.000000 -0.317500 -1.481700
      Read-in Center 2249 is at   0.000000 -0.264600 -1.481700
      Read-in Center 2250 is at   0.000000 -0.211700 -1.481700
      Read-in Center 2251 is at   0.000000 -0.158800 -1.481700
      Read-in Center 2252 is at   0.000000 -0.105800 -1.481700
      Read-in Center 2253 is at   0.000000 -0.052900 -1.481700
      Read-in Center 2254 is at   0.000000  0.000000 -1.481700
      Read-in Center 2255 is at   0.000000  0.052900 -1.481700
      Read-in Center 2256 is at   0.000000  0.105800 -1.481700
      Read-in Center 2257 is at   0.000000  0.158800 -1.481700
      Read-in Center 2258 is at   0.000000  0.211700 -1.481700
      Read-in Center 2259 is at   0.000000  0.264600 -1.481700
      Read-in Center 2260 is at   0.000000  0.317500 -1.481700
      Read-in Center 2261 is at   0.000000  0.370400 -1.481700
      Read-in Center 2262 is at   0.000000  0.423300 -1.481700
      Read-in Center 2263 is at   0.000000  0.476300 -1.481700
      Read-in Center 2264 is at   0.000000  0.529200 -1.481700
      Read-in Center 2265 is at   0.000000  0.582100 -1.481700
      Read-in Center 2266 is at   0.000000  0.635000 -1.481700
      Read-in Center 2267 is at   0.000000  0.687900 -1.481700
      Read-in Center 2268 is at   0.000000  0.740800 -1.481700
      Read-in Center 2269 is at   0.000000  0.793800 -1.481700
      Read-in Center 2270 is at   0.000000  0.846700 -1.481700
      Read-in Center 2271 is at   0.000000  0.899600 -1.481700
      Read-in Center 2272 is at   0.000000  0.952500 -1.481700
      Read-in Center 2273 is at   0.000000  1.005400 -1.481700
      Read-in Center 2274 is at   0.000000  1.058400 -1.481700
      Read-in Center 2275 is at   0.000000  1.111300 -1.481700
      Read-in Center 2276 is at   0.000000  1.164200 -1.481700
      Read-in Center 2277 is at   0.000000  1.217100 -1.481700
      Read-in Center 2278 is at   0.000000  1.270000 -1.481700
      Read-in Center 2279 is at   0.000000  1.322900 -1.481700
      Read-in Center 2280 is at   0.000000  1.375900 -1.481700
      Read-in Center 2281 is at   0.000000  1.428800 -1.481700
      Read-in Center 2282 is at   0.000000  1.481700 -1.481700
      Read-in Center 2283 is at   0.000000  1.534600 -1.481700
      Read-in Center 2284 is at   0.000000  1.587500 -1.481700
      Read-in Center 2285 is at   0.000000  1.640400 -1.481700
      Read-in Center 2286 is at   0.000000  1.693400 -1.481700
      Read-in Center 2287 is at   0.000000  1.746300 -1.481700
      Read-in Center 2288 is at   0.000000  1.799200 -1.481700
      Read-in Center 2289 is at   0.000000  1.852100 -1.481700
      Read-in Center 2290 is at   0.000000  1.905000 -1.481700
      Read-in Center 2291 is at   0.000000  1.958000 -1.481700
      Read-in Center 2292 is at   0.000000  2.010900 -1.481700
      Read-in Center 2293 is at   0.000000  2.063800 -1.481700
      Read-in Center 2294 is at   0.000000  2.116700 -1.481700
      Read-in Center 2295 is at   0.000000  2.169600 -1.481700
      Read-in Center 2296 is at   0.000000  2.222500 -1.481700
      Read-in Center 2297 is at   0.000000  2.275500 -1.481700
      Read-in Center 2298 is at   0.000000  2.328400 -1.481700
      Read-in Center 2299 is at   0.000000  2.381300 -1.481700
      Read-in Center 2300 is at   0.000000  2.434200 -1.481700
      Read-in Center 2301 is at   0.000000  2.487100 -1.481700
      Read-in Center 2302 is at   0.000000  2.540000 -1.481700
      Read-in Center 2303 is at   0.000000  2.593000 -1.481700
      Read-in Center 2304 is at   0.000000 -2.645900 -1.428800
      Read-in Center 2305 is at   0.000000 -2.593000 -1.428800
      Read-in Center 2306 is at   0.000000 -2.540000 -1.428800
      Read-in Center 2307 is at   0.000000 -2.487100 -1.428800
      Read-in Center 2308 is at   0.000000 -2.434200 -1.428800
      Read-in Center 2309 is at   0.000000 -2.381300 -1.428800
      Read-in Center 2310 is at   0.000000 -2.328400 -1.428800
      Read-in Center 2311 is at   0.000000 -2.275500 -1.428800
      Read-in Center 2312 is at   0.000000 -2.222500 -1.428800
      Read-in Center 2313 is at   0.000000 -2.169600 -1.428800
      Read-in Center 2314 is at   0.000000 -2.116700 -1.428800
      Read-in Center 2315 is at   0.000000 -2.063800 -1.428800
      Read-in Center 2316 is at   0.000000 -2.010900 -1.428800
      Read-in Center 2317 is at   0.000000 -1.958000 -1.428800
      Read-in Center 2318 is at   0.000000 -1.905000 -1.428800
      Read-in Center 2319 is at   0.000000 -1.852100 -1.428800
      Read-in Center 2320 is at   0.000000 -1.799200 -1.428800
      Read-in Center 2321 is at   0.000000 -1.746300 -1.428800
      Read-in Center 2322 is at   0.000000 -1.693400 -1.428800
      Read-in Center 2323 is at   0.000000 -1.640400 -1.428800
      Read-in Center 2324 is at   0.000000 -1.587500 -1.428800
      Read-in Center 2325 is at   0.000000 -1.534600 -1.428800
      Read-in Center 2326 is at   0.000000 -1.481700 -1.428800
      Read-in Center 2327 is at   0.000000 -1.428800 -1.428800
      Read-in Center 2328 is at   0.000000 -1.375900 -1.428800
      Read-in Center 2329 is at   0.000000 -1.322900 -1.428800
      Read-in Center 2330 is at   0.000000 -1.270000 -1.428800
      Read-in Center 2331 is at   0.000000 -1.217100 -1.428800
      Read-in Center 2332 is at   0.000000 -1.164200 -1.428800
      Read-in Center 2333 is at   0.000000 -1.111300 -1.428800
      Read-in Center 2334 is at   0.000000 -1.058400 -1.428800
      Read-in Center 2335 is at   0.000000 -1.005400 -1.428800
      Read-in Center 2336 is at   0.000000 -0.952500 -1.428800
      Read-in Center 2337 is at   0.000000 -0.899600 -1.428800
      Read-in Center 2338 is at   0.000000 -0.846700 -1.428800
      Read-in Center 2339 is at   0.000000 -0.793800 -1.428800
      Read-in Center 2340 is at   0.000000 -0.740800 -1.428800
      Read-in Center 2341 is at   0.000000 -0.687900 -1.428800
      Read-in Center 2342 is at   0.000000 -0.635000 -1.428800
      Read-in Center 2343 is at   0.000000 -0.582100 -1.428800
      Read-in Center 2344 is at   0.000000 -0.529200 -1.428800
      Read-in Center 2345 is at   0.000000 -0.476300 -1.428800
      Read-in Center 2346 is at   0.000000 -0.423300 -1.428800
      Read-in Center 2347 is at   0.000000 -0.370400 -1.428800
      Read-in Center 2348 is at   0.000000 -0.317500 -1.428800
      Read-in Center 2349 is at   0.000000 -0.264600 -1.428800
      Read-in Center 2350 is at   0.000000 -0.211700 -1.428800
      Read-in Center 2351 is at   0.000000 -0.158800 -1.428800
      Read-in Center 2352 is at   0.000000 -0.105800 -1.428800
      Read-in Center 2353 is at   0.000000 -0.052900 -1.428800
      Read-in Center 2354 is at   0.000000  0.000000 -1.428800
      Read-in Center 2355 is at   0.000000  0.052900 -1.428800
      Read-in Center 2356 is at   0.000000  0.105800 -1.428800
      Read-in Center 2357 is at   0.000000  0.158800 -1.428800
      Read-in Center 2358 is at   0.000000  0.211700 -1.428800
      Read-in Center 2359 is at   0.000000  0.264600 -1.428800
      Read-in Center 2360 is at   0.000000  0.317500 -1.428800
      Read-in Center 2361 is at   0.000000  0.370400 -1.428800
      Read-in Center 2362 is at   0.000000  0.423300 -1.428800
      Read-in Center 2363 is at   0.000000  0.476300 -1.428800
      Read-in Center 2364 is at   0.000000  0.529200 -1.428800
      Read-in Center 2365 is at   0.000000  0.582100 -1.428800
      Read-in Center 2366 is at   0.000000  0.635000 -1.428800
      Read-in Center 2367 is at   0.000000  0.687900 -1.428800
      Read-in Center 2368 is at   0.000000  0.740800 -1.428800
      Read-in Center 2369 is at   0.000000  0.793800 -1.428800
      Read-in Center 2370 is at   0.000000  0.846700 -1.428800
      Read-in Center 2371 is at   0.000000  0.899600 -1.428800
      Read-in Center 2372 is at   0.000000  0.952500 -1.428800
      Read-in Center 2373 is at   0.000000  1.005400 -1.428800
      Read-in Center 2374 is at   0.000000  1.058400 -1.428800
      Read-in Center 2375 is at   0.000000  1.111300 -1.428800
      Read-in Center 2376 is at   0.000000  1.164200 -1.428800
      Read-in Center 2377 is at   0.000000  1.217100 -1.428800
      Read-in Center 2378 is at   0.000000  1.270000 -1.428800
      Read-in Center 2379 is at   0.000000  1.322900 -1.428800
      Read-in Center 2380 is at   0.000000  1.375900 -1.428800
      Read-in Center 2381 is at   0.000000  1.428800 -1.428800
      Read-in Center 2382 is at   0.000000  1.481700 -1.428800
      Read-in Center 2383 is at   0.000000  1.534600 -1.428800
      Read-in Center 2384 is at   0.000000  1.587500 -1.428800
      Read-in Center 2385 is at   0.000000  1.640400 -1.428800
      Read-in Center 2386 is at   0.000000  1.693400 -1.428800
      Read-in Center 2387 is at   0.000000  1.746300 -1.428800
      Read-in Center 2388 is at   0.000000  1.799200 -1.428800
      Read-in Center 2389 is at   0.000000  1.852100 -1.428800
      Read-in Center 2390 is at   0.000000  1.905000 -1.428800
      Read-in Center 2391 is at   0.000000  1.958000 -1.428800
      Read-in Center 2392 is at   0.000000  2.010900 -1.428800
      Read-in Center 2393 is at   0.000000  2.063800 -1.428800
      Read-in Center 2394 is at   0.000000  2.116700 -1.428800
      Read-in Center 2395 is at   0.000000  2.169600 -1.428800
      Read-in Center 2396 is at   0.000000  2.222500 -1.428800
      Read-in Center 2397 is at   0.000000  2.275500 -1.428800
      Read-in Center 2398 is at   0.000000  2.328400 -1.428800
      Read-in Center 2399 is at   0.000000  2.381300 -1.428800
      Read-in Center 2400 is at   0.000000  2.434200 -1.428800
      Read-in Center 2401 is at   0.000000  2.487100 -1.428800
      Read-in Center 2402 is at   0.000000  2.540000 -1.428800
      Read-in Center 2403 is at   0.000000  2.593000 -1.428800
      Read-in Center 2404 is at   0.000000 -2.645900 -1.375900
      Read-in Center 2405 is at   0.000000 -2.593000 -1.375900
      Read-in Center 2406 is at   0.000000 -2.540000 -1.375900
      Read-in Center 2407 is at   0.000000 -2.487100 -1.375900
      Read-in Center 2408 is at   0.000000 -2.434200 -1.375900
      Read-in Center 2409 is at   0.000000 -2.381300 -1.375900
      Read-in Center 2410 is at   0.000000 -2.328400 -1.375900
      Read-in Center 2411 is at   0.000000 -2.275500 -1.375900
      Read-in Center 2412 is at   0.000000 -2.222500 -1.375900
      Read-in Center 2413 is at   0.000000 -2.169600 -1.375900
      Read-in Center 2414 is at   0.000000 -2.116700 -1.375900
      Read-in Center 2415 is at   0.000000 -2.063800 -1.375900
      Read-in Center 2416 is at   0.000000 -2.010900 -1.375900
      Read-in Center 2417 is at   0.000000 -1.958000 -1.375900
      Read-in Center 2418 is at   0.000000 -1.905000 -1.375900
      Read-in Center 2419 is at   0.000000 -1.852100 -1.375900
      Read-in Center 2420 is at   0.000000 -1.799200 -1.375900
      Read-in Center 2421 is at   0.000000 -1.746300 -1.375900
      Read-in Center 2422 is at   0.000000 -1.693400 -1.375900
      Read-in Center 2423 is at   0.000000 -1.640400 -1.375900
      Read-in Center 2424 is at   0.000000 -1.587500 -1.375900
      Read-in Center 2425 is at   0.000000 -1.534600 -1.375900
      Read-in Center 2426 is at   0.000000 -1.481700 -1.375900
      Read-in Center 2427 is at   0.000000 -1.428800 -1.375900
      Read-in Center 2428 is at   0.000000 -1.375900 -1.375900
      Read-in Center 2429 is at   0.000000 -1.322900 -1.375900
      Read-in Center 2430 is at   0.000000 -1.270000 -1.375900
      Read-in Center 2431 is at   0.000000 -1.217100 -1.375900
      Read-in Center 2432 is at   0.000000 -1.164200 -1.375900
      Read-in Center 2433 is at   0.000000 -1.111300 -1.375900
      Read-in Center 2434 is at   0.000000 -1.058400 -1.375900
      Read-in Center 2435 is at   0.000000 -1.005400 -1.375900
      Read-in Center 2436 is at   0.000000 -0.952500 -1.375900
      Read-in Center 2437 is at   0.000000 -0.899600 -1.375900
      Read-in Center 2438 is at   0.000000 -0.846700 -1.375900
      Read-in Center 2439 is at   0.000000 -0.793800 -1.375900
      Read-in Center 2440 is at   0.000000 -0.740800 -1.375900
      Read-in Center 2441 is at   0.000000 -0.687900 -1.375900
      Read-in Center 2442 is at   0.000000 -0.635000 -1.375900
      Read-in Center 2443 is at   0.000000 -0.582100 -1.375900
      Read-in Center 2444 is at   0.000000 -0.529200 -1.375900
      Read-in Center 2445 is at   0.000000 -0.476300 -1.375900
      Read-in Center 2446 is at   0.000000 -0.423300 -1.375900
      Read-in Center 2447 is at   0.000000 -0.370400 -1.375900
      Read-in Center 2448 is at   0.000000 -0.317500 -1.375900
      Read-in Center 2449 is at   0.000000 -0.264600 -1.375900
      Read-in Center 2450 is at   0.000000 -0.211700 -1.375900
      Read-in Center 2451 is at   0.000000 -0.158800 -1.375900
      Read-in Center 2452 is at   0.000000 -0.105800 -1.375900
      Read-in Center 2453 is at   0.000000 -0.052900 -1.375900
      Read-in Center 2454 is at   0.000000  0.000000 -1.375900
      Read-in Center 2455 is at   0.000000  0.052900 -1.375900
      Read-in Center 2456 is at   0.000000  0.105800 -1.375900
      Read-in Center 2457 is at   0.000000  0.158800 -1.375900
      Read-in Center 2458 is at   0.000000  0.211700 -1.375900
      Read-in Center 2459 is at   0.000000  0.264600 -1.375900
      Read-in Center 2460 is at   0.000000  0.317500 -1.375900
      Read-in Center 2461 is at   0.000000  0.370400 -1.375900
      Read-in Center 2462 is at   0.000000  0.423300 -1.375900
      Read-in Center 2463 is at   0.000000  0.476300 -1.375900
      Read-in Center 2464 is at   0.000000  0.529200 -1.375900
      Read-in Center 2465 is at   0.000000  0.582100 -1.375900
      Read-in Center 2466 is at   0.000000  0.635000 -1.375900
      Read-in Center 2467 is at   0.000000  0.687900 -1.375900
      Read-in Center 2468 is at   0.000000  0.740800 -1.375900
      Read-in Center 2469 is at   0.000000  0.793800 -1.375900
      Read-in Center 2470 is at   0.000000  0.846700 -1.375900
      Read-in Center 2471 is at   0.000000  0.899600 -1.375900
      Read-in Center 2472 is at   0.000000  0.952500 -1.375900
      Read-in Center 2473 is at   0.000000  1.005400 -1.375900
      Read-in Center 2474 is at   0.000000  1.058400 -1.375900
      Read-in Center 2475 is at   0.000000  1.111300 -1.375900
      Read-in Center 2476 is at   0.000000  1.164200 -1.375900
      Read-in Center 2477 is at   0.000000  1.217100 -1.375900
      Read-in Center 2478 is at   0.000000  1.270000 -1.375900
      Read-in Center 2479 is at   0.000000  1.322900 -1.375900
      Read-in Center 2480 is at   0.000000  1.375900 -1.375900
      Read-in Center 2481 is at   0.000000  1.428800 -1.375900
      Read-in Center 2482 is at   0.000000  1.481700 -1.375900
      Read-in Center 2483 is at   0.000000  1.534600 -1.375900
      Read-in Center 2484 is at   0.000000  1.587500 -1.375900
      Read-in Center 2485 is at   0.000000  1.640400 -1.375900
      Read-in Center 2486 is at   0.000000  1.693400 -1.375900
      Read-in Center 2487 is at   0.000000  1.746300 -1.375900
      Read-in Center 2488 is at   0.000000  1.799200 -1.375900
      Read-in Center 2489 is at   0.000000  1.852100 -1.375900
      Read-in Center 2490 is at   0.000000  1.905000 -1.375900
      Read-in Center 2491 is at   0.000000  1.958000 -1.375900
      Read-in Center 2492 is at   0.000000  2.010900 -1.375900
      Read-in Center 2493 is at   0.000000  2.063800 -1.375900
      Read-in Center 2494 is at   0.000000  2.116700 -1.375900
      Read-in Center 2495 is at   0.000000  2.169600 -1.375900
      Read-in Center 2496 is at   0.000000  2.222500 -1.375900
      Read-in Center 2497 is at   0.000000  2.275500 -1.375900
      Read-in Center 2498 is at   0.000000  2.328400 -1.375900
      Read-in Center 2499 is at   0.000000  2.381300 -1.375900
      Read-in Center 2500 is at   0.000000  2.434200 -1.375900
      Read-in Center 2501 is at   0.000000  2.487100 -1.375900
      Read-in Center 2502 is at   0.000000  2.540000 -1.375900
      Read-in Center 2503 is at   0.000000  2.593000 -1.375900
      Read-in Center 2504 is at   0.000000 -2.645900 -1.322900
      Read-in Center 2505 is at   0.000000 -2.593000 -1.322900
      Read-in Center 2506 is at   0.000000 -2.540000 -1.322900
      Read-in Center 2507 is at   0.000000 -2.487100 -1.322900
      Read-in Center 2508 is at   0.000000 -2.434200 -1.322900
      Read-in Center 2509 is at   0.000000 -2.381300 -1.322900
      Read-in Center 2510 is at   0.000000 -2.328400 -1.322900
      Read-in Center 2511 is at   0.000000 -2.275500 -1.322900
      Read-in Center 2512 is at   0.000000 -2.222500 -1.322900
      Read-in Center 2513 is at   0.000000 -2.169600 -1.322900
      Read-in Center 2514 is at   0.000000 -2.116700 -1.322900
      Read-in Center 2515 is at   0.000000 -2.063800 -1.322900
      Read-in Center 2516 is at   0.000000 -2.010900 -1.322900
      Read-in Center 2517 is at   0.000000 -1.958000 -1.322900
      Read-in Center 2518 is at   0.000000 -1.905000 -1.322900
      Read-in Center 2519 is at   0.000000 -1.852100 -1.322900
      Read-in Center 2520 is at   0.000000 -1.799200 -1.322900
      Read-in Center 2521 is at   0.000000 -1.746300 -1.322900
      Read-in Center 2522 is at   0.000000 -1.693400 -1.322900
      Read-in Center 2523 is at   0.000000 -1.640400 -1.322900
      Read-in Center 2524 is at   0.000000 -1.587500 -1.322900
      Read-in Center 2525 is at   0.000000 -1.534600 -1.322900
      Read-in Center 2526 is at   0.000000 -1.481700 -1.322900
      Read-in Center 2527 is at   0.000000 -1.428800 -1.322900
      Read-in Center 2528 is at   0.000000 -1.375900 -1.322900
      Read-in Center 2529 is at   0.000000 -1.322900 -1.322900
      Read-in Center 2530 is at   0.000000 -1.270000 -1.322900
      Read-in Center 2531 is at   0.000000 -1.217100 -1.322900
      Read-in Center 2532 is at   0.000000 -1.164200 -1.322900
      Read-in Center 2533 is at   0.000000 -1.111300 -1.322900
      Read-in Center 2534 is at   0.000000 -1.058400 -1.322900
      Read-in Center 2535 is at   0.000000 -1.005400 -1.322900
      Read-in Center 2536 is at   0.000000 -0.952500 -1.322900
      Read-in Center 2537 is at   0.000000 -0.899600 -1.322900
      Read-in Center 2538 is at   0.000000 -0.846700 -1.322900
      Read-in Center 2539 is at   0.000000 -0.793800 -1.322900
      Read-in Center 2540 is at   0.000000 -0.740800 -1.322900
      Read-in Center 2541 is at   0.000000 -0.687900 -1.322900
      Read-in Center 2542 is at   0.000000 -0.635000 -1.322900
      Read-in Center 2543 is at   0.000000 -0.582100 -1.322900
      Read-in Center 2544 is at   0.000000 -0.529200 -1.322900
      Read-in Center 2545 is at   0.000000 -0.476300 -1.322900
      Read-in Center 2546 is at   0.000000 -0.423300 -1.322900
      Read-in Center 2547 is at   0.000000 -0.370400 -1.322900
      Read-in Center 2548 is at   0.000000 -0.317500 -1.322900
      Read-in Center 2549 is at   0.000000 -0.264600 -1.322900
      Read-in Center 2550 is at   0.000000 -0.211700 -1.322900
      Read-in Center 2551 is at   0.000000 -0.158800 -1.322900
      Read-in Center 2552 is at   0.000000 -0.105800 -1.322900
      Read-in Center 2553 is at   0.000000 -0.052900 -1.322900
      Read-in Center 2554 is at   0.000000  0.000000 -1.322900
      Read-in Center 2555 is at   0.000000  0.052900 -1.322900
      Read-in Center 2556 is at   0.000000  0.105800 -1.322900
      Read-in Center 2557 is at   0.000000  0.158800 -1.322900
      Read-in Center 2558 is at   0.000000  0.211700 -1.322900
      Read-in Center 2559 is at   0.000000  0.264600 -1.322900
      Read-in Center 2560 is at   0.000000  0.317500 -1.322900
      Read-in Center 2561 is at   0.000000  0.370400 -1.322900
      Read-in Center 2562 is at   0.000000  0.423300 -1.322900
      Read-in Center 2563 is at   0.000000  0.476300 -1.322900
      Read-in Center 2564 is at   0.000000  0.529200 -1.322900
      Read-in Center 2565 is at   0.000000  0.582100 -1.322900
      Read-in Center 2566 is at   0.000000  0.635000 -1.322900
      Read-in Center 2567 is at   0.000000  0.687900 -1.322900
      Read-in Center 2568 is at   0.000000  0.740800 -1.322900
      Read-in Center 2569 is at   0.000000  0.793800 -1.322900
      Read-in Center 2570 is at   0.000000  0.846700 -1.322900
      Read-in Center 2571 is at   0.000000  0.899600 -1.322900
      Read-in Center 2572 is at   0.000000  0.952500 -1.322900
      Read-in Center 2573 is at   0.000000  1.005400 -1.322900
      Read-in Center 2574 is at   0.000000  1.058400 -1.322900
      Read-in Center 2575 is at   0.000000  1.111300 -1.322900
      Read-in Center 2576 is at   0.000000  1.164200 -1.322900
      Read-in Center 2577 is at   0.000000  1.217100 -1.322900
      Read-in Center 2578 is at   0.000000  1.270000 -1.322900
      Read-in Center 2579 is at   0.000000  1.322900 -1.322900
      Read-in Center 2580 is at   0.000000  1.375900 -1.322900
      Read-in Center 2581 is at   0.000000  1.428800 -1.322900
      Read-in Center 2582 is at   0.000000  1.481700 -1.322900
      Read-in Center 2583 is at   0.000000  1.534600 -1.322900
      Read-in Center 2584 is at   0.000000  1.587500 -1.322900
      Read-in Center 2585 is at   0.000000  1.640400 -1.322900
      Read-in Center 2586 is at   0.000000  1.693400 -1.322900
      Read-in Center 2587 is at   0.000000  1.746300 -1.322900
      Read-in Center 2588 is at   0.000000  1.799200 -1.322900
      Read-in Center 2589 is at   0.000000  1.852100 -1.322900
      Read-in Center 2590 is at   0.000000  1.905000 -1.322900
      Read-in Center 2591 is at   0.000000  1.958000 -1.322900
      Read-in Center 2592 is at   0.000000  2.010900 -1.322900
      Read-in Center 2593 is at   0.000000  2.063800 -1.322900
      Read-in Center 2594 is at   0.000000  2.116700 -1.322900
      Read-in Center 2595 is at   0.000000  2.169600 -1.322900
      Read-in Center 2596 is at   0.000000  2.222500 -1.322900
      Read-in Center 2597 is at   0.000000  2.275500 -1.322900
      Read-in Center 2598 is at   0.000000  2.328400 -1.322900
      Read-in Center 2599 is at   0.000000  2.381300 -1.322900
      Read-in Center 2600 is at   0.000000  2.434200 -1.322900
      Read-in Center 2601 is at   0.000000  2.487100 -1.322900
      Read-in Center 2602 is at   0.000000  2.540000 -1.322900
      Read-in Center 2603 is at   0.000000  2.593000 -1.322900
      Read-in Center 2604 is at   0.000000 -2.645900 -1.270000
      Read-in Center 2605 is at   0.000000 -2.593000 -1.270000
      Read-in Center 2606 is at   0.000000 -2.540000 -1.270000
      Read-in Center 2607 is at   0.000000 -2.487100 -1.270000
      Read-in Center 2608 is at   0.000000 -2.434200 -1.270000
      Read-in Center 2609 is at   0.000000 -2.381300 -1.270000
      Read-in Center 2610 is at   0.000000 -2.328400 -1.270000
      Read-in Center 2611 is at   0.000000 -2.275500 -1.270000
      Read-in Center 2612 is at   0.000000 -2.222500 -1.270000
      Read-in Center 2613 is at   0.000000 -2.169600 -1.270000
      Read-in Center 2614 is at   0.000000 -2.116700 -1.270000
      Read-in Center 2615 is at   0.000000 -2.063800 -1.270000
      Read-in Center 2616 is at   0.000000 -2.010900 -1.270000
      Read-in Center 2617 is at   0.000000 -1.958000 -1.270000
      Read-in Center 2618 is at   0.000000 -1.905000 -1.270000
      Read-in Center 2619 is at   0.000000 -1.852100 -1.270000
      Read-in Center 2620 is at   0.000000 -1.799200 -1.270000
      Read-in Center 2621 is at   0.000000 -1.746300 -1.270000
      Read-in Center 2622 is at   0.000000 -1.693400 -1.270000
      Read-in Center 2623 is at   0.000000 -1.640400 -1.270000
      Read-in Center 2624 is at   0.000000 -1.587500 -1.270000
      Read-in Center 2625 is at   0.000000 -1.534600 -1.270000
      Read-in Center 2626 is at   0.000000 -1.481700 -1.270000
      Read-in Center 2627 is at   0.000000 -1.428800 -1.270000
      Read-in Center 2628 is at   0.000000 -1.375900 -1.270000
      Read-in Center 2629 is at   0.000000 -1.322900 -1.270000
      Read-in Center 2630 is at   0.000000 -1.270000 -1.270000
      Read-in Center 2631 is at   0.000000 -1.217100 -1.270000
      Read-in Center 2632 is at   0.000000 -1.164200 -1.270000
      Read-in Center 2633 is at   0.000000 -1.111300 -1.270000
      Read-in Center 2634 is at   0.000000 -1.058400 -1.270000
      Read-in Center 2635 is at   0.000000 -1.005400 -1.270000
      Read-in Center 2636 is at   0.000000 -0.952500 -1.270000
      Read-in Center 2637 is at   0.000000 -0.899600 -1.270000
      Read-in Center 2638 is at   0.000000 -0.846700 -1.270000
      Read-in Center 2639 is at   0.000000 -0.793800 -1.270000
      Read-in Center 2640 is at   0.000000 -0.740800 -1.270000
      Read-in Center 2641 is at   0.000000 -0.687900 -1.270000
      Read-in Center 2642 is at   0.000000 -0.635000 -1.270000
      Read-in Center 2643 is at   0.000000 -0.582100 -1.270000
      Read-in Center 2644 is at   0.000000 -0.529200 -1.270000
      Read-in Center 2645 is at   0.000000 -0.476300 -1.270000
      Read-in Center 2646 is at   0.000000 -0.423300 -1.270000
      Read-in Center 2647 is at   0.000000 -0.370400 -1.270000
      Read-in Center 2648 is at   0.000000 -0.317500 -1.270000
      Read-in Center 2649 is at   0.000000 -0.264600 -1.270000
      Read-in Center 2650 is at   0.000000 -0.211700 -1.270000
      Read-in Center 2651 is at   0.000000 -0.158800 -1.270000
      Read-in Center 2652 is at   0.000000 -0.105800 -1.270000
      Read-in Center 2653 is at   0.000000 -0.052900 -1.270000
      Read-in Center 2654 is at   0.000000  0.000000 -1.270000
      Read-in Center 2655 is at   0.000000  0.052900 -1.270000
      Read-in Center 2656 is at   0.000000  0.105800 -1.270000
      Read-in Center 2657 is at   0.000000  0.158800 -1.270000
      Read-in Center 2658 is at   0.000000  0.211700 -1.270000
      Read-in Center 2659 is at   0.000000  0.264600 -1.270000
      Read-in Center 2660 is at   0.000000  0.317500 -1.270000
      Read-in Center 2661 is at   0.000000  0.370400 -1.270000
      Read-in Center 2662 is at   0.000000  0.423300 -1.270000
      Read-in Center 2663 is at   0.000000  0.476300 -1.270000
      Read-in Center 2664 is at   0.000000  0.529200 -1.270000
      Read-in Center 2665 is at   0.000000  0.582100 -1.270000
      Read-in Center 2666 is at   0.000000  0.635000 -1.270000
      Read-in Center 2667 is at   0.000000  0.687900 -1.270000
      Read-in Center 2668 is at   0.000000  0.740800 -1.270000
      Read-in Center 2669 is at   0.000000  0.793800 -1.270000
      Read-in Center 2670 is at   0.000000  0.846700 -1.270000
      Read-in Center 2671 is at   0.000000  0.899600 -1.270000
      Read-in Center 2672 is at   0.000000  0.952500 -1.270000
      Read-in Center 2673 is at   0.000000  1.005400 -1.270000
      Read-in Center 2674 is at   0.000000  1.058400 -1.270000
      Read-in Center 2675 is at   0.000000  1.111300 -1.270000
      Read-in Center 2676 is at   0.000000  1.164200 -1.270000
      Read-in Center 2677 is at   0.000000  1.217100 -1.270000
      Read-in Center 2678 is at   0.000000  1.270000 -1.270000
      Read-in Center 2679 is at   0.000000  1.322900 -1.270000
      Read-in Center 2680 is at   0.000000  1.375900 -1.270000
      Read-in Center 2681 is at   0.000000  1.428800 -1.270000
      Read-in Center 2682 is at   0.000000  1.481700 -1.270000
      Read-in Center 2683 is at   0.000000  1.534600 -1.270000
      Read-in Center 2684 is at   0.000000  1.587500 -1.270000
      Read-in Center 2685 is at   0.000000  1.640400 -1.270000
      Read-in Center 2686 is at   0.000000  1.693400 -1.270000
      Read-in Center 2687 is at   0.000000  1.746300 -1.270000
      Read-in Center 2688 is at   0.000000  1.799200 -1.270000
      Read-in Center 2689 is at   0.000000  1.852100 -1.270000
      Read-in Center 2690 is at   0.000000  1.905000 -1.270000
      Read-in Center 2691 is at   0.000000  1.958000 -1.270000
      Read-in Center 2692 is at   0.000000  2.010900 -1.270000
      Read-in Center 2693 is at   0.000000  2.063800 -1.270000
      Read-in Center 2694 is at   0.000000  2.116700 -1.270000
      Read-in Center 2695 is at   0.000000  2.169600 -1.270000
      Read-in Center 2696 is at   0.000000  2.222500 -1.270000
      Read-in Center 2697 is at   0.000000  2.275500 -1.270000
      Read-in Center 2698 is at   0.000000  2.328400 -1.270000
      Read-in Center 2699 is at   0.000000  2.381300 -1.270000
      Read-in Center 2700 is at   0.000000  2.434200 -1.270000
      Read-in Center 2701 is at   0.000000  2.487100 -1.270000
      Read-in Center 2702 is at   0.000000  2.540000 -1.270000
      Read-in Center 2703 is at   0.000000  2.593000 -1.270000
      Read-in Center 2704 is at   0.000000 -2.645900 -1.217100
      Read-in Center 2705 is at   0.000000 -2.593000 -1.217100
      Read-in Center 2706 is at   0.000000 -2.540000 -1.217100
      Read-in Center 2707 is at   0.000000 -2.487100 -1.217100
      Read-in Center 2708 is at   0.000000 -2.434200 -1.217100
      Read-in Center 2709 is at   0.000000 -2.381300 -1.217100
      Read-in Center 2710 is at   0.000000 -2.328400 -1.217100
      Read-in Center 2711 is at   0.000000 -2.275500 -1.217100
      Read-in Center 2712 is at   0.000000 -2.222500 -1.217100
      Read-in Center 2713 is at   0.000000 -2.169600 -1.217100
      Read-in Center 2714 is at   0.000000 -2.116700 -1.217100
      Read-in Center 2715 is at   0.000000 -2.063800 -1.217100
      Read-in Center 2716 is at   0.000000 -2.010900 -1.217100
      Read-in Center 2717 is at   0.000000 -1.958000 -1.217100
      Read-in Center 2718 is at   0.000000 -1.905000 -1.217100
      Read-in Center 2719 is at   0.000000 -1.852100 -1.217100
      Read-in Center 2720 is at   0.000000 -1.799200 -1.217100
      Read-in Center 2721 is at   0.000000 -1.746300 -1.217100
      Read-in Center 2722 is at   0.000000 -1.693400 -1.217100
      Read-in Center 2723 is at   0.000000 -1.640400 -1.217100
      Read-in Center 2724 is at   0.000000 -1.587500 -1.217100
      Read-in Center 2725 is at   0.000000 -1.534600 -1.217100
      Read-in Center 2726 is at   0.000000 -1.481700 -1.217100
      Read-in Center 2727 is at   0.000000 -1.428800 -1.217100
      Read-in Center 2728 is at   0.000000 -1.375900 -1.217100
      Read-in Center 2729 is at   0.000000 -1.322900 -1.217100
      Read-in Center 2730 is at   0.000000 -1.270000 -1.217100
      Read-in Center 2731 is at   0.000000 -1.217100 -1.217100
      Read-in Center 2732 is at   0.000000 -1.164200 -1.217100
      Read-in Center 2733 is at   0.000000 -1.111300 -1.217100
      Read-in Center 2734 is at   0.000000 -1.058400 -1.217100
      Read-in Center 2735 is at   0.000000 -1.005400 -1.217100
      Read-in Center 2736 is at   0.000000 -0.952500 -1.217100
      Read-in Center 2737 is at   0.000000 -0.899600 -1.217100
      Read-in Center 2738 is at   0.000000 -0.846700 -1.217100
      Read-in Center 2739 is at   0.000000 -0.793800 -1.217100
      Read-in Center 2740 is at   0.000000 -0.740800 -1.217100
      Read-in Center 2741 is at   0.000000 -0.687900 -1.217100
      Read-in Center 2742 is at   0.000000 -0.635000 -1.217100
      Read-in Center 2743 is at   0.000000 -0.582100 -1.217100
      Read-in Center 2744 is at   0.000000 -0.529200 -1.217100
      Read-in Center 2745 is at   0.000000 -0.476300 -1.217100
      Read-in Center 2746 is at   0.000000 -0.423300 -1.217100
      Read-in Center 2747 is at   0.000000 -0.370400 -1.217100
      Read-in Center 2748 is at   0.000000 -0.317500 -1.217100
      Read-in Center 2749 is at   0.000000 -0.264600 -1.217100
      Read-in Center 2750 is at   0.000000 -0.211700 -1.217100
      Read-in Center 2751 is at   0.000000 -0.158800 -1.217100
      Read-in Center 2752 is at   0.000000 -0.105800 -1.217100
      Read-in Center 2753 is at   0.000000 -0.052900 -1.217100
      Read-in Center 2754 is at   0.000000  0.000000 -1.217100
      Read-in Center 2755 is at   0.000000  0.052900 -1.217100
      Read-in Center 2756 is at   0.000000  0.105800 -1.217100
      Read-in Center 2757 is at   0.000000  0.158800 -1.217100
      Read-in Center 2758 is at   0.000000  0.211700 -1.217100
      Read-in Center 2759 is at   0.000000  0.264600 -1.217100
      Read-in Center 2760 is at   0.000000  0.317500 -1.217100
      Read-in Center 2761 is at   0.000000  0.370400 -1.217100
      Read-in Center 2762 is at   0.000000  0.423300 -1.217100
      Read-in Center 2763 is at   0.000000  0.476300 -1.217100
      Read-in Center 2764 is at   0.000000  0.529200 -1.217100
      Read-in Center 2765 is at   0.000000  0.582100 -1.217100
      Read-in Center 2766 is at   0.000000  0.635000 -1.217100
      Read-in Center 2767 is at   0.000000  0.687900 -1.217100
      Read-in Center 2768 is at   0.000000  0.740800 -1.217100
      Read-in Center 2769 is at   0.000000  0.793800 -1.217100
      Read-in Center 2770 is at   0.000000  0.846700 -1.217100
      Read-in Center 2771 is at   0.000000  0.899600 -1.217100
      Read-in Center 2772 is at   0.000000  0.952500 -1.217100
      Read-in Center 2773 is at   0.000000  1.005400 -1.217100
      Read-in Center 2774 is at   0.000000  1.058400 -1.217100
      Read-in Center 2775 is at   0.000000  1.111300 -1.217100
      Read-in Center 2776 is at   0.000000  1.164200 -1.217100
      Read-in Center 2777 is at   0.000000  1.217100 -1.217100
      Read-in Center 2778 is at   0.000000  1.270000 -1.217100
      Read-in Center 2779 is at   0.000000  1.322900 -1.217100
      Read-in Center 2780 is at   0.000000  1.375900 -1.217100
      Read-in Center 2781 is at   0.000000  1.428800 -1.217100
      Read-in Center 2782 is at   0.000000  1.481700 -1.217100
      Read-in Center 2783 is at   0.000000  1.534600 -1.217100
      Read-in Center 2784 is at   0.000000  1.587500 -1.217100
      Read-in Center 2785 is at   0.000000  1.640400 -1.217100
      Read-in Center 2786 is at   0.000000  1.693400 -1.217100
      Read-in Center 2787 is at   0.000000  1.746300 -1.217100
      Read-in Center 2788 is at   0.000000  1.799200 -1.217100
      Read-in Center 2789 is at   0.000000  1.852100 -1.217100
      Read-in Center 2790 is at   0.000000  1.905000 -1.217100
      Read-in Center 2791 is at   0.000000  1.958000 -1.217100
      Read-in Center 2792 is at   0.000000  2.010900 -1.217100
      Read-in Center 2793 is at   0.000000  2.063800 -1.217100
      Read-in Center 2794 is at   0.000000  2.116700 -1.217100
      Read-in Center 2795 is at   0.000000  2.169600 -1.217100
      Read-in Center 2796 is at   0.000000  2.222500 -1.217100
      Read-in Center 2797 is at   0.000000  2.275500 -1.217100
      Read-in Center 2798 is at   0.000000  2.328400 -1.217100
      Read-in Center 2799 is at   0.000000  2.381300 -1.217100
      Read-in Center 2800 is at   0.000000  2.434200 -1.217100
      Read-in Center 2801 is at   0.000000  2.487100 -1.217100
      Read-in Center 2802 is at   0.000000  2.540000 -1.217100
      Read-in Center 2803 is at   0.000000  2.593000 -1.217100
      Read-in Center 2804 is at   0.000000 -2.645900 -1.164200
      Read-in Center 2805 is at   0.000000 -2.593000 -1.164200
      Read-in Center 2806 is at   0.000000 -2.540000 -1.164200
      Read-in Center 2807 is at   0.000000 -2.487100 -1.164200
      Read-in Center 2808 is at   0.000000 -2.434200 -1.164200
      Read-in Center 2809 is at   0.000000 -2.381300 -1.164200
      Read-in Center 2810 is at   0.000000 -2.328400 -1.164200
      Read-in Center 2811 is at   0.000000 -2.275500 -1.164200
      Read-in Center 2812 is at   0.000000 -2.222500 -1.164200
      Read-in Center 2813 is at   0.000000 -2.169600 -1.164200
      Read-in Center 2814 is at   0.000000 -2.116700 -1.164200
      Read-in Center 2815 is at   0.000000 -2.063800 -1.164200
      Read-in Center 2816 is at   0.000000 -2.010900 -1.164200
      Read-in Center 2817 is at   0.000000 -1.958000 -1.164200
      Read-in Center 2818 is at   0.000000 -1.905000 -1.164200
      Read-in Center 2819 is at   0.000000 -1.852100 -1.164200
      Read-in Center 2820 is at   0.000000 -1.799200 -1.164200
      Read-in Center 2821 is at   0.000000 -1.746300 -1.164200
      Read-in Center 2822 is at   0.000000 -1.693400 -1.164200
      Read-in Center 2823 is at   0.000000 -1.640400 -1.164200
      Read-in Center 2824 is at   0.000000 -1.587500 -1.164200
      Read-in Center 2825 is at   0.000000 -1.534600 -1.164200
      Read-in Center 2826 is at   0.000000 -1.481700 -1.164200
      Read-in Center 2827 is at   0.000000 -1.428800 -1.164200
      Read-in Center 2828 is at   0.000000 -1.375900 -1.164200
      Read-in Center 2829 is at   0.000000 -1.322900 -1.164200
      Read-in Center 2830 is at   0.000000 -1.270000 -1.164200
      Read-in Center 2831 is at   0.000000 -1.217100 -1.164200
      Read-in Center 2832 is at   0.000000 -1.164200 -1.164200
      Read-in Center 2833 is at   0.000000 -1.111300 -1.164200
      Read-in Center 2834 is at   0.000000 -1.058400 -1.164200
      Read-in Center 2835 is at   0.000000 -1.005400 -1.164200
      Read-in Center 2836 is at   0.000000 -0.952500 -1.164200
      Read-in Center 2837 is at   0.000000 -0.899600 -1.164200
      Read-in Center 2838 is at   0.000000 -0.846700 -1.164200
      Read-in Center 2839 is at   0.000000 -0.793800 -1.164200
      Read-in Center 2840 is at   0.000000 -0.740800 -1.164200
      Read-in Center 2841 is at   0.000000 -0.687900 -1.164200
      Read-in Center 2842 is at   0.000000 -0.635000 -1.164200
      Read-in Center 2843 is at   0.000000 -0.582100 -1.164200
      Read-in Center 2844 is at   0.000000 -0.529200 -1.164200
      Read-in Center 2845 is at   0.000000 -0.476300 -1.164200
      Read-in Center 2846 is at   0.000000 -0.423300 -1.164200
      Read-in Center 2847 is at   0.000000 -0.370400 -1.164200
      Read-in Center 2848 is at   0.000000 -0.317500 -1.164200
      Read-in Center 2849 is at   0.000000 -0.264600 -1.164200
      Read-in Center 2850 is at   0.000000 -0.211700 -1.164200
      Read-in Center 2851 is at   0.000000 -0.158800 -1.164200
      Read-in Center 2852 is at   0.000000 -0.105800 -1.164200
      Read-in Center 2853 is at   0.000000 -0.052900 -1.164200
      Read-in Center 2854 is at   0.000000  0.000000 -1.164200
      Read-in Center 2855 is at   0.000000  0.052900 -1.164200
      Read-in Center 2856 is at   0.000000  0.105800 -1.164200
      Read-in Center 2857 is at   0.000000  0.158800 -1.164200
      Read-in Center 2858 is at   0.000000  0.211700 -1.164200
      Read-in Center 2859 is at   0.000000  0.264600 -1.164200
      Read-in Center 2860 is at   0.000000  0.317500 -1.164200
      Read-in Center 2861 is at   0.000000  0.370400 -1.164200
      Read-in Center 2862 is at   0.000000  0.423300 -1.164200
      Read-in Center 2863 is at   0.000000  0.476300 -1.164200
      Read-in Center 2864 is at   0.000000  0.529200 -1.164200
      Read-in Center 2865 is at   0.000000  0.582100 -1.164200
      Read-in Center 2866 is at   0.000000  0.635000 -1.164200
      Read-in Center 2867 is at   0.000000  0.687900 -1.164200
      Read-in Center 2868 is at   0.000000  0.740800 -1.164200
      Read-in Center 2869 is at   0.000000  0.793800 -1.164200
      Read-in Center 2870 is at   0.000000  0.846700 -1.164200
      Read-in Center 2871 is at   0.000000  0.899600 -1.164200
      Read-in Center 2872 is at   0.000000  0.952500 -1.164200
      Read-in Center 2873 is at   0.000000  1.005400 -1.164200
      Read-in Center 2874 is at   0.000000  1.058400 -1.164200
      Read-in Center 2875 is at   0.000000  1.111300 -1.164200
      Read-in Center 2876 is at   0.000000  1.164200 -1.164200
      Read-in Center 2877 is at   0.000000  1.217100 -1.164200
      Read-in Center 2878 is at   0.000000  1.270000 -1.164200
      Read-in Center 2879 is at   0.000000  1.322900 -1.164200
      Read-in Center 2880 is at   0.000000  1.375900 -1.164200
      Read-in Center 2881 is at   0.000000  1.428800 -1.164200
      Read-in Center 2882 is at   0.000000  1.481700 -1.164200
      Read-in Center 2883 is at   0.000000  1.534600 -1.164200
      Read-in Center 2884 is at   0.000000  1.587500 -1.164200
      Read-in Center 2885 is at   0.000000  1.640400 -1.164200
      Read-in Center 2886 is at   0.000000  1.693400 -1.164200
      Read-in Center 2887 is at   0.000000  1.746300 -1.164200
      Read-in Center 2888 is at   0.000000  1.799200 -1.164200
      Read-in Center 2889 is at   0.000000  1.852100 -1.164200
      Read-in Center 2890 is at   0.000000  1.905000 -1.164200
      Read-in Center 2891 is at   0.000000  1.958000 -1.164200
      Read-in Center 2892 is at   0.000000  2.010900 -1.164200
      Read-in Center 2893 is at   0.000000  2.063800 -1.164200
      Read-in Center 2894 is at   0.000000  2.116700 -1.164200
      Read-in Center 2895 is at   0.000000  2.169600 -1.164200
      Read-in Center 2896 is at   0.000000  2.222500 -1.164200
      Read-in Center 2897 is at   0.000000  2.275500 -1.164200
      Read-in Center 2898 is at   0.000000  2.328400 -1.164200
      Read-in Center 2899 is at   0.000000  2.381300 -1.164200
      Read-in Center 2900 is at   0.000000  2.434200 -1.164200
      Read-in Center 2901 is at   0.000000  2.487100 -1.164200
      Read-in Center 2902 is at   0.000000  2.540000 -1.164200
      Read-in Center 2903 is at   0.000000  2.593000 -1.164200
      Read-in Center 2904 is at   0.000000 -2.645900 -1.111300
      Read-in Center 2905 is at   0.000000 -2.593000 -1.111300
      Read-in Center 2906 is at   0.000000 -2.540000 -1.111300
      Read-in Center 2907 is at   0.000000 -2.487100 -1.111300
      Read-in Center 2908 is at   0.000000 -2.434200 -1.111300
      Read-in Center 2909 is at   0.000000 -2.381300 -1.111300
      Read-in Center 2910 is at   0.000000 -2.328400 -1.111300
      Read-in Center 2911 is at   0.000000 -2.275500 -1.111300
      Read-in Center 2912 is at   0.000000 -2.222500 -1.111300
      Read-in Center 2913 is at   0.000000 -2.169600 -1.111300
      Read-in Center 2914 is at   0.000000 -2.116700 -1.111300
      Read-in Center 2915 is at   0.000000 -2.063800 -1.111300
      Read-in Center 2916 is at   0.000000 -2.010900 -1.111300
      Read-in Center 2917 is at   0.000000 -1.958000 -1.111300
      Read-in Center 2918 is at   0.000000 -1.905000 -1.111300
      Read-in Center 2919 is at   0.000000 -1.852100 -1.111300
      Read-in Center 2920 is at   0.000000 -1.799200 -1.111300
      Read-in Center 2921 is at   0.000000 -1.746300 -1.111300
      Read-in Center 2922 is at   0.000000 -1.693400 -1.111300
      Read-in Center 2923 is at   0.000000 -1.640400 -1.111300
      Read-in Center 2924 is at   0.000000 -1.587500 -1.111300
      Read-in Center 2925 is at   0.000000 -1.534600 -1.111300
      Read-in Center 2926 is at   0.000000 -1.481700 -1.111300
      Read-in Center 2927 is at   0.000000 -1.428800 -1.111300
      Read-in Center 2928 is at   0.000000 -1.375900 -1.111300
      Read-in Center 2929 is at   0.000000 -1.322900 -1.111300
      Read-in Center 2930 is at   0.000000 -1.270000 -1.111300
      Read-in Center 2931 is at   0.000000 -1.217100 -1.111300
      Read-in Center 2932 is at   0.000000 -1.164200 -1.111300
      Read-in Center 2933 is at   0.000000 -1.111300 -1.111300
      Read-in Center 2934 is at   0.000000 -1.058400 -1.111300
      Read-in Center 2935 is at   0.000000 -1.005400 -1.111300
      Read-in Center 2936 is at   0.000000 -0.952500 -1.111300
      Read-in Center 2937 is at   0.000000 -0.899600 -1.111300
      Read-in Center 2938 is at   0.000000 -0.846700 -1.111300
      Read-in Center 2939 is at   0.000000 -0.793800 -1.111300
      Read-in Center 2940 is at   0.000000 -0.740800 -1.111300
      Read-in Center 2941 is at   0.000000 -0.687900 -1.111300
      Read-in Center 2942 is at   0.000000 -0.635000 -1.111300
      Read-in Center 2943 is at   0.000000 -0.582100 -1.111300
      Read-in Center 2944 is at   0.000000 -0.529200 -1.111300
      Read-in Center 2945 is at   0.000000 -0.476300 -1.111300
      Read-in Center 2946 is at   0.000000 -0.423300 -1.111300
      Read-in Center 2947 is at   0.000000 -0.370400 -1.111300
      Read-in Center 2948 is at   0.000000 -0.317500 -1.111300
      Read-in Center 2949 is at   0.000000 -0.264600 -1.111300
      Read-in Center 2950 is at   0.000000 -0.211700 -1.111300
      Read-in Center 2951 is at   0.000000 -0.158800 -1.111300
      Read-in Center 2952 is at   0.000000 -0.105800 -1.111300
      Read-in Center 2953 is at   0.000000 -0.052900 -1.111300
      Read-in Center 2954 is at   0.000000  0.000000 -1.111300
      Read-in Center 2955 is at   0.000000  0.052900 -1.111300
      Read-in Center 2956 is at   0.000000  0.105800 -1.111300
      Read-in Center 2957 is at   0.000000  0.158800 -1.111300
      Read-in Center 2958 is at   0.000000  0.211700 -1.111300
      Read-in Center 2959 is at   0.000000  0.264600 -1.111300
      Read-in Center 2960 is at   0.000000  0.317500 -1.111300
      Read-in Center 2961 is at   0.000000  0.370400 -1.111300
      Read-in Center 2962 is at   0.000000  0.423300 -1.111300
      Read-in Center 2963 is at   0.000000  0.476300 -1.111300
      Read-in Center 2964 is at   0.000000  0.529200 -1.111300
      Read-in Center 2965 is at   0.000000  0.582100 -1.111300
      Read-in Center 2966 is at   0.000000  0.635000 -1.111300
      Read-in Center 2967 is at   0.000000  0.687900 -1.111300
      Read-in Center 2968 is at   0.000000  0.740800 -1.111300
      Read-in Center 2969 is at   0.000000  0.793800 -1.111300
      Read-in Center 2970 is at   0.000000  0.846700 -1.111300
      Read-in Center 2971 is at   0.000000  0.899600 -1.111300
      Read-in Center 2972 is at   0.000000  0.952500 -1.111300
      Read-in Center 2973 is at   0.000000  1.005400 -1.111300
      Read-in Center 2974 is at   0.000000  1.058400 -1.111300
      Read-in Center 2975 is at   0.000000  1.111300 -1.111300
      Read-in Center 2976 is at   0.000000  1.164200 -1.111300
      Read-in Center 2977 is at   0.000000  1.217100 -1.111300
      Read-in Center 2978 is at   0.000000  1.270000 -1.111300
      Read-in Center 2979 is at   0.000000  1.322900 -1.111300
      Read-in Center 2980 is at   0.000000  1.375900 -1.111300
      Read-in Center 2981 is at   0.000000  1.428800 -1.111300
      Read-in Center 2982 is at   0.000000  1.481700 -1.111300
      Read-in Center 2983 is at   0.000000  1.534600 -1.111300
      Read-in Center 2984 is at   0.000000  1.587500 -1.111300
      Read-in Center 2985 is at   0.000000  1.640400 -1.111300
      Read-in Center 2986 is at   0.000000  1.693400 -1.111300
      Read-in Center 2987 is at   0.000000  1.746300 -1.111300
      Read-in Center 2988 is at   0.000000  1.799200 -1.111300
      Read-in Center 2989 is at   0.000000  1.852100 -1.111300
      Read-in Center 2990 is at   0.000000  1.905000 -1.111300
      Read-in Center 2991 is at   0.000000  1.958000 -1.111300
      Read-in Center 2992 is at   0.000000  2.010900 -1.111300
      Read-in Center 2993 is at   0.000000  2.063800 -1.111300
      Read-in Center 2994 is at   0.000000  2.116700 -1.111300
      Read-in Center 2995 is at   0.000000  2.169600 -1.111300
      Read-in Center 2996 is at   0.000000  2.222500 -1.111300
      Read-in Center 2997 is at   0.000000  2.275500 -1.111300
      Read-in Center 2998 is at   0.000000  2.328400 -1.111300
      Read-in Center 2999 is at   0.000000  2.381300 -1.111300
      Read-in Center 3000 is at   0.000000  2.434200 -1.111300
      Read-in Center 3001 is at   0.000000  2.487100 -1.111300
      Read-in Center 3002 is at   0.000000  2.540000 -1.111300
      Read-in Center 3003 is at   0.000000  2.593000 -1.111300
      Read-in Center 3004 is at   0.000000 -2.645900 -1.058400
      Read-in Center 3005 is at   0.000000 -2.593000 -1.058400
      Read-in Center 3006 is at   0.000000 -2.540000 -1.058400
      Read-in Center 3007 is at   0.000000 -2.487100 -1.058400
      Read-in Center 3008 is at   0.000000 -2.434200 -1.058400
      Read-in Center 3009 is at   0.000000 -2.381300 -1.058400
      Read-in Center 3010 is at   0.000000 -2.328400 -1.058400
      Read-in Center 3011 is at   0.000000 -2.275500 -1.058400
      Read-in Center 3012 is at   0.000000 -2.222500 -1.058400
      Read-in Center 3013 is at   0.000000 -2.169600 -1.058400
      Read-in Center 3014 is at   0.000000 -2.116700 -1.058400
      Read-in Center 3015 is at   0.000000 -2.063800 -1.058400
      Read-in Center 3016 is at   0.000000 -2.010900 -1.058400
      Read-in Center 3017 is at   0.000000 -1.958000 -1.058400
      Read-in Center 3018 is at   0.000000 -1.905000 -1.058400
      Read-in Center 3019 is at   0.000000 -1.852100 -1.058400
      Read-in Center 3020 is at   0.000000 -1.799200 -1.058400
      Read-in Center 3021 is at   0.000000 -1.746300 -1.058400
      Read-in Center 3022 is at   0.000000 -1.693400 -1.058400
      Read-in Center 3023 is at   0.000000 -1.640400 -1.058400
      Read-in Center 3024 is at   0.000000 -1.587500 -1.058400
      Read-in Center 3025 is at   0.000000 -1.534600 -1.058400
      Read-in Center 3026 is at   0.000000 -1.481700 -1.058400
      Read-in Center 3027 is at   0.000000 -1.428800 -1.058400
      Read-in Center 3028 is at   0.000000 -1.375900 -1.058400
      Read-in Center 3029 is at   0.000000 -1.322900 -1.058400
      Read-in Center 3030 is at   0.000000 -1.270000 -1.058400
      Read-in Center 3031 is at   0.000000 -1.217100 -1.058400
      Read-in Center 3032 is at   0.000000 -1.164200 -1.058400
      Read-in Center 3033 is at   0.000000 -1.111300 -1.058400
      Read-in Center 3034 is at   0.000000 -1.058400 -1.058400
      Read-in Center 3035 is at   0.000000 -1.005400 -1.058400
      Read-in Center 3036 is at   0.000000 -0.952500 -1.058400
      Read-in Center 3037 is at   0.000000 -0.899600 -1.058400
      Read-in Center 3038 is at   0.000000 -0.846700 -1.058400
      Read-in Center 3039 is at   0.000000 -0.793800 -1.058400
      Read-in Center 3040 is at   0.000000 -0.740800 -1.058400
      Read-in Center 3041 is at   0.000000 -0.687900 -1.058400
      Read-in Center 3042 is at   0.000000 -0.635000 -1.058400
      Read-in Center 3043 is at   0.000000 -0.582100 -1.058400
      Read-in Center 3044 is at   0.000000 -0.529200 -1.058400
      Read-in Center 3045 is at   0.000000 -0.476300 -1.058400
      Read-in Center 3046 is at   0.000000 -0.423300 -1.058400
      Read-in Center 3047 is at   0.000000 -0.370400 -1.058400
      Read-in Center 3048 is at   0.000000 -0.317500 -1.058400
      Read-in Center 3049 is at   0.000000 -0.264600 -1.058400
      Read-in Center 3050 is at   0.000000 -0.211700 -1.058400
      Read-in Center 3051 is at   0.000000 -0.158800 -1.058400
      Read-in Center 3052 is at   0.000000 -0.105800 -1.058400
      Read-in Center 3053 is at   0.000000 -0.052900 -1.058400
      Read-in Center 3054 is at   0.000000  0.000000 -1.058400
      Read-in Center 3055 is at   0.000000  0.052900 -1.058400
      Read-in Center 3056 is at   0.000000  0.105800 -1.058400
      Read-in Center 3057 is at   0.000000  0.158800 -1.058400
      Read-in Center 3058 is at   0.000000  0.211700 -1.058400
      Read-in Center 3059 is at   0.000000  0.264600 -1.058400
      Read-in Center 3060 is at   0.000000  0.317500 -1.058400
      Read-in Center 3061 is at   0.000000  0.370400 -1.058400
      Read-in Center 3062 is at   0.000000  0.423300 -1.058400
      Read-in Center 3063 is at   0.000000  0.476300 -1.058400
      Read-in Center 3064 is at   0.000000  0.529200 -1.058400
      Read-in Center 3065 is at   0.000000  0.582100 -1.058400
      Read-in Center 3066 is at   0.000000  0.635000 -1.058400
      Read-in Center 3067 is at   0.000000  0.687900 -1.058400
      Read-in Center 3068 is at   0.000000  0.740800 -1.058400
      Read-in Center 3069 is at   0.000000  0.793800 -1.058400
      Read-in Center 3070 is at   0.000000  0.846700 -1.058400
      Read-in Center 3071 is at   0.000000  0.899600 -1.058400
      Read-in Center 3072 is at   0.000000  0.952500 -1.058400
      Read-in Center 3073 is at   0.000000  1.005400 -1.058400
      Read-in Center 3074 is at   0.000000  1.058400 -1.058400
      Read-in Center 3075 is at   0.000000  1.111300 -1.058400
      Read-in Center 3076 is at   0.000000  1.164200 -1.058400
      Read-in Center 3077 is at   0.000000  1.217100 -1.058400
      Read-in Center 3078 is at   0.000000  1.270000 -1.058400
      Read-in Center 3079 is at   0.000000  1.322900 -1.058400
      Read-in Center 3080 is at   0.000000  1.375900 -1.058400
      Read-in Center 3081 is at   0.000000  1.428800 -1.058400
      Read-in Center 3082 is at   0.000000  1.481700 -1.058400
      Read-in Center 3083 is at   0.000000  1.534600 -1.058400
      Read-in Center 3084 is at   0.000000  1.587500 -1.058400
      Read-in Center 3085 is at   0.000000  1.640400 -1.058400
      Read-in Center 3086 is at   0.000000  1.693400 -1.058400
      Read-in Center 3087 is at   0.000000  1.746300 -1.058400
      Read-in Center 3088 is at   0.000000  1.799200 -1.058400
      Read-in Center 3089 is at   0.000000  1.852100 -1.058400
      Read-in Center 3090 is at   0.000000  1.905000 -1.058400
      Read-in Center 3091 is at   0.000000  1.958000 -1.058400
      Read-in Center 3092 is at   0.000000  2.010900 -1.058400
      Read-in Center 3093 is at   0.000000  2.063800 -1.058400
      Read-in Center 3094 is at   0.000000  2.116700 -1.058400
      Read-in Center 3095 is at   0.000000  2.169600 -1.058400
      Read-in Center 3096 is at   0.000000  2.222500 -1.058400
      Read-in Center 3097 is at   0.000000  2.275500 -1.058400
      Read-in Center 3098 is at   0.000000  2.328400 -1.058400
      Read-in Center 3099 is at   0.000000  2.381300 -1.058400
      Read-in Center 3100 is at   0.000000  2.434200 -1.058400
      Read-in Center 3101 is at   0.000000  2.487100 -1.058400
      Read-in Center 3102 is at   0.000000  2.540000 -1.058400
      Read-in Center 3103 is at   0.000000  2.593000 -1.058400
      Read-in Center 3104 is at   0.000000 -2.645900 -1.005400
      Read-in Center 3105 is at   0.000000 -2.593000 -1.005400
      Read-in Center 3106 is at   0.000000 -2.540000 -1.005400
      Read-in Center 3107 is at   0.000000 -2.487100 -1.005400
      Read-in Center 3108 is at   0.000000 -2.434200 -1.005400
      Read-in Center 3109 is at   0.000000 -2.381300 -1.005400
      Read-in Center 3110 is at   0.000000 -2.328400 -1.005400
      Read-in Center 3111 is at   0.000000 -2.275500 -1.005400
      Read-in Center 3112 is at   0.000000 -2.222500 -1.005400
      Read-in Center 3113 is at   0.000000 -2.169600 -1.005400
      Read-in Center 3114 is at   0.000000 -2.116700 -1.005400
      Read-in Center 3115 is at   0.000000 -2.063800 -1.005400
      Read-in Center 3116 is at   0.000000 -2.010900 -1.005400
      Read-in Center 3117 is at   0.000000 -1.958000 -1.005400
      Read-in Center 3118 is at   0.000000 -1.905000 -1.005400
      Read-in Center 3119 is at   0.000000 -1.852100 -1.005400
      Read-in Center 3120 is at   0.000000 -1.799200 -1.005400
      Read-in Center 3121 is at   0.000000 -1.746300 -1.005400
      Read-in Center 3122 is at   0.000000 -1.693400 -1.005400
      Read-in Center 3123 is at   0.000000 -1.640400 -1.005400
      Read-in Center 3124 is at   0.000000 -1.587500 -1.005400
      Read-in Center 3125 is at   0.000000 -1.534600 -1.005400
      Read-in Center 3126 is at   0.000000 -1.481700 -1.005400
      Read-in Center 3127 is at   0.000000 -1.428800 -1.005400
      Read-in Center 3128 is at   0.000000 -1.375900 -1.005400
      Read-in Center 3129 is at   0.000000 -1.322900 -1.005400
      Read-in Center 3130 is at   0.000000 -1.270000 -1.005400
      Read-in Center 3131 is at   0.000000 -1.217100 -1.005400
      Read-in Center 3132 is at   0.000000 -1.164200 -1.005400
      Read-in Center 3133 is at   0.000000 -1.111300 -1.005400
      Read-in Center 3134 is at   0.000000 -1.058400 -1.005400
      Read-in Center 3135 is at   0.000000 -1.005400 -1.005400
      Read-in Center 3136 is at   0.000000 -0.952500 -1.005400
      Read-in Center 3137 is at   0.000000 -0.899600 -1.005400
      Read-in Center 3138 is at   0.000000 -0.846700 -1.005400
      Read-in Center 3139 is at   0.000000 -0.793800 -1.005400
      Read-in Center 3140 is at   0.000000 -0.740800 -1.005400
      Read-in Center 3141 is at   0.000000 -0.687900 -1.005400
      Read-in Center 3142 is at   0.000000 -0.635000 -1.005400
      Read-in Center 3143 is at   0.000000 -0.582100 -1.005400
      Read-in Center 3144 is at   0.000000 -0.529200 -1.005400
      Read-in Center 3145 is at   0.000000 -0.476300 -1.005400
      Read-in Center 3146 is at   0.000000 -0.423300 -1.005400
      Read-in Center 3147 is at   0.000000 -0.370400 -1.005400
      Read-in Center 3148 is at   0.000000 -0.317500 -1.005400
      Read-in Center 3149 is at   0.000000 -0.264600 -1.005400
      Read-in Center 3150 is at   0.000000 -0.211700 -1.005400
      Read-in Center 3151 is at   0.000000 -0.158800 -1.005400
      Read-in Center 3152 is at   0.000000 -0.105800 -1.005400
      Read-in Center 3153 is at   0.000000 -0.052900 -1.005400
      Read-in Center 3154 is at   0.000000  0.000000 -1.005400
      Read-in Center 3155 is at   0.000000  0.052900 -1.005400
      Read-in Center 3156 is at   0.000000  0.105800 -1.005400
      Read-in Center 3157 is at   0.000000  0.158800 -1.005400
      Read-in Center 3158 is at   0.000000  0.211700 -1.005400
      Read-in Center 3159 is at   0.000000  0.264600 -1.005400
      Read-in Center 3160 is at   0.000000  0.317500 -1.005400
      Read-in Center 3161 is at   0.000000  0.370400 -1.005400
      Read-in Center 3162 is at   0.000000  0.423300 -1.005400
      Read-in Center 3163 is at   0.000000  0.476300 -1.005400
      Read-in Center 3164 is at   0.000000  0.529200 -1.005400
      Read-in Center 3165 is at   0.000000  0.582100 -1.005400
      Read-in Center 3166 is at   0.000000  0.635000 -1.005400
      Read-in Center 3167 is at   0.000000  0.687900 -1.005400
      Read-in Center 3168 is at   0.000000  0.740800 -1.005400
      Read-in Center 3169 is at   0.000000  0.793800 -1.005400
      Read-in Center 3170 is at   0.000000  0.846700 -1.005400
      Read-in Center 3171 is at   0.000000  0.899600 -1.005400
      Read-in Center 3172 is at   0.000000  0.952500 -1.005400
      Read-in Center 3173 is at   0.000000  1.005400 -1.005400
      Read-in Center 3174 is at   0.000000  1.058400 -1.005400
      Read-in Center 3175 is at   0.000000  1.111300 -1.005400
      Read-in Center 3176 is at   0.000000  1.164200 -1.005400
      Read-in Center 3177 is at   0.000000  1.217100 -1.005400
      Read-in Center 3178 is at   0.000000  1.270000 -1.005400
      Read-in Center 3179 is at   0.000000  1.322900 -1.005400
      Read-in Center 3180 is at   0.000000  1.375900 -1.005400
      Read-in Center 3181 is at   0.000000  1.428800 -1.005400
      Read-in Center 3182 is at   0.000000  1.481700 -1.005400
      Read-in Center 3183 is at   0.000000  1.534600 -1.005400
      Read-in Center 3184 is at   0.000000  1.587500 -1.005400
      Read-in Center 3185 is at   0.000000  1.640400 -1.005400
      Read-in Center 3186 is at   0.000000  1.693400 -1.005400
      Read-in Center 3187 is at   0.000000  1.746300 -1.005400
      Read-in Center 3188 is at   0.000000  1.799200 -1.005400
      Read-in Center 3189 is at   0.000000  1.852100 -1.005400
      Read-in Center 3190 is at   0.000000  1.905000 -1.005400
      Read-in Center 3191 is at   0.000000  1.958000 -1.005400
      Read-in Center 3192 is at   0.000000  2.010900 -1.005400
      Read-in Center 3193 is at   0.000000  2.063800 -1.005400
      Read-in Center 3194 is at   0.000000  2.116700 -1.005400
      Read-in Center 3195 is at   0.000000  2.169600 -1.005400
      Read-in Center 3196 is at   0.000000  2.222500 -1.005400
      Read-in Center 3197 is at   0.000000  2.275500 -1.005400
      Read-in Center 3198 is at   0.000000  2.328400 -1.005400
      Read-in Center 3199 is at   0.000000  2.381300 -1.005400
      Read-in Center 3200 is at   0.000000  2.434200 -1.005400
      Read-in Center 3201 is at   0.000000  2.487100 -1.005400
      Read-in Center 3202 is at   0.000000  2.540000 -1.005400
      Read-in Center 3203 is at   0.000000  2.593000 -1.005400
      Read-in Center 3204 is at   0.000000 -2.645900 -0.952500
      Read-in Center 3205 is at   0.000000 -2.593000 -0.952500
      Read-in Center 3206 is at   0.000000 -2.540000 -0.952500
      Read-in Center 3207 is at   0.000000 -2.487100 -0.952500
      Read-in Center 3208 is at   0.000000 -2.434200 -0.952500
      Read-in Center 3209 is at   0.000000 -2.381300 -0.952500
      Read-in Center 3210 is at   0.000000 -2.328400 -0.952500
      Read-in Center 3211 is at   0.000000 -2.275500 -0.952500
      Read-in Center 3212 is at   0.000000 -2.222500 -0.952500
      Read-in Center 3213 is at   0.000000 -2.169600 -0.952500
      Read-in Center 3214 is at   0.000000 -2.116700 -0.952500
      Read-in Center 3215 is at   0.000000 -2.063800 -0.952500
      Read-in Center 3216 is at   0.000000 -2.010900 -0.952500
      Read-in Center 3217 is at   0.000000 -1.958000 -0.952500
      Read-in Center 3218 is at   0.000000 -1.905000 -0.952500
      Read-in Center 3219 is at   0.000000 -1.852100 -0.952500
      Read-in Center 3220 is at   0.000000 -1.799200 -0.952500
      Read-in Center 3221 is at   0.000000 -1.746300 -0.952500
      Read-in Center 3222 is at   0.000000 -1.693400 -0.952500
      Read-in Center 3223 is at   0.000000 -1.640400 -0.952500
      Read-in Center 3224 is at   0.000000 -1.587500 -0.952500
      Read-in Center 3225 is at   0.000000 -1.534600 -0.952500
      Read-in Center 3226 is at   0.000000 -1.481700 -0.952500
      Read-in Center 3227 is at   0.000000 -1.428800 -0.952500
      Read-in Center 3228 is at   0.000000 -1.375900 -0.952500
      Read-in Center 3229 is at   0.000000 -1.322900 -0.952500
      Read-in Center 3230 is at   0.000000 -1.270000 -0.952500
      Read-in Center 3231 is at   0.000000 -1.217100 -0.952500
      Read-in Center 3232 is at   0.000000 -1.164200 -0.952500
      Read-in Center 3233 is at   0.000000 -1.111300 -0.952500
      Read-in Center 3234 is at   0.000000 -1.058400 -0.952500
      Read-in Center 3235 is at   0.000000 -1.005400 -0.952500
      Read-in Center 3236 is at   0.000000 -0.952500 -0.952500
      Read-in Center 3237 is at   0.000000 -0.899600 -0.952500
      Read-in Center 3238 is at   0.000000 -0.846700 -0.952500
      Read-in Center 3239 is at   0.000000 -0.793800 -0.952500
      Read-in Center 3240 is at   0.000000 -0.740800 -0.952500
      Read-in Center 3241 is at   0.000000 -0.687900 -0.952500
      Read-in Center 3242 is at   0.000000 -0.635000 -0.952500
      Read-in Center 3243 is at   0.000000 -0.582100 -0.952500
      Read-in Center 3244 is at   0.000000 -0.529200 -0.952500
      Read-in Center 3245 is at   0.000000 -0.476300 -0.952500
      Read-in Center 3246 is at   0.000000 -0.423300 -0.952500
      Read-in Center 3247 is at   0.000000 -0.370400 -0.952500
      Read-in Center 3248 is at   0.000000 -0.317500 -0.952500
      Read-in Center 3249 is at   0.000000 -0.264600 -0.952500
      Read-in Center 3250 is at   0.000000 -0.211700 -0.952500
      Read-in Center 3251 is at   0.000000 -0.158800 -0.952500
      Read-in Center 3252 is at   0.000000 -0.105800 -0.952500
      Read-in Center 3253 is at   0.000000 -0.052900 -0.952500
      Read-in Center 3254 is at   0.000000  0.000000 -0.952500
      Read-in Center 3255 is at   0.000000  0.052900 -0.952500
      Read-in Center 3256 is at   0.000000  0.105800 -0.952500
      Read-in Center 3257 is at   0.000000  0.158800 -0.952500
      Read-in Center 3258 is at   0.000000  0.211700 -0.952500
      Read-in Center 3259 is at   0.000000  0.264600 -0.952500
      Read-in Center 3260 is at   0.000000  0.317500 -0.952500
      Read-in Center 3261 is at   0.000000  0.370400 -0.952500
      Read-in Center 3262 is at   0.000000  0.423300 -0.952500
      Read-in Center 3263 is at   0.000000  0.476300 -0.952500
      Read-in Center 3264 is at   0.000000  0.529200 -0.952500
      Read-in Center 3265 is at   0.000000  0.582100 -0.952500
      Read-in Center 3266 is at   0.000000  0.635000 -0.952500
      Read-in Center 3267 is at   0.000000  0.687900 -0.952500
      Read-in Center 3268 is at   0.000000  0.740800 -0.952500
      Read-in Center 3269 is at   0.000000  0.793800 -0.952500
      Read-in Center 3270 is at   0.000000  0.846700 -0.952500
      Read-in Center 3271 is at   0.000000  0.899600 -0.952500
      Read-in Center 3272 is at   0.000000  0.952500 -0.952500
      Read-in Center 3273 is at   0.000000  1.005400 -0.952500
      Read-in Center 3274 is at   0.000000  1.058400 -0.952500
      Read-in Center 3275 is at   0.000000  1.111300 -0.952500
      Read-in Center 3276 is at   0.000000  1.164200 -0.952500
      Read-in Center 3277 is at   0.000000  1.217100 -0.952500
      Read-in Center 3278 is at   0.000000  1.270000 -0.952500
      Read-in Center 3279 is at   0.000000  1.322900 -0.952500
      Read-in Center 3280 is at   0.000000  1.375900 -0.952500
      Read-in Center 3281 is at   0.000000  1.428800 -0.952500
      Read-in Center 3282 is at   0.000000  1.481700 -0.952500
      Read-in Center 3283 is at   0.000000  1.534600 -0.952500
      Read-in Center 3284 is at   0.000000  1.587500 -0.952500
      Read-in Center 3285 is at   0.000000  1.640400 -0.952500
      Read-in Center 3286 is at   0.000000  1.693400 -0.952500
      Read-in Center 3287 is at   0.000000  1.746300 -0.952500
      Read-in Center 3288 is at   0.000000  1.799200 -0.952500
      Read-in Center 3289 is at   0.000000  1.852100 -0.952500
      Read-in Center 3290 is at   0.000000  1.905000 -0.952500
      Read-in Center 3291 is at   0.000000  1.958000 -0.952500
      Read-in Center 3292 is at   0.000000  2.010900 -0.952500
      Read-in Center 3293 is at   0.000000  2.063800 -0.952500
      Read-in Center 3294 is at   0.000000  2.116700 -0.952500
      Read-in Center 3295 is at   0.000000  2.169600 -0.952500
      Read-in Center 3296 is at   0.000000  2.222500 -0.952500
      Read-in Center 3297 is at   0.000000  2.275500 -0.952500
      Read-in Center 3298 is at   0.000000  2.328400 -0.952500
      Read-in Center 3299 is at   0.000000  2.381300 -0.952500
      Read-in Center 3300 is at   0.000000  2.434200 -0.952500
      Read-in Center 3301 is at   0.000000  2.487100 -0.952500
      Read-in Center 3302 is at   0.000000  2.540000 -0.952500
      Read-in Center 3303 is at   0.000000  2.593000 -0.952500
      Read-in Center 3304 is at   0.000000 -2.645900 -0.899600
      Read-in Center 3305 is at   0.000000 -2.593000 -0.899600
      Read-in Center 3306 is at   0.000000 -2.540000 -0.899600
      Read-in Center 3307 is at   0.000000 -2.487100 -0.899600
      Read-in Center 3308 is at   0.000000 -2.434200 -0.899600
      Read-in Center 3309 is at   0.000000 -2.381300 -0.899600
      Read-in Center 3310 is at   0.000000 -2.328400 -0.899600
      Read-in Center 3311 is at   0.000000 -2.275500 -0.899600
      Read-in Center 3312 is at   0.000000 -2.222500 -0.899600
      Read-in Center 3313 is at   0.000000 -2.169600 -0.899600
      Read-in Center 3314 is at   0.000000 -2.116700 -0.899600
      Read-in Center 3315 is at   0.000000 -2.063800 -0.899600
      Read-in Center 3316 is at   0.000000 -2.010900 -0.899600
      Read-in Center 3317 is at   0.000000 -1.958000 -0.899600
      Read-in Center 3318 is at   0.000000 -1.905000 -0.899600
      Read-in Center 3319 is at   0.000000 -1.852100 -0.899600
      Read-in Center 3320 is at   0.000000 -1.799200 -0.899600
      Read-in Center 3321 is at   0.000000 -1.746300 -0.899600
      Read-in Center 3322 is at   0.000000 -1.693400 -0.899600
      Read-in Center 3323 is at   0.000000 -1.640400 -0.899600
      Read-in Center 3324 is at   0.000000 -1.587500 -0.899600
      Read-in Center 3325 is at   0.000000 -1.534600 -0.899600
      Read-in Center 3326 is at   0.000000 -1.481700 -0.899600
      Read-in Center 3327 is at   0.000000 -1.428800 -0.899600
      Read-in Center 3328 is at   0.000000 -1.375900 -0.899600
      Read-in Center 3329 is at   0.000000 -1.322900 -0.899600
      Read-in Center 3330 is at   0.000000 -1.270000 -0.899600
      Read-in Center 3331 is at   0.000000 -1.217100 -0.899600
      Read-in Center 3332 is at   0.000000 -1.164200 -0.899600
      Read-in Center 3333 is at   0.000000 -1.111300 -0.899600
      Read-in Center 3334 is at   0.000000 -1.058400 -0.899600
      Read-in Center 3335 is at   0.000000 -1.005400 -0.899600
      Read-in Center 3336 is at   0.000000 -0.952500 -0.899600
      Read-in Center 3337 is at   0.000000 -0.899600 -0.899600
      Read-in Center 3338 is at   0.000000 -0.846700 -0.899600
      Read-in Center 3339 is at   0.000000 -0.793800 -0.899600
      Read-in Center 3340 is at   0.000000 -0.740800 -0.899600
      Read-in Center 3341 is at   0.000000 -0.687900 -0.899600
      Read-in Center 3342 is at   0.000000 -0.635000 -0.899600
      Read-in Center 3343 is at   0.000000 -0.582100 -0.899600
      Read-in Center 3344 is at   0.000000 -0.529200 -0.899600
      Read-in Center 3345 is at   0.000000 -0.476300 -0.899600
      Read-in Center 3346 is at   0.000000 -0.423300 -0.899600
      Read-in Center 3347 is at   0.000000 -0.370400 -0.899600
      Read-in Center 3348 is at   0.000000 -0.317500 -0.899600
      Read-in Center 3349 is at   0.000000 -0.264600 -0.899600
      Read-in Center 3350 is at   0.000000 -0.211700 -0.899600
      Read-in Center 3351 is at   0.000000 -0.158800 -0.899600
      Read-in Center 3352 is at   0.000000 -0.105800 -0.899600
      Read-in Center 3353 is at   0.000000 -0.052900 -0.899600
      Read-in Center 3354 is at   0.000000  0.000000 -0.899600
      Read-in Center 3355 is at   0.000000  0.052900 -0.899600
      Read-in Center 3356 is at   0.000000  0.105800 -0.899600
      Read-in Center 3357 is at   0.000000  0.158800 -0.899600
      Read-in Center 3358 is at   0.000000  0.211700 -0.899600
      Read-in Center 3359 is at   0.000000  0.264600 -0.899600
      Read-in Center 3360 is at   0.000000  0.317500 -0.899600
      Read-in Center 3361 is at   0.000000  0.370400 -0.899600
      Read-in Center 3362 is at   0.000000  0.423300 -0.899600
      Read-in Center 3363 is at   0.000000  0.476300 -0.899600
      Read-in Center 3364 is at   0.000000  0.529200 -0.899600
      Read-in Center 3365 is at   0.000000  0.582100 -0.899600
      Read-in Center 3366 is at   0.000000  0.635000 -0.899600
      Read-in Center 3367 is at   0.000000  0.687900 -0.899600
      Read-in Center 3368 is at   0.000000  0.740800 -0.899600
      Read-in Center 3369 is at   0.000000  0.793800 -0.899600
      Read-in Center 3370 is at   0.000000  0.846700 -0.899600
      Read-in Center 3371 is at   0.000000  0.899600 -0.899600
      Read-in Center 3372 is at   0.000000  0.952500 -0.899600
      Read-in Center 3373 is at   0.000000  1.005400 -0.899600
      Read-in Center 3374 is at   0.000000  1.058400 -0.899600
      Read-in Center 3375 is at   0.000000  1.111300 -0.899600
      Read-in Center 3376 is at   0.000000  1.164200 -0.899600
      Read-in Center 3377 is at   0.000000  1.217100 -0.899600
      Read-in Center 3378 is at   0.000000  1.270000 -0.899600
      Read-in Center 3379 is at   0.000000  1.322900 -0.899600
      Read-in Center 3380 is at   0.000000  1.375900 -0.899600
      Read-in Center 3381 is at   0.000000  1.428800 -0.899600
      Read-in Center 3382 is at   0.000000  1.481700 -0.899600
      Read-in Center 3383 is at   0.000000  1.534600 -0.899600
      Read-in Center 3384 is at   0.000000  1.587500 -0.899600
      Read-in Center 3385 is at   0.000000  1.640400 -0.899600
      Read-in Center 3386 is at   0.000000  1.693400 -0.899600
      Read-in Center 3387 is at   0.000000  1.746300 -0.899600
      Read-in Center 3388 is at   0.000000  1.799200 -0.899600
      Read-in Center 3389 is at   0.000000  1.852100 -0.899600
      Read-in Center 3390 is at   0.000000  1.905000 -0.899600
      Read-in Center 3391 is at   0.000000  1.958000 -0.899600
      Read-in Center 3392 is at   0.000000  2.010900 -0.899600
      Read-in Center 3393 is at   0.000000  2.063800 -0.899600
      Read-in Center 3394 is at   0.000000  2.116700 -0.899600
      Read-in Center 3395 is at   0.000000  2.169600 -0.899600
      Read-in Center 3396 is at   0.000000  2.222500 -0.899600
      Read-in Center 3397 is at   0.000000  2.275500 -0.899600
      Read-in Center 3398 is at   0.000000  2.328400 -0.899600
      Read-in Center 3399 is at   0.000000  2.381300 -0.899600
      Read-in Center 3400 is at   0.000000  2.434200 -0.899600
      Read-in Center 3401 is at   0.000000  2.487100 -0.899600
      Read-in Center 3402 is at   0.000000  2.540000 -0.899600
      Read-in Center 3403 is at   0.000000  2.593000 -0.899600
      Read-in Center 3404 is at   0.000000 -2.645900 -0.846700
      Read-in Center 3405 is at   0.000000 -2.593000 -0.846700
      Read-in Center 3406 is at   0.000000 -2.540000 -0.846700
      Read-in Center 3407 is at   0.000000 -2.487100 -0.846700
      Read-in Center 3408 is at   0.000000 -2.434200 -0.846700
      Read-in Center 3409 is at   0.000000 -2.381300 -0.846700
      Read-in Center 3410 is at   0.000000 -2.328400 -0.846700
      Read-in Center 3411 is at   0.000000 -2.275500 -0.846700
      Read-in Center 3412 is at   0.000000 -2.222500 -0.846700
      Read-in Center 3413 is at   0.000000 -2.169600 -0.846700
      Read-in Center 3414 is at   0.000000 -2.116700 -0.846700
      Read-in Center 3415 is at   0.000000 -2.063800 -0.846700
      Read-in Center 3416 is at   0.000000 -2.010900 -0.846700
      Read-in Center 3417 is at   0.000000 -1.958000 -0.846700
      Read-in Center 3418 is at   0.000000 -1.905000 -0.846700
      Read-in Center 3419 is at   0.000000 -1.852100 -0.846700
      Read-in Center 3420 is at   0.000000 -1.799200 -0.846700
      Read-in Center 3421 is at   0.000000 -1.746300 -0.846700
      Read-in Center 3422 is at   0.000000 -1.693400 -0.846700
      Read-in Center 3423 is at   0.000000 -1.640400 -0.846700
      Read-in Center 3424 is at   0.000000 -1.587500 -0.846700
      Read-in Center 3425 is at   0.000000 -1.534600 -0.846700
      Read-in Center 3426 is at   0.000000 -1.481700 -0.846700
      Read-in Center 3427 is at   0.000000 -1.428800 -0.846700
      Read-in Center 3428 is at   0.000000 -1.375900 -0.846700
      Read-in Center 3429 is at   0.000000 -1.322900 -0.846700
      Read-in Center 3430 is at   0.000000 -1.270000 -0.846700
      Read-in Center 3431 is at   0.000000 -1.217100 -0.846700
      Read-in Center 3432 is at   0.000000 -1.164200 -0.846700
      Read-in Center 3433 is at   0.000000 -1.111300 -0.846700
      Read-in Center 3434 is at   0.000000 -1.058400 -0.846700
      Read-in Center 3435 is at   0.000000 -1.005400 -0.846700
      Read-in Center 3436 is at   0.000000 -0.952500 -0.846700
      Read-in Center 3437 is at   0.000000 -0.899600 -0.846700
      Read-in Center 3438 is at   0.000000 -0.846700 -0.846700
      Read-in Center 3439 is at   0.000000 -0.793800 -0.846700
      Read-in Center 3440 is at   0.000000 -0.740800 -0.846700
      Read-in Center 3441 is at   0.000000 -0.687900 -0.846700
      Read-in Center 3442 is at   0.000000 -0.635000 -0.846700
      Read-in Center 3443 is at   0.000000 -0.582100 -0.846700
      Read-in Center 3444 is at   0.000000 -0.529200 -0.846700
      Read-in Center 3445 is at   0.000000 -0.476300 -0.846700
      Read-in Center 3446 is at   0.000000 -0.423300 -0.846700
      Read-in Center 3447 is at   0.000000 -0.370400 -0.846700
      Read-in Center 3448 is at   0.000000 -0.317500 -0.846700
      Read-in Center 3449 is at   0.000000 -0.264600 -0.846700
      Read-in Center 3450 is at   0.000000 -0.211700 -0.846700
      Read-in Center 3451 is at   0.000000 -0.158800 -0.846700
      Read-in Center 3452 is at   0.000000 -0.105800 -0.846700
      Read-in Center 3453 is at   0.000000 -0.052900 -0.846700
      Read-in Center 3454 is at   0.000000  0.000000 -0.846700
      Read-in Center 3455 is at   0.000000  0.052900 -0.846700
      Read-in Center 3456 is at   0.000000  0.105800 -0.846700
      Read-in Center 3457 is at   0.000000  0.158800 -0.846700
      Read-in Center 3458 is at   0.000000  0.211700 -0.846700
      Read-in Center 3459 is at   0.000000  0.264600 -0.846700
      Read-in Center 3460 is at   0.000000  0.317500 -0.846700
      Read-in Center 3461 is at   0.000000  0.370400 -0.846700
      Read-in Center 3462 is at   0.000000  0.423300 -0.846700
      Read-in Center 3463 is at   0.000000  0.476300 -0.846700
      Read-in Center 3464 is at   0.000000  0.529200 -0.846700
      Read-in Center 3465 is at   0.000000  0.582100 -0.846700
      Read-in Center 3466 is at   0.000000  0.635000 -0.846700
      Read-in Center 3467 is at   0.000000  0.687900 -0.846700
      Read-in Center 3468 is at   0.000000  0.740800 -0.846700
      Read-in Center 3469 is at   0.000000  0.793800 -0.846700
      Read-in Center 3470 is at   0.000000  0.846700 -0.846700
      Read-in Center 3471 is at   0.000000  0.899600 -0.846700
      Read-in Center 3472 is at   0.000000  0.952500 -0.846700
      Read-in Center 3473 is at   0.000000  1.005400 -0.846700
      Read-in Center 3474 is at   0.000000  1.058400 -0.846700
      Read-in Center 3475 is at   0.000000  1.111300 -0.846700
      Read-in Center 3476 is at   0.000000  1.164200 -0.846700
      Read-in Center 3477 is at   0.000000  1.217100 -0.846700
      Read-in Center 3478 is at   0.000000  1.270000 -0.846700
      Read-in Center 3479 is at   0.000000  1.322900 -0.846700
      Read-in Center 3480 is at   0.000000  1.375900 -0.846700
      Read-in Center 3481 is at   0.000000  1.428800 -0.846700
      Read-in Center 3482 is at   0.000000  1.481700 -0.846700
      Read-in Center 3483 is at   0.000000  1.534600 -0.846700
      Read-in Center 3484 is at   0.000000  1.587500 -0.846700
      Read-in Center 3485 is at   0.000000  1.640400 -0.846700
      Read-in Center 3486 is at   0.000000  1.693400 -0.846700
      Read-in Center 3487 is at   0.000000  1.746300 -0.846700
      Read-in Center 3488 is at   0.000000  1.799200 -0.846700
      Read-in Center 3489 is at   0.000000  1.852100 -0.846700
      Read-in Center 3490 is at   0.000000  1.905000 -0.846700
      Read-in Center 3491 is at   0.000000  1.958000 -0.846700
      Read-in Center 3492 is at   0.000000  2.010900 -0.846700
      Read-in Center 3493 is at   0.000000  2.063800 -0.846700
      Read-in Center 3494 is at   0.000000  2.116700 -0.846700
      Read-in Center 3495 is at   0.000000  2.169600 -0.846700
      Read-in Center 3496 is at   0.000000  2.222500 -0.846700
      Read-in Center 3497 is at   0.000000  2.275500 -0.846700
      Read-in Center 3498 is at   0.000000  2.328400 -0.846700
      Read-in Center 3499 is at   0.000000  2.381300 -0.846700
      Read-in Center 3500 is at   0.000000  2.434200 -0.846700
      Read-in Center 3501 is at   0.000000  2.487100 -0.846700
      Read-in Center 3502 is at   0.000000  2.540000 -0.846700
      Read-in Center 3503 is at   0.000000  2.593000 -0.846700
      Read-in Center 3504 is at   0.000000 -2.645900 -0.793800
      Read-in Center 3505 is at   0.000000 -2.593000 -0.793800
      Read-in Center 3506 is at   0.000000 -2.540000 -0.793800
      Read-in Center 3507 is at   0.000000 -2.487100 -0.793800
      Read-in Center 3508 is at   0.000000 -2.434200 -0.793800
      Read-in Center 3509 is at   0.000000 -2.381300 -0.793800
      Read-in Center 3510 is at   0.000000 -2.328400 -0.793800
      Read-in Center 3511 is at   0.000000 -2.275500 -0.793800
      Read-in Center 3512 is at   0.000000 -2.222500 -0.793800
      Read-in Center 3513 is at   0.000000 -2.169600 -0.793800
      Read-in Center 3514 is at   0.000000 -2.116700 -0.793800
      Read-in Center 3515 is at   0.000000 -2.063800 -0.793800
      Read-in Center 3516 is at   0.000000 -2.010900 -0.793800
      Read-in Center 3517 is at   0.000000 -1.958000 -0.793800
      Read-in Center 3518 is at   0.000000 -1.905000 -0.793800
      Read-in Center 3519 is at   0.000000 -1.852100 -0.793800
      Read-in Center 3520 is at   0.000000 -1.799200 -0.793800
      Read-in Center 3521 is at   0.000000 -1.746300 -0.793800
      Read-in Center 3522 is at   0.000000 -1.693400 -0.793800
      Read-in Center 3523 is at   0.000000 -1.640400 -0.793800
      Read-in Center 3524 is at   0.000000 -1.587500 -0.793800
      Read-in Center 3525 is at   0.000000 -1.534600 -0.793800
      Read-in Center 3526 is at   0.000000 -1.481700 -0.793800
      Read-in Center 3527 is at   0.000000 -1.428800 -0.793800
      Read-in Center 3528 is at   0.000000 -1.375900 -0.793800
      Read-in Center 3529 is at   0.000000 -1.322900 -0.793800
      Read-in Center 3530 is at   0.000000 -1.270000 -0.793800
      Read-in Center 3531 is at   0.000000 -1.217100 -0.793800
      Read-in Center 3532 is at   0.000000 -1.164200 -0.793800
      Read-in Center 3533 is at   0.000000 -1.111300 -0.793800
      Read-in Center 3534 is at   0.000000 -1.058400 -0.793800
      Read-in Center 3535 is at   0.000000 -1.005400 -0.793800
      Read-in Center 3536 is at   0.000000 -0.952500 -0.793800
      Read-in Center 3537 is at   0.000000 -0.899600 -0.793800
      Read-in Center 3538 is at   0.000000 -0.846700 -0.793800
      Read-in Center 3539 is at   0.000000 -0.793800 -0.793800
      Read-in Center 3540 is at   0.000000 -0.740800 -0.793800
      Read-in Center 3541 is at   0.000000 -0.687900 -0.793800
      Read-in Center 3542 is at   0.000000 -0.635000 -0.793800
      Read-in Center 3543 is at   0.000000 -0.582100 -0.793800
      Read-in Center 3544 is at   0.000000 -0.529200 -0.793800
      Read-in Center 3545 is at   0.000000 -0.476300 -0.793800
      Read-in Center 3546 is at   0.000000 -0.423300 -0.793800
      Read-in Center 3547 is at   0.000000 -0.370400 -0.793800
      Read-in Center 3548 is at   0.000000 -0.317500 -0.793800
      Read-in Center 3549 is at   0.000000 -0.264600 -0.793800
      Read-in Center 3550 is at   0.000000 -0.211700 -0.793800
      Read-in Center 3551 is at   0.000000 -0.158800 -0.793800
      Read-in Center 3552 is at   0.000000 -0.105800 -0.793800
      Read-in Center 3553 is at   0.000000 -0.052900 -0.793800
      Read-in Center 3554 is at   0.000000  0.000000 -0.793800
      Read-in Center 3555 is at   0.000000  0.052900 -0.793800
      Read-in Center 3556 is at   0.000000  0.105800 -0.793800
      Read-in Center 3557 is at   0.000000  0.158800 -0.793800
      Read-in Center 3558 is at   0.000000  0.211700 -0.793800
      Read-in Center 3559 is at   0.000000  0.264600 -0.793800
      Read-in Center 3560 is at   0.000000  0.317500 -0.793800
      Read-in Center 3561 is at   0.000000  0.370400 -0.793800
      Read-in Center 3562 is at   0.000000  0.423300 -0.793800
      Read-in Center 3563 is at   0.000000  0.476300 -0.793800
      Read-in Center 3564 is at   0.000000  0.529200 -0.793800
      Read-in Center 3565 is at   0.000000  0.582100 -0.793800
      Read-in Center 3566 is at   0.000000  0.635000 -0.793800
      Read-in Center 3567 is at   0.000000  0.687900 -0.793800
      Read-in Center 3568 is at   0.000000  0.740800 -0.793800
      Read-in Center 3569 is at   0.000000  0.793800 -0.793800
      Read-in Center 3570 is at   0.000000  0.846700 -0.793800
      Read-in Center 3571 is at   0.000000  0.899600 -0.793800
      Read-in Center 3572 is at   0.000000  0.952500 -0.793800
      Read-in Center 3573 is at   0.000000  1.005400 -0.793800
      Read-in Center 3574 is at   0.000000  1.058400 -0.793800
      Read-in Center 3575 is at   0.000000  1.111300 -0.793800
      Read-in Center 3576 is at   0.000000  1.164200 -0.793800
      Read-in Center 3577 is at   0.000000  1.217100 -0.793800
      Read-in Center 3578 is at   0.000000  1.270000 -0.793800
      Read-in Center 3579 is at   0.000000  1.322900 -0.793800
      Read-in Center 3580 is at   0.000000  1.375900 -0.793800
      Read-in Center 3581 is at   0.000000  1.428800 -0.793800
      Read-in Center 3582 is at   0.000000  1.481700 -0.793800
      Read-in Center 3583 is at   0.000000  1.534600 -0.793800
      Read-in Center 3584 is at   0.000000  1.587500 -0.793800
      Read-in Center 3585 is at   0.000000  1.640400 -0.793800
      Read-in Center 3586 is at   0.000000  1.693400 -0.793800
      Read-in Center 3587 is at   0.000000  1.746300 -0.793800
      Read-in Center 3588 is at   0.000000  1.799200 -0.793800
      Read-in Center 3589 is at   0.000000  1.852100 -0.793800
      Read-in Center 3590 is at   0.000000  1.905000 -0.793800
      Read-in Center 3591 is at   0.000000  1.958000 -0.793800
      Read-in Center 3592 is at   0.000000  2.010900 -0.793800
      Read-in Center 3593 is at   0.000000  2.063800 -0.793800
      Read-in Center 3594 is at   0.000000  2.116700 -0.793800
      Read-in Center 3595 is at   0.000000  2.169600 -0.793800
      Read-in Center 3596 is at   0.000000  2.222500 -0.793800
      Read-in Center 3597 is at   0.000000  2.275500 -0.793800
      Read-in Center 3598 is at   0.000000  2.328400 -0.793800
      Read-in Center 3599 is at   0.000000  2.381300 -0.793800
      Read-in Center 3600 is at   0.000000  2.434200 -0.793800
      Read-in Center 3601 is at   0.000000  2.487100 -0.793800
      Read-in Center 3602 is at   0.000000  2.540000 -0.793800
      Read-in Center 3603 is at   0.000000  2.593000 -0.793800
      Read-in Center 3604 is at   0.000000 -2.645900 -0.740800
      Read-in Center 3605 is at   0.000000 -2.593000 -0.740800
      Read-in Center 3606 is at   0.000000 -2.540000 -0.740800
      Read-in Center 3607 is at   0.000000 -2.487100 -0.740800
      Read-in Center 3608 is at   0.000000 -2.434200 -0.740800
      Read-in Center 3609 is at   0.000000 -2.381300 -0.740800
      Read-in Center 3610 is at   0.000000 -2.328400 -0.740800
      Read-in Center 3611 is at   0.000000 -2.275500 -0.740800
      Read-in Center 3612 is at   0.000000 -2.222500 -0.740800
      Read-in Center 3613 is at   0.000000 -2.169600 -0.740800
      Read-in Center 3614 is at   0.000000 -2.116700 -0.740800
      Read-in Center 3615 is at   0.000000 -2.063800 -0.740800
      Read-in Center 3616 is at   0.000000 -2.010900 -0.740800
      Read-in Center 3617 is at   0.000000 -1.958000 -0.740800
      Read-in Center 3618 is at   0.000000 -1.905000 -0.740800
      Read-in Center 3619 is at   0.000000 -1.852100 -0.740800
      Read-in Center 3620 is at   0.000000 -1.799200 -0.740800
      Read-in Center 3621 is at   0.000000 -1.746300 -0.740800
      Read-in Center 3622 is at   0.000000 -1.693400 -0.740800
      Read-in Center 3623 is at   0.000000 -1.640400 -0.740800
      Read-in Center 3624 is at   0.000000 -1.587500 -0.740800
      Read-in Center 3625 is at   0.000000 -1.534600 -0.740800
      Read-in Center 3626 is at   0.000000 -1.481700 -0.740800
      Read-in Center 3627 is at   0.000000 -1.428800 -0.740800
      Read-in Center 3628 is at   0.000000 -1.375900 -0.740800
      Read-in Center 3629 is at   0.000000 -1.322900 -0.740800
      Read-in Center 3630 is at   0.000000 -1.270000 -0.740800
      Read-in Center 3631 is at   0.000000 -1.217100 -0.740800
      Read-in Center 3632 is at   0.000000 -1.164200 -0.740800
      Read-in Center 3633 is at   0.000000 -1.111300 -0.740800
      Read-in Center 3634 is at   0.000000 -1.058400 -0.740800
      Read-in Center 3635 is at   0.000000 -1.005400 -0.740800
      Read-in Center 3636 is at   0.000000 -0.952500 -0.740800
      Read-in Center 3637 is at   0.000000 -0.899600 -0.740800
      Read-in Center 3638 is at   0.000000 -0.846700 -0.740800
      Read-in Center 3639 is at   0.000000 -0.793800 -0.740800
      Read-in Center 3640 is at   0.000000 -0.740800 -0.740800
      Read-in Center 3641 is at   0.000000 -0.687900 -0.740800
      Read-in Center 3642 is at   0.000000 -0.635000 -0.740800
      Read-in Center 3643 is at   0.000000 -0.582100 -0.740800
      Read-in Center 3644 is at   0.000000 -0.529200 -0.740800
      Read-in Center 3645 is at   0.000000 -0.476300 -0.740800
      Read-in Center 3646 is at   0.000000 -0.423300 -0.740800
      Read-in Center 3647 is at   0.000000 -0.370400 -0.740800
      Read-in Center 3648 is at   0.000000 -0.317500 -0.740800
      Read-in Center 3649 is at   0.000000 -0.264600 -0.740800
      Read-in Center 3650 is at   0.000000 -0.211700 -0.740800
      Read-in Center 3651 is at   0.000000 -0.158800 -0.740800
      Read-in Center 3652 is at   0.000000 -0.105800 -0.740800
      Read-in Center 3653 is at   0.000000 -0.052900 -0.740800
      Read-in Center 3654 is at   0.000000  0.000000 -0.740800
      Read-in Center 3655 is at   0.000000  0.052900 -0.740800
      Read-in Center 3656 is at   0.000000  0.105800 -0.740800
      Read-in Center 3657 is at   0.000000  0.158800 -0.740800
      Read-in Center 3658 is at   0.000000  0.211700 -0.740800
      Read-in Center 3659 is at   0.000000  0.264600 -0.740800
      Read-in Center 3660 is at   0.000000  0.317500 -0.740800
      Read-in Center 3661 is at   0.000000  0.370400 -0.740800
      Read-in Center 3662 is at   0.000000  0.423300 -0.740800
      Read-in Center 3663 is at   0.000000  0.476300 -0.740800
      Read-in Center 3664 is at   0.000000  0.529200 -0.740800
      Read-in Center 3665 is at   0.000000  0.582100 -0.740800
      Read-in Center 3666 is at   0.000000  0.635000 -0.740800
      Read-in Center 3667 is at   0.000000  0.687900 -0.740800
      Read-in Center 3668 is at   0.000000  0.740800 -0.740800
      Read-in Center 3669 is at   0.000000  0.793800 -0.740800
      Read-in Center 3670 is at   0.000000  0.846700 -0.740800
      Read-in Center 3671 is at   0.000000  0.899600 -0.740800
      Read-in Center 3672 is at   0.000000  0.952500 -0.740800
      Read-in Center 3673 is at   0.000000  1.005400 -0.740800
      Read-in Center 3674 is at   0.000000  1.058400 -0.740800
      Read-in Center 3675 is at   0.000000  1.111300 -0.740800
      Read-in Center 3676 is at   0.000000  1.164200 -0.740800
      Read-in Center 3677 is at   0.000000  1.217100 -0.740800
      Read-in Center 3678 is at   0.000000  1.270000 -0.740800
      Read-in Center 3679 is at   0.000000  1.322900 -0.740800
      Read-in Center 3680 is at   0.000000  1.375900 -0.740800
      Read-in Center 3681 is at   0.000000  1.428800 -0.740800
      Read-in Center 3682 is at   0.000000  1.481700 -0.740800
      Read-in Center 3683 is at   0.000000  1.534600 -0.740800
      Read-in Center 3684 is at   0.000000  1.587500 -0.740800
      Read-in Center 3685 is at   0.000000  1.640400 -0.740800
      Read-in Center 3686 is at   0.000000  1.693400 -0.740800
      Read-in Center 3687 is at   0.000000  1.746300 -0.740800
      Read-in Center 3688 is at   0.000000  1.799200 -0.740800
      Read-in Center 3689 is at   0.000000  1.852100 -0.740800
      Read-in Center 3690 is at   0.000000  1.905000 -0.740800
      Read-in Center 3691 is at   0.000000  1.958000 -0.740800
      Read-in Center 3692 is at   0.000000  2.010900 -0.740800
      Read-in Center 3693 is at   0.000000  2.063800 -0.740800
      Read-in Center 3694 is at   0.000000  2.116700 -0.740800
      Read-in Center 3695 is at   0.000000  2.169600 -0.740800
      Read-in Center 3696 is at   0.000000  2.222500 -0.740800
      Read-in Center 3697 is at   0.000000  2.275500 -0.740800
      Read-in Center 3698 is at   0.000000  2.328400 -0.740800
      Read-in Center 3699 is at   0.000000  2.381300 -0.740800
      Read-in Center 3700 is at   0.000000  2.434200 -0.740800
      Read-in Center 3701 is at   0.000000  2.487100 -0.740800
      Read-in Center 3702 is at   0.000000  2.540000 -0.740800
      Read-in Center 3703 is at   0.000000  2.593000 -0.740800
      Read-in Center 3704 is at   0.000000 -2.645900 -0.687900
      Read-in Center 3705 is at   0.000000 -2.593000 -0.687900
      Read-in Center 3706 is at   0.000000 -2.540000 -0.687900
      Read-in Center 3707 is at   0.000000 -2.487100 -0.687900
      Read-in Center 3708 is at   0.000000 -2.434200 -0.687900
      Read-in Center 3709 is at   0.000000 -2.381300 -0.687900
      Read-in Center 3710 is at   0.000000 -2.328400 -0.687900
      Read-in Center 3711 is at   0.000000 -2.275500 -0.687900
      Read-in Center 3712 is at   0.000000 -2.222500 -0.687900
      Read-in Center 3713 is at   0.000000 -2.169600 -0.687900
      Read-in Center 3714 is at   0.000000 -2.116700 -0.687900
      Read-in Center 3715 is at   0.000000 -2.063800 -0.687900
      Read-in Center 3716 is at   0.000000 -2.010900 -0.687900
      Read-in Center 3717 is at   0.000000 -1.958000 -0.687900
      Read-in Center 3718 is at   0.000000 -1.905000 -0.687900
      Read-in Center 3719 is at   0.000000 -1.852100 -0.687900
      Read-in Center 3720 is at   0.000000 -1.799200 -0.687900
      Read-in Center 3721 is at   0.000000 -1.746300 -0.687900
      Read-in Center 3722 is at   0.000000 -1.693400 -0.687900
      Read-in Center 3723 is at   0.000000 -1.640400 -0.687900
      Read-in Center 3724 is at   0.000000 -1.587500 -0.687900
      Read-in Center 3725 is at   0.000000 -1.534600 -0.687900
      Read-in Center 3726 is at   0.000000 -1.481700 -0.687900
      Read-in Center 3727 is at   0.000000 -1.428800 -0.687900
      Read-in Center 3728 is at   0.000000 -1.375900 -0.687900
      Read-in Center 3729 is at   0.000000 -1.322900 -0.687900
      Read-in Center 3730 is at   0.000000 -1.270000 -0.687900
      Read-in Center 3731 is at   0.000000 -1.217100 -0.687900
      Read-in Center 3732 is at   0.000000 -1.164200 -0.687900
      Read-in Center 3733 is at   0.000000 -1.111300 -0.687900
      Read-in Center 3734 is at   0.000000 -1.058400 -0.687900
      Read-in Center 3735 is at   0.000000 -1.005400 -0.687900
      Read-in Center 3736 is at   0.000000 -0.952500 -0.687900
      Read-in Center 3737 is at   0.000000 -0.899600 -0.687900
      Read-in Center 3738 is at   0.000000 -0.846700 -0.687900
      Read-in Center 3739 is at   0.000000 -0.793800 -0.687900
      Read-in Center 3740 is at   0.000000 -0.740800 -0.687900
      Read-in Center 3741 is at   0.000000 -0.687900 -0.687900
      Read-in Center 3742 is at   0.000000 -0.635000 -0.687900
      Read-in Center 3743 is at   0.000000 -0.582100 -0.687900
      Read-in Center 3744 is at   0.000000 -0.529200 -0.687900
      Read-in Center 3745 is at   0.000000 -0.476300 -0.687900
      Read-in Center 3746 is at   0.000000 -0.423300 -0.687900
      Read-in Center 3747 is at   0.000000 -0.370400 -0.687900
      Read-in Center 3748 is at   0.000000 -0.317500 -0.687900
      Read-in Center 3749 is at   0.000000 -0.264600 -0.687900
      Read-in Center 3750 is at   0.000000 -0.211700 -0.687900
      Read-in Center 3751 is at   0.000000 -0.158800 -0.687900
      Read-in Center 3752 is at   0.000000 -0.105800 -0.687900
      Read-in Center 3753 is at   0.000000 -0.052900 -0.687900
      Read-in Center 3754 is at   0.000000  0.000000 -0.687900
      Read-in Center 3755 is at   0.000000  0.052900 -0.687900
      Read-in Center 3756 is at   0.000000  0.105800 -0.687900
      Read-in Center 3757 is at   0.000000  0.158800 -0.687900
      Read-in Center 3758 is at   0.000000  0.211700 -0.687900
      Read-in Center 3759 is at   0.000000  0.264600 -0.687900
      Read-in Center 3760 is at   0.000000  0.317500 -0.687900
      Read-in Center 3761 is at   0.000000  0.370400 -0.687900
      Read-in Center 3762 is at   0.000000  0.423300 -0.687900
      Read-in Center 3763 is at   0.000000  0.476300 -0.687900
      Read-in Center 3764 is at   0.000000  0.529200 -0.687900
      Read-in Center 3765 is at   0.000000  0.582100 -0.687900
      Read-in Center 3766 is at   0.000000  0.635000 -0.687900
      Read-in Center 3767 is at   0.000000  0.687900 -0.687900
      Read-in Center 3768 is at   0.000000  0.740800 -0.687900
      Read-in Center 3769 is at   0.000000  0.793800 -0.687900
      Read-in Center 3770 is at   0.000000  0.846700 -0.687900
      Read-in Center 3771 is at   0.000000  0.899600 -0.687900
      Read-in Center 3772 is at   0.000000  0.952500 -0.687900
      Read-in Center 3773 is at   0.000000  1.005400 -0.687900
      Read-in Center 3774 is at   0.000000  1.058400 -0.687900
      Read-in Center 3775 is at   0.000000  1.111300 -0.687900
      Read-in Center 3776 is at   0.000000  1.164200 -0.687900
      Read-in Center 3777 is at   0.000000  1.217100 -0.687900
      Read-in Center 3778 is at   0.000000  1.270000 -0.687900
      Read-in Center 3779 is at   0.000000  1.322900 -0.687900
      Read-in Center 3780 is at   0.000000  1.375900 -0.687900
      Read-in Center 3781 is at   0.000000  1.428800 -0.687900
      Read-in Center 3782 is at   0.000000  1.481700 -0.687900
      Read-in Center 3783 is at   0.000000  1.534600 -0.687900
      Read-in Center 3784 is at   0.000000  1.587500 -0.687900
      Read-in Center 3785 is at   0.000000  1.640400 -0.687900
      Read-in Center 3786 is at   0.000000  1.693400 -0.687900
      Read-in Center 3787 is at   0.000000  1.746300 -0.687900
      Read-in Center 3788 is at   0.000000  1.799200 -0.687900
      Read-in Center 3789 is at   0.000000  1.852100 -0.687900
      Read-in Center 3790 is at   0.000000  1.905000 -0.687900
      Read-in Center 3791 is at   0.000000  1.958000 -0.687900
      Read-in Center 3792 is at   0.000000  2.010900 -0.687900
      Read-in Center 3793 is at   0.000000  2.063800 -0.687900
      Read-in Center 3794 is at   0.000000  2.116700 -0.687900
      Read-in Center 3795 is at   0.000000  2.169600 -0.687900
      Read-in Center 3796 is at   0.000000  2.222500 -0.687900
      Read-in Center 3797 is at   0.000000  2.275500 -0.687900
      Read-in Center 3798 is at   0.000000  2.328400 -0.687900
      Read-in Center 3799 is at   0.000000  2.381300 -0.687900
      Read-in Center 3800 is at   0.000000  2.434200 -0.687900
      Read-in Center 3801 is at   0.000000  2.487100 -0.687900
      Read-in Center 3802 is at   0.000000  2.540000 -0.687900
      Read-in Center 3803 is at   0.000000  2.593000 -0.687900
      Read-in Center 3804 is at   0.000000 -2.645900 -0.635000
      Read-in Center 3805 is at   0.000000 -2.593000 -0.635000
      Read-in Center 3806 is at   0.000000 -2.540000 -0.635000
      Read-in Center 3807 is at   0.000000 -2.487100 -0.635000
      Read-in Center 3808 is at   0.000000 -2.434200 -0.635000
      Read-in Center 3809 is at   0.000000 -2.381300 -0.635000
      Read-in Center 3810 is at   0.000000 -2.328400 -0.635000
      Read-in Center 3811 is at   0.000000 -2.275500 -0.635000
      Read-in Center 3812 is at   0.000000 -2.222500 -0.635000
      Read-in Center 3813 is at   0.000000 -2.169600 -0.635000
      Read-in Center 3814 is at   0.000000 -2.116700 -0.635000
      Read-in Center 3815 is at   0.000000 -2.063800 -0.635000
      Read-in Center 3816 is at   0.000000 -2.010900 -0.635000
      Read-in Center 3817 is at   0.000000 -1.958000 -0.635000
      Read-in Center 3818 is at   0.000000 -1.905000 -0.635000
      Read-in Center 3819 is at   0.000000 -1.852100 -0.635000
      Read-in Center 3820 is at   0.000000 -1.799200 -0.635000
      Read-in Center 3821 is at   0.000000 -1.746300 -0.635000
      Read-in Center 3822 is at   0.000000 -1.693400 -0.635000
      Read-in Center 3823 is at   0.000000 -1.640400 -0.635000
      Read-in Center 3824 is at   0.000000 -1.587500 -0.635000
      Read-in Center 3825 is at   0.000000 -1.534600 -0.635000
      Read-in Center 3826 is at   0.000000 -1.481700 -0.635000
      Read-in Center 3827 is at   0.000000 -1.428800 -0.635000
      Read-in Center 3828 is at   0.000000 -1.375900 -0.635000
      Read-in Center 3829 is at   0.000000 -1.322900 -0.635000
      Read-in Center 3830 is at   0.000000 -1.270000 -0.635000
      Read-in Center 3831 is at   0.000000 -1.217100 -0.635000
      Read-in Center 3832 is at   0.000000 -1.164200 -0.635000
      Read-in Center 3833 is at   0.000000 -1.111300 -0.635000
      Read-in Center 3834 is at   0.000000 -1.058400 -0.635000
      Read-in Center 3835 is at   0.000000 -1.005400 -0.635000
      Read-in Center 3836 is at   0.000000 -0.952500 -0.635000
      Read-in Center 3837 is at   0.000000 -0.899600 -0.635000
      Read-in Center 3838 is at   0.000000 -0.846700 -0.635000
      Read-in Center 3839 is at   0.000000 -0.793800 -0.635000
      Read-in Center 3840 is at   0.000000 -0.740800 -0.635000
      Read-in Center 3841 is at   0.000000 -0.687900 -0.635000
      Read-in Center 3842 is at   0.000000 -0.635000 -0.635000
      Read-in Center 3843 is at   0.000000 -0.582100 -0.635000
      Read-in Center 3844 is at   0.000000 -0.529200 -0.635000
      Read-in Center 3845 is at   0.000000 -0.476300 -0.635000
      Read-in Center 3846 is at   0.000000 -0.423300 -0.635000
      Read-in Center 3847 is at   0.000000 -0.370400 -0.635000
      Read-in Center 3848 is at   0.000000 -0.317500 -0.635000
      Read-in Center 3849 is at   0.000000 -0.264600 -0.635000
      Read-in Center 3850 is at   0.000000 -0.211700 -0.635000
      Read-in Center 3851 is at   0.000000 -0.158800 -0.635000
      Read-in Center 3852 is at   0.000000 -0.105800 -0.635000
      Read-in Center 3853 is at   0.000000 -0.052900 -0.635000
      Read-in Center 3854 is at   0.000000  0.000000 -0.635000
      Read-in Center 3855 is at   0.000000  0.052900 -0.635000
      Read-in Center 3856 is at   0.000000  0.105800 -0.635000
      Read-in Center 3857 is at   0.000000  0.158800 -0.635000
      Read-in Center 3858 is at   0.000000  0.211700 -0.635000
      Read-in Center 3859 is at   0.000000  0.264600 -0.635000
      Read-in Center 3860 is at   0.000000  0.317500 -0.635000
      Read-in Center 3861 is at   0.000000  0.370400 -0.635000
      Read-in Center 3862 is at   0.000000  0.423300 -0.635000
      Read-in Center 3863 is at   0.000000  0.476300 -0.635000
      Read-in Center 3864 is at   0.000000  0.529200 -0.635000
      Read-in Center 3865 is at   0.000000  0.582100 -0.635000
      Read-in Center 3866 is at   0.000000  0.635000 -0.635000
      Read-in Center 3867 is at   0.000000  0.687900 -0.635000
      Read-in Center 3868 is at   0.000000  0.740800 -0.635000
      Read-in Center 3869 is at   0.000000  0.793800 -0.635000
      Read-in Center 3870 is at   0.000000  0.846700 -0.635000
      Read-in Center 3871 is at   0.000000  0.899600 -0.635000
      Read-in Center 3872 is at   0.000000  0.952500 -0.635000
      Read-in Center 3873 is at   0.000000  1.005400 -0.635000
      Read-in Center 3874 is at   0.000000  1.058400 -0.635000
      Read-in Center 3875 is at   0.000000  1.111300 -0.635000
      Read-in Center 3876 is at   0.000000  1.164200 -0.635000
      Read-in Center 3877 is at   0.000000  1.217100 -0.635000
      Read-in Center 3878 is at   0.000000  1.270000 -0.635000
      Read-in Center 3879 is at   0.000000  1.322900 -0.635000
      Read-in Center 3880 is at   0.000000  1.375900 -0.635000
      Read-in Center 3881 is at   0.000000  1.428800 -0.635000
      Read-in Center 3882 is at   0.000000  1.481700 -0.635000
      Read-in Center 3883 is at   0.000000  1.534600 -0.635000
      Read-in Center 3884 is at   0.000000  1.587500 -0.635000
      Read-in Center 3885 is at   0.000000  1.640400 -0.635000
      Read-in Center 3886 is at   0.000000  1.693400 -0.635000
      Read-in Center 3887 is at   0.000000  1.746300 -0.635000
      Read-in Center 3888 is at   0.000000  1.799200 -0.635000
      Read-in Center 3889 is at   0.000000  1.852100 -0.635000
      Read-in Center 3890 is at   0.000000  1.905000 -0.635000
      Read-in Center 3891 is at   0.000000  1.958000 -0.635000
      Read-in Center 3892 is at   0.000000  2.010900 -0.635000
      Read-in Center 3893 is at   0.000000  2.063800 -0.635000
      Read-in Center 3894 is at   0.000000  2.116700 -0.635000
      Read-in Center 3895 is at   0.000000  2.169600 -0.635000
      Read-in Center 3896 is at   0.000000  2.222500 -0.635000
      Read-in Center 3897 is at   0.000000  2.275500 -0.635000
      Read-in Center 3898 is at   0.000000  2.328400 -0.635000
      Read-in Center 3899 is at   0.000000  2.381300 -0.635000
      Read-in Center 3900 is at   0.000000  2.434200 -0.635000
      Read-in Center 3901 is at   0.000000  2.487100 -0.635000
      Read-in Center 3902 is at   0.000000  2.540000 -0.635000
      Read-in Center 3903 is at   0.000000  2.593000 -0.635000
      Read-in Center 3904 is at   0.000000 -2.645900 -0.582100
      Read-in Center 3905 is at   0.000000 -2.593000 -0.582100
      Read-in Center 3906 is at   0.000000 -2.540000 -0.582100
      Read-in Center 3907 is at   0.000000 -2.487100 -0.582100
      Read-in Center 3908 is at   0.000000 -2.434200 -0.582100
      Read-in Center 3909 is at   0.000000 -2.381300 -0.582100
      Read-in Center 3910 is at   0.000000 -2.328400 -0.582100
      Read-in Center 3911 is at   0.000000 -2.275500 -0.582100
      Read-in Center 3912 is at   0.000000 -2.222500 -0.582100
      Read-in Center 3913 is at   0.000000 -2.169600 -0.582100
      Read-in Center 3914 is at   0.000000 -2.116700 -0.582100
      Read-in Center 3915 is at   0.000000 -2.063800 -0.582100
      Read-in Center 3916 is at   0.000000 -2.010900 -0.582100
      Read-in Center 3917 is at   0.000000 -1.958000 -0.582100
      Read-in Center 3918 is at   0.000000 -1.905000 -0.582100
      Read-in Center 3919 is at   0.000000 -1.852100 -0.582100
      Read-in Center 3920 is at   0.000000 -1.799200 -0.582100
      Read-in Center 3921 is at   0.000000 -1.746300 -0.582100
      Read-in Center 3922 is at   0.000000 -1.693400 -0.582100
      Read-in Center 3923 is at   0.000000 -1.640400 -0.582100
      Read-in Center 3924 is at   0.000000 -1.587500 -0.582100
      Read-in Center 3925 is at   0.000000 -1.534600 -0.582100
      Read-in Center 3926 is at   0.000000 -1.481700 -0.582100
      Read-in Center 3927 is at   0.000000 -1.428800 -0.582100
      Read-in Center 3928 is at   0.000000 -1.375900 -0.582100
      Read-in Center 3929 is at   0.000000 -1.322900 -0.582100
      Read-in Center 3930 is at   0.000000 -1.270000 -0.582100
      Read-in Center 3931 is at   0.000000 -1.217100 -0.582100
      Read-in Center 3932 is at   0.000000 -1.164200 -0.582100
      Read-in Center 3933 is at   0.000000 -1.111300 -0.582100
      Read-in Center 3934 is at   0.000000 -1.058400 -0.582100
      Read-in Center 3935 is at   0.000000 -1.005400 -0.582100
      Read-in Center 3936 is at   0.000000 -0.952500 -0.582100
      Read-in Center 3937 is at   0.000000 -0.899600 -0.582100
      Read-in Center 3938 is at   0.000000 -0.846700 -0.582100
      Read-in Center 3939 is at   0.000000 -0.793800 -0.582100
      Read-in Center 3940 is at   0.000000 -0.740800 -0.582100
      Read-in Center 3941 is at   0.000000 -0.687900 -0.582100
      Read-in Center 3942 is at   0.000000 -0.635000 -0.582100
      Read-in Center 3943 is at   0.000000 -0.582100 -0.582100
      Read-in Center 3944 is at   0.000000 -0.529200 -0.582100
      Read-in Center 3945 is at   0.000000 -0.476300 -0.582100
      Read-in Center 3946 is at   0.000000 -0.423300 -0.582100
      Read-in Center 3947 is at   0.000000 -0.370400 -0.582100
      Read-in Center 3948 is at   0.000000 -0.317500 -0.582100
      Read-in Center 3949 is at   0.000000 -0.264600 -0.582100
      Read-in Center 3950 is at   0.000000 -0.211700 -0.582100
      Read-in Center 3951 is at   0.000000 -0.158800 -0.582100
      Read-in Center 3952 is at   0.000000 -0.105800 -0.582100
      Read-in Center 3953 is at   0.000000 -0.052900 -0.582100
      Read-in Center 3954 is at   0.000000  0.000000 -0.582100
      Read-in Center 3955 is at   0.000000  0.052900 -0.582100
      Read-in Center 3956 is at   0.000000  0.105800 -0.582100
      Read-in Center 3957 is at   0.000000  0.158800 -0.582100
      Read-in Center 3958 is at   0.000000  0.211700 -0.582100
      Read-in Center 3959 is at   0.000000  0.264600 -0.582100
      Read-in Center 3960 is at   0.000000  0.317500 -0.582100
      Read-in Center 3961 is at   0.000000  0.370400 -0.582100
      Read-in Center 3962 is at   0.000000  0.423300 -0.582100
      Read-in Center 3963 is at   0.000000  0.476300 -0.582100
      Read-in Center 3964 is at   0.000000  0.529200 -0.582100
      Read-in Center 3965 is at   0.000000  0.582100 -0.582100
      Read-in Center 3966 is at   0.000000  0.635000 -0.582100
      Read-in Center 3967 is at   0.000000  0.687900 -0.582100
      Read-in Center 3968 is at   0.000000  0.740800 -0.582100
      Read-in Center 3969 is at   0.000000  0.793800 -0.582100
      Read-in Center 3970 is at   0.000000  0.846700 -0.582100
      Read-in Center 3971 is at   0.000000  0.899600 -0.582100
      Read-in Center 3972 is at   0.000000  0.952500 -0.582100
      Read-in Center 3973 is at   0.000000  1.005400 -0.582100
      Read-in Center 3974 is at   0.000000  1.058400 -0.582100
      Read-in Center 3975 is at   0.000000  1.111300 -0.582100
      Read-in Center 3976 is at   0.000000  1.164200 -0.582100
      Read-in Center 3977 is at   0.000000  1.217100 -0.582100
      Read-in Center 3978 is at   0.000000  1.270000 -0.582100
      Read-in Center 3979 is at   0.000000  1.322900 -0.582100
      Read-in Center 3980 is at   0.000000  1.375900 -0.582100
      Read-in Center 3981 is at   0.000000  1.428800 -0.582100
      Read-in Center 3982 is at   0.000000  1.481700 -0.582100
      Read-in Center 3983 is at   0.000000  1.534600 -0.582100
      Read-in Center 3984 is at   0.000000  1.587500 -0.582100
      Read-in Center 3985 is at   0.000000  1.640400 -0.582100
      Read-in Center 3986 is at   0.000000  1.693400 -0.582100
      Read-in Center 3987 is at   0.000000  1.746300 -0.582100
      Read-in Center 3988 is at   0.000000  1.799200 -0.582100
      Read-in Center 3989 is at   0.000000  1.852100 -0.582100
      Read-in Center 3990 is at   0.000000  1.905000 -0.582100
      Read-in Center 3991 is at   0.000000  1.958000 -0.582100
      Read-in Center 3992 is at   0.000000  2.010900 -0.582100
      Read-in Center 3993 is at   0.000000  2.063800 -0.582100
      Read-in Center 3994 is at   0.000000  2.116700 -0.582100
      Read-in Center 3995 is at   0.000000  2.169600 -0.582100
      Read-in Center 3996 is at   0.000000  2.222500 -0.582100
      Read-in Center 3997 is at   0.000000  2.275500 -0.582100
      Read-in Center 3998 is at   0.000000  2.328400 -0.582100
      Read-in Center 3999 is at   0.000000  2.381300 -0.582100
      Read-in Center 4000 is at   0.000000  2.434200 -0.582100
      Read-in Center 4001 is at   0.000000  2.487100 -0.582100
      Read-in Center 4002 is at   0.000000  2.540000 -0.582100
      Read-in Center 4003 is at   0.000000  2.593000 -0.582100
      Read-in Center 4004 is at   0.000000 -2.645900 -0.529200
      Read-in Center 4005 is at   0.000000 -2.593000 -0.529200
      Read-in Center 4006 is at   0.000000 -2.540000 -0.529200
      Read-in Center 4007 is at   0.000000 -2.487100 -0.529200
      Read-in Center 4008 is at   0.000000 -2.434200 -0.529200
      Read-in Center 4009 is at   0.000000 -2.381300 -0.529200
      Read-in Center 4010 is at   0.000000 -2.328400 -0.529200
      Read-in Center 4011 is at   0.000000 -2.275500 -0.529200
      Read-in Center 4012 is at   0.000000 -2.222500 -0.529200
      Read-in Center 4013 is at   0.000000 -2.169600 -0.529200
      Read-in Center 4014 is at   0.000000 -2.116700 -0.529200
      Read-in Center 4015 is at   0.000000 -2.063800 -0.529200
      Read-in Center 4016 is at   0.000000 -2.010900 -0.529200
      Read-in Center 4017 is at   0.000000 -1.958000 -0.529200
      Read-in Center 4018 is at   0.000000 -1.905000 -0.529200
      Read-in Center 4019 is at   0.000000 -1.852100 -0.529200
      Read-in Center 4020 is at   0.000000 -1.799200 -0.529200
      Read-in Center 4021 is at   0.000000 -1.746300 -0.529200
      Read-in Center 4022 is at   0.000000 -1.693400 -0.529200
      Read-in Center 4023 is at   0.000000 -1.640400 -0.529200
      Read-in Center 4024 is at   0.000000 -1.587500 -0.529200
      Read-in Center 4025 is at   0.000000 -1.534600 -0.529200
      Read-in Center 4026 is at   0.000000 -1.481700 -0.529200
      Read-in Center 4027 is at   0.000000 -1.428800 -0.529200
      Read-in Center 4028 is at   0.000000 -1.375900 -0.529200
      Read-in Center 4029 is at   0.000000 -1.322900 -0.529200
      Read-in Center 4030 is at   0.000000 -1.270000 -0.529200
      Read-in Center 4031 is at   0.000000 -1.217100 -0.529200
      Read-in Center 4032 is at   0.000000 -1.164200 -0.529200
      Read-in Center 4033 is at   0.000000 -1.111300 -0.529200
      Read-in Center 4034 is at   0.000000 -1.058400 -0.529200
      Read-in Center 4035 is at   0.000000 -1.005400 -0.529200
      Read-in Center 4036 is at   0.000000 -0.952500 -0.529200
      Read-in Center 4037 is at   0.000000 -0.899600 -0.529200
      Read-in Center 4038 is at   0.000000 -0.846700 -0.529200
      Read-in Center 4039 is at   0.000000 -0.793800 -0.529200
      Read-in Center 4040 is at   0.000000 -0.740800 -0.529200
      Read-in Center 4041 is at   0.000000 -0.687900 -0.529200
      Read-in Center 4042 is at   0.000000 -0.635000 -0.529200
      Read-in Center 4043 is at   0.000000 -0.582100 -0.529200
      Read-in Center 4044 is at   0.000000 -0.529200 -0.529200
      Read-in Center 4045 is at   0.000000 -0.476300 -0.529200
      Read-in Center 4046 is at   0.000000 -0.423300 -0.529200
      Read-in Center 4047 is at   0.000000 -0.370400 -0.529200
      Read-in Center 4048 is at   0.000000 -0.317500 -0.529200
      Read-in Center 4049 is at   0.000000 -0.264600 -0.529200
      Read-in Center 4050 is at   0.000000 -0.211700 -0.529200
      Read-in Center 4051 is at   0.000000 -0.158800 -0.529200
      Read-in Center 4052 is at   0.000000 -0.105800 -0.529200
      Read-in Center 4053 is at   0.000000 -0.052900 -0.529200
      Read-in Center 4054 is at   0.000000  0.000000 -0.529200
      Read-in Center 4055 is at   0.000000  0.052900 -0.529200
      Read-in Center 4056 is at   0.000000  0.105800 -0.529200
      Read-in Center 4057 is at   0.000000  0.158800 -0.529200
      Read-in Center 4058 is at   0.000000  0.211700 -0.529200
      Read-in Center 4059 is at   0.000000  0.264600 -0.529200
      Read-in Center 4060 is at   0.000000  0.317500 -0.529200
      Read-in Center 4061 is at   0.000000  0.370400 -0.529200
      Read-in Center 4062 is at   0.000000  0.423300 -0.529200
      Read-in Center 4063 is at   0.000000  0.476300 -0.529200
      Read-in Center 4064 is at   0.000000  0.529200 -0.529200
      Read-in Center 4065 is at   0.000000  0.582100 -0.529200
      Read-in Center 4066 is at   0.000000  0.635000 -0.529200
      Read-in Center 4067 is at   0.000000  0.687900 -0.529200
      Read-in Center 4068 is at   0.000000  0.740800 -0.529200
      Read-in Center 4069 is at   0.000000  0.793800 -0.529200
      Read-in Center 4070 is at   0.000000  0.846700 -0.529200
      Read-in Center 4071 is at   0.000000  0.899600 -0.529200
      Read-in Center 4072 is at   0.000000  0.952500 -0.529200
      Read-in Center 4073 is at   0.000000  1.005400 -0.529200
      Read-in Center 4074 is at   0.000000  1.058400 -0.529200
      Read-in Center 4075 is at   0.000000  1.111300 -0.529200
      Read-in Center 4076 is at   0.000000  1.164200 -0.529200
      Read-in Center 4077 is at   0.000000  1.217100 -0.529200
      Read-in Center 4078 is at   0.000000  1.270000 -0.529200
      Read-in Center 4079 is at   0.000000  1.322900 -0.529200
      Read-in Center 4080 is at   0.000000  1.375900 -0.529200
      Read-in Center 4081 is at   0.000000  1.428800 -0.529200
      Read-in Center 4082 is at   0.000000  1.481700 -0.529200
      Read-in Center 4083 is at   0.000000  1.534600 -0.529200
      Read-in Center 4084 is at   0.000000  1.587500 -0.529200
      Read-in Center 4085 is at   0.000000  1.640400 -0.529200
      Read-in Center 4086 is at   0.000000  1.693400 -0.529200
      Read-in Center 4087 is at   0.000000  1.746300 -0.529200
      Read-in Center 4088 is at   0.000000  1.799200 -0.529200
      Read-in Center 4089 is at   0.000000  1.852100 -0.529200
      Read-in Center 4090 is at   0.000000  1.905000 -0.529200
      Read-in Center 4091 is at   0.000000  1.958000 -0.529200
      Read-in Center 4092 is at   0.000000  2.010900 -0.529200
      Read-in Center 4093 is at   0.000000  2.063800 -0.529200
      Read-in Center 4094 is at   0.000000  2.116700 -0.529200
      Read-in Center 4095 is at   0.000000  2.169600 -0.529200
      Read-in Center 4096 is at   0.000000  2.222500 -0.529200
      Read-in Center 4097 is at   0.000000  2.275500 -0.529200
      Read-in Center 4098 is at   0.000000  2.328400 -0.529200
      Read-in Center 4099 is at   0.000000  2.381300 -0.529200
      Read-in Center 4100 is at   0.000000  2.434200 -0.529200
      Read-in Center 4101 is at   0.000000  2.487100 -0.529200
      Read-in Center 4102 is at   0.000000  2.540000 -0.529200
      Read-in Center 4103 is at   0.000000  2.593000 -0.529200
      Read-in Center 4104 is at   0.000000 -2.645900 -0.476300
      Read-in Center 4105 is at   0.000000 -2.593000 -0.476300
      Read-in Center 4106 is at   0.000000 -2.540000 -0.476300
      Read-in Center 4107 is at   0.000000 -2.487100 -0.476300
      Read-in Center 4108 is at   0.000000 -2.434200 -0.476300
      Read-in Center 4109 is at   0.000000 -2.381300 -0.476300
      Read-in Center 4110 is at   0.000000 -2.328400 -0.476300
      Read-in Center 4111 is at   0.000000 -2.275500 -0.476300
      Read-in Center 4112 is at   0.000000 -2.222500 -0.476300
      Read-in Center 4113 is at   0.000000 -2.169600 -0.476300
      Read-in Center 4114 is at   0.000000 -2.116700 -0.476300
      Read-in Center 4115 is at   0.000000 -2.063800 -0.476300
      Read-in Center 4116 is at   0.000000 -2.010900 -0.476300
      Read-in Center 4117 is at   0.000000 -1.958000 -0.476300
      Read-in Center 4118 is at   0.000000 -1.905000 -0.476300
      Read-in Center 4119 is at   0.000000 -1.852100 -0.476300
      Read-in Center 4120 is at   0.000000 -1.799200 -0.476300
      Read-in Center 4121 is at   0.000000 -1.746300 -0.476300
      Read-in Center 4122 is at   0.000000 -1.693400 -0.476300
      Read-in Center 4123 is at   0.000000 -1.640400 -0.476300
      Read-in Center 4124 is at   0.000000 -1.587500 -0.476300
      Read-in Center 4125 is at   0.000000 -1.534600 -0.476300
      Read-in Center 4126 is at   0.000000 -1.481700 -0.476300
      Read-in Center 4127 is at   0.000000 -1.428800 -0.476300
      Read-in Center 4128 is at   0.000000 -1.375900 -0.476300
      Read-in Center 4129 is at   0.000000 -1.322900 -0.476300
      Read-in Center 4130 is at   0.000000 -1.270000 -0.476300
      Read-in Center 4131 is at   0.000000 -1.217100 -0.476300
      Read-in Center 4132 is at   0.000000 -1.164200 -0.476300
      Read-in Center 4133 is at   0.000000 -1.111300 -0.476300
      Read-in Center 4134 is at   0.000000 -1.058400 -0.476300
      Read-in Center 4135 is at   0.000000 -1.005400 -0.476300
      Read-in Center 4136 is at   0.000000 -0.952500 -0.476300
      Read-in Center 4137 is at   0.000000 -0.899600 -0.476300
      Read-in Center 4138 is at   0.000000 -0.846700 -0.476300
      Read-in Center 4139 is at   0.000000 -0.793800 -0.476300
      Read-in Center 4140 is at   0.000000 -0.740800 -0.476300
      Read-in Center 4141 is at   0.000000 -0.687900 -0.476300
      Read-in Center 4142 is at   0.000000 -0.635000 -0.476300
      Read-in Center 4143 is at   0.000000 -0.582100 -0.476300
      Read-in Center 4144 is at   0.000000 -0.529200 -0.476300
      Read-in Center 4145 is at   0.000000 -0.476300 -0.476300
      Read-in Center 4146 is at   0.000000 -0.423300 -0.476300
      Read-in Center 4147 is at   0.000000 -0.370400 -0.476300
      Read-in Center 4148 is at   0.000000 -0.317500 -0.476300
      Read-in Center 4149 is at   0.000000 -0.264600 -0.476300
      Read-in Center 4150 is at   0.000000 -0.211700 -0.476300
      Read-in Center 4151 is at   0.000000 -0.158800 -0.476300
      Read-in Center 4152 is at   0.000000 -0.105800 -0.476300
      Read-in Center 4153 is at   0.000000 -0.052900 -0.476300
      Read-in Center 4154 is at   0.000000  0.000000 -0.476300
      Read-in Center 4155 is at   0.000000  0.052900 -0.476300
      Read-in Center 4156 is at   0.000000  0.105800 -0.476300
      Read-in Center 4157 is at   0.000000  0.158800 -0.476300
      Read-in Center 4158 is at   0.000000  0.211700 -0.476300
      Read-in Center 4159 is at   0.000000  0.264600 -0.476300
      Read-in Center 4160 is at   0.000000  0.317500 -0.476300
      Read-in Center 4161 is at   0.000000  0.370400 -0.476300
      Read-in Center 4162 is at   0.000000  0.423300 -0.476300
      Read-in Center 4163 is at   0.000000  0.476300 -0.476300
      Read-in Center 4164 is at   0.000000  0.529200 -0.476300
      Read-in Center 4165 is at   0.000000  0.582100 -0.476300
      Read-in Center 4166 is at   0.000000  0.635000 -0.476300
      Read-in Center 4167 is at   0.000000  0.687900 -0.476300
      Read-in Center 4168 is at   0.000000  0.740800 -0.476300
      Read-in Center 4169 is at   0.000000  0.793800 -0.476300
      Read-in Center 4170 is at   0.000000  0.846700 -0.476300
      Read-in Center 4171 is at   0.000000  0.899600 -0.476300
      Read-in Center 4172 is at   0.000000  0.952500 -0.476300
      Read-in Center 4173 is at   0.000000  1.005400 -0.476300
      Read-in Center 4174 is at   0.000000  1.058400 -0.476300
      Read-in Center 4175 is at   0.000000  1.111300 -0.476300
      Read-in Center 4176 is at   0.000000  1.164200 -0.476300
      Read-in Center 4177 is at   0.000000  1.217100 -0.476300
      Read-in Center 4178 is at   0.000000  1.270000 -0.476300
      Read-in Center 4179 is at   0.000000  1.322900 -0.476300
      Read-in Center 4180 is at   0.000000  1.375900 -0.476300
      Read-in Center 4181 is at   0.000000  1.428800 -0.476300
      Read-in Center 4182 is at   0.000000  1.481700 -0.476300
      Read-in Center 4183 is at   0.000000  1.534600 -0.476300
      Read-in Center 4184 is at   0.000000  1.587500 -0.476300
      Read-in Center 4185 is at   0.000000  1.640400 -0.476300
      Read-in Center 4186 is at   0.000000  1.693400 -0.476300
      Read-in Center 4187 is at   0.000000  1.746300 -0.476300
      Read-in Center 4188 is at   0.000000  1.799200 -0.476300
      Read-in Center 4189 is at   0.000000  1.852100 -0.476300
      Read-in Center 4190 is at   0.000000  1.905000 -0.476300
      Read-in Center 4191 is at   0.000000  1.958000 -0.476300
      Read-in Center 4192 is at   0.000000  2.010900 -0.476300
      Read-in Center 4193 is at   0.000000  2.063800 -0.476300
      Read-in Center 4194 is at   0.000000  2.116700 -0.476300
      Read-in Center 4195 is at   0.000000  2.169600 -0.476300
      Read-in Center 4196 is at   0.000000  2.222500 -0.476300
      Read-in Center 4197 is at   0.000000  2.275500 -0.476300
      Read-in Center 4198 is at   0.000000  2.328400 -0.476300
      Read-in Center 4199 is at   0.000000  2.381300 -0.476300
      Read-in Center 4200 is at   0.000000  2.434200 -0.476300
      Read-in Center 4201 is at   0.000000  2.487100 -0.476300
      Read-in Center 4202 is at   0.000000  2.540000 -0.476300
      Read-in Center 4203 is at   0.000000  2.593000 -0.476300
      Read-in Center 4204 is at   0.000000 -2.645900 -0.423300
      Read-in Center 4205 is at   0.000000 -2.593000 -0.423300
      Read-in Center 4206 is at   0.000000 -2.540000 -0.423300
      Read-in Center 4207 is at   0.000000 -2.487100 -0.423300
      Read-in Center 4208 is at   0.000000 -2.434200 -0.423300
      Read-in Center 4209 is at   0.000000 -2.381300 -0.423300
      Read-in Center 4210 is at   0.000000 -2.328400 -0.423300
      Read-in Center 4211 is at   0.000000 -2.275500 -0.423300
      Read-in Center 4212 is at   0.000000 -2.222500 -0.423300
      Read-in Center 4213 is at   0.000000 -2.169600 -0.423300
      Read-in Center 4214 is at   0.000000 -2.116700 -0.423300
      Read-in Center 4215 is at   0.000000 -2.063800 -0.423300
      Read-in Center 4216 is at   0.000000 -2.010900 -0.423300
      Read-in Center 4217 is at   0.000000 -1.958000 -0.423300
      Read-in Center 4218 is at   0.000000 -1.905000 -0.423300
      Read-in Center 4219 is at   0.000000 -1.852100 -0.423300
      Read-in Center 4220 is at   0.000000 -1.799200 -0.423300
      Read-in Center 4221 is at   0.000000 -1.746300 -0.423300
      Read-in Center 4222 is at   0.000000 -1.693400 -0.423300
      Read-in Center 4223 is at   0.000000 -1.640400 -0.423300
      Read-in Center 4224 is at   0.000000 -1.587500 -0.423300
      Read-in Center 4225 is at   0.000000 -1.534600 -0.423300
      Read-in Center 4226 is at   0.000000 -1.481700 -0.423300
      Read-in Center 4227 is at   0.000000 -1.428800 -0.423300
      Read-in Center 4228 is at   0.000000 -1.375900 -0.423300
      Read-in Center 4229 is at   0.000000 -1.322900 -0.423300
      Read-in Center 4230 is at   0.000000 -1.270000 -0.423300
      Read-in Center 4231 is at   0.000000 -1.217100 -0.423300
      Read-in Center 4232 is at   0.000000 -1.164200 -0.423300
      Read-in Center 4233 is at   0.000000 -1.111300 -0.423300
      Read-in Center 4234 is at   0.000000 -1.058400 -0.423300
      Read-in Center 4235 is at   0.000000 -1.005400 -0.423300
      Read-in Center 4236 is at   0.000000 -0.952500 -0.423300
      Read-in Center 4237 is at   0.000000 -0.899600 -0.423300
      Read-in Center 4238 is at   0.000000 -0.846700 -0.423300
      Read-in Center 4239 is at   0.000000 -0.793800 -0.423300
      Read-in Center 4240 is at   0.000000 -0.740800 -0.423300
      Read-in Center 4241 is at   0.000000 -0.687900 -0.423300
      Read-in Center 4242 is at   0.000000 -0.635000 -0.423300
      Read-in Center 4243 is at   0.000000 -0.582100 -0.423300
      Read-in Center 4244 is at   0.000000 -0.529200 -0.423300
      Read-in Center 4245 is at   0.000000 -0.476300 -0.423300
      Read-in Center 4246 is at   0.000000 -0.423300 -0.423300
      Read-in Center 4247 is at   0.000000 -0.370400 -0.423300
      Read-in Center 4248 is at   0.000000 -0.317500 -0.423300
      Read-in Center 4249 is at   0.000000 -0.264600 -0.423300
      Read-in Center 4250 is at   0.000000 -0.211700 -0.423300
      Read-in Center 4251 is at   0.000000 -0.158800 -0.423300
      Read-in Center 4252 is at   0.000000 -0.105800 -0.423300
      Read-in Center 4253 is at   0.000000 -0.052900 -0.423300
      Read-in Center 4254 is at   0.000000  0.000000 -0.423300
      Read-in Center 4255 is at   0.000000  0.052900 -0.423300
      Read-in Center 4256 is at   0.000000  0.105800 -0.423300
      Read-in Center 4257 is at   0.000000  0.158800 -0.423300
      Read-in Center 4258 is at   0.000000  0.211700 -0.423300
      Read-in Center 4259 is at   0.000000  0.264600 -0.423300
      Read-in Center 4260 is at   0.000000  0.317500 -0.423300
      Read-in Center 4261 is at   0.000000  0.370400 -0.423300
      Read-in Center 4262 is at   0.000000  0.423300 -0.423300
      Read-in Center 4263 is at   0.000000  0.476300 -0.423300
      Read-in Center 4264 is at   0.000000  0.529200 -0.423300
      Read-in Center 4265 is at   0.000000  0.582100 -0.423300
      Read-in Center 4266 is at   0.000000  0.635000 -0.423300
      Read-in Center 4267 is at   0.000000  0.687900 -0.423300
      Read-in Center 4268 is at   0.000000  0.740800 -0.423300
      Read-in Center 4269 is at   0.000000  0.793800 -0.423300
      Read-in Center 4270 is at   0.000000  0.846700 -0.423300
      Read-in Center 4271 is at   0.000000  0.899600 -0.423300
      Read-in Center 4272 is at   0.000000  0.952500 -0.423300
      Read-in Center 4273 is at   0.000000  1.005400 -0.423300
      Read-in Center 4274 is at   0.000000  1.058400 -0.423300
      Read-in Center 4275 is at   0.000000  1.111300 -0.423300
      Read-in Center 4276 is at   0.000000  1.164200 -0.423300
      Read-in Center 4277 is at   0.000000  1.217100 -0.423300
      Read-in Center 4278 is at   0.000000  1.270000 -0.423300
      Read-in Center 4279 is at   0.000000  1.322900 -0.423300
      Read-in Center 4280 is at   0.000000  1.375900 -0.423300
      Read-in Center 4281 is at   0.000000  1.428800 -0.423300
      Read-in Center 4282 is at   0.000000  1.481700 -0.423300
      Read-in Center 4283 is at   0.000000  1.534600 -0.423300
      Read-in Center 4284 is at   0.000000  1.587500 -0.423300
      Read-in Center 4285 is at   0.000000  1.640400 -0.423300
      Read-in Center 4286 is at   0.000000  1.693400 -0.423300
      Read-in Center 4287 is at   0.000000  1.746300 -0.423300
      Read-in Center 4288 is at   0.000000  1.799200 -0.423300
      Read-in Center 4289 is at   0.000000  1.852100 -0.423300
      Read-in Center 4290 is at   0.000000  1.905000 -0.423300
      Read-in Center 4291 is at   0.000000  1.958000 -0.423300
      Read-in Center 4292 is at   0.000000  2.010900 -0.423300
      Read-in Center 4293 is at   0.000000  2.063800 -0.423300
      Read-in Center 4294 is at   0.000000  2.116700 -0.423300
      Read-in Center 4295 is at   0.000000  2.169600 -0.423300
      Read-in Center 4296 is at   0.000000  2.222500 -0.423300
      Read-in Center 4297 is at   0.000000  2.275500 -0.423300
      Read-in Center 4298 is at   0.000000  2.328400 -0.423300
      Read-in Center 4299 is at   0.000000  2.381300 -0.423300
      Read-in Center 4300 is at   0.000000  2.434200 -0.423300
      Read-in Center 4301 is at   0.000000  2.487100 -0.423300
      Read-in Center 4302 is at   0.000000  2.540000 -0.423300
      Read-in Center 4303 is at   0.000000  2.593000 -0.423300
      Read-in Center 4304 is at   0.000000 -2.645900 -0.370400
      Read-in Center 4305 is at   0.000000 -2.593000 -0.370400
      Read-in Center 4306 is at   0.000000 -2.540000 -0.370400
      Read-in Center 4307 is at   0.000000 -2.487100 -0.370400
      Read-in Center 4308 is at   0.000000 -2.434200 -0.370400
      Read-in Center 4309 is at   0.000000 -2.381300 -0.370400
      Read-in Center 4310 is at   0.000000 -2.328400 -0.370400
      Read-in Center 4311 is at   0.000000 -2.275500 -0.370400
      Read-in Center 4312 is at   0.000000 -2.222500 -0.370400
      Read-in Center 4313 is at   0.000000 -2.169600 -0.370400
      Read-in Center 4314 is at   0.000000 -2.116700 -0.370400
      Read-in Center 4315 is at   0.000000 -2.063800 -0.370400
      Read-in Center 4316 is at   0.000000 -2.010900 -0.370400
      Read-in Center 4317 is at   0.000000 -1.958000 -0.370400
      Read-in Center 4318 is at   0.000000 -1.905000 -0.370400
      Read-in Center 4319 is at   0.000000 -1.852100 -0.370400
      Read-in Center 4320 is at   0.000000 -1.799200 -0.370400
      Read-in Center 4321 is at   0.000000 -1.746300 -0.370400
      Read-in Center 4322 is at   0.000000 -1.693400 -0.370400
      Read-in Center 4323 is at   0.000000 -1.640400 -0.370400
      Read-in Center 4324 is at   0.000000 -1.587500 -0.370400
      Read-in Center 4325 is at   0.000000 -1.534600 -0.370400
      Read-in Center 4326 is at   0.000000 -1.481700 -0.370400
      Read-in Center 4327 is at   0.000000 -1.428800 -0.370400
      Read-in Center 4328 is at   0.000000 -1.375900 -0.370400
      Read-in Center 4329 is at   0.000000 -1.322900 -0.370400
      Read-in Center 4330 is at   0.000000 -1.270000 -0.370400
      Read-in Center 4331 is at   0.000000 -1.217100 -0.370400
      Read-in Center 4332 is at   0.000000 -1.164200 -0.370400
      Read-in Center 4333 is at   0.000000 -1.111300 -0.370400
      Read-in Center 4334 is at   0.000000 -1.058400 -0.370400
      Read-in Center 4335 is at   0.000000 -1.005400 -0.370400
      Read-in Center 4336 is at   0.000000 -0.952500 -0.370400
      Read-in Center 4337 is at   0.000000 -0.899600 -0.370400
      Read-in Center 4338 is at   0.000000 -0.846700 -0.370400
      Read-in Center 4339 is at   0.000000 -0.793800 -0.370400
      Read-in Center 4340 is at   0.000000 -0.740800 -0.370400
      Read-in Center 4341 is at   0.000000 -0.687900 -0.370400
      Read-in Center 4342 is at   0.000000 -0.635000 -0.370400
      Read-in Center 4343 is at   0.000000 -0.582100 -0.370400
      Read-in Center 4344 is at   0.000000 -0.529200 -0.370400
      Read-in Center 4345 is at   0.000000 -0.476300 -0.370400
      Read-in Center 4346 is at   0.000000 -0.423300 -0.370400
      Read-in Center 4347 is at   0.000000 -0.370400 -0.370400
      Read-in Center 4348 is at   0.000000 -0.317500 -0.370400
      Read-in Center 4349 is at   0.000000 -0.264600 -0.370400
      Read-in Center 4350 is at   0.000000 -0.211700 -0.370400
      Read-in Center 4351 is at   0.000000 -0.158800 -0.370400
      Read-in Center 4352 is at   0.000000 -0.105800 -0.370400
      Read-in Center 4353 is at   0.000000 -0.052900 -0.370400
      Read-in Center 4354 is at   0.000000  0.000000 -0.370400
      Read-in Center 4355 is at   0.000000  0.052900 -0.370400
      Read-in Center 4356 is at   0.000000  0.105800 -0.370400
      Read-in Center 4357 is at   0.000000  0.158800 -0.370400
      Read-in Center 4358 is at   0.000000  0.211700 -0.370400
      Read-in Center 4359 is at   0.000000  0.264600 -0.370400
      Read-in Center 4360 is at   0.000000  0.317500 -0.370400
      Read-in Center 4361 is at   0.000000  0.370400 -0.370400
      Read-in Center 4362 is at   0.000000  0.423300 -0.370400
      Read-in Center 4363 is at   0.000000  0.476300 -0.370400
      Read-in Center 4364 is at   0.000000  0.529200 -0.370400
      Read-in Center 4365 is at   0.000000  0.582100 -0.370400
      Read-in Center 4366 is at   0.000000  0.635000 -0.370400
      Read-in Center 4367 is at   0.000000  0.687900 -0.370400
      Read-in Center 4368 is at   0.000000  0.740800 -0.370400
      Read-in Center 4369 is at   0.000000  0.793800 -0.370400
      Read-in Center 4370 is at   0.000000  0.846700 -0.370400
      Read-in Center 4371 is at   0.000000  0.899600 -0.370400
      Read-in Center 4372 is at   0.000000  0.952500 -0.370400
      Read-in Center 4373 is at   0.000000  1.005400 -0.370400
      Read-in Center 4374 is at   0.000000  1.058400 -0.370400
      Read-in Center 4375 is at   0.000000  1.111300 -0.370400
      Read-in Center 4376 is at   0.000000  1.164200 -0.370400
      Read-in Center 4377 is at   0.000000  1.217100 -0.370400
      Read-in Center 4378 is at   0.000000  1.270000 -0.370400
      Read-in Center 4379 is at   0.000000  1.322900 -0.370400
      Read-in Center 4380 is at   0.000000  1.375900 -0.370400
      Read-in Center 4381 is at   0.000000  1.428800 -0.370400
      Read-in Center 4382 is at   0.000000  1.481700 -0.370400
      Read-in Center 4383 is at   0.000000  1.534600 -0.370400
      Read-in Center 4384 is at   0.000000  1.587500 -0.370400
      Read-in Center 4385 is at   0.000000  1.640400 -0.370400
      Read-in Center 4386 is at   0.000000  1.693400 -0.370400
      Read-in Center 4387 is at   0.000000  1.746300 -0.370400
      Read-in Center 4388 is at   0.000000  1.799200 -0.370400
      Read-in Center 4389 is at   0.000000  1.852100 -0.370400
      Read-in Center 4390 is at   0.000000  1.905000 -0.370400
      Read-in Center 4391 is at   0.000000  1.958000 -0.370400
      Read-in Center 4392 is at   0.000000  2.010900 -0.370400
      Read-in Center 4393 is at   0.000000  2.063800 -0.370400
      Read-in Center 4394 is at   0.000000  2.116700 -0.370400
      Read-in Center 4395 is at   0.000000  2.169600 -0.370400
      Read-in Center 4396 is at   0.000000  2.222500 -0.370400
      Read-in Center 4397 is at   0.000000  2.275500 -0.370400
      Read-in Center 4398 is at   0.000000  2.328400 -0.370400
      Read-in Center 4399 is at   0.000000  2.381300 -0.370400
      Read-in Center 4400 is at   0.000000  2.434200 -0.370400
      Read-in Center 4401 is at   0.000000  2.487100 -0.370400
      Read-in Center 4402 is at   0.000000  2.540000 -0.370400
      Read-in Center 4403 is at   0.000000  2.593000 -0.370400
      Read-in Center 4404 is at   0.000000 -2.645900 -0.317500
      Read-in Center 4405 is at   0.000000 -2.593000 -0.317500
      Read-in Center 4406 is at   0.000000 -2.540000 -0.317500
      Read-in Center 4407 is at   0.000000 -2.487100 -0.317500
      Read-in Center 4408 is at   0.000000 -2.434200 -0.317500
      Read-in Center 4409 is at   0.000000 -2.381300 -0.317500
      Read-in Center 4410 is at   0.000000 -2.328400 -0.317500
      Read-in Center 4411 is at   0.000000 -2.275500 -0.317500
      Read-in Center 4412 is at   0.000000 -2.222500 -0.317500
      Read-in Center 4413 is at   0.000000 -2.169600 -0.317500
      Read-in Center 4414 is at   0.000000 -2.116700 -0.317500
      Read-in Center 4415 is at   0.000000 -2.063800 -0.317500
      Read-in Center 4416 is at   0.000000 -2.010900 -0.317500
      Read-in Center 4417 is at   0.000000 -1.958000 -0.317500
      Read-in Center 4418 is at   0.000000 -1.905000 -0.317500
      Read-in Center 4419 is at   0.000000 -1.852100 -0.317500
      Read-in Center 4420 is at   0.000000 -1.799200 -0.317500
      Read-in Center 4421 is at   0.000000 -1.746300 -0.317500
      Read-in Center 4422 is at   0.000000 -1.693400 -0.317500
      Read-in Center 4423 is at   0.000000 -1.640400 -0.317500
      Read-in Center 4424 is at   0.000000 -1.587500 -0.317500
      Read-in Center 4425 is at   0.000000 -1.534600 -0.317500
      Read-in Center 4426 is at   0.000000 -1.481700 -0.317500
      Read-in Center 4427 is at   0.000000 -1.428800 -0.317500
      Read-in Center 4428 is at   0.000000 -1.375900 -0.317500
      Read-in Center 4429 is at   0.000000 -1.322900 -0.317500
      Read-in Center 4430 is at   0.000000 -1.270000 -0.317500
      Read-in Center 4431 is at   0.000000 -1.217100 -0.317500
      Read-in Center 4432 is at   0.000000 -1.164200 -0.317500
      Read-in Center 4433 is at   0.000000 -1.111300 -0.317500
      Read-in Center 4434 is at   0.000000 -1.058400 -0.317500
      Read-in Center 4435 is at   0.000000 -1.005400 -0.317500
      Read-in Center 4436 is at   0.000000 -0.952500 -0.317500
      Read-in Center 4437 is at   0.000000 -0.899600 -0.317500
      Read-in Center 4438 is at   0.000000 -0.846700 -0.317500
      Read-in Center 4439 is at   0.000000 -0.793800 -0.317500
      Read-in Center 4440 is at   0.000000 -0.740800 -0.317500
      Read-in Center 4441 is at   0.000000 -0.687900 -0.317500
      Read-in Center 4442 is at   0.000000 -0.635000 -0.317500
      Read-in Center 4443 is at   0.000000 -0.582100 -0.317500
      Read-in Center 4444 is at   0.000000 -0.529200 -0.317500
      Read-in Center 4445 is at   0.000000 -0.476300 -0.317500
      Read-in Center 4446 is at   0.000000 -0.423300 -0.317500
      Read-in Center 4447 is at   0.000000 -0.370400 -0.317500
      Read-in Center 4448 is at   0.000000 -0.317500 -0.317500
      Read-in Center 4449 is at   0.000000 -0.264600 -0.317500
      Read-in Center 4450 is at   0.000000 -0.211700 -0.317500
      Read-in Center 4451 is at   0.000000 -0.158800 -0.317500
      Read-in Center 4452 is at   0.000000 -0.105800 -0.317500
      Read-in Center 4453 is at   0.000000 -0.052900 -0.317500
      Read-in Center 4454 is at   0.000000  0.000000 -0.317500
      Read-in Center 4455 is at   0.000000  0.052900 -0.317500
      Read-in Center 4456 is at   0.000000  0.105800 -0.317500
      Read-in Center 4457 is at   0.000000  0.158800 -0.317500
      Read-in Center 4458 is at   0.000000  0.211700 -0.317500
      Read-in Center 4459 is at   0.000000  0.264600 -0.317500
      Read-in Center 4460 is at   0.000000  0.317500 -0.317500
      Read-in Center 4461 is at   0.000000  0.370400 -0.317500
      Read-in Center 4462 is at   0.000000  0.423300 -0.317500
      Read-in Center 4463 is at   0.000000  0.476300 -0.317500
      Read-in Center 4464 is at   0.000000  0.529200 -0.317500
      Read-in Center 4465 is at   0.000000  0.582100 -0.317500
      Read-in Center 4466 is at   0.000000  0.635000 -0.317500
      Read-in Center 4467 is at   0.000000  0.687900 -0.317500
      Read-in Center 4468 is at   0.000000  0.740800 -0.317500
      Read-in Center 4469 is at   0.000000  0.793800 -0.317500
      Read-in Center 4470 is at   0.000000  0.846700 -0.317500
      Read-in Center 4471 is at   0.000000  0.899600 -0.317500
      Read-in Center 4472 is at   0.000000  0.952500 -0.317500
      Read-in Center 4473 is at   0.000000  1.005400 -0.317500
      Read-in Center 4474 is at   0.000000  1.058400 -0.317500
      Read-in Center 4475 is at   0.000000  1.111300 -0.317500
      Read-in Center 4476 is at   0.000000  1.164200 -0.317500
      Read-in Center 4477 is at   0.000000  1.217100 -0.317500
      Read-in Center 4478 is at   0.000000  1.270000 -0.317500
      Read-in Center 4479 is at   0.000000  1.322900 -0.317500
      Read-in Center 4480 is at   0.000000  1.375900 -0.317500
      Read-in Center 4481 is at   0.000000  1.428800 -0.317500
      Read-in Center 4482 is at   0.000000  1.481700 -0.317500
      Read-in Center 4483 is at   0.000000  1.534600 -0.317500
      Read-in Center 4484 is at   0.000000  1.587500 -0.317500
      Read-in Center 4485 is at   0.000000  1.640400 -0.317500
      Read-in Center 4486 is at   0.000000  1.693400 -0.317500
      Read-in Center 4487 is at   0.000000  1.746300 -0.317500
      Read-in Center 4488 is at   0.000000  1.799200 -0.317500
      Read-in Center 4489 is at   0.000000  1.852100 -0.317500
      Read-in Center 4490 is at   0.000000  1.905000 -0.317500
      Read-in Center 4491 is at   0.000000  1.958000 -0.317500
      Read-in Center 4492 is at   0.000000  2.010900 -0.317500
      Read-in Center 4493 is at   0.000000  2.063800 -0.317500
      Read-in Center 4494 is at   0.000000  2.116700 -0.317500
      Read-in Center 4495 is at   0.000000  2.169600 -0.317500
      Read-in Center 4496 is at   0.000000  2.222500 -0.317500
      Read-in Center 4497 is at   0.000000  2.275500 -0.317500
      Read-in Center 4498 is at   0.000000  2.328400 -0.317500
      Read-in Center 4499 is at   0.000000  2.381300 -0.317500
      Read-in Center 4500 is at   0.000000  2.434200 -0.317500
      Read-in Center 4501 is at   0.000000  2.487100 -0.317500
      Read-in Center 4502 is at   0.000000  2.540000 -0.317500
      Read-in Center 4503 is at   0.000000  2.593000 -0.317500
      Read-in Center 4504 is at   0.000000 -2.645900 -0.264600
      Read-in Center 4505 is at   0.000000 -2.593000 -0.264600
      Read-in Center 4506 is at   0.000000 -2.540000 -0.264600
      Read-in Center 4507 is at   0.000000 -2.487100 -0.264600
      Read-in Center 4508 is at   0.000000 -2.434200 -0.264600
      Read-in Center 4509 is at   0.000000 -2.381300 -0.264600
      Read-in Center 4510 is at   0.000000 -2.328400 -0.264600
      Read-in Center 4511 is at   0.000000 -2.275500 -0.264600
      Read-in Center 4512 is at   0.000000 -2.222500 -0.264600
      Read-in Center 4513 is at   0.000000 -2.169600 -0.264600
      Read-in Center 4514 is at   0.000000 -2.116700 -0.264600
      Read-in Center 4515 is at   0.000000 -2.063800 -0.264600
      Read-in Center 4516 is at   0.000000 -2.010900 -0.264600
      Read-in Center 4517 is at   0.000000 -1.958000 -0.264600
      Read-in Center 4518 is at   0.000000 -1.905000 -0.264600
      Read-in Center 4519 is at   0.000000 -1.852100 -0.264600
      Read-in Center 4520 is at   0.000000 -1.799200 -0.264600
      Read-in Center 4521 is at   0.000000 -1.746300 -0.264600
      Read-in Center 4522 is at   0.000000 -1.693400 -0.264600
      Read-in Center 4523 is at   0.000000 -1.640400 -0.264600
      Read-in Center 4524 is at   0.000000 -1.587500 -0.264600
      Read-in Center 4525 is at   0.000000 -1.534600 -0.264600
      Read-in Center 4526 is at   0.000000 -1.481700 -0.264600
      Read-in Center 4527 is at   0.000000 -1.428800 -0.264600
      Read-in Center 4528 is at   0.000000 -1.375900 -0.264600
      Read-in Center 4529 is at   0.000000 -1.322900 -0.264600
      Read-in Center 4530 is at   0.000000 -1.270000 -0.264600
      Read-in Center 4531 is at   0.000000 -1.217100 -0.264600
      Read-in Center 4532 is at   0.000000 -1.164200 -0.264600
      Read-in Center 4533 is at   0.000000 -1.111300 -0.264600
      Read-in Center 4534 is at   0.000000 -1.058400 -0.264600
      Read-in Center 4535 is at   0.000000 -1.005400 -0.264600
      Read-in Center 4536 is at   0.000000 -0.952500 -0.264600
      Read-in Center 4537 is at   0.000000 -0.899600 -0.264600
      Read-in Center 4538 is at   0.000000 -0.846700 -0.264600
      Read-in Center 4539 is at   0.000000 -0.793800 -0.264600
      Read-in Center 4540 is at   0.000000 -0.740800 -0.264600
      Read-in Center 4541 is at   0.000000 -0.687900 -0.264600
      Read-in Center 4542 is at   0.000000 -0.635000 -0.264600
      Read-in Center 4543 is at   0.000000 -0.582100 -0.264600
      Read-in Center 4544 is at   0.000000 -0.529200 -0.264600
      Read-in Center 4545 is at   0.000000 -0.476300 -0.264600
      Read-in Center 4546 is at   0.000000 -0.423300 -0.264600
      Read-in Center 4547 is at   0.000000 -0.370400 -0.264600
      Read-in Center 4548 is at   0.000000 -0.317500 -0.264600
      Read-in Center 4549 is at   0.000000 -0.264600 -0.264600
      Read-in Center 4550 is at   0.000000 -0.211700 -0.264600
      Read-in Center 4551 is at   0.000000 -0.158800 -0.264600
      Read-in Center 4552 is at   0.000000 -0.105800 -0.264600
      Read-in Center 4553 is at   0.000000 -0.052900 -0.264600
      Read-in Center 4554 is at   0.000000  0.000000 -0.264600
      Read-in Center 4555 is at   0.000000  0.052900 -0.264600
      Read-in Center 4556 is at   0.000000  0.105800 -0.264600
      Read-in Center 4557 is at   0.000000  0.158800 -0.264600
      Read-in Center 4558 is at   0.000000  0.211700 -0.264600
      Read-in Center 4559 is at   0.000000  0.264600 -0.264600
      Read-in Center 4560 is at   0.000000  0.317500 -0.264600
      Read-in Center 4561 is at   0.000000  0.370400 -0.264600
      Read-in Center 4562 is at   0.000000  0.423300 -0.264600
      Read-in Center 4563 is at   0.000000  0.476300 -0.264600
      Read-in Center 4564 is at   0.000000  0.529200 -0.264600
      Read-in Center 4565 is at   0.000000  0.582100 -0.264600
      Read-in Center 4566 is at   0.000000  0.635000 -0.264600
      Read-in Center 4567 is at   0.000000  0.687900 -0.264600
      Read-in Center 4568 is at   0.000000  0.740800 -0.264600
      Read-in Center 4569 is at   0.000000  0.793800 -0.264600
      Read-in Center 4570 is at   0.000000  0.846700 -0.264600
      Read-in Center 4571 is at   0.000000  0.899600 -0.264600
      Read-in Center 4572 is at   0.000000  0.952500 -0.264600
      Read-in Center 4573 is at   0.000000  1.005400 -0.264600
      Read-in Center 4574 is at   0.000000  1.058400 -0.264600
      Read-in Center 4575 is at   0.000000  1.111300 -0.264600
      Read-in Center 4576 is at   0.000000  1.164200 -0.264600
      Read-in Center 4577 is at   0.000000  1.217100 -0.264600
      Read-in Center 4578 is at   0.000000  1.270000 -0.264600
      Read-in Center 4579 is at   0.000000  1.322900 -0.264600
      Read-in Center 4580 is at   0.000000  1.375900 -0.264600
      Read-in Center 4581 is at   0.000000  1.428800 -0.264600
      Read-in Center 4582 is at   0.000000  1.481700 -0.264600
      Read-in Center 4583 is at   0.000000  1.534600 -0.264600
      Read-in Center 4584 is at   0.000000  1.587500 -0.264600
      Read-in Center 4585 is at   0.000000  1.640400 -0.264600
      Read-in Center 4586 is at   0.000000  1.693400 -0.264600
      Read-in Center 4587 is at   0.000000  1.746300 -0.264600
      Read-in Center 4588 is at   0.000000  1.799200 -0.264600
      Read-in Center 4589 is at   0.000000  1.852100 -0.264600
      Read-in Center 4590 is at   0.000000  1.905000 -0.264600
      Read-in Center 4591 is at   0.000000  1.958000 -0.264600
      Read-in Center 4592 is at   0.000000  2.010900 -0.264600
      Read-in Center 4593 is at   0.000000  2.063800 -0.264600
      Read-in Center 4594 is at   0.000000  2.116700 -0.264600
      Read-in Center 4595 is at   0.000000  2.169600 -0.264600
      Read-in Center 4596 is at   0.000000  2.222500 -0.264600
      Read-in Center 4597 is at   0.000000  2.275500 -0.264600
      Read-in Center 4598 is at   0.000000  2.328400 -0.264600
      Read-in Center 4599 is at   0.000000  2.381300 -0.264600
      Read-in Center 4600 is at   0.000000  2.434200 -0.264600
      Read-in Center 4601 is at   0.000000  2.487100 -0.264600
      Read-in Center 4602 is at   0.000000  2.540000 -0.264600
      Read-in Center 4603 is at   0.000000  2.593000 -0.264600
      Read-in Center 4604 is at   0.000000 -2.645900 -0.211700
      Read-in Center 4605 is at   0.000000 -2.593000 -0.211700
      Read-in Center 4606 is at   0.000000 -2.540000 -0.211700
      Read-in Center 4607 is at   0.000000 -2.487100 -0.211700
      Read-in Center 4608 is at   0.000000 -2.434200 -0.211700
      Read-in Center 4609 is at   0.000000 -2.381300 -0.211700
      Read-in Center 4610 is at   0.000000 -2.328400 -0.211700
      Read-in Center 4611 is at   0.000000 -2.275500 -0.211700
      Read-in Center 4612 is at   0.000000 -2.222500 -0.211700
      Read-in Center 4613 is at   0.000000 -2.169600 -0.211700
      Read-in Center 4614 is at   0.000000 -2.116700 -0.211700
      Read-in Center 4615 is at   0.000000 -2.063800 -0.211700
      Read-in Center 4616 is at   0.000000 -2.010900 -0.211700
      Read-in Center 4617 is at   0.000000 -1.958000 -0.211700
      Read-in Center 4618 is at   0.000000 -1.905000 -0.211700
      Read-in Center 4619 is at   0.000000 -1.852100 -0.211700
      Read-in Center 4620 is at   0.000000 -1.799200 -0.211700
      Read-in Center 4621 is at   0.000000 -1.746300 -0.211700
      Read-in Center 4622 is at   0.000000 -1.693400 -0.211700
      Read-in Center 4623 is at   0.000000 -1.640400 -0.211700
      Read-in Center 4624 is at   0.000000 -1.587500 -0.211700
      Read-in Center 4625 is at   0.000000 -1.534600 -0.211700
      Read-in Center 4626 is at   0.000000 -1.481700 -0.211700
      Read-in Center 4627 is at   0.000000 -1.428800 -0.211700
      Read-in Center 4628 is at   0.000000 -1.375900 -0.211700
      Read-in Center 4629 is at   0.000000 -1.322900 -0.211700
      Read-in Center 4630 is at   0.000000 -1.270000 -0.211700
      Read-in Center 4631 is at   0.000000 -1.217100 -0.211700
      Read-in Center 4632 is at   0.000000 -1.164200 -0.211700
      Read-in Center 4633 is at   0.000000 -1.111300 -0.211700
      Read-in Center 4634 is at   0.000000 -1.058400 -0.211700
      Read-in Center 4635 is at   0.000000 -1.005400 -0.211700
      Read-in Center 4636 is at   0.000000 -0.952500 -0.211700
      Read-in Center 4637 is at   0.000000 -0.899600 -0.211700
      Read-in Center 4638 is at   0.000000 -0.846700 -0.211700
      Read-in Center 4639 is at   0.000000 -0.793800 -0.211700
      Read-in Center 4640 is at   0.000000 -0.740800 -0.211700
      Read-in Center 4641 is at   0.000000 -0.687900 -0.211700
      Read-in Center 4642 is at   0.000000 -0.635000 -0.211700
      Read-in Center 4643 is at   0.000000 -0.582100 -0.211700
      Read-in Center 4644 is at   0.000000 -0.529200 -0.211700
      Read-in Center 4645 is at   0.000000 -0.476300 -0.211700
      Read-in Center 4646 is at   0.000000 -0.423300 -0.211700
      Read-in Center 4647 is at   0.000000 -0.370400 -0.211700
      Read-in Center 4648 is at   0.000000 -0.317500 -0.211700
      Read-in Center 4649 is at   0.000000 -0.264600 -0.211700
      Read-in Center 4650 is at   0.000000 -0.211700 -0.211700
      Read-in Center 4651 is at   0.000000 -0.158800 -0.211700
      Read-in Center 4652 is at   0.000000 -0.105800 -0.211700
      Read-in Center 4653 is at   0.000000 -0.052900 -0.211700
      Read-in Center 4654 is at   0.000000  0.000000 -0.211700
      Read-in Center 4655 is at   0.000000  0.052900 -0.211700
      Read-in Center 4656 is at   0.000000  0.105800 -0.211700
      Read-in Center 4657 is at   0.000000  0.158800 -0.211700
      Read-in Center 4658 is at   0.000000  0.211700 -0.211700
      Read-in Center 4659 is at   0.000000  0.264600 -0.211700
      Read-in Center 4660 is at   0.000000  0.317500 -0.211700
      Read-in Center 4661 is at   0.000000  0.370400 -0.211700
      Read-in Center 4662 is at   0.000000  0.423300 -0.211700
      Read-in Center 4663 is at   0.000000  0.476300 -0.211700
      Read-in Center 4664 is at   0.000000  0.529200 -0.211700
      Read-in Center 4665 is at   0.000000  0.582100 -0.211700
      Read-in Center 4666 is at   0.000000  0.635000 -0.211700
      Read-in Center 4667 is at   0.000000  0.687900 -0.211700
      Read-in Center 4668 is at   0.000000  0.740800 -0.211700
      Read-in Center 4669 is at   0.000000  0.793800 -0.211700
      Read-in Center 4670 is at   0.000000  0.846700 -0.211700
      Read-in Center 4671 is at   0.000000  0.899600 -0.211700
      Read-in Center 4672 is at   0.000000  0.952500 -0.211700
      Read-in Center 4673 is at   0.000000  1.005400 -0.211700
      Read-in Center 4674 is at   0.000000  1.058400 -0.211700
      Read-in Center 4675 is at   0.000000  1.111300 -0.211700
      Read-in Center 4676 is at   0.000000  1.164200 -0.211700
      Read-in Center 4677 is at   0.000000  1.217100 -0.211700
      Read-in Center 4678 is at   0.000000  1.270000 -0.211700
      Read-in Center 4679 is at   0.000000  1.322900 -0.211700
      Read-in Center 4680 is at   0.000000  1.375900 -0.211700
      Read-in Center 4681 is at   0.000000  1.428800 -0.211700
      Read-in Center 4682 is at   0.000000  1.481700 -0.211700
      Read-in Center 4683 is at   0.000000  1.534600 -0.211700
      Read-in Center 4684 is at   0.000000  1.587500 -0.211700
      Read-in Center 4685 is at   0.000000  1.640400 -0.211700
      Read-in Center 4686 is at   0.000000  1.693400 -0.211700
      Read-in Center 4687 is at   0.000000  1.746300 -0.211700
      Read-in Center 4688 is at   0.000000  1.799200 -0.211700
      Read-in Center 4689 is at   0.000000  1.852100 -0.211700
      Read-in Center 4690 is at   0.000000  1.905000 -0.211700
      Read-in Center 4691 is at   0.000000  1.958000 -0.211700
      Read-in Center 4692 is at   0.000000  2.010900 -0.211700
      Read-in Center 4693 is at   0.000000  2.063800 -0.211700
      Read-in Center 4694 is at   0.000000  2.116700 -0.211700
      Read-in Center 4695 is at   0.000000  2.169600 -0.211700
      Read-in Center 4696 is at   0.000000  2.222500 -0.211700
      Read-in Center 4697 is at   0.000000  2.275500 -0.211700
      Read-in Center 4698 is at   0.000000  2.328400 -0.211700
      Read-in Center 4699 is at   0.000000  2.381300 -0.211700
      Read-in Center 4700 is at   0.000000  2.434200 -0.211700
      Read-in Center 4701 is at   0.000000  2.487100 -0.211700
      Read-in Center 4702 is at   0.000000  2.540000 -0.211700
      Read-in Center 4703 is at   0.000000  2.593000 -0.211700
      Read-in Center 4704 is at   0.000000 -2.645900 -0.158800
      Read-in Center 4705 is at   0.000000 -2.593000 -0.158800
      Read-in Center 4706 is at   0.000000 -2.540000 -0.158800
      Read-in Center 4707 is at   0.000000 -2.487100 -0.158800
      Read-in Center 4708 is at   0.000000 -2.434200 -0.158800
      Read-in Center 4709 is at   0.000000 -2.381300 -0.158800
      Read-in Center 4710 is at   0.000000 -2.328400 -0.158800
      Read-in Center 4711 is at   0.000000 -2.275500 -0.158800
      Read-in Center 4712 is at   0.000000 -2.222500 -0.158800
      Read-in Center 4713 is at   0.000000 -2.169600 -0.158800
      Read-in Center 4714 is at   0.000000 -2.116700 -0.158800
      Read-in Center 4715 is at   0.000000 -2.063800 -0.158800
      Read-in Center 4716 is at   0.000000 -2.010900 -0.158800
      Read-in Center 4717 is at   0.000000 -1.958000 -0.158800
      Read-in Center 4718 is at   0.000000 -1.905000 -0.158800
      Read-in Center 4719 is at   0.000000 -1.852100 -0.158800
      Read-in Center 4720 is at   0.000000 -1.799200 -0.158800
      Read-in Center 4721 is at   0.000000 -1.746300 -0.158800
      Read-in Center 4722 is at   0.000000 -1.693400 -0.158800
      Read-in Center 4723 is at   0.000000 -1.640400 -0.158800
      Read-in Center 4724 is at   0.000000 -1.587500 -0.158800
      Read-in Center 4725 is at   0.000000 -1.534600 -0.158800
      Read-in Center 4726 is at   0.000000 -1.481700 -0.158800
      Read-in Center 4727 is at   0.000000 -1.428800 -0.158800
      Read-in Center 4728 is at   0.000000 -1.375900 -0.158800
      Read-in Center 4729 is at   0.000000 -1.322900 -0.158800
      Read-in Center 4730 is at   0.000000 -1.270000 -0.158800
      Read-in Center 4731 is at   0.000000 -1.217100 -0.158800
      Read-in Center 4732 is at   0.000000 -1.164200 -0.158800
      Read-in Center 4733 is at   0.000000 -1.111300 -0.158800
      Read-in Center 4734 is at   0.000000 -1.058400 -0.158800
      Read-in Center 4735 is at   0.000000 -1.005400 -0.158800
      Read-in Center 4736 is at   0.000000 -0.952500 -0.158800
      Read-in Center 4737 is at   0.000000 -0.899600 -0.158800
      Read-in Center 4738 is at   0.000000 -0.846700 -0.158800
      Read-in Center 4739 is at   0.000000 -0.793800 -0.158800
      Read-in Center 4740 is at   0.000000 -0.740800 -0.158800
      Read-in Center 4741 is at   0.000000 -0.687900 -0.158800
      Read-in Center 4742 is at   0.000000 -0.635000 -0.158800
      Read-in Center 4743 is at   0.000000 -0.582100 -0.158800
      Read-in Center 4744 is at   0.000000 -0.529200 -0.158800
      Read-in Center 4745 is at   0.000000 -0.476300 -0.158800
      Read-in Center 4746 is at   0.000000 -0.423300 -0.158800
      Read-in Center 4747 is at   0.000000 -0.370400 -0.158800
      Read-in Center 4748 is at   0.000000 -0.317500 -0.158800
      Read-in Center 4749 is at   0.000000 -0.264600 -0.158800
      Read-in Center 4750 is at   0.000000 -0.211700 -0.158800
      Read-in Center 4751 is at   0.000000 -0.158800 -0.158800
      Read-in Center 4752 is at   0.000000 -0.105800 -0.158800
      Read-in Center 4753 is at   0.000000 -0.052900 -0.158800
      Read-in Center 4754 is at   0.000000  0.000000 -0.158800
      Read-in Center 4755 is at   0.000000  0.052900 -0.158800
      Read-in Center 4756 is at   0.000000  0.105800 -0.158800
      Read-in Center 4757 is at   0.000000  0.158800 -0.158800
      Read-in Center 4758 is at   0.000000  0.211700 -0.158800
      Read-in Center 4759 is at   0.000000  0.264600 -0.158800
      Read-in Center 4760 is at   0.000000  0.317500 -0.158800
      Read-in Center 4761 is at   0.000000  0.370400 -0.158800
      Read-in Center 4762 is at   0.000000  0.423300 -0.158800
      Read-in Center 4763 is at   0.000000  0.476300 -0.158800
      Read-in Center 4764 is at   0.000000  0.529200 -0.158800
      Read-in Center 4765 is at   0.000000  0.582100 -0.158800
      Read-in Center 4766 is at   0.000000  0.635000 -0.158800
      Read-in Center 4767 is at   0.000000  0.687900 -0.158800
      Read-in Center 4768 is at   0.000000  0.740800 -0.158800
      Read-in Center 4769 is at   0.000000  0.793800 -0.158800
      Read-in Center 4770 is at   0.000000  0.846700 -0.158800
      Read-in Center 4771 is at   0.000000  0.899600 -0.158800
      Read-in Center 4772 is at   0.000000  0.952500 -0.158800
      Read-in Center 4773 is at   0.000000  1.005400 -0.158800
      Read-in Center 4774 is at   0.000000  1.058400 -0.158800
      Read-in Center 4775 is at   0.000000  1.111300 -0.158800
      Read-in Center 4776 is at   0.000000  1.164200 -0.158800
      Read-in Center 4777 is at   0.000000  1.217100 -0.158800
      Read-in Center 4778 is at   0.000000  1.270000 -0.158800
      Read-in Center 4779 is at   0.000000  1.322900 -0.158800
      Read-in Center 4780 is at   0.000000  1.375900 -0.158800
      Read-in Center 4781 is at   0.000000  1.428800 -0.158800
      Read-in Center 4782 is at   0.000000  1.481700 -0.158800
      Read-in Center 4783 is at   0.000000  1.534600 -0.158800
      Read-in Center 4784 is at   0.000000  1.587500 -0.158800
      Read-in Center 4785 is at   0.000000  1.640400 -0.158800
      Read-in Center 4786 is at   0.000000  1.693400 -0.158800
      Read-in Center 4787 is at   0.000000  1.746300 -0.158800
      Read-in Center 4788 is at   0.000000  1.799200 -0.158800
      Read-in Center 4789 is at   0.000000  1.852100 -0.158800
      Read-in Center 4790 is at   0.000000  1.905000 -0.158800
      Read-in Center 4791 is at   0.000000  1.958000 -0.158800
      Read-in Center 4792 is at   0.000000  2.010900 -0.158800
      Read-in Center 4793 is at   0.000000  2.063800 -0.158800
      Read-in Center 4794 is at   0.000000  2.116700 -0.158800
      Read-in Center 4795 is at   0.000000  2.169600 -0.158800
      Read-in Center 4796 is at   0.000000  2.222500 -0.158800
      Read-in Center 4797 is at   0.000000  2.275500 -0.158800
      Read-in Center 4798 is at   0.000000  2.328400 -0.158800
      Read-in Center 4799 is at   0.000000  2.381300 -0.158800
      Read-in Center 4800 is at   0.000000  2.434200 -0.158800
      Read-in Center 4801 is at   0.000000  2.487100 -0.158800
      Read-in Center 4802 is at   0.000000  2.540000 -0.158800
      Read-in Center 4803 is at   0.000000  2.593000 -0.158800
      Read-in Center 4804 is at   0.000000 -2.645900 -0.105800
      Read-in Center 4805 is at   0.000000 -2.593000 -0.105800
      Read-in Center 4806 is at   0.000000 -2.540000 -0.105800
      Read-in Center 4807 is at   0.000000 -2.487100 -0.105800
      Read-in Center 4808 is at   0.000000 -2.434200 -0.105800
      Read-in Center 4809 is at   0.000000 -2.381300 -0.105800
      Read-in Center 4810 is at   0.000000 -2.328400 -0.105800
      Read-in Center 4811 is at   0.000000 -2.275500 -0.105800
      Read-in Center 4812 is at   0.000000 -2.222500 -0.105800
      Read-in Center 4813 is at   0.000000 -2.169600 -0.105800
      Read-in Center 4814 is at   0.000000 -2.116700 -0.105800
      Read-in Center 4815 is at   0.000000 -2.063800 -0.105800
      Read-in Center 4816 is at   0.000000 -2.010900 -0.105800
      Read-in Center 4817 is at   0.000000 -1.958000 -0.105800
      Read-in Center 4818 is at   0.000000 -1.905000 -0.105800
      Read-in Center 4819 is at   0.000000 -1.852100 -0.105800
      Read-in Center 4820 is at   0.000000 -1.799200 -0.105800
      Read-in Center 4821 is at   0.000000 -1.746300 -0.105800
      Read-in Center 4822 is at   0.000000 -1.693400 -0.105800
      Read-in Center 4823 is at   0.000000 -1.640400 -0.105800
      Read-in Center 4824 is at   0.000000 -1.587500 -0.105800
      Read-in Center 4825 is at   0.000000 -1.534600 -0.105800
      Read-in Center 4826 is at   0.000000 -1.481700 -0.105800
      Read-in Center 4827 is at   0.000000 -1.428800 -0.105800
      Read-in Center 4828 is at   0.000000 -1.375900 -0.105800
      Read-in Center 4829 is at   0.000000 -1.322900 -0.105800
      Read-in Center 4830 is at   0.000000 -1.270000 -0.105800
      Read-in Center 4831 is at   0.000000 -1.217100 -0.105800
      Read-in Center 4832 is at   0.000000 -1.164200 -0.105800
      Read-in Center 4833 is at   0.000000 -1.111300 -0.105800
      Read-in Center 4834 is at   0.000000 -1.058400 -0.105800
      Read-in Center 4835 is at   0.000000 -1.005400 -0.105800
      Read-in Center 4836 is at   0.000000 -0.952500 -0.105800
      Read-in Center 4837 is at   0.000000 -0.899600 -0.105800
      Read-in Center 4838 is at   0.000000 -0.846700 -0.105800
      Read-in Center 4839 is at   0.000000 -0.793800 -0.105800
      Read-in Center 4840 is at   0.000000 -0.740800 -0.105800
      Read-in Center 4841 is at   0.000000 -0.687900 -0.105800
      Read-in Center 4842 is at   0.000000 -0.635000 -0.105800
      Read-in Center 4843 is at   0.000000 -0.582100 -0.105800
      Read-in Center 4844 is at   0.000000 -0.529200 -0.105800
      Read-in Center 4845 is at   0.000000 -0.476300 -0.105800
      Read-in Center 4846 is at   0.000000 -0.423300 -0.105800
      Read-in Center 4847 is at   0.000000 -0.370400 -0.105800
      Read-in Center 4848 is at   0.000000 -0.317500 -0.105800
      Read-in Center 4849 is at   0.000000 -0.264600 -0.105800
      Read-in Center 4850 is at   0.000000 -0.211700 -0.105800
      Read-in Center 4851 is at   0.000000 -0.158800 -0.105800
      Read-in Center 4852 is at   0.000000 -0.105800 -0.105800
      Read-in Center 4853 is at   0.000000 -0.052900 -0.105800
      Read-in Center 4854 is at   0.000000  0.000000 -0.105800
      Read-in Center 4855 is at   0.000000  0.052900 -0.105800
      Read-in Center 4856 is at   0.000000  0.105800 -0.105800
      Read-in Center 4857 is at   0.000000  0.158800 -0.105800
      Read-in Center 4858 is at   0.000000  0.211700 -0.105800
      Read-in Center 4859 is at   0.000000  0.264600 -0.105800
      Read-in Center 4860 is at   0.000000  0.317500 -0.105800
      Read-in Center 4861 is at   0.000000  0.370400 -0.105800
      Read-in Center 4862 is at   0.000000  0.423300 -0.105800
      Read-in Center 4863 is at   0.000000  0.476300 -0.105800
      Read-in Center 4864 is at   0.000000  0.529200 -0.105800
      Read-in Center 4865 is at   0.000000  0.582100 -0.105800
      Read-in Center 4866 is at   0.000000  0.635000 -0.105800
      Read-in Center 4867 is at   0.000000  0.687900 -0.105800
      Read-in Center 4868 is at   0.000000  0.740800 -0.105800
      Read-in Center 4869 is at   0.000000  0.793800 -0.105800
      Read-in Center 4870 is at   0.000000  0.846700 -0.105800
      Read-in Center 4871 is at   0.000000  0.899600 -0.105800
      Read-in Center 4872 is at   0.000000  0.952500 -0.105800
      Read-in Center 4873 is at   0.000000  1.005400 -0.105800
      Read-in Center 4874 is at   0.000000  1.058400 -0.105800
      Read-in Center 4875 is at   0.000000  1.111300 -0.105800
      Read-in Center 4876 is at   0.000000  1.164200 -0.105800
      Read-in Center 4877 is at   0.000000  1.217100 -0.105800
      Read-in Center 4878 is at   0.000000  1.270000 -0.105800
      Read-in Center 4879 is at   0.000000  1.322900 -0.105800
      Read-in Center 4880 is at   0.000000  1.375900 -0.105800
      Read-in Center 4881 is at   0.000000  1.428800 -0.105800
      Read-in Center 4882 is at   0.000000  1.481700 -0.105800
      Read-in Center 4883 is at   0.000000  1.534600 -0.105800
      Read-in Center 4884 is at   0.000000  1.587500 -0.105800
      Read-in Center 4885 is at   0.000000  1.640400 -0.105800
      Read-in Center 4886 is at   0.000000  1.693400 -0.105800
      Read-in Center 4887 is at   0.000000  1.746300 -0.105800
      Read-in Center 4888 is at   0.000000  1.799200 -0.105800
      Read-in Center 4889 is at   0.000000  1.852100 -0.105800
      Read-in Center 4890 is at   0.000000  1.905000 -0.105800
      Read-in Center 4891 is at   0.000000  1.958000 -0.105800
      Read-in Center 4892 is at   0.000000  2.010900 -0.105800
      Read-in Center 4893 is at   0.000000  2.063800 -0.105800
      Read-in Center 4894 is at   0.000000  2.116700 -0.105800
      Read-in Center 4895 is at   0.000000  2.169600 -0.105800
      Read-in Center 4896 is at   0.000000  2.222500 -0.105800
      Read-in Center 4897 is at   0.000000  2.275500 -0.105800
      Read-in Center 4898 is at   0.000000  2.328400 -0.105800
      Read-in Center 4899 is at   0.000000  2.381300 -0.105800
      Read-in Center 4900 is at   0.000000  2.434200 -0.105800
      Read-in Center 4901 is at   0.000000  2.487100 -0.105800
      Read-in Center 4902 is at   0.000000  2.540000 -0.105800
      Read-in Center 4903 is at   0.000000  2.593000 -0.105800
      Read-in Center 4904 is at   0.000000 -2.645900 -0.052900
      Read-in Center 4905 is at   0.000000 -2.593000 -0.052900
      Read-in Center 4906 is at   0.000000 -2.540000 -0.052900
      Read-in Center 4907 is at   0.000000 -2.487100 -0.052900
      Read-in Center 4908 is at   0.000000 -2.434200 -0.052900
      Read-in Center 4909 is at   0.000000 -2.381300 -0.052900
      Read-in Center 4910 is at   0.000000 -2.328400 -0.052900
      Read-in Center 4911 is at   0.000000 -2.275500 -0.052900
      Read-in Center 4912 is at   0.000000 -2.222500 -0.052900
      Read-in Center 4913 is at   0.000000 -2.169600 -0.052900
      Read-in Center 4914 is at   0.000000 -2.116700 -0.052900
      Read-in Center 4915 is at   0.000000 -2.063800 -0.052900
      Read-in Center 4916 is at   0.000000 -2.010900 -0.052900
      Read-in Center 4917 is at   0.000000 -1.958000 -0.052900
      Read-in Center 4918 is at   0.000000 -1.905000 -0.052900
      Read-in Center 4919 is at   0.000000 -1.852100 -0.052900
      Read-in Center 4920 is at   0.000000 -1.799200 -0.052900
      Read-in Center 4921 is at   0.000000 -1.746300 -0.052900
      Read-in Center 4922 is at   0.000000 -1.693400 -0.052900
      Read-in Center 4923 is at   0.000000 -1.640400 -0.052900
      Read-in Center 4924 is at   0.000000 -1.587500 -0.052900
      Read-in Center 4925 is at   0.000000 -1.534600 -0.052900
      Read-in Center 4926 is at   0.000000 -1.481700 -0.052900
      Read-in Center 4927 is at   0.000000 -1.428800 -0.052900
      Read-in Center 4928 is at   0.000000 -1.375900 -0.052900
      Read-in Center 4929 is at   0.000000 -1.322900 -0.052900
      Read-in Center 4930 is at   0.000000 -1.270000 -0.052900
      Read-in Center 4931 is at   0.000000 -1.217100 -0.052900
      Read-in Center 4932 is at   0.000000 -1.164200 -0.052900
      Read-in Center 4933 is at   0.000000 -1.111300 -0.052900
      Read-in Center 4934 is at   0.000000 -1.058400 -0.052900
      Read-in Center 4935 is at   0.000000 -1.005400 -0.052900
      Read-in Center 4936 is at   0.000000 -0.952500 -0.052900
      Read-in Center 4937 is at   0.000000 -0.899600 -0.052900
      Read-in Center 4938 is at   0.000000 -0.846700 -0.052900
      Read-in Center 4939 is at   0.000000 -0.793800 -0.052900
      Read-in Center 4940 is at   0.000000 -0.740800 -0.052900
      Read-in Center 4941 is at   0.000000 -0.687900 -0.052900
      Read-in Center 4942 is at   0.000000 -0.635000 -0.052900
      Read-in Center 4943 is at   0.000000 -0.582100 -0.052900
      Read-in Center 4944 is at   0.000000 -0.529200 -0.052900
      Read-in Center 4945 is at   0.000000 -0.476300 -0.052900
      Read-in Center 4946 is at   0.000000 -0.423300 -0.052900
      Read-in Center 4947 is at   0.000000 -0.370400 -0.052900
      Read-in Center 4948 is at   0.000000 -0.317500 -0.052900
      Read-in Center 4949 is at   0.000000 -0.264600 -0.052900
      Read-in Center 4950 is at   0.000000 -0.211700 -0.052900
      Read-in Center 4951 is at   0.000000 -0.158800 -0.052900
      Read-in Center 4952 is at   0.000000 -0.105800 -0.052900
      Read-in Center 4953 is at   0.000000 -0.052900 -0.052900
      Read-in Center 4954 is at   0.000000  0.000000 -0.052900
      Read-in Center 4955 is at   0.000000  0.052900 -0.052900
      Read-in Center 4956 is at   0.000000  0.105800 -0.052900
      Read-in Center 4957 is at   0.000000  0.158800 -0.052900
      Read-in Center 4958 is at   0.000000  0.211700 -0.052900
      Read-in Center 4959 is at   0.000000  0.264600 -0.052900
      Read-in Center 4960 is at   0.000000  0.317500 -0.052900
      Read-in Center 4961 is at   0.000000  0.370400 -0.052900
      Read-in Center 4962 is at   0.000000  0.423300 -0.052900
      Read-in Center 4963 is at   0.000000  0.476300 -0.052900
      Read-in Center 4964 is at   0.000000  0.529200 -0.052900
      Read-in Center 4965 is at   0.000000  0.582100 -0.052900
      Read-in Center 4966 is at   0.000000  0.635000 -0.052900
      Read-in Center 4967 is at   0.000000  0.687900 -0.052900
      Read-in Center 4968 is at   0.000000  0.740800 -0.052900
      Read-in Center 4969 is at   0.000000  0.793800 -0.052900
      Read-in Center 4970 is at   0.000000  0.846700 -0.052900
      Read-in Center 4971 is at   0.000000  0.899600 -0.052900
      Read-in Center 4972 is at   0.000000  0.952500 -0.052900
      Read-in Center 4973 is at   0.000000  1.005400 -0.052900
      Read-in Center 4974 is at   0.000000  1.058400 -0.052900
      Read-in Center 4975 is at   0.000000  1.111300 -0.052900
      Read-in Center 4976 is at   0.000000  1.164200 -0.052900
      Read-in Center 4977 is at   0.000000  1.217100 -0.052900
      Read-in Center 4978 is at   0.000000  1.270000 -0.052900
      Read-in Center 4979 is at   0.000000  1.322900 -0.052900
      Read-in Center 4980 is at   0.000000  1.375900 -0.052900
      Read-in Center 4981 is at   0.000000  1.428800 -0.052900
      Read-in Center 4982 is at   0.000000  1.481700 -0.052900
      Read-in Center 4983 is at   0.000000  1.534600 -0.052900
      Read-in Center 4984 is at   0.000000  1.587500 -0.052900
      Read-in Center 4985 is at   0.000000  1.640400 -0.052900
      Read-in Center 4986 is at   0.000000  1.693400 -0.052900
      Read-in Center 4987 is at   0.000000  1.746300 -0.052900
      Read-in Center 4988 is at   0.000000  1.799200 -0.052900
      Read-in Center 4989 is at   0.000000  1.852100 -0.052900
      Read-in Center 4990 is at   0.000000  1.905000 -0.052900
      Read-in Center 4991 is at   0.000000  1.958000 -0.052900
      Read-in Center 4992 is at   0.000000  2.010900 -0.052900
      Read-in Center 4993 is at   0.000000  2.063800 -0.052900
      Read-in Center 4994 is at   0.000000  2.116700 -0.052900
      Read-in Center 4995 is at   0.000000  2.169600 -0.052900
      Read-in Center 4996 is at   0.000000  2.222500 -0.052900
      Read-in Center 4997 is at   0.000000  2.275500 -0.052900
      Read-in Center 4998 is at   0.000000  2.328400 -0.052900
      Read-in Center 4999 is at   0.000000  2.381300 -0.052900
      Read-in Center 5000 is at   0.000000  2.434200 -0.052900
      Read-in Center 5001 is at   0.000000  2.487100 -0.052900
      Read-in Center 5002 is at   0.000000  2.540000 -0.052900
      Read-in Center 5003 is at   0.000000  2.593000 -0.052900
      Read-in Center 5004 is at   0.000000 -2.645900  0.000000
      Read-in Center 5005 is at   0.000000 -2.593000  0.000000
      Read-in Center 5006 is at   0.000000 -2.540000  0.000000
      Read-in Center 5007 is at   0.000000 -2.487100  0.000000
      Read-in Center 5008 is at   0.000000 -2.434200  0.000000
      Read-in Center 5009 is at   0.000000 -2.381300  0.000000
      Read-in Center 5010 is at   0.000000 -2.328400  0.000000
      Read-in Center 5011 is at   0.000000 -2.275500  0.000000
      Read-in Center 5012 is at   0.000000 -2.222500  0.000000
      Read-in Center 5013 is at   0.000000 -2.169600  0.000000
      Read-in Center 5014 is at   0.000000 -2.116700  0.000000
      Read-in Center 5015 is at   0.000000 -2.063800  0.000000
      Read-in Center 5016 is at   0.000000 -2.010900  0.000000
      Read-in Center 5017 is at   0.000000 -1.958000  0.000000
      Read-in Center 5018 is at   0.000000 -1.905000  0.000000
      Read-in Center 5019 is at   0.000000 -1.852100  0.000000
      Read-in Center 5020 is at   0.000000 -1.799200  0.000000
      Read-in Center 5021 is at   0.000000 -1.746300  0.000000
      Read-in Center 5022 is at   0.000000 -1.693400  0.000000
      Read-in Center 5023 is at   0.000000 -1.640400  0.000000
      Read-in Center 5024 is at   0.000000 -1.587500  0.000000
      Read-in Center 5025 is at   0.000000 -1.534600  0.000000
      Read-in Center 5026 is at   0.000000 -1.481700  0.000000
      Read-in Center 5027 is at   0.000000 -1.428800  0.000000
      Read-in Center 5028 is at   0.000000 -1.375900  0.000000
      Read-in Center 5029 is at   0.000000 -1.322900  0.000000
      Read-in Center 5030 is at   0.000000 -1.270000  0.000000
      Read-in Center 5031 is at   0.000000 -1.217100  0.000000
      Read-in Center 5032 is at   0.000000 -1.164200  0.000000
      Read-in Center 5033 is at   0.000000 -1.111300  0.000000
      Read-in Center 5034 is at   0.000000 -1.058400  0.000000
      Read-in Center 5035 is at   0.000000 -1.005400  0.000000
      Read-in Center 5036 is at   0.000000 -0.952500  0.000000
      Read-in Center 5037 is at   0.000000 -0.899600  0.000000
      Read-in Center 5038 is at   0.000000 -0.846700  0.000000
      Read-in Center 5039 is at   0.000000 -0.793800  0.000000
      Read-in Center 5040 is at   0.000000 -0.740800  0.000000
      Read-in Center 5041 is at   0.000000 -0.687900  0.000000
      Read-in Center 5042 is at   0.000000 -0.635000  0.000000
      Read-in Center 5043 is at   0.000000 -0.582100  0.000000
      Read-in Center 5044 is at   0.000000 -0.529200  0.000000
      Read-in Center 5045 is at   0.000000 -0.476300  0.000000
      Read-in Center 5046 is at   0.000000 -0.423300  0.000000
      Read-in Center 5047 is at   0.000000 -0.370400  0.000000
      Read-in Center 5048 is at   0.000000 -0.317500  0.000000
      Read-in Center 5049 is at   0.000000 -0.264600  0.000000
      Read-in Center 5050 is at   0.000000 -0.211700  0.000000
      Read-in Center 5051 is at   0.000000 -0.158800  0.000000
      Read-in Center 5052 is at   0.000000 -0.105800  0.000000
      Read-in Center 5053 is at   0.000000 -0.052900  0.000000
      Read-in Center 5054 is at   0.000000  0.000000  0.000000
      Read-in Center 5055 is at   0.000000  0.052900  0.000000
      Read-in Center 5056 is at   0.000000  0.105800  0.000000
      Read-in Center 5057 is at   0.000000  0.158800  0.000000
      Read-in Center 5058 is at   0.000000  0.211700  0.000000
      Read-in Center 5059 is at   0.000000  0.264600  0.000000
      Read-in Center 5060 is at   0.000000  0.317500  0.000000
      Read-in Center 5061 is at   0.000000  0.370400  0.000000
      Read-in Center 5062 is at   0.000000  0.423300  0.000000
      Read-in Center 5063 is at   0.000000  0.476300  0.000000
      Read-in Center 5064 is at   0.000000  0.529200  0.000000
      Read-in Center 5065 is at   0.000000  0.582100  0.000000
      Read-in Center 5066 is at   0.000000  0.635000  0.000000
      Read-in Center 5067 is at   0.000000  0.687900  0.000000
      Read-in Center 5068 is at   0.000000  0.740800  0.000000
      Read-in Center 5069 is at   0.000000  0.793800  0.000000
      Read-in Center 5070 is at   0.000000  0.846700  0.000000
      Read-in Center 5071 is at   0.000000  0.899600  0.000000
      Read-in Center 5072 is at   0.000000  0.952500  0.000000
      Read-in Center 5073 is at   0.000000  1.005400  0.000000
      Read-in Center 5074 is at   0.000000  1.058400  0.000000
      Read-in Center 5075 is at   0.000000  1.111300  0.000000
      Read-in Center 5076 is at   0.000000  1.164200  0.000000
      Read-in Center 5077 is at   0.000000  1.217100  0.000000
      Read-in Center 5078 is at   0.000000  1.270000  0.000000
      Read-in Center 5079 is at   0.000000  1.322900  0.000000
      Read-in Center 5080 is at   0.000000  1.375900  0.000000
      Read-in Center 5081 is at   0.000000  1.428800  0.000000
      Read-in Center 5082 is at   0.000000  1.481700  0.000000
      Read-in Center 5083 is at   0.000000  1.534600  0.000000
      Read-in Center 5084 is at   0.000000  1.587500  0.000000
      Read-in Center 5085 is at   0.000000  1.640400  0.000000
      Read-in Center 5086 is at   0.000000  1.693400  0.000000
      Read-in Center 5087 is at   0.000000  1.746300  0.000000
      Read-in Center 5088 is at   0.000000  1.799200  0.000000
      Read-in Center 5089 is at   0.000000  1.852100  0.000000
      Read-in Center 5090 is at   0.000000  1.905000  0.000000
      Read-in Center 5091 is at   0.000000  1.958000  0.000000
      Read-in Center 5092 is at   0.000000  2.010900  0.000000
      Read-in Center 5093 is at   0.000000  2.063800  0.000000
      Read-in Center 5094 is at   0.000000  2.116700  0.000000
      Read-in Center 5095 is at   0.000000  2.169600  0.000000
      Read-in Center 5096 is at   0.000000  2.222500  0.000000
      Read-in Center 5097 is at   0.000000  2.275500  0.000000
      Read-in Center 5098 is at   0.000000  2.328400  0.000000
      Read-in Center 5099 is at   0.000000  2.381300  0.000000
      Read-in Center 5100 is at   0.000000  2.434200  0.000000
      Read-in Center 5101 is at   0.000000  2.487100  0.000000
      Read-in Center 5102 is at   0.000000  2.540000  0.000000
      Read-in Center 5103 is at   0.000000  2.593000  0.000000
      Read-in Center 5104 is at   0.000000 -2.645900  0.052900
      Read-in Center 5105 is at   0.000000 -2.593000  0.052900
      Read-in Center 5106 is at   0.000000 -2.540000  0.052900
      Read-in Center 5107 is at   0.000000 -2.487100  0.052900
      Read-in Center 5108 is at   0.000000 -2.434200  0.052900
      Read-in Center 5109 is at   0.000000 -2.381300  0.052900
      Read-in Center 5110 is at   0.000000 -2.328400  0.052900
      Read-in Center 5111 is at   0.000000 -2.275500  0.052900
      Read-in Center 5112 is at   0.000000 -2.222500  0.052900
      Read-in Center 5113 is at   0.000000 -2.169600  0.052900
      Read-in Center 5114 is at   0.000000 -2.116700  0.052900
      Read-in Center 5115 is at   0.000000 -2.063800  0.052900
      Read-in Center 5116 is at   0.000000 -2.010900  0.052900
      Read-in Center 5117 is at   0.000000 -1.958000  0.052900
      Read-in Center 5118 is at   0.000000 -1.905000  0.052900
      Read-in Center 5119 is at   0.000000 -1.852100  0.052900
      Read-in Center 5120 is at   0.000000 -1.799200  0.052900
      Read-in Center 5121 is at   0.000000 -1.746300  0.052900
      Read-in Center 5122 is at   0.000000 -1.693400  0.052900
      Read-in Center 5123 is at   0.000000 -1.640400  0.052900
      Read-in Center 5124 is at   0.000000 -1.587500  0.052900
      Read-in Center 5125 is at   0.000000 -1.534600  0.052900
      Read-in Center 5126 is at   0.000000 -1.481700  0.052900
      Read-in Center 5127 is at   0.000000 -1.428800  0.052900
      Read-in Center 5128 is at   0.000000 -1.375900  0.052900
      Read-in Center 5129 is at   0.000000 -1.322900  0.052900
      Read-in Center 5130 is at   0.000000 -1.270000  0.052900
      Read-in Center 5131 is at   0.000000 -1.217100  0.052900
      Read-in Center 5132 is at   0.000000 -1.164200  0.052900
      Read-in Center 5133 is at   0.000000 -1.111300  0.052900
      Read-in Center 5134 is at   0.000000 -1.058400  0.052900
      Read-in Center 5135 is at   0.000000 -1.005400  0.052900
      Read-in Center 5136 is at   0.000000 -0.952500  0.052900
      Read-in Center 5137 is at   0.000000 -0.899600  0.052900
      Read-in Center 5138 is at   0.000000 -0.846700  0.052900
      Read-in Center 5139 is at   0.000000 -0.793800  0.052900
      Read-in Center 5140 is at   0.000000 -0.740800  0.052900
      Read-in Center 5141 is at   0.000000 -0.687900  0.052900
      Read-in Center 5142 is at   0.000000 -0.635000  0.052900
      Read-in Center 5143 is at   0.000000 -0.582100  0.052900
      Read-in Center 5144 is at   0.000000 -0.529200  0.052900
      Read-in Center 5145 is at   0.000000 -0.476300  0.052900
      Read-in Center 5146 is at   0.000000 -0.423300  0.052900
      Read-in Center 5147 is at   0.000000 -0.370400  0.052900
      Read-in Center 5148 is at   0.000000 -0.317500  0.052900
      Read-in Center 5149 is at   0.000000 -0.264600  0.052900
      Read-in Center 5150 is at   0.000000 -0.211700  0.052900
      Read-in Center 5151 is at   0.000000 -0.158800  0.052900
      Read-in Center 5152 is at   0.000000 -0.105800  0.052900
      Read-in Center 5153 is at   0.000000 -0.052900  0.052900
      Read-in Center 5154 is at   0.000000  0.000000  0.052900
      Read-in Center 5155 is at   0.000000  0.052900  0.052900
      Read-in Center 5156 is at   0.000000  0.105800  0.052900
      Read-in Center 5157 is at   0.000000  0.158800  0.052900
      Read-in Center 5158 is at   0.000000  0.211700  0.052900
      Read-in Center 5159 is at   0.000000  0.264600  0.052900
      Read-in Center 5160 is at   0.000000  0.317500  0.052900
      Read-in Center 5161 is at   0.000000  0.370400  0.052900
      Read-in Center 5162 is at   0.000000  0.423300  0.052900
      Read-in Center 5163 is at   0.000000  0.476300  0.052900
      Read-in Center 5164 is at   0.000000  0.529200  0.052900
      Read-in Center 5165 is at   0.000000  0.582100  0.052900
      Read-in Center 5166 is at   0.000000  0.635000  0.052900
      Read-in Center 5167 is at   0.000000  0.687900  0.052900
      Read-in Center 5168 is at   0.000000  0.740800  0.052900
      Read-in Center 5169 is at   0.000000  0.793800  0.052900
      Read-in Center 5170 is at   0.000000  0.846700  0.052900
      Read-in Center 5171 is at   0.000000  0.899600  0.052900
      Read-in Center 5172 is at   0.000000  0.952500  0.052900
      Read-in Center 5173 is at   0.000000  1.005400  0.052900
      Read-in Center 5174 is at   0.000000  1.058400  0.052900
      Read-in Center 5175 is at   0.000000  1.111300  0.052900
      Read-in Center 5176 is at   0.000000  1.164200  0.052900
      Read-in Center 5177 is at   0.000000  1.217100  0.052900
      Read-in Center 5178 is at   0.000000  1.270000  0.052900
      Read-in Center 5179 is at   0.000000  1.322900  0.052900
      Read-in Center 5180 is at   0.000000  1.375900  0.052900
      Read-in Center 5181 is at   0.000000  1.428800  0.052900
      Read-in Center 5182 is at   0.000000  1.481700  0.052900
      Read-in Center 5183 is at   0.000000  1.534600  0.052900
      Read-in Center 5184 is at   0.000000  1.587500  0.052900
      Read-in Center 5185 is at   0.000000  1.640400  0.052900
      Read-in Center 5186 is at   0.000000  1.693400  0.052900
      Read-in Center 5187 is at   0.000000  1.746300  0.052900
      Read-in Center 5188 is at   0.000000  1.799200  0.052900
      Read-in Center 5189 is at   0.000000  1.852100  0.052900
      Read-in Center 5190 is at   0.000000  1.905000  0.052900
      Read-in Center 5191 is at   0.000000  1.958000  0.052900
      Read-in Center 5192 is at   0.000000  2.010900  0.052900
      Read-in Center 5193 is at   0.000000  2.063800  0.052900
      Read-in Center 5194 is at   0.000000  2.116700  0.052900
      Read-in Center 5195 is at   0.000000  2.169600  0.052900
      Read-in Center 5196 is at   0.000000  2.222500  0.052900
      Read-in Center 5197 is at   0.000000  2.275500  0.052900
      Read-in Center 5198 is at   0.000000  2.328400  0.052900
      Read-in Center 5199 is at   0.000000  2.381300  0.052900
      Read-in Center 5200 is at   0.000000  2.434200  0.052900
      Read-in Center 5201 is at   0.000000  2.487100  0.052900
      Read-in Center 5202 is at   0.000000  2.540000  0.052900
      Read-in Center 5203 is at   0.000000  2.593000  0.052900
      Read-in Center 5204 is at   0.000000 -2.645900  0.105800
      Read-in Center 5205 is at   0.000000 -2.593000  0.105800
      Read-in Center 5206 is at   0.000000 -2.540000  0.105800
      Read-in Center 5207 is at   0.000000 -2.487100  0.105800
      Read-in Center 5208 is at   0.000000 -2.434200  0.105800
      Read-in Center 5209 is at   0.000000 -2.381300  0.105800
      Read-in Center 5210 is at   0.000000 -2.328400  0.105800
      Read-in Center 5211 is at   0.000000 -2.275500  0.105800
      Read-in Center 5212 is at   0.000000 -2.222500  0.105800
      Read-in Center 5213 is at   0.000000 -2.169600  0.105800
      Read-in Center 5214 is at   0.000000 -2.116700  0.105800
      Read-in Center 5215 is at   0.000000 -2.063800  0.105800
      Read-in Center 5216 is at   0.000000 -2.010900  0.105800
      Read-in Center 5217 is at   0.000000 -1.958000  0.105800
      Read-in Center 5218 is at   0.000000 -1.905000  0.105800
      Read-in Center 5219 is at   0.000000 -1.852100  0.105800
      Read-in Center 5220 is at   0.000000 -1.799200  0.105800
      Read-in Center 5221 is at   0.000000 -1.746300  0.105800
      Read-in Center 5222 is at   0.000000 -1.693400  0.105800
      Read-in Center 5223 is at   0.000000 -1.640400  0.105800
      Read-in Center 5224 is at   0.000000 -1.587500  0.105800
      Read-in Center 5225 is at   0.000000 -1.534600  0.105800
      Read-in Center 5226 is at   0.000000 -1.481700  0.105800
      Read-in Center 5227 is at   0.000000 -1.428800  0.105800
      Read-in Center 5228 is at   0.000000 -1.375900  0.105800
      Read-in Center 5229 is at   0.000000 -1.322900  0.105800
      Read-in Center 5230 is at   0.000000 -1.270000  0.105800
      Read-in Center 5231 is at   0.000000 -1.217100  0.105800
      Read-in Center 5232 is at   0.000000 -1.164200  0.105800
      Read-in Center 5233 is at   0.000000 -1.111300  0.105800
      Read-in Center 5234 is at   0.000000 -1.058400  0.105800
      Read-in Center 5235 is at   0.000000 -1.005400  0.105800
      Read-in Center 5236 is at   0.000000 -0.952500  0.105800
      Read-in Center 5237 is at   0.000000 -0.899600  0.105800
      Read-in Center 5238 is at   0.000000 -0.846700  0.105800
      Read-in Center 5239 is at   0.000000 -0.793800  0.105800
      Read-in Center 5240 is at   0.000000 -0.740800  0.105800
      Read-in Center 5241 is at   0.000000 -0.687900  0.105800
      Read-in Center 5242 is at   0.000000 -0.635000  0.105800
      Read-in Center 5243 is at   0.000000 -0.582100  0.105800
      Read-in Center 5244 is at   0.000000 -0.529200  0.105800
      Read-in Center 5245 is at   0.000000 -0.476300  0.105800
      Read-in Center 5246 is at   0.000000 -0.423300  0.105800
      Read-in Center 5247 is at   0.000000 -0.370400  0.105800
      Read-in Center 5248 is at   0.000000 -0.317500  0.105800
      Read-in Center 5249 is at   0.000000 -0.264600  0.105800
      Read-in Center 5250 is at   0.000000 -0.211700  0.105800
      Read-in Center 5251 is at   0.000000 -0.158800  0.105800
      Read-in Center 5252 is at   0.000000 -0.105800  0.105800
      Read-in Center 5253 is at   0.000000 -0.052900  0.105800
      Read-in Center 5254 is at   0.000000  0.000000  0.105800
      Read-in Center 5255 is at   0.000000  0.052900  0.105800
      Read-in Center 5256 is at   0.000000  0.105800  0.105800
      Read-in Center 5257 is at   0.000000  0.158800  0.105800
      Read-in Center 5258 is at   0.000000  0.211700  0.105800
      Read-in Center 5259 is at   0.000000  0.264600  0.105800
      Read-in Center 5260 is at   0.000000  0.317500  0.105800
      Read-in Center 5261 is at   0.000000  0.370400  0.105800
      Read-in Center 5262 is at   0.000000  0.423300  0.105800
      Read-in Center 5263 is at   0.000000  0.476300  0.105800
      Read-in Center 5264 is at   0.000000  0.529200  0.105800
      Read-in Center 5265 is at   0.000000  0.582100  0.105800
      Read-in Center 5266 is at   0.000000  0.635000  0.105800
      Read-in Center 5267 is at   0.000000  0.687900  0.105800
      Read-in Center 5268 is at   0.000000  0.740800  0.105800
      Read-in Center 5269 is at   0.000000  0.793800  0.105800
      Read-in Center 5270 is at   0.000000  0.846700  0.105800
      Read-in Center 5271 is at   0.000000  0.899600  0.105800
      Read-in Center 5272 is at   0.000000  0.952500  0.105800
      Read-in Center 5273 is at   0.000000  1.005400  0.105800
      Read-in Center 5274 is at   0.000000  1.058400  0.105800
      Read-in Center 5275 is at   0.000000  1.111300  0.105800
      Read-in Center 5276 is at   0.000000  1.164200  0.105800
      Read-in Center 5277 is at   0.000000  1.217100  0.105800
      Read-in Center 5278 is at   0.000000  1.270000  0.105800
      Read-in Center 5279 is at   0.000000  1.322900  0.105800
      Read-in Center 5280 is at   0.000000  1.375900  0.105800
      Read-in Center 5281 is at   0.000000  1.428800  0.105800
      Read-in Center 5282 is at   0.000000  1.481700  0.105800
      Read-in Center 5283 is at   0.000000  1.534600  0.105800
      Read-in Center 5284 is at   0.000000  1.587500  0.105800
      Read-in Center 5285 is at   0.000000  1.640400  0.105800
      Read-in Center 5286 is at   0.000000  1.693400  0.105800
      Read-in Center 5287 is at   0.000000  1.746300  0.105800
      Read-in Center 5288 is at   0.000000  1.799200  0.105800
      Read-in Center 5289 is at   0.000000  1.852100  0.105800
      Read-in Center 5290 is at   0.000000  1.905000  0.105800
      Read-in Center 5291 is at   0.000000  1.958000  0.105800
      Read-in Center 5292 is at   0.000000  2.010900  0.105800
      Read-in Center 5293 is at   0.000000  2.063800  0.105800
      Read-in Center 5294 is at   0.000000  2.116700  0.105800
      Read-in Center 5295 is at   0.000000  2.169600  0.105800
      Read-in Center 5296 is at   0.000000  2.222500  0.105800
      Read-in Center 5297 is at   0.000000  2.275500  0.105800
      Read-in Center 5298 is at   0.000000  2.328400  0.105800
      Read-in Center 5299 is at   0.000000  2.381300  0.105800
      Read-in Center 5300 is at   0.000000  2.434200  0.105800
      Read-in Center 5301 is at   0.000000  2.487100  0.105800
      Read-in Center 5302 is at   0.000000  2.540000  0.105800
      Read-in Center 5303 is at   0.000000  2.593000  0.105800
      Read-in Center 5304 is at   0.000000 -2.645900  0.158800
      Read-in Center 5305 is at   0.000000 -2.593000  0.158800
      Read-in Center 5306 is at   0.000000 -2.540000  0.158800
      Read-in Center 5307 is at   0.000000 -2.487100  0.158800
      Read-in Center 5308 is at   0.000000 -2.434200  0.158800
      Read-in Center 5309 is at   0.000000 -2.381300  0.158800
      Read-in Center 5310 is at   0.000000 -2.328400  0.158800
      Read-in Center 5311 is at   0.000000 -2.275500  0.158800
      Read-in Center 5312 is at   0.000000 -2.222500  0.158800
      Read-in Center 5313 is at   0.000000 -2.169600  0.158800
      Read-in Center 5314 is at   0.000000 -2.116700  0.158800
      Read-in Center 5315 is at   0.000000 -2.063800  0.158800
      Read-in Center 5316 is at   0.000000 -2.010900  0.158800
      Read-in Center 5317 is at   0.000000 -1.958000  0.158800
      Read-in Center 5318 is at   0.000000 -1.905000  0.158800
      Read-in Center 5319 is at   0.000000 -1.852100  0.158800
      Read-in Center 5320 is at   0.000000 -1.799200  0.158800
      Read-in Center 5321 is at   0.000000 -1.746300  0.158800
      Read-in Center 5322 is at   0.000000 -1.693400  0.158800
      Read-in Center 5323 is at   0.000000 -1.640400  0.158800
      Read-in Center 5324 is at   0.000000 -1.587500  0.158800
      Read-in Center 5325 is at   0.000000 -1.534600  0.158800
      Read-in Center 5326 is at   0.000000 -1.481700  0.158800
      Read-in Center 5327 is at   0.000000 -1.428800  0.158800
      Read-in Center 5328 is at   0.000000 -1.375900  0.158800
      Read-in Center 5329 is at   0.000000 -1.322900  0.158800
      Read-in Center 5330 is at   0.000000 -1.270000  0.158800
      Read-in Center 5331 is at   0.000000 -1.217100  0.158800
      Read-in Center 5332 is at   0.000000 -1.164200  0.158800
      Read-in Center 5333 is at   0.000000 -1.111300  0.158800
      Read-in Center 5334 is at   0.000000 -1.058400  0.158800
      Read-in Center 5335 is at   0.000000 -1.005400  0.158800
      Read-in Center 5336 is at   0.000000 -0.952500  0.158800
      Read-in Center 5337 is at   0.000000 -0.899600  0.158800
      Read-in Center 5338 is at   0.000000 -0.846700  0.158800
      Read-in Center 5339 is at   0.000000 -0.793800  0.158800
      Read-in Center 5340 is at   0.000000 -0.740800  0.158800
      Read-in Center 5341 is at   0.000000 -0.687900  0.158800
      Read-in Center 5342 is at   0.000000 -0.635000  0.158800
      Read-in Center 5343 is at   0.000000 -0.582100  0.158800
      Read-in Center 5344 is at   0.000000 -0.529200  0.158800
      Read-in Center 5345 is at   0.000000 -0.476300  0.158800
      Read-in Center 5346 is at   0.000000 -0.423300  0.158800
      Read-in Center 5347 is at   0.000000 -0.370400  0.158800
      Read-in Center 5348 is at   0.000000 -0.317500  0.158800
      Read-in Center 5349 is at   0.000000 -0.264600  0.158800
      Read-in Center 5350 is at   0.000000 -0.211700  0.158800
      Read-in Center 5351 is at   0.000000 -0.158800  0.158800
      Read-in Center 5352 is at   0.000000 -0.105800  0.158800
      Read-in Center 5353 is at   0.000000 -0.052900  0.158800
      Read-in Center 5354 is at   0.000000  0.000000  0.158800
      Read-in Center 5355 is at   0.000000  0.052900  0.158800
      Read-in Center 5356 is at   0.000000  0.105800  0.158800
      Read-in Center 5357 is at   0.000000  0.158800  0.158800
      Read-in Center 5358 is at   0.000000  0.211700  0.158800
      Read-in Center 5359 is at   0.000000  0.264600  0.158800
      Read-in Center 5360 is at   0.000000  0.317500  0.158800
      Read-in Center 5361 is at   0.000000  0.370400  0.158800
      Read-in Center 5362 is at   0.000000  0.423300  0.158800
      Read-in Center 5363 is at   0.000000  0.476300  0.158800
      Read-in Center 5364 is at   0.000000  0.529200  0.158800
      Read-in Center 5365 is at   0.000000  0.582100  0.158800
      Read-in Center 5366 is at   0.000000  0.635000  0.158800
      Read-in Center 5367 is at   0.000000  0.687900  0.158800
      Read-in Center 5368 is at   0.000000  0.740800  0.158800
      Read-in Center 5369 is at   0.000000  0.793800  0.158800
      Read-in Center 5370 is at   0.000000  0.846700  0.158800
      Read-in Center 5371 is at   0.000000  0.899600  0.158800
      Read-in Center 5372 is at   0.000000  0.952500  0.158800
      Read-in Center 5373 is at   0.000000  1.005400  0.158800
      Read-in Center 5374 is at   0.000000  1.058400  0.158800
      Read-in Center 5375 is at   0.000000  1.111300  0.158800
      Read-in Center 5376 is at   0.000000  1.164200  0.158800
      Read-in Center 5377 is at   0.000000  1.217100  0.158800
      Read-in Center 5378 is at   0.000000  1.270000  0.158800
      Read-in Center 5379 is at   0.000000  1.322900  0.158800
      Read-in Center 5380 is at   0.000000  1.375900  0.158800
      Read-in Center 5381 is at   0.000000  1.428800  0.158800
      Read-in Center 5382 is at   0.000000  1.481700  0.158800
      Read-in Center 5383 is at   0.000000  1.534600  0.158800
      Read-in Center 5384 is at   0.000000  1.587500  0.158800
      Read-in Center 5385 is at   0.000000  1.640400  0.158800
      Read-in Center 5386 is at   0.000000  1.693400  0.158800
      Read-in Center 5387 is at   0.000000  1.746300  0.158800
      Read-in Center 5388 is at   0.000000  1.799200  0.158800
      Read-in Center 5389 is at   0.000000  1.852100  0.158800
      Read-in Center 5390 is at   0.000000  1.905000  0.158800
      Read-in Center 5391 is at   0.000000  1.958000  0.158800
      Read-in Center 5392 is at   0.000000  2.010900  0.158800
      Read-in Center 5393 is at   0.000000  2.063800  0.158800
      Read-in Center 5394 is at   0.000000  2.116700  0.158800
      Read-in Center 5395 is at   0.000000  2.169600  0.158800
      Read-in Center 5396 is at   0.000000  2.222500  0.158800
      Read-in Center 5397 is at   0.000000  2.275500  0.158800
      Read-in Center 5398 is at   0.000000  2.328400  0.158800
      Read-in Center 5399 is at   0.000000  2.381300  0.158800
      Read-in Center 5400 is at   0.000000  2.434200  0.158800
      Read-in Center 5401 is at   0.000000  2.487100  0.158800
      Read-in Center 5402 is at   0.000000  2.540000  0.158800
      Read-in Center 5403 is at   0.000000  2.593000  0.158800
      Read-in Center 5404 is at   0.000000 -2.645900  0.211700
      Read-in Center 5405 is at   0.000000 -2.593000  0.211700
      Read-in Center 5406 is at   0.000000 -2.540000  0.211700
      Read-in Center 5407 is at   0.000000 -2.487100  0.211700
      Read-in Center 5408 is at   0.000000 -2.434200  0.211700
      Read-in Center 5409 is at   0.000000 -2.381300  0.211700
      Read-in Center 5410 is at   0.000000 -2.328400  0.211700
      Read-in Center 5411 is at   0.000000 -2.275500  0.211700
      Read-in Center 5412 is at   0.000000 -2.222500  0.211700
      Read-in Center 5413 is at   0.000000 -2.169600  0.211700
      Read-in Center 5414 is at   0.000000 -2.116700  0.211700
      Read-in Center 5415 is at   0.000000 -2.063800  0.211700
      Read-in Center 5416 is at   0.000000 -2.010900  0.211700
      Read-in Center 5417 is at   0.000000 -1.958000  0.211700
      Read-in Center 5418 is at   0.000000 -1.905000  0.211700
      Read-in Center 5419 is at   0.000000 -1.852100  0.211700
      Read-in Center 5420 is at   0.000000 -1.799200  0.211700
      Read-in Center 5421 is at   0.000000 -1.746300  0.211700
      Read-in Center 5422 is at   0.000000 -1.693400  0.211700
      Read-in Center 5423 is at   0.000000 -1.640400  0.211700
      Read-in Center 5424 is at   0.000000 -1.587500  0.211700
      Read-in Center 5425 is at   0.000000 -1.534600  0.211700
      Read-in Center 5426 is at   0.000000 -1.481700  0.211700
      Read-in Center 5427 is at   0.000000 -1.428800  0.211700
      Read-in Center 5428 is at   0.000000 -1.375900  0.211700
      Read-in Center 5429 is at   0.000000 -1.322900  0.211700
      Read-in Center 5430 is at   0.000000 -1.270000  0.211700
      Read-in Center 5431 is at   0.000000 -1.217100  0.211700
      Read-in Center 5432 is at   0.000000 -1.164200  0.211700
      Read-in Center 5433 is at   0.000000 -1.111300  0.211700
      Read-in Center 5434 is at   0.000000 -1.058400  0.211700
      Read-in Center 5435 is at   0.000000 -1.005400  0.211700
      Read-in Center 5436 is at   0.000000 -0.952500  0.211700
      Read-in Center 5437 is at   0.000000 -0.899600  0.211700
      Read-in Center 5438 is at   0.000000 -0.846700  0.211700
      Read-in Center 5439 is at   0.000000 -0.793800  0.211700
      Read-in Center 5440 is at   0.000000 -0.740800  0.211700
      Read-in Center 5441 is at   0.000000 -0.687900  0.211700
      Read-in Center 5442 is at   0.000000 -0.635000  0.211700
      Read-in Center 5443 is at   0.000000 -0.582100  0.211700
      Read-in Center 5444 is at   0.000000 -0.529200  0.211700
      Read-in Center 5445 is at   0.000000 -0.476300  0.211700
      Read-in Center 5446 is at   0.000000 -0.423300  0.211700
      Read-in Center 5447 is at   0.000000 -0.370400  0.211700
      Read-in Center 5448 is at   0.000000 -0.317500  0.211700
      Read-in Center 5449 is at   0.000000 -0.264600  0.211700
      Read-in Center 5450 is at   0.000000 -0.211700  0.211700
      Read-in Center 5451 is at   0.000000 -0.158800  0.211700
      Read-in Center 5452 is at   0.000000 -0.105800  0.211700
      Read-in Center 5453 is at   0.000000 -0.052900  0.211700
      Read-in Center 5454 is at   0.000000  0.000000  0.211700
      Read-in Center 5455 is at   0.000000  0.052900  0.211700
      Read-in Center 5456 is at   0.000000  0.105800  0.211700
      Read-in Center 5457 is at   0.000000  0.158800  0.211700
      Read-in Center 5458 is at   0.000000  0.211700  0.211700
      Read-in Center 5459 is at   0.000000  0.264600  0.211700
      Read-in Center 5460 is at   0.000000  0.317500  0.211700
      Read-in Center 5461 is at   0.000000  0.370400  0.211700
      Read-in Center 5462 is at   0.000000  0.423300  0.211700
      Read-in Center 5463 is at   0.000000  0.476300  0.211700
      Read-in Center 5464 is at   0.000000  0.529200  0.211700
      Read-in Center 5465 is at   0.000000  0.582100  0.211700
      Read-in Center 5466 is at   0.000000  0.635000  0.211700
      Read-in Center 5467 is at   0.000000  0.687900  0.211700
      Read-in Center 5468 is at   0.000000  0.740800  0.211700
      Read-in Center 5469 is at   0.000000  0.793800  0.211700
      Read-in Center 5470 is at   0.000000  0.846700  0.211700
      Read-in Center 5471 is at   0.000000  0.899600  0.211700
      Read-in Center 5472 is at   0.000000  0.952500  0.211700
      Read-in Center 5473 is at   0.000000  1.005400  0.211700
      Read-in Center 5474 is at   0.000000  1.058400  0.211700
      Read-in Center 5475 is at   0.000000  1.111300  0.211700
      Read-in Center 5476 is at   0.000000  1.164200  0.211700
      Read-in Center 5477 is at   0.000000  1.217100  0.211700
      Read-in Center 5478 is at   0.000000  1.270000  0.211700
      Read-in Center 5479 is at   0.000000  1.322900  0.211700
      Read-in Center 5480 is at   0.000000  1.375900  0.211700
      Read-in Center 5481 is at   0.000000  1.428800  0.211700
      Read-in Center 5482 is at   0.000000  1.481700  0.211700
      Read-in Center 5483 is at   0.000000  1.534600  0.211700
      Read-in Center 5484 is at   0.000000  1.587500  0.211700
      Read-in Center 5485 is at   0.000000  1.640400  0.211700
      Read-in Center 5486 is at   0.000000  1.693400  0.211700
      Read-in Center 5487 is at   0.000000  1.746300  0.211700
      Read-in Center 5488 is at   0.000000  1.799200  0.211700
      Read-in Center 5489 is at   0.000000  1.852100  0.211700
      Read-in Center 5490 is at   0.000000  1.905000  0.211700
      Read-in Center 5491 is at   0.000000  1.958000  0.211700
      Read-in Center 5492 is at   0.000000  2.010900  0.211700
      Read-in Center 5493 is at   0.000000  2.063800  0.211700
      Read-in Center 5494 is at   0.000000  2.116700  0.211700
      Read-in Center 5495 is at   0.000000  2.169600  0.211700
      Read-in Center 5496 is at   0.000000  2.222500  0.211700
      Read-in Center 5497 is at   0.000000  2.275500  0.211700
      Read-in Center 5498 is at   0.000000  2.328400  0.211700
      Read-in Center 5499 is at   0.000000  2.381300  0.211700
      Read-in Center 5500 is at   0.000000  2.434200  0.211700
      Read-in Center 5501 is at   0.000000  2.487100  0.211700
      Read-in Center 5502 is at   0.000000  2.540000  0.211700
      Read-in Center 5503 is at   0.000000  2.593000  0.211700
      Read-in Center 5504 is at   0.000000 -2.645900  0.264600
      Read-in Center 5505 is at   0.000000 -2.593000  0.264600
      Read-in Center 5506 is at   0.000000 -2.540000  0.264600
      Read-in Center 5507 is at   0.000000 -2.487100  0.264600
      Read-in Center 5508 is at   0.000000 -2.434200  0.264600
      Read-in Center 5509 is at   0.000000 -2.381300  0.264600
      Read-in Center 5510 is at   0.000000 -2.328400  0.264600
      Read-in Center 5511 is at   0.000000 -2.275500  0.264600
      Read-in Center 5512 is at   0.000000 -2.222500  0.264600
      Read-in Center 5513 is at   0.000000 -2.169600  0.264600
      Read-in Center 5514 is at   0.000000 -2.116700  0.264600
      Read-in Center 5515 is at   0.000000 -2.063800  0.264600
      Read-in Center 5516 is at   0.000000 -2.010900  0.264600
      Read-in Center 5517 is at   0.000000 -1.958000  0.264600
      Read-in Center 5518 is at   0.000000 -1.905000  0.264600
      Read-in Center 5519 is at   0.000000 -1.852100  0.264600
      Read-in Center 5520 is at   0.000000 -1.799200  0.264600
      Read-in Center 5521 is at   0.000000 -1.746300  0.264600
      Read-in Center 5522 is at   0.000000 -1.693400  0.264600
      Read-in Center 5523 is at   0.000000 -1.640400  0.264600
      Read-in Center 5524 is at   0.000000 -1.587500  0.264600
      Read-in Center 5525 is at   0.000000 -1.534600  0.264600
      Read-in Center 5526 is at   0.000000 -1.481700  0.264600
      Read-in Center 5527 is at   0.000000 -1.428800  0.264600
      Read-in Center 5528 is at   0.000000 -1.375900  0.264600
      Read-in Center 5529 is at   0.000000 -1.322900  0.264600
      Read-in Center 5530 is at   0.000000 -1.270000  0.264600
      Read-in Center 5531 is at   0.000000 -1.217100  0.264600
      Read-in Center 5532 is at   0.000000 -1.164200  0.264600
      Read-in Center 5533 is at   0.000000 -1.111300  0.264600
      Read-in Center 5534 is at   0.000000 -1.058400  0.264600
      Read-in Center 5535 is at   0.000000 -1.005400  0.264600
      Read-in Center 5536 is at   0.000000 -0.952500  0.264600
      Read-in Center 5537 is at   0.000000 -0.899600  0.264600
      Read-in Center 5538 is at   0.000000 -0.846700  0.264600
      Read-in Center 5539 is at   0.000000 -0.793800  0.264600
      Read-in Center 5540 is at   0.000000 -0.740800  0.264600
      Read-in Center 5541 is at   0.000000 -0.687900  0.264600
      Read-in Center 5542 is at   0.000000 -0.635000  0.264600
      Read-in Center 5543 is at   0.000000 -0.582100  0.264600
      Read-in Center 5544 is at   0.000000 -0.529200  0.264600
      Read-in Center 5545 is at   0.000000 -0.476300  0.264600
      Read-in Center 5546 is at   0.000000 -0.423300  0.264600
      Read-in Center 5547 is at   0.000000 -0.370400  0.264600
      Read-in Center 5548 is at   0.000000 -0.317500  0.264600
      Read-in Center 5549 is at   0.000000 -0.264600  0.264600
      Read-in Center 5550 is at   0.000000 -0.211700  0.264600
      Read-in Center 5551 is at   0.000000 -0.158800  0.264600
      Read-in Center 5552 is at   0.000000 -0.105800  0.264600
      Read-in Center 5553 is at   0.000000 -0.052900  0.264600
      Read-in Center 5554 is at   0.000000  0.000000  0.264600
      Read-in Center 5555 is at   0.000000  0.052900  0.264600
      Read-in Center 5556 is at   0.000000  0.105800  0.264600
      Read-in Center 5557 is at   0.000000  0.158800  0.264600
      Read-in Center 5558 is at   0.000000  0.211700  0.264600
      Read-in Center 5559 is at   0.000000  0.264600  0.264600
      Read-in Center 5560 is at   0.000000  0.317500  0.264600
      Read-in Center 5561 is at   0.000000  0.370400  0.264600
      Read-in Center 5562 is at   0.000000  0.423300  0.264600
      Read-in Center 5563 is at   0.000000  0.476300  0.264600
      Read-in Center 5564 is at   0.000000  0.529200  0.264600
      Read-in Center 5565 is at   0.000000  0.582100  0.264600
      Read-in Center 5566 is at   0.000000  0.635000  0.264600
      Read-in Center 5567 is at   0.000000  0.687900  0.264600
      Read-in Center 5568 is at   0.000000  0.740800  0.264600
      Read-in Center 5569 is at   0.000000  0.793800  0.264600
      Read-in Center 5570 is at   0.000000  0.846700  0.264600
      Read-in Center 5571 is at   0.000000  0.899600  0.264600
      Read-in Center 5572 is at   0.000000  0.952500  0.264600
      Read-in Center 5573 is at   0.000000  1.005400  0.264600
      Read-in Center 5574 is at   0.000000  1.058400  0.264600
      Read-in Center 5575 is at   0.000000  1.111300  0.264600
      Read-in Center 5576 is at   0.000000  1.164200  0.264600
      Read-in Center 5577 is at   0.000000  1.217100  0.264600
      Read-in Center 5578 is at   0.000000  1.270000  0.264600
      Read-in Center 5579 is at   0.000000  1.322900  0.264600
      Read-in Center 5580 is at   0.000000  1.375900  0.264600
      Read-in Center 5581 is at   0.000000  1.428800  0.264600
      Read-in Center 5582 is at   0.000000  1.481700  0.264600
      Read-in Center 5583 is at   0.000000  1.534600  0.264600
      Read-in Center 5584 is at   0.000000  1.587500  0.264600
      Read-in Center 5585 is at   0.000000  1.640400  0.264600
      Read-in Center 5586 is at   0.000000  1.693400  0.264600
      Read-in Center 5587 is at   0.000000  1.746300  0.264600
      Read-in Center 5588 is at   0.000000  1.799200  0.264600
      Read-in Center 5589 is at   0.000000  1.852100  0.264600
      Read-in Center 5590 is at   0.000000  1.905000  0.264600
      Read-in Center 5591 is at   0.000000  1.958000  0.264600
      Read-in Center 5592 is at   0.000000  2.010900  0.264600
      Read-in Center 5593 is at   0.000000  2.063800  0.264600
      Read-in Center 5594 is at   0.000000  2.116700  0.264600
      Read-in Center 5595 is at   0.000000  2.169600  0.264600
      Read-in Center 5596 is at   0.000000  2.222500  0.264600
      Read-in Center 5597 is at   0.000000  2.275500  0.264600
      Read-in Center 5598 is at   0.000000  2.328400  0.264600
      Read-in Center 5599 is at   0.000000  2.381300  0.264600
      Read-in Center 5600 is at   0.000000  2.434200  0.264600
      Read-in Center 5601 is at   0.000000  2.487100  0.264600
      Read-in Center 5602 is at   0.000000  2.540000  0.264600
      Read-in Center 5603 is at   0.000000  2.593000  0.264600
      Read-in Center 5604 is at   0.000000 -2.645900  0.317500
      Read-in Center 5605 is at   0.000000 -2.593000  0.317500
      Read-in Center 5606 is at   0.000000 -2.540000  0.317500
      Read-in Center 5607 is at   0.000000 -2.487100  0.317500
      Read-in Center 5608 is at   0.000000 -2.434200  0.317500
      Read-in Center 5609 is at   0.000000 -2.381300  0.317500
      Read-in Center 5610 is at   0.000000 -2.328400  0.317500
      Read-in Center 5611 is at   0.000000 -2.275500  0.317500
      Read-in Center 5612 is at   0.000000 -2.222500  0.317500
      Read-in Center 5613 is at   0.000000 -2.169600  0.317500
      Read-in Center 5614 is at   0.000000 -2.116700  0.317500
      Read-in Center 5615 is at   0.000000 -2.063800  0.317500
      Read-in Center 5616 is at   0.000000 -2.010900  0.317500
      Read-in Center 5617 is at   0.000000 -1.958000  0.317500
      Read-in Center 5618 is at   0.000000 -1.905000  0.317500
      Read-in Center 5619 is at   0.000000 -1.852100  0.317500
      Read-in Center 5620 is at   0.000000 -1.799200  0.317500
      Read-in Center 5621 is at   0.000000 -1.746300  0.317500
      Read-in Center 5622 is at   0.000000 -1.693400  0.317500
      Read-in Center 5623 is at   0.000000 -1.640400  0.317500
      Read-in Center 5624 is at   0.000000 -1.587500  0.317500
      Read-in Center 5625 is at   0.000000 -1.534600  0.317500
      Read-in Center 5626 is at   0.000000 -1.481700  0.317500
      Read-in Center 5627 is at   0.000000 -1.428800  0.317500
      Read-in Center 5628 is at   0.000000 -1.375900  0.317500
      Read-in Center 5629 is at   0.000000 -1.322900  0.317500
      Read-in Center 5630 is at   0.000000 -1.270000  0.317500
      Read-in Center 5631 is at   0.000000 -1.217100  0.317500
      Read-in Center 5632 is at   0.000000 -1.164200  0.317500
      Read-in Center 5633 is at   0.000000 -1.111300  0.317500
      Read-in Center 5634 is at   0.000000 -1.058400  0.317500
      Read-in Center 5635 is at   0.000000 -1.005400  0.317500
      Read-in Center 5636 is at   0.000000 -0.952500  0.317500
      Read-in Center 5637 is at   0.000000 -0.899600  0.317500
      Read-in Center 5638 is at   0.000000 -0.846700  0.317500
      Read-in Center 5639 is at   0.000000 -0.793800  0.317500
      Read-in Center 5640 is at   0.000000 -0.740800  0.317500
      Read-in Center 5641 is at   0.000000 -0.687900  0.317500
      Read-in Center 5642 is at   0.000000 -0.635000  0.317500
      Read-in Center 5643 is at   0.000000 -0.582100  0.317500
      Read-in Center 5644 is at   0.000000 -0.529200  0.317500
      Read-in Center 5645 is at   0.000000 -0.476300  0.317500
      Read-in Center 5646 is at   0.000000 -0.423300  0.317500
      Read-in Center 5647 is at   0.000000 -0.370400  0.317500
      Read-in Center 5648 is at   0.000000 -0.317500  0.317500
      Read-in Center 5649 is at   0.000000 -0.264600  0.317500
      Read-in Center 5650 is at   0.000000 -0.211700  0.317500
      Read-in Center 5651 is at   0.000000 -0.158800  0.317500
      Read-in Center 5652 is at   0.000000 -0.105800  0.317500
      Read-in Center 5653 is at   0.000000 -0.052900  0.317500
      Read-in Center 5654 is at   0.000000  0.000000  0.317500
      Read-in Center 5655 is at   0.000000  0.052900  0.317500
      Read-in Center 5656 is at   0.000000  0.105800  0.317500
      Read-in Center 5657 is at   0.000000  0.158800  0.317500
      Read-in Center 5658 is at   0.000000  0.211700  0.317500
      Read-in Center 5659 is at   0.000000  0.264600  0.317500
      Read-in Center 5660 is at   0.000000  0.317500  0.317500
      Read-in Center 5661 is at   0.000000  0.370400  0.317500
      Read-in Center 5662 is at   0.000000  0.423300  0.317500
      Read-in Center 5663 is at   0.000000  0.476300  0.317500
      Read-in Center 5664 is at   0.000000  0.529200  0.317500
      Read-in Center 5665 is at   0.000000  0.582100  0.317500
      Read-in Center 5666 is at   0.000000  0.635000  0.317500
      Read-in Center 5667 is at   0.000000  0.687900  0.317500
      Read-in Center 5668 is at   0.000000  0.740800  0.317500
      Read-in Center 5669 is at   0.000000  0.793800  0.317500
      Read-in Center 5670 is at   0.000000  0.846700  0.317500
      Read-in Center 5671 is at   0.000000  0.899600  0.317500
      Read-in Center 5672 is at   0.000000  0.952500  0.317500
      Read-in Center 5673 is at   0.000000  1.005400  0.317500
      Read-in Center 5674 is at   0.000000  1.058400  0.317500
      Read-in Center 5675 is at   0.000000  1.111300  0.317500
      Read-in Center 5676 is at   0.000000  1.164200  0.317500
      Read-in Center 5677 is at   0.000000  1.217100  0.317500
      Read-in Center 5678 is at   0.000000  1.270000  0.317500
      Read-in Center 5679 is at   0.000000  1.322900  0.317500
      Read-in Center 5680 is at   0.000000  1.375900  0.317500
      Read-in Center 5681 is at   0.000000  1.428800  0.317500
      Read-in Center 5682 is at   0.000000  1.481700  0.317500
      Read-in Center 5683 is at   0.000000  1.534600  0.317500
      Read-in Center 5684 is at   0.000000  1.587500  0.317500
      Read-in Center 5685 is at   0.000000  1.640400  0.317500
      Read-in Center 5686 is at   0.000000  1.693400  0.317500
      Read-in Center 5687 is at   0.000000  1.746300  0.317500
      Read-in Center 5688 is at   0.000000  1.799200  0.317500
      Read-in Center 5689 is at   0.000000  1.852100  0.317500
      Read-in Center 5690 is at   0.000000  1.905000  0.317500
      Read-in Center 5691 is at   0.000000  1.958000  0.317500
      Read-in Center 5692 is at   0.000000  2.010900  0.317500
      Read-in Center 5693 is at   0.000000  2.063800  0.317500
      Read-in Center 5694 is at   0.000000  2.116700  0.317500
      Read-in Center 5695 is at   0.000000  2.169600  0.317500
      Read-in Center 5696 is at   0.000000  2.222500  0.317500
      Read-in Center 5697 is at   0.000000  2.275500  0.317500
      Read-in Center 5698 is at   0.000000  2.328400  0.317500
      Read-in Center 5699 is at   0.000000  2.381300  0.317500
      Read-in Center 5700 is at   0.000000  2.434200  0.317500
      Read-in Center 5701 is at   0.000000  2.487100  0.317500
      Read-in Center 5702 is at   0.000000  2.540000  0.317500
      Read-in Center 5703 is at   0.000000  2.593000  0.317500
      Read-in Center 5704 is at   0.000000 -2.645900  0.370400
      Read-in Center 5705 is at   0.000000 -2.593000  0.370400
      Read-in Center 5706 is at   0.000000 -2.540000  0.370400
      Read-in Center 5707 is at   0.000000 -2.487100  0.370400
      Read-in Center 5708 is at   0.000000 -2.434200  0.370400
      Read-in Center 5709 is at   0.000000 -2.381300  0.370400
      Read-in Center 5710 is at   0.000000 -2.328400  0.370400
      Read-in Center 5711 is at   0.000000 -2.275500  0.370400
      Read-in Center 5712 is at   0.000000 -2.222500  0.370400
      Read-in Center 5713 is at   0.000000 -2.169600  0.370400
      Read-in Center 5714 is at   0.000000 -2.116700  0.370400
      Read-in Center 5715 is at   0.000000 -2.063800  0.370400
      Read-in Center 5716 is at   0.000000 -2.010900  0.370400
      Read-in Center 5717 is at   0.000000 -1.958000  0.370400
      Read-in Center 5718 is at   0.000000 -1.905000  0.370400
      Read-in Center 5719 is at   0.000000 -1.852100  0.370400
      Read-in Center 5720 is at   0.000000 -1.799200  0.370400
      Read-in Center 5721 is at   0.000000 -1.746300  0.370400
      Read-in Center 5722 is at   0.000000 -1.693400  0.370400
      Read-in Center 5723 is at   0.000000 -1.640400  0.370400
      Read-in Center 5724 is at   0.000000 -1.587500  0.370400
      Read-in Center 5725 is at   0.000000 -1.534600  0.370400
      Read-in Center 5726 is at   0.000000 -1.481700  0.370400
      Read-in Center 5727 is at   0.000000 -1.428800  0.370400
      Read-in Center 5728 is at   0.000000 -1.375900  0.370400
      Read-in Center 5729 is at   0.000000 -1.322900  0.370400
      Read-in Center 5730 is at   0.000000 -1.270000  0.370400
      Read-in Center 5731 is at   0.000000 -1.217100  0.370400
      Read-in Center 5732 is at   0.000000 -1.164200  0.370400
      Read-in Center 5733 is at   0.000000 -1.111300  0.370400
      Read-in Center 5734 is at   0.000000 -1.058400  0.370400
      Read-in Center 5735 is at   0.000000 -1.005400  0.370400
      Read-in Center 5736 is at   0.000000 -0.952500  0.370400
      Read-in Center 5737 is at   0.000000 -0.899600  0.370400
      Read-in Center 5738 is at   0.000000 -0.846700  0.370400
      Read-in Center 5739 is at   0.000000 -0.793800  0.370400
      Read-in Center 5740 is at   0.000000 -0.740800  0.370400
      Read-in Center 5741 is at   0.000000 -0.687900  0.370400
      Read-in Center 5742 is at   0.000000 -0.635000  0.370400
      Read-in Center 5743 is at   0.000000 -0.582100  0.370400
      Read-in Center 5744 is at   0.000000 -0.529200  0.370400
      Read-in Center 5745 is at   0.000000 -0.476300  0.370400
      Read-in Center 5746 is at   0.000000 -0.423300  0.370400
      Read-in Center 5747 is at   0.000000 -0.370400  0.370400
      Read-in Center 5748 is at   0.000000 -0.317500  0.370400
      Read-in Center 5749 is at   0.000000 -0.264600  0.370400
      Read-in Center 5750 is at   0.000000 -0.211700  0.370400
      Read-in Center 5751 is at   0.000000 -0.158800  0.370400
      Read-in Center 5752 is at   0.000000 -0.105800  0.370400
      Read-in Center 5753 is at   0.000000 -0.052900  0.370400
      Read-in Center 5754 is at   0.000000  0.000000  0.370400
      Read-in Center 5755 is at   0.000000  0.052900  0.370400
      Read-in Center 5756 is at   0.000000  0.105800  0.370400
      Read-in Center 5757 is at   0.000000  0.158800  0.370400
      Read-in Center 5758 is at   0.000000  0.211700  0.370400
      Read-in Center 5759 is at   0.000000  0.264600  0.370400
      Read-in Center 5760 is at   0.000000  0.317500  0.370400
      Read-in Center 5761 is at   0.000000  0.370400  0.370400
      Read-in Center 5762 is at   0.000000  0.423300  0.370400
      Read-in Center 5763 is at   0.000000  0.476300  0.370400
      Read-in Center 5764 is at   0.000000  0.529200  0.370400
      Read-in Center 5765 is at   0.000000  0.582100  0.370400
      Read-in Center 5766 is at   0.000000  0.635000  0.370400
      Read-in Center 5767 is at   0.000000  0.687900  0.370400
      Read-in Center 5768 is at   0.000000  0.740800  0.370400
      Read-in Center 5769 is at   0.000000  0.793800  0.370400
      Read-in Center 5770 is at   0.000000  0.846700  0.370400
      Read-in Center 5771 is at   0.000000  0.899600  0.370400
      Read-in Center 5772 is at   0.000000  0.952500  0.370400
      Read-in Center 5773 is at   0.000000  1.005400  0.370400
      Read-in Center 5774 is at   0.000000  1.058400  0.370400
      Read-in Center 5775 is at   0.000000  1.111300  0.370400
      Read-in Center 5776 is at   0.000000  1.164200  0.370400
      Read-in Center 5777 is at   0.000000  1.217100  0.370400
      Read-in Center 5778 is at   0.000000  1.270000  0.370400
      Read-in Center 5779 is at   0.000000  1.322900  0.370400
      Read-in Center 5780 is at   0.000000  1.375900  0.370400
      Read-in Center 5781 is at   0.000000  1.428800  0.370400
      Read-in Center 5782 is at   0.000000  1.481700  0.370400
      Read-in Center 5783 is at   0.000000  1.534600  0.370400
      Read-in Center 5784 is at   0.000000  1.587500  0.370400
      Read-in Center 5785 is at   0.000000  1.640400  0.370400
      Read-in Center 5786 is at   0.000000  1.693400  0.370400
      Read-in Center 5787 is at   0.000000  1.746300  0.370400
      Read-in Center 5788 is at   0.000000  1.799200  0.370400
      Read-in Center 5789 is at   0.000000  1.852100  0.370400
      Read-in Center 5790 is at   0.000000  1.905000  0.370400
      Read-in Center 5791 is at   0.000000  1.958000  0.370400
      Read-in Center 5792 is at   0.000000  2.010900  0.370400
      Read-in Center 5793 is at   0.000000  2.063800  0.370400
      Read-in Center 5794 is at   0.000000  2.116700  0.370400
      Read-in Center 5795 is at   0.000000  2.169600  0.370400
      Read-in Center 5796 is at   0.000000  2.222500  0.370400
      Read-in Center 5797 is at   0.000000  2.275500  0.370400
      Read-in Center 5798 is at   0.000000  2.328400  0.370400
      Read-in Center 5799 is at   0.000000  2.381300  0.370400
      Read-in Center 5800 is at   0.000000  2.434200  0.370400
      Read-in Center 5801 is at   0.000000  2.487100  0.370400
      Read-in Center 5802 is at   0.000000  2.540000  0.370400
      Read-in Center 5803 is at   0.000000  2.593000  0.370400
      Read-in Center 5804 is at   0.000000 -2.645900  0.423300
      Read-in Center 5805 is at   0.000000 -2.593000  0.423300
      Read-in Center 5806 is at   0.000000 -2.540000  0.423300
      Read-in Center 5807 is at   0.000000 -2.487100  0.423300
      Read-in Center 5808 is at   0.000000 -2.434200  0.423300
      Read-in Center 5809 is at   0.000000 -2.381300  0.423300
      Read-in Center 5810 is at   0.000000 -2.328400  0.423300
      Read-in Center 5811 is at   0.000000 -2.275500  0.423300
      Read-in Center 5812 is at   0.000000 -2.222500  0.423300
      Read-in Center 5813 is at   0.000000 -2.169600  0.423300
      Read-in Center 5814 is at   0.000000 -2.116700  0.423300
      Read-in Center 5815 is at   0.000000 -2.063800  0.423300
      Read-in Center 5816 is at   0.000000 -2.010900  0.423300
      Read-in Center 5817 is at   0.000000 -1.958000  0.423300
      Read-in Center 5818 is at   0.000000 -1.905000  0.423300
      Read-in Center 5819 is at   0.000000 -1.852100  0.423300
      Read-in Center 5820 is at   0.000000 -1.799200  0.423300
      Read-in Center 5821 is at   0.000000 -1.746300  0.423300
      Read-in Center 5822 is at   0.000000 -1.693400  0.423300
      Read-in Center 5823 is at   0.000000 -1.640400  0.423300
      Read-in Center 5824 is at   0.000000 -1.587500  0.423300
      Read-in Center 5825 is at   0.000000 -1.534600  0.423300
      Read-in Center 5826 is at   0.000000 -1.481700  0.423300
      Read-in Center 5827 is at   0.000000 -1.428800  0.423300
      Read-in Center 5828 is at   0.000000 -1.375900  0.423300
      Read-in Center 5829 is at   0.000000 -1.322900  0.423300
      Read-in Center 5830 is at   0.000000 -1.270000  0.423300
      Read-in Center 5831 is at   0.000000 -1.217100  0.423300
      Read-in Center 5832 is at   0.000000 -1.164200  0.423300
      Read-in Center 5833 is at   0.000000 -1.111300  0.423300
      Read-in Center 5834 is at   0.000000 -1.058400  0.423300
      Read-in Center 5835 is at   0.000000 -1.005400  0.423300
      Read-in Center 5836 is at   0.000000 -0.952500  0.423300
      Read-in Center 5837 is at   0.000000 -0.899600  0.423300
      Read-in Center 5838 is at   0.000000 -0.846700  0.423300
      Read-in Center 5839 is at   0.000000 -0.793800  0.423300
      Read-in Center 5840 is at   0.000000 -0.740800  0.423300
      Read-in Center 5841 is at   0.000000 -0.687900  0.423300
      Read-in Center 5842 is at   0.000000 -0.635000  0.423300
      Read-in Center 5843 is at   0.000000 -0.582100  0.423300
      Read-in Center 5844 is at   0.000000 -0.529200  0.423300
      Read-in Center 5845 is at   0.000000 -0.476300  0.423300
      Read-in Center 5846 is at   0.000000 -0.423300  0.423300
      Read-in Center 5847 is at   0.000000 -0.370400  0.423300
      Read-in Center 5848 is at   0.000000 -0.317500  0.423300
      Read-in Center 5849 is at   0.000000 -0.264600  0.423300
      Read-in Center 5850 is at   0.000000 -0.211700  0.423300
      Read-in Center 5851 is at   0.000000 -0.158800  0.423300
      Read-in Center 5852 is at   0.000000 -0.105800  0.423300
      Read-in Center 5853 is at   0.000000 -0.052900  0.423300
      Read-in Center 5854 is at   0.000000  0.000000  0.423300
      Read-in Center 5855 is at   0.000000  0.052900  0.423300
      Read-in Center 5856 is at   0.000000  0.105800  0.423300
      Read-in Center 5857 is at   0.000000  0.158800  0.423300
      Read-in Center 5858 is at   0.000000  0.211700  0.423300
      Read-in Center 5859 is at   0.000000  0.264600  0.423300
      Read-in Center 5860 is at   0.000000  0.317500  0.423300
      Read-in Center 5861 is at   0.000000  0.370400  0.423300
      Read-in Center 5862 is at   0.000000  0.423300  0.423300
      Read-in Center 5863 is at   0.000000  0.476300  0.423300
      Read-in Center 5864 is at   0.000000  0.529200  0.423300
      Read-in Center 5865 is at   0.000000  0.582100  0.423300
      Read-in Center 5866 is at   0.000000  0.635000  0.423300
      Read-in Center 5867 is at   0.000000  0.687900  0.423300
      Read-in Center 5868 is at   0.000000  0.740800  0.423300
      Read-in Center 5869 is at   0.000000  0.793800  0.423300
      Read-in Center 5870 is at   0.000000  0.846700  0.423300
      Read-in Center 5871 is at   0.000000  0.899600  0.423300
      Read-in Center 5872 is at   0.000000  0.952500  0.423300
      Read-in Center 5873 is at   0.000000  1.005400  0.423300
      Read-in Center 5874 is at   0.000000  1.058400  0.423300
      Read-in Center 5875 is at   0.000000  1.111300  0.423300
      Read-in Center 5876 is at   0.000000  1.164200  0.423300
      Read-in Center 5877 is at   0.000000  1.217100  0.423300
      Read-in Center 5878 is at   0.000000  1.270000  0.423300
      Read-in Center 5879 is at   0.000000  1.322900  0.423300
      Read-in Center 5880 is at   0.000000  1.375900  0.423300
      Read-in Center 5881 is at   0.000000  1.428800  0.423300
      Read-in Center 5882 is at   0.000000  1.481700  0.423300
      Read-in Center 5883 is at   0.000000  1.534600  0.423300
      Read-in Center 5884 is at   0.000000  1.587500  0.423300
      Read-in Center 5885 is at   0.000000  1.640400  0.423300
      Read-in Center 5886 is at   0.000000  1.693400  0.423300
      Read-in Center 5887 is at   0.000000  1.746300  0.423300
      Read-in Center 5888 is at   0.000000  1.799200  0.423300
      Read-in Center 5889 is at   0.000000  1.852100  0.423300
      Read-in Center 5890 is at   0.000000  1.905000  0.423300
      Read-in Center 5891 is at   0.000000  1.958000  0.423300
      Read-in Center 5892 is at   0.000000  2.010900  0.423300
      Read-in Center 5893 is at   0.000000  2.063800  0.423300
      Read-in Center 5894 is at   0.000000  2.116700  0.423300
      Read-in Center 5895 is at   0.000000  2.169600  0.423300
      Read-in Center 5896 is at   0.000000  2.222500  0.423300
      Read-in Center 5897 is at   0.000000  2.275500  0.423300
      Read-in Center 5898 is at   0.000000  2.328400  0.423300
      Read-in Center 5899 is at   0.000000  2.381300  0.423300
      Read-in Center 5900 is at   0.000000  2.434200  0.423300
      Read-in Center 5901 is at   0.000000  2.487100  0.423300
      Read-in Center 5902 is at   0.000000  2.540000  0.423300
      Read-in Center 5903 is at   0.000000  2.593000  0.423300
      Read-in Center 5904 is at   0.000000 -2.645900  0.476300
      Read-in Center 5905 is at   0.000000 -2.593000  0.476300
      Read-in Center 5906 is at   0.000000 -2.540000  0.476300
      Read-in Center 5907 is at   0.000000 -2.487100  0.476300
      Read-in Center 5908 is at   0.000000 -2.434200  0.476300
      Read-in Center 5909 is at   0.000000 -2.381300  0.476300
      Read-in Center 5910 is at   0.000000 -2.328400  0.476300
      Read-in Center 5911 is at   0.000000 -2.275500  0.476300
      Read-in Center 5912 is at   0.000000 -2.222500  0.476300
      Read-in Center 5913 is at   0.000000 -2.169600  0.476300
      Read-in Center 5914 is at   0.000000 -2.116700  0.476300
      Read-in Center 5915 is at   0.000000 -2.063800  0.476300
      Read-in Center 5916 is at   0.000000 -2.010900  0.476300
      Read-in Center 5917 is at   0.000000 -1.958000  0.476300
      Read-in Center 5918 is at   0.000000 -1.905000  0.476300
      Read-in Center 5919 is at   0.000000 -1.852100  0.476300
      Read-in Center 5920 is at   0.000000 -1.799200  0.476300
      Read-in Center 5921 is at   0.000000 -1.746300  0.476300
      Read-in Center 5922 is at   0.000000 -1.693400  0.476300
      Read-in Center 5923 is at   0.000000 -1.640400  0.476300
      Read-in Center 5924 is at   0.000000 -1.587500  0.476300
      Read-in Center 5925 is at   0.000000 -1.534600  0.476300
      Read-in Center 5926 is at   0.000000 -1.481700  0.476300
      Read-in Center 5927 is at   0.000000 -1.428800  0.476300
      Read-in Center 5928 is at   0.000000 -1.375900  0.476300
      Read-in Center 5929 is at   0.000000 -1.322900  0.476300
      Read-in Center 5930 is at   0.000000 -1.270000  0.476300
      Read-in Center 5931 is at   0.000000 -1.217100  0.476300
      Read-in Center 5932 is at   0.000000 -1.164200  0.476300
      Read-in Center 5933 is at   0.000000 -1.111300  0.476300
      Read-in Center 5934 is at   0.000000 -1.058400  0.476300
      Read-in Center 5935 is at   0.000000 -1.005400  0.476300
      Read-in Center 5936 is at   0.000000 -0.952500  0.476300
      Read-in Center 5937 is at   0.000000 -0.899600  0.476300
      Read-in Center 5938 is at   0.000000 -0.846700  0.476300
      Read-in Center 5939 is at   0.000000 -0.793800  0.476300
      Read-in Center 5940 is at   0.000000 -0.740800  0.476300
      Read-in Center 5941 is at   0.000000 -0.687900  0.476300
      Read-in Center 5942 is at   0.000000 -0.635000  0.476300
      Read-in Center 5943 is at   0.000000 -0.582100  0.476300
      Read-in Center 5944 is at   0.000000 -0.529200  0.476300
      Read-in Center 5945 is at   0.000000 -0.476300  0.476300
      Read-in Center 5946 is at   0.000000 -0.423300  0.476300
      Read-in Center 5947 is at   0.000000 -0.370400  0.476300
      Read-in Center 5948 is at   0.000000 -0.317500  0.476300
      Read-in Center 5949 is at   0.000000 -0.264600  0.476300
      Read-in Center 5950 is at   0.000000 -0.211700  0.476300
      Read-in Center 5951 is at   0.000000 -0.158800  0.476300
      Read-in Center 5952 is at   0.000000 -0.105800  0.476300
      Read-in Center 5953 is at   0.000000 -0.052900  0.476300
      Read-in Center 5954 is at   0.000000  0.000000  0.476300
      Read-in Center 5955 is at   0.000000  0.052900  0.476300
      Read-in Center 5956 is at   0.000000  0.105800  0.476300
      Read-in Center 5957 is at   0.000000  0.158800  0.476300
      Read-in Center 5958 is at   0.000000  0.211700  0.476300
      Read-in Center 5959 is at   0.000000  0.264600  0.476300
      Read-in Center 5960 is at   0.000000  0.317500  0.476300
      Read-in Center 5961 is at   0.000000  0.370400  0.476300
      Read-in Center 5962 is at   0.000000  0.423300  0.476300
      Read-in Center 5963 is at   0.000000  0.476300  0.476300
      Read-in Center 5964 is at   0.000000  0.529200  0.476300
      Read-in Center 5965 is at   0.000000  0.582100  0.476300
      Read-in Center 5966 is at   0.000000  0.635000  0.476300
      Read-in Center 5967 is at   0.000000  0.687900  0.476300
      Read-in Center 5968 is at   0.000000  0.740800  0.476300
      Read-in Center 5969 is at   0.000000  0.793800  0.476300
      Read-in Center 5970 is at   0.000000  0.846700  0.476300
      Read-in Center 5971 is at   0.000000  0.899600  0.476300
      Read-in Center 5972 is at   0.000000  0.952500  0.476300
      Read-in Center 5973 is at   0.000000  1.005400  0.476300
      Read-in Center 5974 is at   0.000000  1.058400  0.476300
      Read-in Center 5975 is at   0.000000  1.111300  0.476300
      Read-in Center 5976 is at   0.000000  1.164200  0.476300
      Read-in Center 5977 is at   0.000000  1.217100  0.476300
      Read-in Center 5978 is at   0.000000  1.270000  0.476300
      Read-in Center 5979 is at   0.000000  1.322900  0.476300
      Read-in Center 5980 is at   0.000000  1.375900  0.476300
      Read-in Center 5981 is at   0.000000  1.428800  0.476300
      Read-in Center 5982 is at   0.000000  1.481700  0.476300
      Read-in Center 5983 is at   0.000000  1.534600  0.476300
      Read-in Center 5984 is at   0.000000  1.587500  0.476300
      Read-in Center 5985 is at   0.000000  1.640400  0.476300
      Read-in Center 5986 is at   0.000000  1.693400  0.476300
      Read-in Center 5987 is at   0.000000  1.746300  0.476300
      Read-in Center 5988 is at   0.000000  1.799200  0.476300
      Read-in Center 5989 is at   0.000000  1.852100  0.476300
      Read-in Center 5990 is at   0.000000  1.905000  0.476300
      Read-in Center 5991 is at   0.000000  1.958000  0.476300
      Read-in Center 5992 is at   0.000000  2.010900  0.476300
      Read-in Center 5993 is at   0.000000  2.063800  0.476300
      Read-in Center 5994 is at   0.000000  2.116700  0.476300
      Read-in Center 5995 is at   0.000000  2.169600  0.476300
      Read-in Center 5996 is at   0.000000  2.222500  0.476300
      Read-in Center 5997 is at   0.000000  2.275500  0.476300
      Read-in Center 5998 is at   0.000000  2.328400  0.476300
      Read-in Center 5999 is at   0.000000  2.381300  0.476300
      Read-in Center 6000 is at   0.000000  2.434200  0.476300
      Read-in Center 6001 is at   0.000000  2.487100  0.476300
      Read-in Center 6002 is at   0.000000  2.540000  0.476300
      Read-in Center 6003 is at   0.000000  2.593000  0.476300
      Read-in Center 6004 is at   0.000000 -2.645900  0.529200
      Read-in Center 6005 is at   0.000000 -2.593000  0.529200
      Read-in Center 6006 is at   0.000000 -2.540000  0.529200
      Read-in Center 6007 is at   0.000000 -2.487100  0.529200
      Read-in Center 6008 is at   0.000000 -2.434200  0.529200
      Read-in Center 6009 is at   0.000000 -2.381300  0.529200
      Read-in Center 6010 is at   0.000000 -2.328400  0.529200
      Read-in Center 6011 is at   0.000000 -2.275500  0.529200
      Read-in Center 6012 is at   0.000000 -2.222500  0.529200
      Read-in Center 6013 is at   0.000000 -2.169600  0.529200
      Read-in Center 6014 is at   0.000000 -2.116700  0.529200
      Read-in Center 6015 is at   0.000000 -2.063800  0.529200
      Read-in Center 6016 is at   0.000000 -2.010900  0.529200
      Read-in Center 6017 is at   0.000000 -1.958000  0.529200
      Read-in Center 6018 is at   0.000000 -1.905000  0.529200
      Read-in Center 6019 is at   0.000000 -1.852100  0.529200
      Read-in Center 6020 is at   0.000000 -1.799200  0.529200
      Read-in Center 6021 is at   0.000000 -1.746300  0.529200
      Read-in Center 6022 is at   0.000000 -1.693400  0.529200
      Read-in Center 6023 is at   0.000000 -1.640400  0.529200
      Read-in Center 6024 is at   0.000000 -1.587500  0.529200
      Read-in Center 6025 is at   0.000000 -1.534600  0.529200
      Read-in Center 6026 is at   0.000000 -1.481700  0.529200
      Read-in Center 6027 is at   0.000000 -1.428800  0.529200
      Read-in Center 6028 is at   0.000000 -1.375900  0.529200
      Read-in Center 6029 is at   0.000000 -1.322900  0.529200
      Read-in Center 6030 is at   0.000000 -1.270000  0.529200
      Read-in Center 6031 is at   0.000000 -1.217100  0.529200
      Read-in Center 6032 is at   0.000000 -1.164200  0.529200
      Read-in Center 6033 is at   0.000000 -1.111300  0.529200
      Read-in Center 6034 is at   0.000000 -1.058400  0.529200
      Read-in Center 6035 is at   0.000000 -1.005400  0.529200
      Read-in Center 6036 is at   0.000000 -0.952500  0.529200
      Read-in Center 6037 is at   0.000000 -0.899600  0.529200
      Read-in Center 6038 is at   0.000000 -0.846700  0.529200
      Read-in Center 6039 is at   0.000000 -0.793800  0.529200
      Read-in Center 6040 is at   0.000000 -0.740800  0.529200
      Read-in Center 6041 is at   0.000000 -0.687900  0.529200
      Read-in Center 6042 is at   0.000000 -0.635000  0.529200
      Read-in Center 6043 is at   0.000000 -0.582100  0.529200
      Read-in Center 6044 is at   0.000000 -0.529200  0.529200
      Read-in Center 6045 is at   0.000000 -0.476300  0.529200
      Read-in Center 6046 is at   0.000000 -0.423300  0.529200
      Read-in Center 6047 is at   0.000000 -0.370400  0.529200
      Read-in Center 6048 is at   0.000000 -0.317500  0.529200
      Read-in Center 6049 is at   0.000000 -0.264600  0.529200
      Read-in Center 6050 is at   0.000000 -0.211700  0.529200
      Read-in Center 6051 is at   0.000000 -0.158800  0.529200
      Read-in Center 6052 is at   0.000000 -0.105800  0.529200
      Read-in Center 6053 is at   0.000000 -0.052900  0.529200
      Read-in Center 6054 is at   0.000000  0.000000  0.529200
      Read-in Center 6055 is at   0.000000  0.052900  0.529200
      Read-in Center 6056 is at   0.000000  0.105800  0.529200
      Read-in Center 6057 is at   0.000000  0.158800  0.529200
      Read-in Center 6058 is at   0.000000  0.211700  0.529200
      Read-in Center 6059 is at   0.000000  0.264600  0.529200
      Read-in Center 6060 is at   0.000000  0.317500  0.529200
      Read-in Center 6061 is at   0.000000  0.370400  0.529200
      Read-in Center 6062 is at   0.000000  0.423300  0.529200
      Read-in Center 6063 is at   0.000000  0.476300  0.529200
      Read-in Center 6064 is at   0.000000  0.529200  0.529200
      Read-in Center 6065 is at   0.000000  0.582100  0.529200
      Read-in Center 6066 is at   0.000000  0.635000  0.529200
      Read-in Center 6067 is at   0.000000  0.687900  0.529200
      Read-in Center 6068 is at   0.000000  0.740800  0.529200
      Read-in Center 6069 is at   0.000000  0.793800  0.529200
      Read-in Center 6070 is at   0.000000  0.846700  0.529200
      Read-in Center 6071 is at   0.000000  0.899600  0.529200
      Read-in Center 6072 is at   0.000000  0.952500  0.529200
      Read-in Center 6073 is at   0.000000  1.005400  0.529200
      Read-in Center 6074 is at   0.000000  1.058400  0.529200
      Read-in Center 6075 is at   0.000000  1.111300  0.529200
      Read-in Center 6076 is at   0.000000  1.164200  0.529200
      Read-in Center 6077 is at   0.000000  1.217100  0.529200
      Read-in Center 6078 is at   0.000000  1.270000  0.529200
      Read-in Center 6079 is at   0.000000  1.322900  0.529200
      Read-in Center 6080 is at   0.000000  1.375900  0.529200
      Read-in Center 6081 is at   0.000000  1.428800  0.529200
      Read-in Center 6082 is at   0.000000  1.481700  0.529200
      Read-in Center 6083 is at   0.000000  1.534600  0.529200
      Read-in Center 6084 is at   0.000000  1.587500  0.529200
      Read-in Center 6085 is at   0.000000  1.640400  0.529200
      Read-in Center 6086 is at   0.000000  1.693400  0.529200
      Read-in Center 6087 is at   0.000000  1.746300  0.529200
      Read-in Center 6088 is at   0.000000  1.799200  0.529200
      Read-in Center 6089 is at   0.000000  1.852100  0.529200
      Read-in Center 6090 is at   0.000000  1.905000  0.529200
      Read-in Center 6091 is at   0.000000  1.958000  0.529200
      Read-in Center 6092 is at   0.000000  2.010900  0.529200
      Read-in Center 6093 is at   0.000000  2.063800  0.529200
      Read-in Center 6094 is at   0.000000  2.116700  0.529200
      Read-in Center 6095 is at   0.000000  2.169600  0.529200
      Read-in Center 6096 is at   0.000000  2.222500  0.529200
      Read-in Center 6097 is at   0.000000  2.275500  0.529200
      Read-in Center 6098 is at   0.000000  2.328400  0.529200
      Read-in Center 6099 is at   0.000000  2.381300  0.529200
      Read-in Center 6100 is at   0.000000  2.434200  0.529200
      Read-in Center 6101 is at   0.000000  2.487100  0.529200
      Read-in Center 6102 is at   0.000000  2.540000  0.529200
      Read-in Center 6103 is at   0.000000  2.593000  0.529200
      Read-in Center 6104 is at   0.000000 -2.645900  0.582100
      Read-in Center 6105 is at   0.000000 -2.593000  0.582100
      Read-in Center 6106 is at   0.000000 -2.540000  0.582100
      Read-in Center 6107 is at   0.000000 -2.487100  0.582100
      Read-in Center 6108 is at   0.000000 -2.434200  0.582100
      Read-in Center 6109 is at   0.000000 -2.381300  0.582100
      Read-in Center 6110 is at   0.000000 -2.328400  0.582100
      Read-in Center 6111 is at   0.000000 -2.275500  0.582100
      Read-in Center 6112 is at   0.000000 -2.222500  0.582100
      Read-in Center 6113 is at   0.000000 -2.169600  0.582100
      Read-in Center 6114 is at   0.000000 -2.116700  0.582100
      Read-in Center 6115 is at   0.000000 -2.063800  0.582100
      Read-in Center 6116 is at   0.000000 -2.010900  0.582100
      Read-in Center 6117 is at   0.000000 -1.958000  0.582100
      Read-in Center 6118 is at   0.000000 -1.905000  0.582100
      Read-in Center 6119 is at   0.000000 -1.852100  0.582100
      Read-in Center 6120 is at   0.000000 -1.799200  0.582100
      Read-in Center 6121 is at   0.000000 -1.746300  0.582100
      Read-in Center 6122 is at   0.000000 -1.693400  0.582100
      Read-in Center 6123 is at   0.000000 -1.640400  0.582100
      Read-in Center 6124 is at   0.000000 -1.587500  0.582100
      Read-in Center 6125 is at   0.000000 -1.534600  0.582100
      Read-in Center 6126 is at   0.000000 -1.481700  0.582100
      Read-in Center 6127 is at   0.000000 -1.428800  0.582100
      Read-in Center 6128 is at   0.000000 -1.375900  0.582100
      Read-in Center 6129 is at   0.000000 -1.322900  0.582100
      Read-in Center 6130 is at   0.000000 -1.270000  0.582100
      Read-in Center 6131 is at   0.000000 -1.217100  0.582100
      Read-in Center 6132 is at   0.000000 -1.164200  0.582100
      Read-in Center 6133 is at   0.000000 -1.111300  0.582100
      Read-in Center 6134 is at   0.000000 -1.058400  0.582100
      Read-in Center 6135 is at   0.000000 -1.005400  0.582100
      Read-in Center 6136 is at   0.000000 -0.952500  0.582100
      Read-in Center 6137 is at   0.000000 -0.899600  0.582100
      Read-in Center 6138 is at   0.000000 -0.846700  0.582100
      Read-in Center 6139 is at   0.000000 -0.793800  0.582100
      Read-in Center 6140 is at   0.000000 -0.740800  0.582100
      Read-in Center 6141 is at   0.000000 -0.687900  0.582100
      Read-in Center 6142 is at   0.000000 -0.635000  0.582100
      Read-in Center 6143 is at   0.000000 -0.582100  0.582100
      Read-in Center 6144 is at   0.000000 -0.529200  0.582100
      Read-in Center 6145 is at   0.000000 -0.476300  0.582100
      Read-in Center 6146 is at   0.000000 -0.423300  0.582100
      Read-in Center 6147 is at   0.000000 -0.370400  0.582100
      Read-in Center 6148 is at   0.000000 -0.317500  0.582100
      Read-in Center 6149 is at   0.000000 -0.264600  0.582100
      Read-in Center 6150 is at   0.000000 -0.211700  0.582100
      Read-in Center 6151 is at   0.000000 -0.158800  0.582100
      Read-in Center 6152 is at   0.000000 -0.105800  0.582100
      Read-in Center 6153 is at   0.000000 -0.052900  0.582100
      Read-in Center 6154 is at   0.000000  0.000000  0.582100
      Read-in Center 6155 is at   0.000000  0.052900  0.582100
      Read-in Center 6156 is at   0.000000  0.105800  0.582100
      Read-in Center 6157 is at   0.000000  0.158800  0.582100
      Read-in Center 6158 is at   0.000000  0.211700  0.582100
      Read-in Center 6159 is at   0.000000  0.264600  0.582100
      Read-in Center 6160 is at   0.000000  0.317500  0.582100
      Read-in Center 6161 is at   0.000000  0.370400  0.582100
      Read-in Center 6162 is at   0.000000  0.423300  0.582100
      Read-in Center 6163 is at   0.000000  0.476300  0.582100
      Read-in Center 6164 is at   0.000000  0.529200  0.582100
      Read-in Center 6165 is at   0.000000  0.582100  0.582100
      Read-in Center 6166 is at   0.000000  0.635000  0.582100
      Read-in Center 6167 is at   0.000000  0.687900  0.582100
      Read-in Center 6168 is at   0.000000  0.740800  0.582100
      Read-in Center 6169 is at   0.000000  0.793800  0.582100
      Read-in Center 6170 is at   0.000000  0.846700  0.582100
      Read-in Center 6171 is at   0.000000  0.899600  0.582100
      Read-in Center 6172 is at   0.000000  0.952500  0.582100
      Read-in Center 6173 is at   0.000000  1.005400  0.582100
      Read-in Center 6174 is at   0.000000  1.058400  0.582100
      Read-in Center 6175 is at   0.000000  1.111300  0.582100
      Read-in Center 6176 is at   0.000000  1.164200  0.582100
      Read-in Center 6177 is at   0.000000  1.217100  0.582100
      Read-in Center 6178 is at   0.000000  1.270000  0.582100
      Read-in Center 6179 is at   0.000000  1.322900  0.582100
      Read-in Center 6180 is at   0.000000  1.375900  0.582100
      Read-in Center 6181 is at   0.000000  1.428800  0.582100
      Read-in Center 6182 is at   0.000000  1.481700  0.582100
      Read-in Center 6183 is at   0.000000  1.534600  0.582100
      Read-in Center 6184 is at   0.000000  1.587500  0.582100
      Read-in Center 6185 is at   0.000000  1.640400  0.582100
      Read-in Center 6186 is at   0.000000  1.693400  0.582100
      Read-in Center 6187 is at   0.000000  1.746300  0.582100
      Read-in Center 6188 is at   0.000000  1.799200  0.582100
      Read-in Center 6189 is at   0.000000  1.852100  0.582100
      Read-in Center 6190 is at   0.000000  1.905000  0.582100
      Read-in Center 6191 is at   0.000000  1.958000  0.582100
      Read-in Center 6192 is at   0.000000  2.010900  0.582100
      Read-in Center 6193 is at   0.000000  2.063800  0.582100
      Read-in Center 6194 is at   0.000000  2.116700  0.582100
      Read-in Center 6195 is at   0.000000  2.169600  0.582100
      Read-in Center 6196 is at   0.000000  2.222500  0.582100
      Read-in Center 6197 is at   0.000000  2.275500  0.582100
      Read-in Center 6198 is at   0.000000  2.328400  0.582100
      Read-in Center 6199 is at   0.000000  2.381300  0.582100
      Read-in Center 6200 is at   0.000000  2.434200  0.582100
      Read-in Center 6201 is at   0.000000  2.487100  0.582100
      Read-in Center 6202 is at   0.000000  2.540000  0.582100
      Read-in Center 6203 is at   0.000000  2.593000  0.582100
      Read-in Center 6204 is at   0.000000 -2.645900  0.635000
      Read-in Center 6205 is at   0.000000 -2.593000  0.635000
      Read-in Center 6206 is at   0.000000 -2.540000  0.635000
      Read-in Center 6207 is at   0.000000 -2.487100  0.635000
      Read-in Center 6208 is at   0.000000 -2.434200  0.635000
      Read-in Center 6209 is at   0.000000 -2.381300  0.635000
      Read-in Center 6210 is at   0.000000 -2.328400  0.635000
      Read-in Center 6211 is at   0.000000 -2.275500  0.635000
      Read-in Center 6212 is at   0.000000 -2.222500  0.635000
      Read-in Center 6213 is at   0.000000 -2.169600  0.635000
      Read-in Center 6214 is at   0.000000 -2.116700  0.635000
      Read-in Center 6215 is at   0.000000 -2.063800  0.635000
      Read-in Center 6216 is at   0.000000 -2.010900  0.635000
      Read-in Center 6217 is at   0.000000 -1.958000  0.635000
      Read-in Center 6218 is at   0.000000 -1.905000  0.635000
      Read-in Center 6219 is at   0.000000 -1.852100  0.635000
      Read-in Center 6220 is at   0.000000 -1.799200  0.635000
      Read-in Center 6221 is at   0.000000 -1.746300  0.635000
      Read-in Center 6222 is at   0.000000 -1.693400  0.635000
      Read-in Center 6223 is at   0.000000 -1.640400  0.635000
      Read-in Center 6224 is at   0.000000 -1.587500  0.635000
      Read-in Center 6225 is at   0.000000 -1.534600  0.635000
      Read-in Center 6226 is at   0.000000 -1.481700  0.635000
      Read-in Center 6227 is at   0.000000 -1.428800  0.635000
      Read-in Center 6228 is at   0.000000 -1.375900  0.635000
      Read-in Center 6229 is at   0.000000 -1.322900  0.635000
      Read-in Center 6230 is at   0.000000 -1.270000  0.635000
      Read-in Center 6231 is at   0.000000 -1.217100  0.635000
      Read-in Center 6232 is at   0.000000 -1.164200  0.635000
      Read-in Center 6233 is at   0.000000 -1.111300  0.635000
      Read-in Center 6234 is at   0.000000 -1.058400  0.635000
      Read-in Center 6235 is at   0.000000 -1.005400  0.635000
      Read-in Center 6236 is at   0.000000 -0.952500  0.635000
      Read-in Center 6237 is at   0.000000 -0.899600  0.635000
      Read-in Center 6238 is at   0.000000 -0.846700  0.635000
      Read-in Center 6239 is at   0.000000 -0.793800  0.635000
      Read-in Center 6240 is at   0.000000 -0.740800  0.635000
      Read-in Center 6241 is at   0.000000 -0.687900  0.635000
      Read-in Center 6242 is at   0.000000 -0.635000  0.635000
      Read-in Center 6243 is at   0.000000 -0.582100  0.635000
      Read-in Center 6244 is at   0.000000 -0.529200  0.635000
      Read-in Center 6245 is at   0.000000 -0.476300  0.635000
      Read-in Center 6246 is at   0.000000 -0.423300  0.635000
      Read-in Center 6247 is at   0.000000 -0.370400  0.635000
      Read-in Center 6248 is at   0.000000 -0.317500  0.635000
      Read-in Center 6249 is at   0.000000 -0.264600  0.635000
      Read-in Center 6250 is at   0.000000 -0.211700  0.635000
      Read-in Center 6251 is at   0.000000 -0.158800  0.635000
      Read-in Center 6252 is at   0.000000 -0.105800  0.635000
      Read-in Center 6253 is at   0.000000 -0.052900  0.635000
      Read-in Center 6254 is at   0.000000  0.000000  0.635000
      Read-in Center 6255 is at   0.000000  0.052900  0.635000
      Read-in Center 6256 is at   0.000000  0.105800  0.635000
      Read-in Center 6257 is at   0.000000  0.158800  0.635000
      Read-in Center 6258 is at   0.000000  0.211700  0.635000
      Read-in Center 6259 is at   0.000000  0.264600  0.635000
      Read-in Center 6260 is at   0.000000  0.317500  0.635000
      Read-in Center 6261 is at   0.000000  0.370400  0.635000
      Read-in Center 6262 is at   0.000000  0.423300  0.635000
      Read-in Center 6263 is at   0.000000  0.476300  0.635000
      Read-in Center 6264 is at   0.000000  0.529200  0.635000
      Read-in Center 6265 is at   0.000000  0.582100  0.635000
      Read-in Center 6266 is at   0.000000  0.635000  0.635000
      Read-in Center 6267 is at   0.000000  0.687900  0.635000
      Read-in Center 6268 is at   0.000000  0.740800  0.635000
      Read-in Center 6269 is at   0.000000  0.793800  0.635000
      Read-in Center 6270 is at   0.000000  0.846700  0.635000
      Read-in Center 6271 is at   0.000000  0.899600  0.635000
      Read-in Center 6272 is at   0.000000  0.952500  0.635000
      Read-in Center 6273 is at   0.000000  1.005400  0.635000
      Read-in Center 6274 is at   0.000000  1.058400  0.635000
      Read-in Center 6275 is at   0.000000  1.111300  0.635000
      Read-in Center 6276 is at   0.000000  1.164200  0.635000
      Read-in Center 6277 is at   0.000000  1.217100  0.635000
      Read-in Center 6278 is at   0.000000  1.270000  0.635000
      Read-in Center 6279 is at   0.000000  1.322900  0.635000
      Read-in Center 6280 is at   0.000000  1.375900  0.635000
      Read-in Center 6281 is at   0.000000  1.428800  0.635000
      Read-in Center 6282 is at   0.000000  1.481700  0.635000
      Read-in Center 6283 is at   0.000000  1.534600  0.635000
      Read-in Center 6284 is at   0.000000  1.587500  0.635000
      Read-in Center 6285 is at   0.000000  1.640400  0.635000
      Read-in Center 6286 is at   0.000000  1.693400  0.635000
      Read-in Center 6287 is at   0.000000  1.746300  0.635000
      Read-in Center 6288 is at   0.000000  1.799200  0.635000
      Read-in Center 6289 is at   0.000000  1.852100  0.635000
      Read-in Center 6290 is at   0.000000  1.905000  0.635000
      Read-in Center 6291 is at   0.000000  1.958000  0.635000
      Read-in Center 6292 is at   0.000000  2.010900  0.635000
      Read-in Center 6293 is at   0.000000  2.063800  0.635000
      Read-in Center 6294 is at   0.000000  2.116700  0.635000
      Read-in Center 6295 is at   0.000000  2.169600  0.635000
      Read-in Center 6296 is at   0.000000  2.222500  0.635000
      Read-in Center 6297 is at   0.000000  2.275500  0.635000
      Read-in Center 6298 is at   0.000000  2.328400  0.635000
      Read-in Center 6299 is at   0.000000  2.381300  0.635000
      Read-in Center 6300 is at   0.000000  2.434200  0.635000
      Read-in Center 6301 is at   0.000000  2.487100  0.635000
      Read-in Center 6302 is at   0.000000  2.540000  0.635000
      Read-in Center 6303 is at   0.000000  2.593000  0.635000
      Read-in Center 6304 is at   0.000000 -2.645900  0.687900
      Read-in Center 6305 is at   0.000000 -2.593000  0.687900
      Read-in Center 6306 is at   0.000000 -2.540000  0.687900
      Read-in Center 6307 is at   0.000000 -2.487100  0.687900
      Read-in Center 6308 is at   0.000000 -2.434200  0.687900
      Read-in Center 6309 is at   0.000000 -2.381300  0.687900
      Read-in Center 6310 is at   0.000000 -2.328400  0.687900
      Read-in Center 6311 is at   0.000000 -2.275500  0.687900
      Read-in Center 6312 is at   0.000000 -2.222500  0.687900
      Read-in Center 6313 is at   0.000000 -2.169600  0.687900
      Read-in Center 6314 is at   0.000000 -2.116700  0.687900
      Read-in Center 6315 is at   0.000000 -2.063800  0.687900
      Read-in Center 6316 is at   0.000000 -2.010900  0.687900
      Read-in Center 6317 is at   0.000000 -1.958000  0.687900
      Read-in Center 6318 is at   0.000000 -1.905000  0.687900
      Read-in Center 6319 is at   0.000000 -1.852100  0.687900
      Read-in Center 6320 is at   0.000000 -1.799200  0.687900
      Read-in Center 6321 is at   0.000000 -1.746300  0.687900
      Read-in Center 6322 is at   0.000000 -1.693400  0.687900
      Read-in Center 6323 is at   0.000000 -1.640400  0.687900
      Read-in Center 6324 is at   0.000000 -1.587500  0.687900
      Read-in Center 6325 is at   0.000000 -1.534600  0.687900
      Read-in Center 6326 is at   0.000000 -1.481700  0.687900
      Read-in Center 6327 is at   0.000000 -1.428800  0.687900
      Read-in Center 6328 is at   0.000000 -1.375900  0.687900
      Read-in Center 6329 is at   0.000000 -1.322900  0.687900
      Read-in Center 6330 is at   0.000000 -1.270000  0.687900
      Read-in Center 6331 is at   0.000000 -1.217100  0.687900
      Read-in Center 6332 is at   0.000000 -1.164200  0.687900
      Read-in Center 6333 is at   0.000000 -1.111300  0.687900
      Read-in Center 6334 is at   0.000000 -1.058400  0.687900
      Read-in Center 6335 is at   0.000000 -1.005400  0.687900
      Read-in Center 6336 is at   0.000000 -0.952500  0.687900
      Read-in Center 6337 is at   0.000000 -0.899600  0.687900
      Read-in Center 6338 is at   0.000000 -0.846700  0.687900
      Read-in Center 6339 is at   0.000000 -0.793800  0.687900
      Read-in Center 6340 is at   0.000000 -0.740800  0.687900
      Read-in Center 6341 is at   0.000000 -0.687900  0.687900
      Read-in Center 6342 is at   0.000000 -0.635000  0.687900
      Read-in Center 6343 is at   0.000000 -0.582100  0.687900
      Read-in Center 6344 is at   0.000000 -0.529200  0.687900
      Read-in Center 6345 is at   0.000000 -0.476300  0.687900
      Read-in Center 6346 is at   0.000000 -0.423300  0.687900
      Read-in Center 6347 is at   0.000000 -0.370400  0.687900
      Read-in Center 6348 is at   0.000000 -0.317500  0.687900
      Read-in Center 6349 is at   0.000000 -0.264600  0.687900
      Read-in Center 6350 is at   0.000000 -0.211700  0.687900
      Read-in Center 6351 is at   0.000000 -0.158800  0.687900
      Read-in Center 6352 is at   0.000000 -0.105800  0.687900
      Read-in Center 6353 is at   0.000000 -0.052900  0.687900
      Read-in Center 6354 is at   0.000000  0.000000  0.687900
      Read-in Center 6355 is at   0.000000  0.052900  0.687900
      Read-in Center 6356 is at   0.000000  0.105800  0.687900
      Read-in Center 6357 is at   0.000000  0.158800  0.687900
      Read-in Center 6358 is at   0.000000  0.211700  0.687900
      Read-in Center 6359 is at   0.000000  0.264600  0.687900
      Read-in Center 6360 is at   0.000000  0.317500  0.687900
      Read-in Center 6361 is at   0.000000  0.370400  0.687900
      Read-in Center 6362 is at   0.000000  0.423300  0.687900
      Read-in Center 6363 is at   0.000000  0.476300  0.687900
      Read-in Center 6364 is at   0.000000  0.529200  0.687900
      Read-in Center 6365 is at   0.000000  0.582100  0.687900
      Read-in Center 6366 is at   0.000000  0.635000  0.687900
      Read-in Center 6367 is at   0.000000  0.687900  0.687900
      Read-in Center 6368 is at   0.000000  0.740800  0.687900
      Read-in Center 6369 is at   0.000000  0.793800  0.687900
      Read-in Center 6370 is at   0.000000  0.846700  0.687900
      Read-in Center 6371 is at   0.000000  0.899600  0.687900
      Read-in Center 6372 is at   0.000000  0.952500  0.687900
      Read-in Center 6373 is at   0.000000  1.005400  0.687900
      Read-in Center 6374 is at   0.000000  1.058400  0.687900
      Read-in Center 6375 is at   0.000000  1.111300  0.687900
      Read-in Center 6376 is at   0.000000  1.164200  0.687900
      Read-in Center 6377 is at   0.000000  1.217100  0.687900
      Read-in Center 6378 is at   0.000000  1.270000  0.687900
      Read-in Center 6379 is at   0.000000  1.322900  0.687900
      Read-in Center 6380 is at   0.000000  1.375900  0.687900
      Read-in Center 6381 is at   0.000000  1.428800  0.687900
      Read-in Center 6382 is at   0.000000  1.481700  0.687900
      Read-in Center 6383 is at   0.000000  1.534600  0.687900
      Read-in Center 6384 is at   0.000000  1.587500  0.687900
      Read-in Center 6385 is at   0.000000  1.640400  0.687900
      Read-in Center 6386 is at   0.000000  1.693400  0.687900
      Read-in Center 6387 is at   0.000000  1.746300  0.687900
      Read-in Center 6388 is at   0.000000  1.799200  0.687900
      Read-in Center 6389 is at   0.000000  1.852100  0.687900
      Read-in Center 6390 is at   0.000000  1.905000  0.687900
      Read-in Center 6391 is at   0.000000  1.958000  0.687900
      Read-in Center 6392 is at   0.000000  2.010900  0.687900
      Read-in Center 6393 is at   0.000000  2.063800  0.687900
      Read-in Center 6394 is at   0.000000  2.116700  0.687900
      Read-in Center 6395 is at   0.000000  2.169600  0.687900
      Read-in Center 6396 is at   0.000000  2.222500  0.687900
      Read-in Center 6397 is at   0.000000  2.275500  0.687900
      Read-in Center 6398 is at   0.000000  2.328400  0.687900
      Read-in Center 6399 is at   0.000000  2.381300  0.687900
      Read-in Center 6400 is at   0.000000  2.434200  0.687900
      Read-in Center 6401 is at   0.000000  2.487100  0.687900
      Read-in Center 6402 is at   0.000000  2.540000  0.687900
      Read-in Center 6403 is at   0.000000  2.593000  0.687900
      Read-in Center 6404 is at   0.000000 -2.645900  0.740800
      Read-in Center 6405 is at   0.000000 -2.593000  0.740800
      Read-in Center 6406 is at   0.000000 -2.540000  0.740800
      Read-in Center 6407 is at   0.000000 -2.487100  0.740800
      Read-in Center 6408 is at   0.000000 -2.434200  0.740800
      Read-in Center 6409 is at   0.000000 -2.381300  0.740800
      Read-in Center 6410 is at   0.000000 -2.328400  0.740800
      Read-in Center 6411 is at   0.000000 -2.275500  0.740800
      Read-in Center 6412 is at   0.000000 -2.222500  0.740800
      Read-in Center 6413 is at   0.000000 -2.169600  0.740800
      Read-in Center 6414 is at   0.000000 -2.116700  0.740800
      Read-in Center 6415 is at   0.000000 -2.063800  0.740800
      Read-in Center 6416 is at   0.000000 -2.010900  0.740800
      Read-in Center 6417 is at   0.000000 -1.958000  0.740800
      Read-in Center 6418 is at   0.000000 -1.905000  0.740800
      Read-in Center 6419 is at   0.000000 -1.852100  0.740800
      Read-in Center 6420 is at   0.000000 -1.799200  0.740800
      Read-in Center 6421 is at   0.000000 -1.746300  0.740800
      Read-in Center 6422 is at   0.000000 -1.693400  0.740800
      Read-in Center 6423 is at   0.000000 -1.640400  0.740800
      Read-in Center 6424 is at   0.000000 -1.587500  0.740800
      Read-in Center 6425 is at   0.000000 -1.534600  0.740800
      Read-in Center 6426 is at   0.000000 -1.481700  0.740800
      Read-in Center 6427 is at   0.000000 -1.428800  0.740800
      Read-in Center 6428 is at   0.000000 -1.375900  0.740800
      Read-in Center 6429 is at   0.000000 -1.322900  0.740800
      Read-in Center 6430 is at   0.000000 -1.270000  0.740800
      Read-in Center 6431 is at   0.000000 -1.217100  0.740800
      Read-in Center 6432 is at   0.000000 -1.164200  0.740800
      Read-in Center 6433 is at   0.000000 -1.111300  0.740800
      Read-in Center 6434 is at   0.000000 -1.058400  0.740800
      Read-in Center 6435 is at   0.000000 -1.005400  0.740800
      Read-in Center 6436 is at   0.000000 -0.952500  0.740800
      Read-in Center 6437 is at   0.000000 -0.899600  0.740800
      Read-in Center 6438 is at   0.000000 -0.846700  0.740800
      Read-in Center 6439 is at   0.000000 -0.793800  0.740800
      Read-in Center 6440 is at   0.000000 -0.740800  0.740800
      Read-in Center 6441 is at   0.000000 -0.687900  0.740800
      Read-in Center 6442 is at   0.000000 -0.635000  0.740800
      Read-in Center 6443 is at   0.000000 -0.582100  0.740800
      Read-in Center 6444 is at   0.000000 -0.529200  0.740800
      Read-in Center 6445 is at   0.000000 -0.476300  0.740800
      Read-in Center 6446 is at   0.000000 -0.423300  0.740800
      Read-in Center 6447 is at   0.000000 -0.370400  0.740800
      Read-in Center 6448 is at   0.000000 -0.317500  0.740800
      Read-in Center 6449 is at   0.000000 -0.264600  0.740800
      Read-in Center 6450 is at   0.000000 -0.211700  0.740800
      Read-in Center 6451 is at   0.000000 -0.158800  0.740800
      Read-in Center 6452 is at   0.000000 -0.105800  0.740800
      Read-in Center 6453 is at   0.000000 -0.052900  0.740800
      Read-in Center 6454 is at   0.000000  0.000000  0.740800
      Read-in Center 6455 is at   0.000000  0.052900  0.740800
      Read-in Center 6456 is at   0.000000  0.105800  0.740800
      Read-in Center 6457 is at   0.000000  0.158800  0.740800
      Read-in Center 6458 is at   0.000000  0.211700  0.740800
      Read-in Center 6459 is at   0.000000  0.264600  0.740800
      Read-in Center 6460 is at   0.000000  0.317500  0.740800
      Read-in Center 6461 is at   0.000000  0.370400  0.740800
      Read-in Center 6462 is at   0.000000  0.423300  0.740800
      Read-in Center 6463 is at   0.000000  0.476300  0.740800
      Read-in Center 6464 is at   0.000000  0.529200  0.740800
      Read-in Center 6465 is at   0.000000  0.582100  0.740800
      Read-in Center 6466 is at   0.000000  0.635000  0.740800
      Read-in Center 6467 is at   0.000000  0.687900  0.740800
      Read-in Center 6468 is at   0.000000  0.740800  0.740800
      Read-in Center 6469 is at   0.000000  0.793800  0.740800
      Read-in Center 6470 is at   0.000000  0.846700  0.740800
      Read-in Center 6471 is at   0.000000  0.899600  0.740800
      Read-in Center 6472 is at   0.000000  0.952500  0.740800
      Read-in Center 6473 is at   0.000000  1.005400  0.740800
      Read-in Center 6474 is at   0.000000  1.058400  0.740800
      Read-in Center 6475 is at   0.000000  1.111300  0.740800
      Read-in Center 6476 is at   0.000000  1.164200  0.740800
      Read-in Center 6477 is at   0.000000  1.217100  0.740800
      Read-in Center 6478 is at   0.000000  1.270000  0.740800
      Read-in Center 6479 is at   0.000000  1.322900  0.740800
      Read-in Center 6480 is at   0.000000  1.375900  0.740800
      Read-in Center 6481 is at   0.000000  1.428800  0.740800
      Read-in Center 6482 is at   0.000000  1.481700  0.740800
      Read-in Center 6483 is at   0.000000  1.534600  0.740800
      Read-in Center 6484 is at   0.000000  1.587500  0.740800
      Read-in Center 6485 is at   0.000000  1.640400  0.740800
      Read-in Center 6486 is at   0.000000  1.693400  0.740800
      Read-in Center 6487 is at   0.000000  1.746300  0.740800
      Read-in Center 6488 is at   0.000000  1.799200  0.740800
      Read-in Center 6489 is at   0.000000  1.852100  0.740800
      Read-in Center 6490 is at   0.000000  1.905000  0.740800
      Read-in Center 6491 is at   0.000000  1.958000  0.740800
      Read-in Center 6492 is at   0.000000  2.010900  0.740800
      Read-in Center 6493 is at   0.000000  2.063800  0.740800
      Read-in Center 6494 is at   0.000000  2.116700  0.740800
      Read-in Center 6495 is at   0.000000  2.169600  0.740800
      Read-in Center 6496 is at   0.000000  2.222500  0.740800
      Read-in Center 6497 is at   0.000000  2.275500  0.740800
      Read-in Center 6498 is at   0.000000  2.328400  0.740800
      Read-in Center 6499 is at   0.000000  2.381300  0.740800
      Read-in Center 6500 is at   0.000000  2.434200  0.740800
      Read-in Center 6501 is at   0.000000  2.487100  0.740800
      Read-in Center 6502 is at   0.000000  2.540000  0.740800
      Read-in Center 6503 is at   0.000000  2.593000  0.740800
      Read-in Center 6504 is at   0.000000 -2.645900  0.793800
      Read-in Center 6505 is at   0.000000 -2.593000  0.793800
      Read-in Center 6506 is at   0.000000 -2.540000  0.793800
      Read-in Center 6507 is at   0.000000 -2.487100  0.793800
      Read-in Center 6508 is at   0.000000 -2.434200  0.793800
      Read-in Center 6509 is at   0.000000 -2.381300  0.793800
      Read-in Center 6510 is at   0.000000 -2.328400  0.793800
      Read-in Center 6511 is at   0.000000 -2.275500  0.793800
      Read-in Center 6512 is at   0.000000 -2.222500  0.793800
      Read-in Center 6513 is at   0.000000 -2.169600  0.793800
      Read-in Center 6514 is at   0.000000 -2.116700  0.793800
      Read-in Center 6515 is at   0.000000 -2.063800  0.793800
      Read-in Center 6516 is at   0.000000 -2.010900  0.793800
      Read-in Center 6517 is at   0.000000 -1.958000  0.793800
      Read-in Center 6518 is at   0.000000 -1.905000  0.793800
      Read-in Center 6519 is at   0.000000 -1.852100  0.793800
      Read-in Center 6520 is at   0.000000 -1.799200  0.793800
      Read-in Center 6521 is at   0.000000 -1.746300  0.793800
      Read-in Center 6522 is at   0.000000 -1.693400  0.793800
      Read-in Center 6523 is at   0.000000 -1.640400  0.793800
      Read-in Center 6524 is at   0.000000 -1.587500  0.793800
      Read-in Center 6525 is at   0.000000 -1.534600  0.793800
      Read-in Center 6526 is at   0.000000 -1.481700  0.793800
      Read-in Center 6527 is at   0.000000 -1.428800  0.793800
      Read-in Center 6528 is at   0.000000 -1.375900  0.793800
      Read-in Center 6529 is at   0.000000 -1.322900  0.793800
      Read-in Center 6530 is at   0.000000 -1.270000  0.793800
      Read-in Center 6531 is at   0.000000 -1.217100  0.793800
      Read-in Center 6532 is at   0.000000 -1.164200  0.793800
      Read-in Center 6533 is at   0.000000 -1.111300  0.793800
      Read-in Center 6534 is at   0.000000 -1.058400  0.793800
      Read-in Center 6535 is at   0.000000 -1.005400  0.793800
      Read-in Center 6536 is at   0.000000 -0.952500  0.793800
      Read-in Center 6537 is at   0.000000 -0.899600  0.793800
      Read-in Center 6538 is at   0.000000 -0.846700  0.793800
      Read-in Center 6539 is at   0.000000 -0.793800  0.793800
      Read-in Center 6540 is at   0.000000 -0.740800  0.793800
      Read-in Center 6541 is at   0.000000 -0.687900  0.793800
      Read-in Center 6542 is at   0.000000 -0.635000  0.793800
      Read-in Center 6543 is at   0.000000 -0.582100  0.793800
      Read-in Center 6544 is at   0.000000 -0.529200  0.793800
      Read-in Center 6545 is at   0.000000 -0.476300  0.793800
      Read-in Center 6546 is at   0.000000 -0.423300  0.793800
      Read-in Center 6547 is at   0.000000 -0.370400  0.793800
      Read-in Center 6548 is at   0.000000 -0.317500  0.793800
      Read-in Center 6549 is at   0.000000 -0.264600  0.793800
      Read-in Center 6550 is at   0.000000 -0.211700  0.793800
      Read-in Center 6551 is at   0.000000 -0.158800  0.793800
      Read-in Center 6552 is at   0.000000 -0.105800  0.793800
      Read-in Center 6553 is at   0.000000 -0.052900  0.793800
      Read-in Center 6554 is at   0.000000  0.000000  0.793800
      Read-in Center 6555 is at   0.000000  0.052900  0.793800
      Read-in Center 6556 is at   0.000000  0.105800  0.793800
      Read-in Center 6557 is at   0.000000  0.158800  0.793800
      Read-in Center 6558 is at   0.000000  0.211700  0.793800
      Read-in Center 6559 is at   0.000000  0.264600  0.793800
      Read-in Center 6560 is at   0.000000  0.317500  0.793800
      Read-in Center 6561 is at   0.000000  0.370400  0.793800
      Read-in Center 6562 is at   0.000000  0.423300  0.793800
      Read-in Center 6563 is at   0.000000  0.476300  0.793800
      Read-in Center 6564 is at   0.000000  0.529200  0.793800
      Read-in Center 6565 is at   0.000000  0.582100  0.793800
      Read-in Center 6566 is at   0.000000  0.635000  0.793800
      Read-in Center 6567 is at   0.000000  0.687900  0.793800
      Read-in Center 6568 is at   0.000000  0.740800  0.793800
      Read-in Center 6569 is at   0.000000  0.793800  0.793800
      Read-in Center 6570 is at   0.000000  0.846700  0.793800
      Read-in Center 6571 is at   0.000000  0.899600  0.793800
      Read-in Center 6572 is at   0.000000  0.952500  0.793800
      Read-in Center 6573 is at   0.000000  1.005400  0.793800
      Read-in Center 6574 is at   0.000000  1.058400  0.793800
      Read-in Center 6575 is at   0.000000  1.111300  0.793800
      Read-in Center 6576 is at   0.000000  1.164200  0.793800
      Read-in Center 6577 is at   0.000000  1.217100  0.793800
      Read-in Center 6578 is at   0.000000  1.270000  0.793800
      Read-in Center 6579 is at   0.000000  1.322900  0.793800
      Read-in Center 6580 is at   0.000000  1.375900  0.793800
      Read-in Center 6581 is at   0.000000  1.428800  0.793800
      Read-in Center 6582 is at   0.000000  1.481700  0.793800
      Read-in Center 6583 is at   0.000000  1.534600  0.793800
      Read-in Center 6584 is at   0.000000  1.587500  0.793800
      Read-in Center 6585 is at   0.000000  1.640400  0.793800
      Read-in Center 6586 is at   0.000000  1.693400  0.793800
      Read-in Center 6587 is at   0.000000  1.746300  0.793800
      Read-in Center 6588 is at   0.000000  1.799200  0.793800
      Read-in Center 6589 is at   0.000000  1.852100  0.793800
      Read-in Center 6590 is at   0.000000  1.905000  0.793800
      Read-in Center 6591 is at   0.000000  1.958000  0.793800
      Read-in Center 6592 is at   0.000000  2.010900  0.793800
      Read-in Center 6593 is at   0.000000  2.063800  0.793800
      Read-in Center 6594 is at   0.000000  2.116700  0.793800
      Read-in Center 6595 is at   0.000000  2.169600  0.793800
      Read-in Center 6596 is at   0.000000  2.222500  0.793800
      Read-in Center 6597 is at   0.000000  2.275500  0.793800
      Read-in Center 6598 is at   0.000000  2.328400  0.793800
      Read-in Center 6599 is at   0.000000  2.381300  0.793800
      Read-in Center 6600 is at   0.000000  2.434200  0.793800
      Read-in Center 6601 is at   0.000000  2.487100  0.793800
      Read-in Center 6602 is at   0.000000  2.540000  0.793800
      Read-in Center 6603 is at   0.000000  2.593000  0.793800
      Read-in Center 6604 is at   0.000000 -2.645900  0.846700
      Read-in Center 6605 is at   0.000000 -2.593000  0.846700
      Read-in Center 6606 is at   0.000000 -2.540000  0.846700
      Read-in Center 6607 is at   0.000000 -2.487100  0.846700
      Read-in Center 6608 is at   0.000000 -2.434200  0.846700
      Read-in Center 6609 is at   0.000000 -2.381300  0.846700
      Read-in Center 6610 is at   0.000000 -2.328400  0.846700
      Read-in Center 6611 is at   0.000000 -2.275500  0.846700
      Read-in Center 6612 is at   0.000000 -2.222500  0.846700
      Read-in Center 6613 is at   0.000000 -2.169600  0.846700
      Read-in Center 6614 is at   0.000000 -2.116700  0.846700
      Read-in Center 6615 is at   0.000000 -2.063800  0.846700
      Read-in Center 6616 is at   0.000000 -2.010900  0.846700
      Read-in Center 6617 is at   0.000000 -1.958000  0.846700
      Read-in Center 6618 is at   0.000000 -1.905000  0.846700
      Read-in Center 6619 is at   0.000000 -1.852100  0.846700
      Read-in Center 6620 is at   0.000000 -1.799200  0.846700
      Read-in Center 6621 is at   0.000000 -1.746300  0.846700
      Read-in Center 6622 is at   0.000000 -1.693400  0.846700
      Read-in Center 6623 is at   0.000000 -1.640400  0.846700
      Read-in Center 6624 is at   0.000000 -1.587500  0.846700
      Read-in Center 6625 is at   0.000000 -1.534600  0.846700
      Read-in Center 6626 is at   0.000000 -1.481700  0.846700
      Read-in Center 6627 is at   0.000000 -1.428800  0.846700
      Read-in Center 6628 is at   0.000000 -1.375900  0.846700
      Read-in Center 6629 is at   0.000000 -1.322900  0.846700
      Read-in Center 6630 is at   0.000000 -1.270000  0.846700
      Read-in Center 6631 is at   0.000000 -1.217100  0.846700
      Read-in Center 6632 is at   0.000000 -1.164200  0.846700
      Read-in Center 6633 is at   0.000000 -1.111300  0.846700
      Read-in Center 6634 is at   0.000000 -1.058400  0.846700
      Read-in Center 6635 is at   0.000000 -1.005400  0.846700
      Read-in Center 6636 is at   0.000000 -0.952500  0.846700
      Read-in Center 6637 is at   0.000000 -0.899600  0.846700
      Read-in Center 6638 is at   0.000000 -0.846700  0.846700
      Read-in Center 6639 is at   0.000000 -0.793800  0.846700
      Read-in Center 6640 is at   0.000000 -0.740800  0.846700
      Read-in Center 6641 is at   0.000000 -0.687900  0.846700
      Read-in Center 6642 is at   0.000000 -0.635000  0.846700
      Read-in Center 6643 is at   0.000000 -0.582100  0.846700
      Read-in Center 6644 is at   0.000000 -0.529200  0.846700
      Read-in Center 6645 is at   0.000000 -0.476300  0.846700
      Read-in Center 6646 is at   0.000000 -0.423300  0.846700
      Read-in Center 6647 is at   0.000000 -0.370400  0.846700
      Read-in Center 6648 is at   0.000000 -0.317500  0.846700
      Read-in Center 6649 is at   0.000000 -0.264600  0.846700
      Read-in Center 6650 is at   0.000000 -0.211700  0.846700
      Read-in Center 6651 is at   0.000000 -0.158800  0.846700
      Read-in Center 6652 is at   0.000000 -0.105800  0.846700
      Read-in Center 6653 is at   0.000000 -0.052900  0.846700
      Read-in Center 6654 is at   0.000000  0.000000  0.846700
      Read-in Center 6655 is at   0.000000  0.052900  0.846700
      Read-in Center 6656 is at   0.000000  0.105800  0.846700
      Read-in Center 6657 is at   0.000000  0.158800  0.846700
      Read-in Center 6658 is at   0.000000  0.211700  0.846700
      Read-in Center 6659 is at   0.000000  0.264600  0.846700
      Read-in Center 6660 is at   0.000000  0.317500  0.846700
      Read-in Center 6661 is at   0.000000  0.370400  0.846700
      Read-in Center 6662 is at   0.000000  0.423300  0.846700
      Read-in Center 6663 is at   0.000000  0.476300  0.846700
      Read-in Center 6664 is at   0.000000  0.529200  0.846700
      Read-in Center 6665 is at   0.000000  0.582100  0.846700
      Read-in Center 6666 is at   0.000000  0.635000  0.846700
      Read-in Center 6667 is at   0.000000  0.687900  0.846700
      Read-in Center 6668 is at   0.000000  0.740800  0.846700
      Read-in Center 6669 is at   0.000000  0.793800  0.846700
      Read-in Center 6670 is at   0.000000  0.846700  0.846700
      Read-in Center 6671 is at   0.000000  0.899600  0.846700
      Read-in Center 6672 is at   0.000000  0.952500  0.846700
      Read-in Center 6673 is at   0.000000  1.005400  0.846700
      Read-in Center 6674 is at   0.000000  1.058400  0.846700
      Read-in Center 6675 is at   0.000000  1.111300  0.846700
      Read-in Center 6676 is at   0.000000  1.164200  0.846700
      Read-in Center 6677 is at   0.000000  1.217100  0.846700
      Read-in Center 6678 is at   0.000000  1.270000  0.846700
      Read-in Center 6679 is at   0.000000  1.322900  0.846700
      Read-in Center 6680 is at   0.000000  1.375900  0.846700
      Read-in Center 6681 is at   0.000000  1.428800  0.846700
      Read-in Center 6682 is at   0.000000  1.481700  0.846700
      Read-in Center 6683 is at   0.000000  1.534600  0.846700
      Read-in Center 6684 is at   0.000000  1.587500  0.846700
      Read-in Center 6685 is at   0.000000  1.640400  0.846700
      Read-in Center 6686 is at   0.000000  1.693400  0.846700
      Read-in Center 6687 is at   0.000000  1.746300  0.846700
      Read-in Center 6688 is at   0.000000  1.799200  0.846700
      Read-in Center 6689 is at   0.000000  1.852100  0.846700
      Read-in Center 6690 is at   0.000000  1.905000  0.846700
      Read-in Center 6691 is at   0.000000  1.958000  0.846700
      Read-in Center 6692 is at   0.000000  2.010900  0.846700
      Read-in Center 6693 is at   0.000000  2.063800  0.846700
      Read-in Center 6694 is at   0.000000  2.116700  0.846700
      Read-in Center 6695 is at   0.000000  2.169600  0.846700
      Read-in Center 6696 is at   0.000000  2.222500  0.846700
      Read-in Center 6697 is at   0.000000  2.275500  0.846700
      Read-in Center 6698 is at   0.000000  2.328400  0.846700
      Read-in Center 6699 is at   0.000000  2.381300  0.846700
      Read-in Center 6700 is at   0.000000  2.434200  0.846700
      Read-in Center 6701 is at   0.000000  2.487100  0.846700
      Read-in Center 6702 is at   0.000000  2.540000  0.846700
      Read-in Center 6703 is at   0.000000  2.593000  0.846700
      Read-in Center 6704 is at   0.000000 -2.645900  0.899600
      Read-in Center 6705 is at   0.000000 -2.593000  0.899600
      Read-in Center 6706 is at   0.000000 -2.540000  0.899600
      Read-in Center 6707 is at   0.000000 -2.487100  0.899600
      Read-in Center 6708 is at   0.000000 -2.434200  0.899600
      Read-in Center 6709 is at   0.000000 -2.381300  0.899600
      Read-in Center 6710 is at   0.000000 -2.328400  0.899600
      Read-in Center 6711 is at   0.000000 -2.275500  0.899600
      Read-in Center 6712 is at   0.000000 -2.222500  0.899600
      Read-in Center 6713 is at   0.000000 -2.169600  0.899600
      Read-in Center 6714 is at   0.000000 -2.116700  0.899600
      Read-in Center 6715 is at   0.000000 -2.063800  0.899600
      Read-in Center 6716 is at   0.000000 -2.010900  0.899600
      Read-in Center 6717 is at   0.000000 -1.958000  0.899600
      Read-in Center 6718 is at   0.000000 -1.905000  0.899600
      Read-in Center 6719 is at   0.000000 -1.852100  0.899600
      Read-in Center 6720 is at   0.000000 -1.799200  0.899600
      Read-in Center 6721 is at   0.000000 -1.746300  0.899600
      Read-in Center 6722 is at   0.000000 -1.693400  0.899600
      Read-in Center 6723 is at   0.000000 -1.640400  0.899600
      Read-in Center 6724 is at   0.000000 -1.587500  0.899600
      Read-in Center 6725 is at   0.000000 -1.534600  0.899600
      Read-in Center 6726 is at   0.000000 -1.481700  0.899600
      Read-in Center 6727 is at   0.000000 -1.428800  0.899600
      Read-in Center 6728 is at   0.000000 -1.375900  0.899600
      Read-in Center 6729 is at   0.000000 -1.322900  0.899600
      Read-in Center 6730 is at   0.000000 -1.270000  0.899600
      Read-in Center 6731 is at   0.000000 -1.217100  0.899600
      Read-in Center 6732 is at   0.000000 -1.164200  0.899600
      Read-in Center 6733 is at   0.000000 -1.111300  0.899600
      Read-in Center 6734 is at   0.000000 -1.058400  0.899600
      Read-in Center 6735 is at   0.000000 -1.005400  0.899600
      Read-in Center 6736 is at   0.000000 -0.952500  0.899600
      Read-in Center 6737 is at   0.000000 -0.899600  0.899600
      Read-in Center 6738 is at   0.000000 -0.846700  0.899600
      Read-in Center 6739 is at   0.000000 -0.793800  0.899600
      Read-in Center 6740 is at   0.000000 -0.740800  0.899600
      Read-in Center 6741 is at   0.000000 -0.687900  0.899600
      Read-in Center 6742 is at   0.000000 -0.635000  0.899600
      Read-in Center 6743 is at   0.000000 -0.582100  0.899600
      Read-in Center 6744 is at   0.000000 -0.529200  0.899600
      Read-in Center 6745 is at   0.000000 -0.476300  0.899600
      Read-in Center 6746 is at   0.000000 -0.423300  0.899600
      Read-in Center 6747 is at   0.000000 -0.370400  0.899600
      Read-in Center 6748 is at   0.000000 -0.317500  0.899600
      Read-in Center 6749 is at   0.000000 -0.264600  0.899600
      Read-in Center 6750 is at   0.000000 -0.211700  0.899600
      Read-in Center 6751 is at   0.000000 -0.158800  0.899600
      Read-in Center 6752 is at   0.000000 -0.105800  0.899600
      Read-in Center 6753 is at   0.000000 -0.052900  0.899600
      Read-in Center 6754 is at   0.000000  0.000000  0.899600
      Read-in Center 6755 is at   0.000000  0.052900  0.899600
      Read-in Center 6756 is at   0.000000  0.105800  0.899600
      Read-in Center 6757 is at   0.000000  0.158800  0.899600
      Read-in Center 6758 is at   0.000000  0.211700  0.899600
      Read-in Center 6759 is at   0.000000  0.264600  0.899600
      Read-in Center 6760 is at   0.000000  0.317500  0.899600
      Read-in Center 6761 is at   0.000000  0.370400  0.899600
      Read-in Center 6762 is at   0.000000  0.423300  0.899600
      Read-in Center 6763 is at   0.000000  0.476300  0.899600
      Read-in Center 6764 is at   0.000000  0.529200  0.899600
      Read-in Center 6765 is at   0.000000  0.582100  0.899600
      Read-in Center 6766 is at   0.000000  0.635000  0.899600
      Read-in Center 6767 is at   0.000000  0.687900  0.899600
      Read-in Center 6768 is at   0.000000  0.740800  0.899600
      Read-in Center 6769 is at   0.000000  0.793800  0.899600
      Read-in Center 6770 is at   0.000000  0.846700  0.899600
      Read-in Center 6771 is at   0.000000  0.899600  0.899600
      Read-in Center 6772 is at   0.000000  0.952500  0.899600
      Read-in Center 6773 is at   0.000000  1.005400  0.899600
      Read-in Center 6774 is at   0.000000  1.058400  0.899600
      Read-in Center 6775 is at   0.000000  1.111300  0.899600
      Read-in Center 6776 is at   0.000000  1.164200  0.899600
      Read-in Center 6777 is at   0.000000  1.217100  0.899600
      Read-in Center 6778 is at   0.000000  1.270000  0.899600
      Read-in Center 6779 is at   0.000000  1.322900  0.899600
      Read-in Center 6780 is at   0.000000  1.375900  0.899600
      Read-in Center 6781 is at   0.000000  1.428800  0.899600
      Read-in Center 6782 is at   0.000000  1.481700  0.899600
      Read-in Center 6783 is at   0.000000  1.534600  0.899600
      Read-in Center 6784 is at   0.000000  1.587500  0.899600
      Read-in Center 6785 is at   0.000000  1.640400  0.899600
      Read-in Center 6786 is at   0.000000  1.693400  0.899600
      Read-in Center 6787 is at   0.000000  1.746300  0.899600
      Read-in Center 6788 is at   0.000000  1.799200  0.899600
      Read-in Center 6789 is at   0.000000  1.852100  0.899600
      Read-in Center 6790 is at   0.000000  1.905000  0.899600
      Read-in Center 6791 is at   0.000000  1.958000  0.899600
      Read-in Center 6792 is at   0.000000  2.010900  0.899600
      Read-in Center 6793 is at   0.000000  2.063800  0.899600
      Read-in Center 6794 is at   0.000000  2.116700  0.899600
      Read-in Center 6795 is at   0.000000  2.169600  0.899600
      Read-in Center 6796 is at   0.000000  2.222500  0.899600
      Read-in Center 6797 is at   0.000000  2.275500  0.899600
      Read-in Center 6798 is at   0.000000  2.328400  0.899600
      Read-in Center 6799 is at   0.000000  2.381300  0.899600
      Read-in Center 6800 is at   0.000000  2.434200  0.899600
      Read-in Center 6801 is at   0.000000  2.487100  0.899600
      Read-in Center 6802 is at   0.000000  2.540000  0.899600
      Read-in Center 6803 is at   0.000000  2.593000  0.899600
      Read-in Center 6804 is at   0.000000 -2.645900  0.952500
      Read-in Center 6805 is at   0.000000 -2.593000  0.952500
      Read-in Center 6806 is at   0.000000 -2.540000  0.952500
      Read-in Center 6807 is at   0.000000 -2.487100  0.952500
      Read-in Center 6808 is at   0.000000 -2.434200  0.952500
      Read-in Center 6809 is at   0.000000 -2.381300  0.952500
      Read-in Center 6810 is at   0.000000 -2.328400  0.952500
      Read-in Center 6811 is at   0.000000 -2.275500  0.952500
      Read-in Center 6812 is at   0.000000 -2.222500  0.952500
      Read-in Center 6813 is at   0.000000 -2.169600  0.952500
      Read-in Center 6814 is at   0.000000 -2.116700  0.952500
      Read-in Center 6815 is at   0.000000 -2.063800  0.952500
      Read-in Center 6816 is at   0.000000 -2.010900  0.952500
      Read-in Center 6817 is at   0.000000 -1.958000  0.952500
      Read-in Center 6818 is at   0.000000 -1.905000  0.952500
      Read-in Center 6819 is at   0.000000 -1.852100  0.952500
      Read-in Center 6820 is at   0.000000 -1.799200  0.952500
      Read-in Center 6821 is at   0.000000 -1.746300  0.952500
      Read-in Center 6822 is at   0.000000 -1.693400  0.952500
      Read-in Center 6823 is at   0.000000 -1.640400  0.952500
      Read-in Center 6824 is at   0.000000 -1.587500  0.952500
      Read-in Center 6825 is at   0.000000 -1.534600  0.952500
      Read-in Center 6826 is at   0.000000 -1.481700  0.952500
      Read-in Center 6827 is at   0.000000 -1.428800  0.952500
      Read-in Center 6828 is at   0.000000 -1.375900  0.952500
      Read-in Center 6829 is at   0.000000 -1.322900  0.952500
      Read-in Center 6830 is at   0.000000 -1.270000  0.952500
      Read-in Center 6831 is at   0.000000 -1.217100  0.952500
      Read-in Center 6832 is at   0.000000 -1.164200  0.952500
      Read-in Center 6833 is at   0.000000 -1.111300  0.952500
      Read-in Center 6834 is at   0.000000 -1.058400  0.952500
      Read-in Center 6835 is at   0.000000 -1.005400  0.952500
      Read-in Center 6836 is at   0.000000 -0.952500  0.952500
      Read-in Center 6837 is at   0.000000 -0.899600  0.952500
      Read-in Center 6838 is at   0.000000 -0.846700  0.952500
      Read-in Center 6839 is at   0.000000 -0.793800  0.952500
      Read-in Center 6840 is at   0.000000 -0.740800  0.952500
      Read-in Center 6841 is at   0.000000 -0.687900  0.952500
      Read-in Center 6842 is at   0.000000 -0.635000  0.952500
      Read-in Center 6843 is at   0.000000 -0.582100  0.952500
      Read-in Center 6844 is at   0.000000 -0.529200  0.952500
      Read-in Center 6845 is at   0.000000 -0.476300  0.952500
      Read-in Center 6846 is at   0.000000 -0.423300  0.952500
      Read-in Center 6847 is at   0.000000 -0.370400  0.952500
      Read-in Center 6848 is at   0.000000 -0.317500  0.952500
      Read-in Center 6849 is at   0.000000 -0.264600  0.952500
      Read-in Center 6850 is at   0.000000 -0.211700  0.952500
      Read-in Center 6851 is at   0.000000 -0.158800  0.952500
      Read-in Center 6852 is at   0.000000 -0.105800  0.952500
      Read-in Center 6853 is at   0.000000 -0.052900  0.952500
      Read-in Center 6854 is at   0.000000  0.000000  0.952500
      Read-in Center 6855 is at   0.000000  0.052900  0.952500
      Read-in Center 6856 is at   0.000000  0.105800  0.952500
      Read-in Center 6857 is at   0.000000  0.158800  0.952500
      Read-in Center 6858 is at   0.000000  0.211700  0.952500
      Read-in Center 6859 is at   0.000000  0.264600  0.952500
      Read-in Center 6860 is at   0.000000  0.317500  0.952500
      Read-in Center 6861 is at   0.000000  0.370400  0.952500
      Read-in Center 6862 is at   0.000000  0.423300  0.952500
      Read-in Center 6863 is at   0.000000  0.476300  0.952500
      Read-in Center 6864 is at   0.000000  0.529200  0.952500
      Read-in Center 6865 is at   0.000000  0.582100  0.952500
      Read-in Center 6866 is at   0.000000  0.635000  0.952500
      Read-in Center 6867 is at   0.000000  0.687900  0.952500
      Read-in Center 6868 is at   0.000000  0.740800  0.952500
      Read-in Center 6869 is at   0.000000  0.793800  0.952500
      Read-in Center 6870 is at   0.000000  0.846700  0.952500
      Read-in Center 6871 is at   0.000000  0.899600  0.952500
      Read-in Center 6872 is at   0.000000  0.952500  0.952500
      Read-in Center 6873 is at   0.000000  1.005400  0.952500
      Read-in Center 6874 is at   0.000000  1.058400  0.952500
      Read-in Center 6875 is at   0.000000  1.111300  0.952500
      Read-in Center 6876 is at   0.000000  1.164200  0.952500
      Read-in Center 6877 is at   0.000000  1.217100  0.952500
      Read-in Center 6878 is at   0.000000  1.270000  0.952500
      Read-in Center 6879 is at   0.000000  1.322900  0.952500
      Read-in Center 6880 is at   0.000000  1.375900  0.952500
      Read-in Center 6881 is at   0.000000  1.428800  0.952500
      Read-in Center 6882 is at   0.000000  1.481700  0.952500
      Read-in Center 6883 is at   0.000000  1.534600  0.952500
      Read-in Center 6884 is at   0.000000  1.587500  0.952500
      Read-in Center 6885 is at   0.000000  1.640400  0.952500
      Read-in Center 6886 is at   0.000000  1.693400  0.952500
      Read-in Center 6887 is at   0.000000  1.746300  0.952500
      Read-in Center 6888 is at   0.000000  1.799200  0.952500
      Read-in Center 6889 is at   0.000000  1.852100  0.952500
      Read-in Center 6890 is at   0.000000  1.905000  0.952500
      Read-in Center 6891 is at   0.000000  1.958000  0.952500
      Read-in Center 6892 is at   0.000000  2.010900  0.952500
      Read-in Center 6893 is at   0.000000  2.063800  0.952500
      Read-in Center 6894 is at   0.000000  2.116700  0.952500
      Read-in Center 6895 is at   0.000000  2.169600  0.952500
      Read-in Center 6896 is at   0.000000  2.222500  0.952500
      Read-in Center 6897 is at   0.000000  2.275500  0.952500
      Read-in Center 6898 is at   0.000000  2.328400  0.952500
      Read-in Center 6899 is at   0.000000  2.381300  0.952500
      Read-in Center 6900 is at   0.000000  2.434200  0.952500
      Read-in Center 6901 is at   0.000000  2.487100  0.952500
      Read-in Center 6902 is at   0.000000  2.540000  0.952500
      Read-in Center 6903 is at   0.000000  2.593000  0.952500
      Read-in Center 6904 is at   0.000000 -2.645900  1.005400
      Read-in Center 6905 is at   0.000000 -2.593000  1.005400
      Read-in Center 6906 is at   0.000000 -2.540000  1.005400
      Read-in Center 6907 is at   0.000000 -2.487100  1.005400
      Read-in Center 6908 is at   0.000000 -2.434200  1.005400
      Read-in Center 6909 is at   0.000000 -2.381300  1.005400
      Read-in Center 6910 is at   0.000000 -2.328400  1.005400
      Read-in Center 6911 is at   0.000000 -2.275500  1.005400
      Read-in Center 6912 is at   0.000000 -2.222500  1.005400
      Read-in Center 6913 is at   0.000000 -2.169600  1.005400
      Read-in Center 6914 is at   0.000000 -2.116700  1.005400
      Read-in Center 6915 is at   0.000000 -2.063800  1.005400
      Read-in Center 6916 is at   0.000000 -2.010900  1.005400
      Read-in Center 6917 is at   0.000000 -1.958000  1.005400
      Read-in Center 6918 is at   0.000000 -1.905000  1.005400
      Read-in Center 6919 is at   0.000000 -1.852100  1.005400
      Read-in Center 6920 is at   0.000000 -1.799200  1.005400
      Read-in Center 6921 is at   0.000000 -1.746300  1.005400
      Read-in Center 6922 is at   0.000000 -1.693400  1.005400
      Read-in Center 6923 is at   0.000000 -1.640400  1.005400
      Read-in Center 6924 is at   0.000000 -1.587500  1.005400
      Read-in Center 6925 is at   0.000000 -1.534600  1.005400
      Read-in Center 6926 is at   0.000000 -1.481700  1.005400
      Read-in Center 6927 is at   0.000000 -1.428800  1.005400
      Read-in Center 6928 is at   0.000000 -1.375900  1.005400
      Read-in Center 6929 is at   0.000000 -1.322900  1.005400
      Read-in Center 6930 is at   0.000000 -1.270000  1.005400
      Read-in Center 6931 is at   0.000000 -1.217100  1.005400
      Read-in Center 6932 is at   0.000000 -1.164200  1.005400
      Read-in Center 6933 is at   0.000000 -1.111300  1.005400
      Read-in Center 6934 is at   0.000000 -1.058400  1.005400
      Read-in Center 6935 is at   0.000000 -1.005400  1.005400
      Read-in Center 6936 is at   0.000000 -0.952500  1.005400
      Read-in Center 6937 is at   0.000000 -0.899600  1.005400
      Read-in Center 6938 is at   0.000000 -0.846700  1.005400
      Read-in Center 6939 is at   0.000000 -0.793800  1.005400
      Read-in Center 6940 is at   0.000000 -0.740800  1.005400
      Read-in Center 6941 is at   0.000000 -0.687900  1.005400
      Read-in Center 6942 is at   0.000000 -0.635000  1.005400
      Read-in Center 6943 is at   0.000000 -0.582100  1.005400
      Read-in Center 6944 is at   0.000000 -0.529200  1.005400
      Read-in Center 6945 is at   0.000000 -0.476300  1.005400
      Read-in Center 6946 is at   0.000000 -0.423300  1.005400
      Read-in Center 6947 is at   0.000000 -0.370400  1.005400
      Read-in Center 6948 is at   0.000000 -0.317500  1.005400
      Read-in Center 6949 is at   0.000000 -0.264600  1.005400
      Read-in Center 6950 is at   0.000000 -0.211700  1.005400
      Read-in Center 6951 is at   0.000000 -0.158800  1.005400
      Read-in Center 6952 is at   0.000000 -0.105800  1.005400
      Read-in Center 6953 is at   0.000000 -0.052900  1.005400
      Read-in Center 6954 is at   0.000000  0.000000  1.005400
      Read-in Center 6955 is at   0.000000  0.052900  1.005400
      Read-in Center 6956 is at   0.000000  0.105800  1.005400
      Read-in Center 6957 is at   0.000000  0.158800  1.005400
      Read-in Center 6958 is at   0.000000  0.211700  1.005400
      Read-in Center 6959 is at   0.000000  0.264600  1.005400
      Read-in Center 6960 is at   0.000000  0.317500  1.005400
      Read-in Center 6961 is at   0.000000  0.370400  1.005400
      Read-in Center 6962 is at   0.000000  0.423300  1.005400
      Read-in Center 6963 is at   0.000000  0.476300  1.005400
      Read-in Center 6964 is at   0.000000  0.529200  1.005400
      Read-in Center 6965 is at   0.000000  0.582100  1.005400
      Read-in Center 6966 is at   0.000000  0.635000  1.005400
      Read-in Center 6967 is at   0.000000  0.687900  1.005400
      Read-in Center 6968 is at   0.000000  0.740800  1.005400
      Read-in Center 6969 is at   0.000000  0.793800  1.005400
      Read-in Center 6970 is at   0.000000  0.846700  1.005400
      Read-in Center 6971 is at   0.000000  0.899600  1.005400
      Read-in Center 6972 is at   0.000000  0.952500  1.005400
      Read-in Center 6973 is at   0.000000  1.005400  1.005400
      Read-in Center 6974 is at   0.000000  1.058400  1.005400
      Read-in Center 6975 is at   0.000000  1.111300  1.005400
      Read-in Center 6976 is at   0.000000  1.164200  1.005400
      Read-in Center 6977 is at   0.000000  1.217100  1.005400
      Read-in Center 6978 is at   0.000000  1.270000  1.005400
      Read-in Center 6979 is at   0.000000  1.322900  1.005400
      Read-in Center 6980 is at   0.000000  1.375900  1.005400
      Read-in Center 6981 is at   0.000000  1.428800  1.005400
      Read-in Center 6982 is at   0.000000  1.481700  1.005400
      Read-in Center 6983 is at   0.000000  1.534600  1.005400
      Read-in Center 6984 is at   0.000000  1.587500  1.005400
      Read-in Center 6985 is at   0.000000  1.640400  1.005400
      Read-in Center 6986 is at   0.000000  1.693400  1.005400
      Read-in Center 6987 is at   0.000000  1.746300  1.005400
      Read-in Center 6988 is at   0.000000  1.799200  1.005400
      Read-in Center 6989 is at   0.000000  1.852100  1.005400
      Read-in Center 6990 is at   0.000000  1.905000  1.005400
      Read-in Center 6991 is at   0.000000  1.958000  1.005400
      Read-in Center 6992 is at   0.000000  2.010900  1.005400
      Read-in Center 6993 is at   0.000000  2.063800  1.005400
      Read-in Center 6994 is at   0.000000  2.116700  1.005400
      Read-in Center 6995 is at   0.000000  2.169600  1.005400
      Read-in Center 6996 is at   0.000000  2.222500  1.005400
      Read-in Center 6997 is at   0.000000  2.275500  1.005400
      Read-in Center 6998 is at   0.000000  2.328400  1.005400
      Read-in Center 6999 is at   0.000000  2.381300  1.005400
      Read-in Center 7000 is at   0.000000  2.434200  1.005400
      Read-in Center 7001 is at   0.000000  2.487100  1.005400
      Read-in Center 7002 is at   0.000000  2.540000  1.005400
      Read-in Center 7003 is at   0.000000  2.593000  1.005400
      Read-in Center 7004 is at   0.000000 -2.645900  1.058400
      Read-in Center 7005 is at   0.000000 -2.593000  1.058400
      Read-in Center 7006 is at   0.000000 -2.540000  1.058400
      Read-in Center 7007 is at   0.000000 -2.487100  1.058400
      Read-in Center 7008 is at   0.000000 -2.434200  1.058400
      Read-in Center 7009 is at   0.000000 -2.381300  1.058400
      Read-in Center 7010 is at   0.000000 -2.328400  1.058400
      Read-in Center 7011 is at   0.000000 -2.275500  1.058400
      Read-in Center 7012 is at   0.000000 -2.222500  1.058400
      Read-in Center 7013 is at   0.000000 -2.169600  1.058400
      Read-in Center 7014 is at   0.000000 -2.116700  1.058400
      Read-in Center 7015 is at   0.000000 -2.063800  1.058400
      Read-in Center 7016 is at   0.000000 -2.010900  1.058400
      Read-in Center 7017 is at   0.000000 -1.958000  1.058400
      Read-in Center 7018 is at   0.000000 -1.905000  1.058400
      Read-in Center 7019 is at   0.000000 -1.852100  1.058400
      Read-in Center 7020 is at   0.000000 -1.799200  1.058400
      Read-in Center 7021 is at   0.000000 -1.746300  1.058400
      Read-in Center 7022 is at   0.000000 -1.693400  1.058400
      Read-in Center 7023 is at   0.000000 -1.640400  1.058400
      Read-in Center 7024 is at   0.000000 -1.587500  1.058400
      Read-in Center 7025 is at   0.000000 -1.534600  1.058400
      Read-in Center 7026 is at   0.000000 -1.481700  1.058400
      Read-in Center 7027 is at   0.000000 -1.428800  1.058400
      Read-in Center 7028 is at   0.000000 -1.375900  1.058400
      Read-in Center 7029 is at   0.000000 -1.322900  1.058400
      Read-in Center 7030 is at   0.000000 -1.270000  1.058400
      Read-in Center 7031 is at   0.000000 -1.217100  1.058400
      Read-in Center 7032 is at   0.000000 -1.164200  1.058400
      Read-in Center 7033 is at   0.000000 -1.111300  1.058400
      Read-in Center 7034 is at   0.000000 -1.058400  1.058400
      Read-in Center 7035 is at   0.000000 -1.005400  1.058400
      Read-in Center 7036 is at   0.000000 -0.952500  1.058400
      Read-in Center 7037 is at   0.000000 -0.899600  1.058400
      Read-in Center 7038 is at   0.000000 -0.846700  1.058400
      Read-in Center 7039 is at   0.000000 -0.793800  1.058400
      Read-in Center 7040 is at   0.000000 -0.740800  1.058400
      Read-in Center 7041 is at   0.000000 -0.687900  1.058400
      Read-in Center 7042 is at   0.000000 -0.635000  1.058400
      Read-in Center 7043 is at   0.000000 -0.582100  1.058400
      Read-in Center 7044 is at   0.000000 -0.529200  1.058400
      Read-in Center 7045 is at   0.000000 -0.476300  1.058400
      Read-in Center 7046 is at   0.000000 -0.423300  1.058400
      Read-in Center 7047 is at   0.000000 -0.370400  1.058400
      Read-in Center 7048 is at   0.000000 -0.317500  1.058400
      Read-in Center 7049 is at   0.000000 -0.264600  1.058400
      Read-in Center 7050 is at   0.000000 -0.211700  1.058400
      Read-in Center 7051 is at   0.000000 -0.158800  1.058400
      Read-in Center 7052 is at   0.000000 -0.105800  1.058400
      Read-in Center 7053 is at   0.000000 -0.052900  1.058400
      Read-in Center 7054 is at   0.000000  0.000000  1.058400
      Read-in Center 7055 is at   0.000000  0.052900  1.058400
      Read-in Center 7056 is at   0.000000  0.105800  1.058400
      Read-in Center 7057 is at   0.000000  0.158800  1.058400
      Read-in Center 7058 is at   0.000000  0.211700  1.058400
      Read-in Center 7059 is at   0.000000  0.264600  1.058400
      Read-in Center 7060 is at   0.000000  0.317500  1.058400
      Read-in Center 7061 is at   0.000000  0.370400  1.058400
      Read-in Center 7062 is at   0.000000  0.423300  1.058400
      Read-in Center 7063 is at   0.000000  0.476300  1.058400
      Read-in Center 7064 is at   0.000000  0.529200  1.058400
      Read-in Center 7065 is at   0.000000  0.582100  1.058400
      Read-in Center 7066 is at   0.000000  0.635000  1.058400
      Read-in Center 7067 is at   0.000000  0.687900  1.058400
      Read-in Center 7068 is at   0.000000  0.740800  1.058400
      Read-in Center 7069 is at   0.000000  0.793800  1.058400
      Read-in Center 7070 is at   0.000000  0.846700  1.058400
      Read-in Center 7071 is at   0.000000  0.899600  1.058400
      Read-in Center 7072 is at   0.000000  0.952500  1.058400
      Read-in Center 7073 is at   0.000000  1.005400  1.058400
      Read-in Center 7074 is at   0.000000  1.058400  1.058400
      Read-in Center 7075 is at   0.000000  1.111300  1.058400
      Read-in Center 7076 is at   0.000000  1.164200  1.058400
      Read-in Center 7077 is at   0.000000  1.217100  1.058400
      Read-in Center 7078 is at   0.000000  1.270000  1.058400
      Read-in Center 7079 is at   0.000000  1.322900  1.058400
      Read-in Center 7080 is at   0.000000  1.375900  1.058400
      Read-in Center 7081 is at   0.000000  1.428800  1.058400
      Read-in Center 7082 is at   0.000000  1.481700  1.058400
      Read-in Center 7083 is at   0.000000  1.534600  1.058400
      Read-in Center 7084 is at   0.000000  1.587500  1.058400
      Read-in Center 7085 is at   0.000000  1.640400  1.058400
      Read-in Center 7086 is at   0.000000  1.693400  1.058400
      Read-in Center 7087 is at   0.000000  1.746300  1.058400
      Read-in Center 7088 is at   0.000000  1.799200  1.058400
      Read-in Center 7089 is at   0.000000  1.852100  1.058400
      Read-in Center 7090 is at   0.000000  1.905000  1.058400
      Read-in Center 7091 is at   0.000000  1.958000  1.058400
      Read-in Center 7092 is at   0.000000  2.010900  1.058400
      Read-in Center 7093 is at   0.000000  2.063800  1.058400
      Read-in Center 7094 is at   0.000000  2.116700  1.058400
      Read-in Center 7095 is at   0.000000  2.169600  1.058400
      Read-in Center 7096 is at   0.000000  2.222500  1.058400
      Read-in Center 7097 is at   0.000000  2.275500  1.058400
      Read-in Center 7098 is at   0.000000  2.328400  1.058400
      Read-in Center 7099 is at   0.000000  2.381300  1.058400
      Read-in Center 7100 is at   0.000000  2.434200  1.058400
      Read-in Center 7101 is at   0.000000  2.487100  1.058400
      Read-in Center 7102 is at   0.000000  2.540000  1.058400
      Read-in Center 7103 is at   0.000000  2.593000  1.058400
      Read-in Center 7104 is at   0.000000 -2.645900  1.111300
      Read-in Center 7105 is at   0.000000 -2.593000  1.111300
      Read-in Center 7106 is at   0.000000 -2.540000  1.111300
      Read-in Center 7107 is at   0.000000 -2.487100  1.111300
      Read-in Center 7108 is at   0.000000 -2.434200  1.111300
      Read-in Center 7109 is at   0.000000 -2.381300  1.111300
      Read-in Center 7110 is at   0.000000 -2.328400  1.111300
      Read-in Center 7111 is at   0.000000 -2.275500  1.111300
      Read-in Center 7112 is at   0.000000 -2.222500  1.111300
      Read-in Center 7113 is at   0.000000 -2.169600  1.111300
      Read-in Center 7114 is at   0.000000 -2.116700  1.111300
      Read-in Center 7115 is at   0.000000 -2.063800  1.111300
      Read-in Center 7116 is at   0.000000 -2.010900  1.111300
      Read-in Center 7117 is at   0.000000 -1.958000  1.111300
      Read-in Center 7118 is at   0.000000 -1.905000  1.111300
      Read-in Center 7119 is at   0.000000 -1.852100  1.111300
      Read-in Center 7120 is at   0.000000 -1.799200  1.111300
      Read-in Center 7121 is at   0.000000 -1.746300  1.111300
      Read-in Center 7122 is at   0.000000 -1.693400  1.111300
      Read-in Center 7123 is at   0.000000 -1.640400  1.111300
      Read-in Center 7124 is at   0.000000 -1.587500  1.111300
      Read-in Center 7125 is at   0.000000 -1.534600  1.111300
      Read-in Center 7126 is at   0.000000 -1.481700  1.111300
      Read-in Center 7127 is at   0.000000 -1.428800  1.111300
      Read-in Center 7128 is at   0.000000 -1.375900  1.111300
      Read-in Center 7129 is at   0.000000 -1.322900  1.111300
      Read-in Center 7130 is at   0.000000 -1.270000  1.111300
      Read-in Center 7131 is at   0.000000 -1.217100  1.111300
      Read-in Center 7132 is at   0.000000 -1.164200  1.111300
      Read-in Center 7133 is at   0.000000 -1.111300  1.111300
      Read-in Center 7134 is at   0.000000 -1.058400  1.111300
      Read-in Center 7135 is at   0.000000 -1.005400  1.111300
      Read-in Center 7136 is at   0.000000 -0.952500  1.111300
      Read-in Center 7137 is at   0.000000 -0.899600  1.111300
      Read-in Center 7138 is at   0.000000 -0.846700  1.111300
      Read-in Center 7139 is at   0.000000 -0.793800  1.111300
      Read-in Center 7140 is at   0.000000 -0.740800  1.111300
      Read-in Center 7141 is at   0.000000 -0.687900  1.111300
      Read-in Center 7142 is at   0.000000 -0.635000  1.111300
      Read-in Center 7143 is at   0.000000 -0.582100  1.111300
      Read-in Center 7144 is at   0.000000 -0.529200  1.111300
      Read-in Center 7145 is at   0.000000 -0.476300  1.111300
      Read-in Center 7146 is at   0.000000 -0.423300  1.111300
      Read-in Center 7147 is at   0.000000 -0.370400  1.111300
      Read-in Center 7148 is at   0.000000 -0.317500  1.111300
      Read-in Center 7149 is at   0.000000 -0.264600  1.111300
      Read-in Center 7150 is at   0.000000 -0.211700  1.111300
      Read-in Center 7151 is at   0.000000 -0.158800  1.111300
      Read-in Center 7152 is at   0.000000 -0.105800  1.111300
      Read-in Center 7153 is at   0.000000 -0.052900  1.111300
      Read-in Center 7154 is at   0.000000  0.000000  1.111300
      Read-in Center 7155 is at   0.000000  0.052900  1.111300
      Read-in Center 7156 is at   0.000000  0.105800  1.111300
      Read-in Center 7157 is at   0.000000  0.158800  1.111300
      Read-in Center 7158 is at   0.000000  0.211700  1.111300
      Read-in Center 7159 is at   0.000000  0.264600  1.111300
      Read-in Center 7160 is at   0.000000  0.317500  1.111300
      Read-in Center 7161 is at   0.000000  0.370400  1.111300
      Read-in Center 7162 is at   0.000000  0.423300  1.111300
      Read-in Center 7163 is at   0.000000  0.476300  1.111300
      Read-in Center 7164 is at   0.000000  0.529200  1.111300
      Read-in Center 7165 is at   0.000000  0.582100  1.111300
      Read-in Center 7166 is at   0.000000  0.635000  1.111300
      Read-in Center 7167 is at   0.000000  0.687900  1.111300
      Read-in Center 7168 is at   0.000000  0.740800  1.111300
      Read-in Center 7169 is at   0.000000  0.793800  1.111300
      Read-in Center 7170 is at   0.000000  0.846700  1.111300
      Read-in Center 7171 is at   0.000000  0.899600  1.111300
      Read-in Center 7172 is at   0.000000  0.952500  1.111300
      Read-in Center 7173 is at   0.000000  1.005400  1.111300
      Read-in Center 7174 is at   0.000000  1.058400  1.111300
      Read-in Center 7175 is at   0.000000  1.111300  1.111300
      Read-in Center 7176 is at   0.000000  1.164200  1.111300
      Read-in Center 7177 is at   0.000000  1.217100  1.111300
      Read-in Center 7178 is at   0.000000  1.270000  1.111300
      Read-in Center 7179 is at   0.000000  1.322900  1.111300
      Read-in Center 7180 is at   0.000000  1.375900  1.111300
      Read-in Center 7181 is at   0.000000  1.428800  1.111300
      Read-in Center 7182 is at   0.000000  1.481700  1.111300
      Read-in Center 7183 is at   0.000000  1.534600  1.111300
      Read-in Center 7184 is at   0.000000  1.587500  1.111300
      Read-in Center 7185 is at   0.000000  1.640400  1.111300
      Read-in Center 7186 is at   0.000000  1.693400  1.111300
      Read-in Center 7187 is at   0.000000  1.746300  1.111300
      Read-in Center 7188 is at   0.000000  1.799200  1.111300
      Read-in Center 7189 is at   0.000000  1.852100  1.111300
      Read-in Center 7190 is at   0.000000  1.905000  1.111300
      Read-in Center 7191 is at   0.000000  1.958000  1.111300
      Read-in Center 7192 is at   0.000000  2.010900  1.111300
      Read-in Center 7193 is at   0.000000  2.063800  1.111300
      Read-in Center 7194 is at   0.000000  2.116700  1.111300
      Read-in Center 7195 is at   0.000000  2.169600  1.111300
      Read-in Center 7196 is at   0.000000  2.222500  1.111300
      Read-in Center 7197 is at   0.000000  2.275500  1.111300
      Read-in Center 7198 is at   0.000000  2.328400  1.111300
      Read-in Center 7199 is at   0.000000  2.381300  1.111300
      Read-in Center 7200 is at   0.000000  2.434200  1.111300
      Read-in Center 7201 is at   0.000000  2.487100  1.111300
      Read-in Center 7202 is at   0.000000  2.540000  1.111300
      Read-in Center 7203 is at   0.000000  2.593000  1.111300
      Read-in Center 7204 is at   0.000000 -2.645900  1.164200
      Read-in Center 7205 is at   0.000000 -2.593000  1.164200
      Read-in Center 7206 is at   0.000000 -2.540000  1.164200
      Read-in Center 7207 is at   0.000000 -2.487100  1.164200
      Read-in Center 7208 is at   0.000000 -2.434200  1.164200
      Read-in Center 7209 is at   0.000000 -2.381300  1.164200
      Read-in Center 7210 is at   0.000000 -2.328400  1.164200
      Read-in Center 7211 is at   0.000000 -2.275500  1.164200
      Read-in Center 7212 is at   0.000000 -2.222500  1.164200
      Read-in Center 7213 is at   0.000000 -2.169600  1.164200
      Read-in Center 7214 is at   0.000000 -2.116700  1.164200
      Read-in Center 7215 is at   0.000000 -2.063800  1.164200
      Read-in Center 7216 is at   0.000000 -2.010900  1.164200
      Read-in Center 7217 is at   0.000000 -1.958000  1.164200
      Read-in Center 7218 is at   0.000000 -1.905000  1.164200
      Read-in Center 7219 is at   0.000000 -1.852100  1.164200
      Read-in Center 7220 is at   0.000000 -1.799200  1.164200
      Read-in Center 7221 is at   0.000000 -1.746300  1.164200
      Read-in Center 7222 is at   0.000000 -1.693400  1.164200
      Read-in Center 7223 is at   0.000000 -1.640400  1.164200
      Read-in Center 7224 is at   0.000000 -1.587500  1.164200
      Read-in Center 7225 is at   0.000000 -1.534600  1.164200
      Read-in Center 7226 is at   0.000000 -1.481700  1.164200
      Read-in Center 7227 is at   0.000000 -1.428800  1.164200
      Read-in Center 7228 is at   0.000000 -1.375900  1.164200
      Read-in Center 7229 is at   0.000000 -1.322900  1.164200
      Read-in Center 7230 is at   0.000000 -1.270000  1.164200
      Read-in Center 7231 is at   0.000000 -1.217100  1.164200
      Read-in Center 7232 is at   0.000000 -1.164200  1.164200
      Read-in Center 7233 is at   0.000000 -1.111300  1.164200
      Read-in Center 7234 is at   0.000000 -1.058400  1.164200
      Read-in Center 7235 is at   0.000000 -1.005400  1.164200
      Read-in Center 7236 is at   0.000000 -0.952500  1.164200
      Read-in Center 7237 is at   0.000000 -0.899600  1.164200
      Read-in Center 7238 is at   0.000000 -0.846700  1.164200
      Read-in Center 7239 is at   0.000000 -0.793800  1.164200
      Read-in Center 7240 is at   0.000000 -0.740800  1.164200
      Read-in Center 7241 is at   0.000000 -0.687900  1.164200
      Read-in Center 7242 is at   0.000000 -0.635000  1.164200
      Read-in Center 7243 is at   0.000000 -0.582100  1.164200
      Read-in Center 7244 is at   0.000000 -0.529200  1.164200
      Read-in Center 7245 is at   0.000000 -0.476300  1.164200
      Read-in Center 7246 is at   0.000000 -0.423300  1.164200
      Read-in Center 7247 is at   0.000000 -0.370400  1.164200
      Read-in Center 7248 is at   0.000000 -0.317500  1.164200
      Read-in Center 7249 is at   0.000000 -0.264600  1.164200
      Read-in Center 7250 is at   0.000000 -0.211700  1.164200
      Read-in Center 7251 is at   0.000000 -0.158800  1.164200
      Read-in Center 7252 is at   0.000000 -0.105800  1.164200
      Read-in Center 7253 is at   0.000000 -0.052900  1.164200
      Read-in Center 7254 is at   0.000000  0.000000  1.164200
      Read-in Center 7255 is at   0.000000  0.052900  1.164200
      Read-in Center 7256 is at   0.000000  0.105800  1.164200
      Read-in Center 7257 is at   0.000000  0.158800  1.164200
      Read-in Center 7258 is at   0.000000  0.211700  1.164200
      Read-in Center 7259 is at   0.000000  0.264600  1.164200
      Read-in Center 7260 is at   0.000000  0.317500  1.164200
      Read-in Center 7261 is at   0.000000  0.370400  1.164200
      Read-in Center 7262 is at   0.000000  0.423300  1.164200
      Read-in Center 7263 is at   0.000000  0.476300  1.164200
      Read-in Center 7264 is at   0.000000  0.529200  1.164200
      Read-in Center 7265 is at   0.000000  0.582100  1.164200
      Read-in Center 7266 is at   0.000000  0.635000  1.164200
      Read-in Center 7267 is at   0.000000  0.687900  1.164200
      Read-in Center 7268 is at   0.000000  0.740800  1.164200
      Read-in Center 7269 is at   0.000000  0.793800  1.164200
      Read-in Center 7270 is at   0.000000  0.846700  1.164200
      Read-in Center 7271 is at   0.000000  0.899600  1.164200
      Read-in Center 7272 is at   0.000000  0.952500  1.164200
      Read-in Center 7273 is at   0.000000  1.005400  1.164200
      Read-in Center 7274 is at   0.000000  1.058400  1.164200
      Read-in Center 7275 is at   0.000000  1.111300  1.164200
      Read-in Center 7276 is at   0.000000  1.164200  1.164200
      Read-in Center 7277 is at   0.000000  1.217100  1.164200
      Read-in Center 7278 is at   0.000000  1.270000  1.164200
      Read-in Center 7279 is at   0.000000  1.322900  1.164200
      Read-in Center 7280 is at   0.000000  1.375900  1.164200
      Read-in Center 7281 is at   0.000000  1.428800  1.164200
      Read-in Center 7282 is at   0.000000  1.481700  1.164200
      Read-in Center 7283 is at   0.000000  1.534600  1.164200
      Read-in Center 7284 is at   0.000000  1.587500  1.164200
      Read-in Center 7285 is at   0.000000  1.640400  1.164200
      Read-in Center 7286 is at   0.000000  1.693400  1.164200
      Read-in Center 7287 is at   0.000000  1.746300  1.164200
      Read-in Center 7288 is at   0.000000  1.799200  1.164200
      Read-in Center 7289 is at   0.000000  1.852100  1.164200
      Read-in Center 7290 is at   0.000000  1.905000  1.164200
      Read-in Center 7291 is at   0.000000  1.958000  1.164200
      Read-in Center 7292 is at   0.000000  2.010900  1.164200
      Read-in Center 7293 is at   0.000000  2.063800  1.164200
      Read-in Center 7294 is at   0.000000  2.116700  1.164200
      Read-in Center 7295 is at   0.000000  2.169600  1.164200
      Read-in Center 7296 is at   0.000000  2.222500  1.164200
      Read-in Center 7297 is at   0.000000  2.275500  1.164200
      Read-in Center 7298 is at   0.000000  2.328400  1.164200
      Read-in Center 7299 is at   0.000000  2.381300  1.164200
      Read-in Center 7300 is at   0.000000  2.434200  1.164200
      Read-in Center 7301 is at   0.000000  2.487100  1.164200
      Read-in Center 7302 is at   0.000000  2.540000  1.164200
      Read-in Center 7303 is at   0.000000  2.593000  1.164200
      Read-in Center 7304 is at   0.000000 -2.645900  1.217100
      Read-in Center 7305 is at   0.000000 -2.593000  1.217100
      Read-in Center 7306 is at   0.000000 -2.540000  1.217100
      Read-in Center 7307 is at   0.000000 -2.487100  1.217100
      Read-in Center 7308 is at   0.000000 -2.434200  1.217100
      Read-in Center 7309 is at   0.000000 -2.381300  1.217100
      Read-in Center 7310 is at   0.000000 -2.328400  1.217100
      Read-in Center 7311 is at   0.000000 -2.275500  1.217100
      Read-in Center 7312 is at   0.000000 -2.222500  1.217100
      Read-in Center 7313 is at   0.000000 -2.169600  1.217100
      Read-in Center 7314 is at   0.000000 -2.116700  1.217100
      Read-in Center 7315 is at   0.000000 -2.063800  1.217100
      Read-in Center 7316 is at   0.000000 -2.010900  1.217100
      Read-in Center 7317 is at   0.000000 -1.958000  1.217100
      Read-in Center 7318 is at   0.000000 -1.905000  1.217100
      Read-in Center 7319 is at   0.000000 -1.852100  1.217100
      Read-in Center 7320 is at   0.000000 -1.799200  1.217100
      Read-in Center 7321 is at   0.000000 -1.746300  1.217100
      Read-in Center 7322 is at   0.000000 -1.693400  1.217100
      Read-in Center 7323 is at   0.000000 -1.640400  1.217100
      Read-in Center 7324 is at   0.000000 -1.587500  1.217100
      Read-in Center 7325 is at   0.000000 -1.534600  1.217100
      Read-in Center 7326 is at   0.000000 -1.481700  1.217100
      Read-in Center 7327 is at   0.000000 -1.428800  1.217100
      Read-in Center 7328 is at   0.000000 -1.375900  1.217100
      Read-in Center 7329 is at   0.000000 -1.322900  1.217100
      Read-in Center 7330 is at   0.000000 -1.270000  1.217100
      Read-in Center 7331 is at   0.000000 -1.217100  1.217100
      Read-in Center 7332 is at   0.000000 -1.164200  1.217100
      Read-in Center 7333 is at   0.000000 -1.111300  1.217100
      Read-in Center 7334 is at   0.000000 -1.058400  1.217100
      Read-in Center 7335 is at   0.000000 -1.005400  1.217100
      Read-in Center 7336 is at   0.000000 -0.952500  1.217100
      Read-in Center 7337 is at   0.000000 -0.899600  1.217100
      Read-in Center 7338 is at   0.000000 -0.846700  1.217100
      Read-in Center 7339 is at   0.000000 -0.793800  1.217100
      Read-in Center 7340 is at   0.000000 -0.740800  1.217100
      Read-in Center 7341 is at   0.000000 -0.687900  1.217100
      Read-in Center 7342 is at   0.000000 -0.635000  1.217100
      Read-in Center 7343 is at   0.000000 -0.582100  1.217100
      Read-in Center 7344 is at   0.000000 -0.529200  1.217100
      Read-in Center 7345 is at   0.000000 -0.476300  1.217100
      Read-in Center 7346 is at   0.000000 -0.423300  1.217100
      Read-in Center 7347 is at   0.000000 -0.370400  1.217100
      Read-in Center 7348 is at   0.000000 -0.317500  1.217100
      Read-in Center 7349 is at   0.000000 -0.264600  1.217100
      Read-in Center 7350 is at   0.000000 -0.211700  1.217100
      Read-in Center 7351 is at   0.000000 -0.158800  1.217100
      Read-in Center 7352 is at   0.000000 -0.105800  1.217100
      Read-in Center 7353 is at   0.000000 -0.052900  1.217100
      Read-in Center 7354 is at   0.000000  0.000000  1.217100
      Read-in Center 7355 is at   0.000000  0.052900  1.217100
      Read-in Center 7356 is at   0.000000  0.105800  1.217100
      Read-in Center 7357 is at   0.000000  0.158800  1.217100
      Read-in Center 7358 is at   0.000000  0.211700  1.217100
      Read-in Center 7359 is at   0.000000  0.264600  1.217100
      Read-in Center 7360 is at   0.000000  0.317500  1.217100
      Read-in Center 7361 is at   0.000000  0.370400  1.217100
      Read-in Center 7362 is at   0.000000  0.423300  1.217100
      Read-in Center 7363 is at   0.000000  0.476300  1.217100
      Read-in Center 7364 is at   0.000000  0.529200  1.217100
      Read-in Center 7365 is at   0.000000  0.582100  1.217100
      Read-in Center 7366 is at   0.000000  0.635000  1.217100
      Read-in Center 7367 is at   0.000000  0.687900  1.217100
      Read-in Center 7368 is at   0.000000  0.740800  1.217100
      Read-in Center 7369 is at   0.000000  0.793800  1.217100
      Read-in Center 7370 is at   0.000000  0.846700  1.217100
      Read-in Center 7371 is at   0.000000  0.899600  1.217100
      Read-in Center 7372 is at   0.000000  0.952500  1.217100
      Read-in Center 7373 is at   0.000000  1.005400  1.217100
      Read-in Center 7374 is at   0.000000  1.058400  1.217100
      Read-in Center 7375 is at   0.000000  1.111300  1.217100
      Read-in Center 7376 is at   0.000000  1.164200  1.217100
      Read-in Center 7377 is at   0.000000  1.217100  1.217100
      Read-in Center 7378 is at   0.000000  1.270000  1.217100
      Read-in Center 7379 is at   0.000000  1.322900  1.217100
      Read-in Center 7380 is at   0.000000  1.375900  1.217100
      Read-in Center 7381 is at   0.000000  1.428800  1.217100
      Read-in Center 7382 is at   0.000000  1.481700  1.217100
      Read-in Center 7383 is at   0.000000  1.534600  1.217100
      Read-in Center 7384 is at   0.000000  1.587500  1.217100
      Read-in Center 7385 is at   0.000000  1.640400  1.217100
      Read-in Center 7386 is at   0.000000  1.693400  1.217100
      Read-in Center 7387 is at   0.000000  1.746300  1.217100
      Read-in Center 7388 is at   0.000000  1.799200  1.217100
      Read-in Center 7389 is at   0.000000  1.852100  1.217100
      Read-in Center 7390 is at   0.000000  1.905000  1.217100
      Read-in Center 7391 is at   0.000000  1.958000  1.217100
      Read-in Center 7392 is at   0.000000  2.010900  1.217100
      Read-in Center 7393 is at   0.000000  2.063800  1.217100
      Read-in Center 7394 is at   0.000000  2.116700  1.217100
      Read-in Center 7395 is at   0.000000  2.169600  1.217100
      Read-in Center 7396 is at   0.000000  2.222500  1.217100
      Read-in Center 7397 is at   0.000000  2.275500  1.217100
      Read-in Center 7398 is at   0.000000  2.328400  1.217100
      Read-in Center 7399 is at   0.000000  2.381300  1.217100
      Read-in Center 7400 is at   0.000000  2.434200  1.217100
      Read-in Center 7401 is at   0.000000  2.487100  1.217100
      Read-in Center 7402 is at   0.000000  2.540000  1.217100
      Read-in Center 7403 is at   0.000000  2.593000  1.217100
      Read-in Center 7404 is at   0.000000 -2.645900  1.270000
      Read-in Center 7405 is at   0.000000 -2.593000  1.270000
      Read-in Center 7406 is at   0.000000 -2.540000  1.270000
      Read-in Center 7407 is at   0.000000 -2.487100  1.270000
      Read-in Center 7408 is at   0.000000 -2.434200  1.270000
      Read-in Center 7409 is at   0.000000 -2.381300  1.270000
      Read-in Center 7410 is at   0.000000 -2.328400  1.270000
      Read-in Center 7411 is at   0.000000 -2.275500  1.270000
      Read-in Center 7412 is at   0.000000 -2.222500  1.270000
      Read-in Center 7413 is at   0.000000 -2.169600  1.270000
      Read-in Center 7414 is at   0.000000 -2.116700  1.270000
      Read-in Center 7415 is at   0.000000 -2.063800  1.270000
      Read-in Center 7416 is at   0.000000 -2.010900  1.270000
      Read-in Center 7417 is at   0.000000 -1.958000  1.270000
      Read-in Center 7418 is at   0.000000 -1.905000  1.270000
      Read-in Center 7419 is at   0.000000 -1.852100  1.270000
      Read-in Center 7420 is at   0.000000 -1.799200  1.270000
      Read-in Center 7421 is at   0.000000 -1.746300  1.270000
      Read-in Center 7422 is at   0.000000 -1.693400  1.270000
      Read-in Center 7423 is at   0.000000 -1.640400  1.270000
      Read-in Center 7424 is at   0.000000 -1.587500  1.270000
      Read-in Center 7425 is at   0.000000 -1.534600  1.270000
      Read-in Center 7426 is at   0.000000 -1.481700  1.270000
      Read-in Center 7427 is at   0.000000 -1.428800  1.270000
      Read-in Center 7428 is at   0.000000 -1.375900  1.270000
      Read-in Center 7429 is at   0.000000 -1.322900  1.270000
      Read-in Center 7430 is at   0.000000 -1.270000  1.270000
      Read-in Center 7431 is at   0.000000 -1.217100  1.270000
      Read-in Center 7432 is at   0.000000 -1.164200  1.270000
      Read-in Center 7433 is at   0.000000 -1.111300  1.270000
      Read-in Center 7434 is at   0.000000 -1.058400  1.270000
      Read-in Center 7435 is at   0.000000 -1.005400  1.270000
      Read-in Center 7436 is at   0.000000 -0.952500  1.270000
      Read-in Center 7437 is at   0.000000 -0.899600  1.270000
      Read-in Center 7438 is at   0.000000 -0.846700  1.270000
      Read-in Center 7439 is at   0.000000 -0.793800  1.270000
      Read-in Center 7440 is at   0.000000 -0.740800  1.270000
      Read-in Center 7441 is at   0.000000 -0.687900  1.270000
      Read-in Center 7442 is at   0.000000 -0.635000  1.270000
      Read-in Center 7443 is at   0.000000 -0.582100  1.270000
      Read-in Center 7444 is at   0.000000 -0.529200  1.270000
      Read-in Center 7445 is at   0.000000 -0.476300  1.270000
      Read-in Center 7446 is at   0.000000 -0.423300  1.270000
      Read-in Center 7447 is at   0.000000 -0.370400  1.270000
      Read-in Center 7448 is at   0.000000 -0.317500  1.270000
      Read-in Center 7449 is at   0.000000 -0.264600  1.270000
      Read-in Center 7450 is at   0.000000 -0.211700  1.270000
      Read-in Center 7451 is at   0.000000 -0.158800  1.270000
      Read-in Center 7452 is at   0.000000 -0.105800  1.270000
      Read-in Center 7453 is at   0.000000 -0.052900  1.270000
      Read-in Center 7454 is at   0.000000  0.000000  1.270000
      Read-in Center 7455 is at   0.000000  0.052900  1.270000
      Read-in Center 7456 is at   0.000000  0.105800  1.270000
      Read-in Center 7457 is at   0.000000  0.158800  1.270000
      Read-in Center 7458 is at   0.000000  0.211700  1.270000
      Read-in Center 7459 is at   0.000000  0.264600  1.270000
      Read-in Center 7460 is at   0.000000  0.317500  1.270000
      Read-in Center 7461 is at   0.000000  0.370400  1.270000
      Read-in Center 7462 is at   0.000000  0.423300  1.270000
      Read-in Center 7463 is at   0.000000  0.476300  1.270000
      Read-in Center 7464 is at   0.000000  0.529200  1.270000
      Read-in Center 7465 is at   0.000000  0.582100  1.270000
      Read-in Center 7466 is at   0.000000  0.635000  1.270000
      Read-in Center 7467 is at   0.000000  0.687900  1.270000
      Read-in Center 7468 is at   0.000000  0.740800  1.270000
      Read-in Center 7469 is at   0.000000  0.793800  1.270000
      Read-in Center 7470 is at   0.000000  0.846700  1.270000
      Read-in Center 7471 is at   0.000000  0.899600  1.270000
      Read-in Center 7472 is at   0.000000  0.952500  1.270000
      Read-in Center 7473 is at   0.000000  1.005400  1.270000
      Read-in Center 7474 is at   0.000000  1.058400  1.270000
      Read-in Center 7475 is at   0.000000  1.111300  1.270000
      Read-in Center 7476 is at   0.000000  1.164200  1.270000
      Read-in Center 7477 is at   0.000000  1.217100  1.270000
      Read-in Center 7478 is at   0.000000  1.270000  1.270000
      Read-in Center 7479 is at   0.000000  1.322900  1.270000
      Read-in Center 7480 is at   0.000000  1.375900  1.270000
      Read-in Center 7481 is at   0.000000  1.428800  1.270000
      Read-in Center 7482 is at   0.000000  1.481700  1.270000
      Read-in Center 7483 is at   0.000000  1.534600  1.270000
      Read-in Center 7484 is at   0.000000  1.587500  1.270000
      Read-in Center 7485 is at   0.000000  1.640400  1.270000
      Read-in Center 7486 is at   0.000000  1.693400  1.270000
      Read-in Center 7487 is at   0.000000  1.746300  1.270000
      Read-in Center 7488 is at   0.000000  1.799200  1.270000
      Read-in Center 7489 is at   0.000000  1.852100  1.270000
      Read-in Center 7490 is at   0.000000  1.905000  1.270000
      Read-in Center 7491 is at   0.000000  1.958000  1.270000
      Read-in Center 7492 is at   0.000000  2.010900  1.270000
      Read-in Center 7493 is at   0.000000  2.063800  1.270000
      Read-in Center 7494 is at   0.000000  2.116700  1.270000
      Read-in Center 7495 is at   0.000000  2.169600  1.270000
      Read-in Center 7496 is at   0.000000  2.222500  1.270000
      Read-in Center 7497 is at   0.000000  2.275500  1.270000
      Read-in Center 7498 is at   0.000000  2.328400  1.270000
      Read-in Center 7499 is at   0.000000  2.381300  1.270000
      Read-in Center 7500 is at   0.000000  2.434200  1.270000
      Read-in Center 7501 is at   0.000000  2.487100  1.270000
      Read-in Center 7502 is at   0.000000  2.540000  1.270000
      Read-in Center 7503 is at   0.000000  2.593000  1.270000
      Read-in Center 7504 is at   0.000000 -2.645900  1.322900
      Read-in Center 7505 is at   0.000000 -2.593000  1.322900
      Read-in Center 7506 is at   0.000000 -2.540000  1.322900
      Read-in Center 7507 is at   0.000000 -2.487100  1.322900
      Read-in Center 7508 is at   0.000000 -2.434200  1.322900
      Read-in Center 7509 is at   0.000000 -2.381300  1.322900
      Read-in Center 7510 is at   0.000000 -2.328400  1.322900
      Read-in Center 7511 is at   0.000000 -2.275500  1.322900
      Read-in Center 7512 is at   0.000000 -2.222500  1.322900
      Read-in Center 7513 is at   0.000000 -2.169600  1.322900
      Read-in Center 7514 is at   0.000000 -2.116700  1.322900
      Read-in Center 7515 is at   0.000000 -2.063800  1.322900
      Read-in Center 7516 is at   0.000000 -2.010900  1.322900
      Read-in Center 7517 is at   0.000000 -1.958000  1.322900
      Read-in Center 7518 is at   0.000000 -1.905000  1.322900
      Read-in Center 7519 is at   0.000000 -1.852100  1.322900
      Read-in Center 7520 is at   0.000000 -1.799200  1.322900
      Read-in Center 7521 is at   0.000000 -1.746300  1.322900
      Read-in Center 7522 is at   0.000000 -1.693400  1.322900
      Read-in Center 7523 is at   0.000000 -1.640400  1.322900
      Read-in Center 7524 is at   0.000000 -1.587500  1.322900
      Read-in Center 7525 is at   0.000000 -1.534600  1.322900
      Read-in Center 7526 is at   0.000000 -1.481700  1.322900
      Read-in Center 7527 is at   0.000000 -1.428800  1.322900
      Read-in Center 7528 is at   0.000000 -1.375900  1.322900
      Read-in Center 7529 is at   0.000000 -1.322900  1.322900
      Read-in Center 7530 is at   0.000000 -1.270000  1.322900
      Read-in Center 7531 is at   0.000000 -1.217100  1.322900
      Read-in Center 7532 is at   0.000000 -1.164200  1.322900
      Read-in Center 7533 is at   0.000000 -1.111300  1.322900
      Read-in Center 7534 is at   0.000000 -1.058400  1.322900
      Read-in Center 7535 is at   0.000000 -1.005400  1.322900
      Read-in Center 7536 is at   0.000000 -0.952500  1.322900
      Read-in Center 7537 is at   0.000000 -0.899600  1.322900
      Read-in Center 7538 is at   0.000000 -0.846700  1.322900
      Read-in Center 7539 is at   0.000000 -0.793800  1.322900
      Read-in Center 7540 is at   0.000000 -0.740800  1.322900
      Read-in Center 7541 is at   0.000000 -0.687900  1.322900
      Read-in Center 7542 is at   0.000000 -0.635000  1.322900
      Read-in Center 7543 is at   0.000000 -0.582100  1.322900
      Read-in Center 7544 is at   0.000000 -0.529200  1.322900
      Read-in Center 7545 is at   0.000000 -0.476300  1.322900
      Read-in Center 7546 is at   0.000000 -0.423300  1.322900
      Read-in Center 7547 is at   0.000000 -0.370400  1.322900
      Read-in Center 7548 is at   0.000000 -0.317500  1.322900
      Read-in Center 7549 is at   0.000000 -0.264600  1.322900
      Read-in Center 7550 is at   0.000000 -0.211700  1.322900
      Read-in Center 7551 is at   0.000000 -0.158800  1.322900
      Read-in Center 7552 is at   0.000000 -0.105800  1.322900
      Read-in Center 7553 is at   0.000000 -0.052900  1.322900
      Read-in Center 7554 is at   0.000000  0.000000  1.322900
      Read-in Center 7555 is at   0.000000  0.052900  1.322900
      Read-in Center 7556 is at   0.000000  0.105800  1.322900
      Read-in Center 7557 is at   0.000000  0.158800  1.322900
      Read-in Center 7558 is at   0.000000  0.211700  1.322900
      Read-in Center 7559 is at   0.000000  0.264600  1.322900
      Read-in Center 7560 is at   0.000000  0.317500  1.322900
      Read-in Center 7561 is at   0.000000  0.370400  1.322900
      Read-in Center 7562 is at   0.000000  0.423300  1.322900
      Read-in Center 7563 is at   0.000000  0.476300  1.322900
      Read-in Center 7564 is at   0.000000  0.529200  1.322900
      Read-in Center 7565 is at   0.000000  0.582100  1.322900
      Read-in Center 7566 is at   0.000000  0.635000  1.322900
      Read-in Center 7567 is at   0.000000  0.687900  1.322900
      Read-in Center 7568 is at   0.000000  0.740800  1.322900
      Read-in Center 7569 is at   0.000000  0.793800  1.322900
      Read-in Center 7570 is at   0.000000  0.846700  1.322900
      Read-in Center 7571 is at   0.000000  0.899600  1.322900
      Read-in Center 7572 is at   0.000000  0.952500  1.322900
      Read-in Center 7573 is at   0.000000  1.005400  1.322900
      Read-in Center 7574 is at   0.000000  1.058400  1.322900
      Read-in Center 7575 is at   0.000000  1.111300  1.322900
      Read-in Center 7576 is at   0.000000  1.164200  1.322900
      Read-in Center 7577 is at   0.000000  1.217100  1.322900
      Read-in Center 7578 is at   0.000000  1.270000  1.322900
      Read-in Center 7579 is at   0.000000  1.322900  1.322900
      Read-in Center 7580 is at   0.000000  1.375900  1.322900
      Read-in Center 7581 is at   0.000000  1.428800  1.322900
      Read-in Center 7582 is at   0.000000  1.481700  1.322900
      Read-in Center 7583 is at   0.000000  1.534600  1.322900
      Read-in Center 7584 is at   0.000000  1.587500  1.322900
      Read-in Center 7585 is at   0.000000  1.640400  1.322900
      Read-in Center 7586 is at   0.000000  1.693400  1.322900
      Read-in Center 7587 is at   0.000000  1.746300  1.322900
      Read-in Center 7588 is at   0.000000  1.799200  1.322900
      Read-in Center 7589 is at   0.000000  1.852100  1.322900
      Read-in Center 7590 is at   0.000000  1.905000  1.322900
      Read-in Center 7591 is at   0.000000  1.958000  1.322900
      Read-in Center 7592 is at   0.000000  2.010900  1.322900
      Read-in Center 7593 is at   0.000000  2.063800  1.322900
      Read-in Center 7594 is at   0.000000  2.116700  1.322900
      Read-in Center 7595 is at   0.000000  2.169600  1.322900
      Read-in Center 7596 is at   0.000000  2.222500  1.322900
      Read-in Center 7597 is at   0.000000  2.275500  1.322900
      Read-in Center 7598 is at   0.000000  2.328400  1.322900
      Read-in Center 7599 is at   0.000000  2.381300  1.322900
      Read-in Center 7600 is at   0.000000  2.434200  1.322900
      Read-in Center 7601 is at   0.000000  2.487100  1.322900
      Read-in Center 7602 is at   0.000000  2.540000  1.322900
      Read-in Center 7603 is at   0.000000  2.593000  1.322900
      Read-in Center 7604 is at   0.000000 -2.645900  1.375900
      Read-in Center 7605 is at   0.000000 -2.593000  1.375900
      Read-in Center 7606 is at   0.000000 -2.540000  1.375900
      Read-in Center 7607 is at   0.000000 -2.487100  1.375900
      Read-in Center 7608 is at   0.000000 -2.434200  1.375900
      Read-in Center 7609 is at   0.000000 -2.381300  1.375900
      Read-in Center 7610 is at   0.000000 -2.328400  1.375900
      Read-in Center 7611 is at   0.000000 -2.275500  1.375900
      Read-in Center 7612 is at   0.000000 -2.222500  1.375900
      Read-in Center 7613 is at   0.000000 -2.169600  1.375900
      Read-in Center 7614 is at   0.000000 -2.116700  1.375900
      Read-in Center 7615 is at   0.000000 -2.063800  1.375900
      Read-in Center 7616 is at   0.000000 -2.010900  1.375900
      Read-in Center 7617 is at   0.000000 -1.958000  1.375900
      Read-in Center 7618 is at   0.000000 -1.905000  1.375900
      Read-in Center 7619 is at   0.000000 -1.852100  1.375900
      Read-in Center 7620 is at   0.000000 -1.799200  1.375900
      Read-in Center 7621 is at   0.000000 -1.746300  1.375900
      Read-in Center 7622 is at   0.000000 -1.693400  1.375900
      Read-in Center 7623 is at   0.000000 -1.640400  1.375900
      Read-in Center 7624 is at   0.000000 -1.587500  1.375900
      Read-in Center 7625 is at   0.000000 -1.534600  1.375900
      Read-in Center 7626 is at   0.000000 -1.481700  1.375900
      Read-in Center 7627 is at   0.000000 -1.428800  1.375900
      Read-in Center 7628 is at   0.000000 -1.375900  1.375900
      Read-in Center 7629 is at   0.000000 -1.322900  1.375900
      Read-in Center 7630 is at   0.000000 -1.270000  1.375900
      Read-in Center 7631 is at   0.000000 -1.217100  1.375900
      Read-in Center 7632 is at   0.000000 -1.164200  1.375900
      Read-in Center 7633 is at   0.000000 -1.111300  1.375900
      Read-in Center 7634 is at   0.000000 -1.058400  1.375900
      Read-in Center 7635 is at   0.000000 -1.005400  1.375900
      Read-in Center 7636 is at   0.000000 -0.952500  1.375900
      Read-in Center 7637 is at   0.000000 -0.899600  1.375900
      Read-in Center 7638 is at   0.000000 -0.846700  1.375900
      Read-in Center 7639 is at   0.000000 -0.793800  1.375900
      Read-in Center 7640 is at   0.000000 -0.740800  1.375900
      Read-in Center 7641 is at   0.000000 -0.687900  1.375900
      Read-in Center 7642 is at   0.000000 -0.635000  1.375900
      Read-in Center 7643 is at   0.000000 -0.582100  1.375900
      Read-in Center 7644 is at   0.000000 -0.529200  1.375900
      Read-in Center 7645 is at   0.000000 -0.476300  1.375900
      Read-in Center 7646 is at   0.000000 -0.423300  1.375900
      Read-in Center 7647 is at   0.000000 -0.370400  1.375900
      Read-in Center 7648 is at   0.000000 -0.317500  1.375900
      Read-in Center 7649 is at   0.000000 -0.264600  1.375900
      Read-in Center 7650 is at   0.000000 -0.211700  1.375900
      Read-in Center 7651 is at   0.000000 -0.158800  1.375900
      Read-in Center 7652 is at   0.000000 -0.105800  1.375900
      Read-in Center 7653 is at   0.000000 -0.052900  1.375900
      Read-in Center 7654 is at   0.000000  0.000000  1.375900
      Read-in Center 7655 is at   0.000000  0.052900  1.375900
      Read-in Center 7656 is at   0.000000  0.105800  1.375900
      Read-in Center 7657 is at   0.000000  0.158800  1.375900
      Read-in Center 7658 is at   0.000000  0.211700  1.375900
      Read-in Center 7659 is at   0.000000  0.264600  1.375900
      Read-in Center 7660 is at   0.000000  0.317500  1.375900
      Read-in Center 7661 is at   0.000000  0.370400  1.375900
      Read-in Center 7662 is at   0.000000  0.423300  1.375900
      Read-in Center 7663 is at   0.000000  0.476300  1.375900
      Read-in Center 7664 is at   0.000000  0.529200  1.375900
      Read-in Center 7665 is at   0.000000  0.582100  1.375900
      Read-in Center 7666 is at   0.000000  0.635000  1.375900
      Read-in Center 7667 is at   0.000000  0.687900  1.375900
      Read-in Center 7668 is at   0.000000  0.740800  1.375900
      Read-in Center 7669 is at   0.000000  0.793800  1.375900
      Read-in Center 7670 is at   0.000000  0.846700  1.375900
      Read-in Center 7671 is at   0.000000  0.899600  1.375900
      Read-in Center 7672 is at   0.000000  0.952500  1.375900
      Read-in Center 7673 is at   0.000000  1.005400  1.375900
      Read-in Center 7674 is at   0.000000  1.058400  1.375900
      Read-in Center 7675 is at   0.000000  1.111300  1.375900
      Read-in Center 7676 is at   0.000000  1.164200  1.375900
      Read-in Center 7677 is at   0.000000  1.217100  1.375900
      Read-in Center 7678 is at   0.000000  1.270000  1.375900
      Read-in Center 7679 is at   0.000000  1.322900  1.375900
      Read-in Center 7680 is at   0.000000  1.375900  1.375900
      Read-in Center 7681 is at   0.000000  1.428800  1.375900
      Read-in Center 7682 is at   0.000000  1.481700  1.375900
      Read-in Center 7683 is at   0.000000  1.534600  1.375900
      Read-in Center 7684 is at   0.000000  1.587500  1.375900
      Read-in Center 7685 is at   0.000000  1.640400  1.375900
      Read-in Center 7686 is at   0.000000  1.693400  1.375900
      Read-in Center 7687 is at   0.000000  1.746300  1.375900
      Read-in Center 7688 is at   0.000000  1.799200  1.375900
      Read-in Center 7689 is at   0.000000  1.852100  1.375900
      Read-in Center 7690 is at   0.000000  1.905000  1.375900
      Read-in Center 7691 is at   0.000000  1.958000  1.375900
      Read-in Center 7692 is at   0.000000  2.010900  1.375900
      Read-in Center 7693 is at   0.000000  2.063800  1.375900
      Read-in Center 7694 is at   0.000000  2.116700  1.375900
      Read-in Center 7695 is at   0.000000  2.169600  1.375900
      Read-in Center 7696 is at   0.000000  2.222500  1.375900
      Read-in Center 7697 is at   0.000000  2.275500  1.375900
      Read-in Center 7698 is at   0.000000  2.328400  1.375900
      Read-in Center 7699 is at   0.000000  2.381300  1.375900
      Read-in Center 7700 is at   0.000000  2.434200  1.375900
      Read-in Center 7701 is at   0.000000  2.487100  1.375900
      Read-in Center 7702 is at   0.000000  2.540000  1.375900
      Read-in Center 7703 is at   0.000000  2.593000  1.375900
      Read-in Center 7704 is at   0.000000 -2.645900  1.428800
      Read-in Center 7705 is at   0.000000 -2.593000  1.428800
      Read-in Center 7706 is at   0.000000 -2.540000  1.428800
      Read-in Center 7707 is at   0.000000 -2.487100  1.428800
      Read-in Center 7708 is at   0.000000 -2.434200  1.428800
      Read-in Center 7709 is at   0.000000 -2.381300  1.428800
      Read-in Center 7710 is at   0.000000 -2.328400  1.428800
      Read-in Center 7711 is at   0.000000 -2.275500  1.428800
      Read-in Center 7712 is at   0.000000 -2.222500  1.428800
      Read-in Center 7713 is at   0.000000 -2.169600  1.428800
      Read-in Center 7714 is at   0.000000 -2.116700  1.428800
      Read-in Center 7715 is at   0.000000 -2.063800  1.428800
      Read-in Center 7716 is at   0.000000 -2.010900  1.428800
      Read-in Center 7717 is at   0.000000 -1.958000  1.428800
      Read-in Center 7718 is at   0.000000 -1.905000  1.428800
      Read-in Center 7719 is at   0.000000 -1.852100  1.428800
      Read-in Center 7720 is at   0.000000 -1.799200  1.428800
      Read-in Center 7721 is at   0.000000 -1.746300  1.428800
      Read-in Center 7722 is at   0.000000 -1.693400  1.428800
      Read-in Center 7723 is at   0.000000 -1.640400  1.428800
      Read-in Center 7724 is at   0.000000 -1.587500  1.428800
      Read-in Center 7725 is at   0.000000 -1.534600  1.428800
      Read-in Center 7726 is at   0.000000 -1.481700  1.428800
      Read-in Center 7727 is at   0.000000 -1.428800  1.428800
      Read-in Center 7728 is at   0.000000 -1.375900  1.428800
      Read-in Center 7729 is at   0.000000 -1.322900  1.428800
      Read-in Center 7730 is at   0.000000 -1.270000  1.428800
      Read-in Center 7731 is at   0.000000 -1.217100  1.428800
      Read-in Center 7732 is at   0.000000 -1.164200  1.428800
      Read-in Center 7733 is at   0.000000 -1.111300  1.428800
      Read-in Center 7734 is at   0.000000 -1.058400  1.428800
      Read-in Center 7735 is at   0.000000 -1.005400  1.428800
      Read-in Center 7736 is at   0.000000 -0.952500  1.428800
      Read-in Center 7737 is at   0.000000 -0.899600  1.428800
      Read-in Center 7738 is at   0.000000 -0.846700  1.428800
      Read-in Center 7739 is at   0.000000 -0.793800  1.428800
      Read-in Center 7740 is at   0.000000 -0.740800  1.428800
      Read-in Center 7741 is at   0.000000 -0.687900  1.428800
      Read-in Center 7742 is at   0.000000 -0.635000  1.428800
      Read-in Center 7743 is at   0.000000 -0.582100  1.428800
      Read-in Center 7744 is at   0.000000 -0.529200  1.428800
      Read-in Center 7745 is at   0.000000 -0.476300  1.428800
      Read-in Center 7746 is at   0.000000 -0.423300  1.428800
      Read-in Center 7747 is at   0.000000 -0.370400  1.428800
      Read-in Center 7748 is at   0.000000 -0.317500  1.428800
      Read-in Center 7749 is at   0.000000 -0.264600  1.428800
      Read-in Center 7750 is at   0.000000 -0.211700  1.428800
      Read-in Center 7751 is at   0.000000 -0.158800  1.428800
      Read-in Center 7752 is at   0.000000 -0.105800  1.428800
      Read-in Center 7753 is at   0.000000 -0.052900  1.428800
      Read-in Center 7754 is at   0.000000  0.000000  1.428800
      Read-in Center 7755 is at   0.000000  0.052900  1.428800
      Read-in Center 7756 is at   0.000000  0.105800  1.428800
      Read-in Center 7757 is at   0.000000  0.158800  1.428800
      Read-in Center 7758 is at   0.000000  0.211700  1.428800
      Read-in Center 7759 is at   0.000000  0.264600  1.428800
      Read-in Center 7760 is at   0.000000  0.317500  1.428800
      Read-in Center 7761 is at   0.000000  0.370400  1.428800
      Read-in Center 7762 is at   0.000000  0.423300  1.428800
      Read-in Center 7763 is at   0.000000  0.476300  1.428800
      Read-in Center 7764 is at   0.000000  0.529200  1.428800
      Read-in Center 7765 is at   0.000000  0.582100  1.428800
      Read-in Center 7766 is at   0.000000  0.635000  1.428800
      Read-in Center 7767 is at   0.000000  0.687900  1.428800
      Read-in Center 7768 is at   0.000000  0.740800  1.428800
      Read-in Center 7769 is at   0.000000  0.793800  1.428800
      Read-in Center 7770 is at   0.000000  0.846700  1.428800
      Read-in Center 7771 is at   0.000000  0.899600  1.428800
      Read-in Center 7772 is at   0.000000  0.952500  1.428800
      Read-in Center 7773 is at   0.000000  1.005400  1.428800
      Read-in Center 7774 is at   0.000000  1.058400  1.428800
      Read-in Center 7775 is at   0.000000  1.111300  1.428800
      Read-in Center 7776 is at   0.000000  1.164200  1.428800
      Read-in Center 7777 is at   0.000000  1.217100  1.428800
      Read-in Center 7778 is at   0.000000  1.270000  1.428800
      Read-in Center 7779 is at   0.000000  1.322900  1.428800
      Read-in Center 7780 is at   0.000000  1.375900  1.428800
      Read-in Center 7781 is at   0.000000  1.428800  1.428800
      Read-in Center 7782 is at   0.000000  1.481700  1.428800
      Read-in Center 7783 is at   0.000000  1.534600  1.428800
      Read-in Center 7784 is at   0.000000  1.587500  1.428800
      Read-in Center 7785 is at   0.000000  1.640400  1.428800
      Read-in Center 7786 is at   0.000000  1.693400  1.428800
      Read-in Center 7787 is at   0.000000  1.746300  1.428800
      Read-in Center 7788 is at   0.000000  1.799200  1.428800
      Read-in Center 7789 is at   0.000000  1.852100  1.428800
      Read-in Center 7790 is at   0.000000  1.905000  1.428800
      Read-in Center 7791 is at   0.000000  1.958000  1.428800
      Read-in Center 7792 is at   0.000000  2.010900  1.428800
      Read-in Center 7793 is at   0.000000  2.063800  1.428800
      Read-in Center 7794 is at   0.000000  2.116700  1.428800
      Read-in Center 7795 is at   0.000000  2.169600  1.428800
      Read-in Center 7796 is at   0.000000  2.222500  1.428800
      Read-in Center 7797 is at   0.000000  2.275500  1.428800
      Read-in Center 7798 is at   0.000000  2.328400  1.428800
      Read-in Center 7799 is at   0.000000  2.381300  1.428800
      Read-in Center 7800 is at   0.000000  2.434200  1.428800
      Read-in Center 7801 is at   0.000000  2.487100  1.428800
      Read-in Center 7802 is at   0.000000  2.540000  1.428800
      Read-in Center 7803 is at   0.000000  2.593000  1.428800
      Read-in Center 7804 is at   0.000000 -2.645900  1.481700
      Read-in Center 7805 is at   0.000000 -2.593000  1.481700
      Read-in Center 7806 is at   0.000000 -2.540000  1.481700
      Read-in Center 7807 is at   0.000000 -2.487100  1.481700
      Read-in Center 7808 is at   0.000000 -2.434200  1.481700
      Read-in Center 7809 is at   0.000000 -2.381300  1.481700
      Read-in Center 7810 is at   0.000000 -2.328400  1.481700
      Read-in Center 7811 is at   0.000000 -2.275500  1.481700
      Read-in Center 7812 is at   0.000000 -2.222500  1.481700
      Read-in Center 7813 is at   0.000000 -2.169600  1.481700
      Read-in Center 7814 is at   0.000000 -2.116700  1.481700
      Read-in Center 7815 is at   0.000000 -2.063800  1.481700
      Read-in Center 7816 is at   0.000000 -2.010900  1.481700
      Read-in Center 7817 is at   0.000000 -1.958000  1.481700
      Read-in Center 7818 is at   0.000000 -1.905000  1.481700
      Read-in Center 7819 is at   0.000000 -1.852100  1.481700
      Read-in Center 7820 is at   0.000000 -1.799200  1.481700
      Read-in Center 7821 is at   0.000000 -1.746300  1.481700
      Read-in Center 7822 is at   0.000000 -1.693400  1.481700
      Read-in Center 7823 is at   0.000000 -1.640400  1.481700
      Read-in Center 7824 is at   0.000000 -1.587500  1.481700
      Read-in Center 7825 is at   0.000000 -1.534600  1.481700
      Read-in Center 7826 is at   0.000000 -1.481700  1.481700
      Read-in Center 7827 is at   0.000000 -1.428800  1.481700
      Read-in Center 7828 is at   0.000000 -1.375900  1.481700
      Read-in Center 7829 is at   0.000000 -1.322900  1.481700
      Read-in Center 7830 is at   0.000000 -1.270000  1.481700
      Read-in Center 7831 is at   0.000000 -1.217100  1.481700
      Read-in Center 7832 is at   0.000000 -1.164200  1.481700
      Read-in Center 7833 is at   0.000000 -1.111300  1.481700
      Read-in Center 7834 is at   0.000000 -1.058400  1.481700
      Read-in Center 7835 is at   0.000000 -1.005400  1.481700
      Read-in Center 7836 is at   0.000000 -0.952500  1.481700
      Read-in Center 7837 is at   0.000000 -0.899600  1.481700
      Read-in Center 7838 is at   0.000000 -0.846700  1.481700
      Read-in Center 7839 is at   0.000000 -0.793800  1.481700
      Read-in Center 7840 is at   0.000000 -0.740800  1.481700
      Read-in Center 7841 is at   0.000000 -0.687900  1.481700
      Read-in Center 7842 is at   0.000000 -0.635000  1.481700
      Read-in Center 7843 is at   0.000000 -0.582100  1.481700
      Read-in Center 7844 is at   0.000000 -0.529200  1.481700
      Read-in Center 7845 is at   0.000000 -0.476300  1.481700
      Read-in Center 7846 is at   0.000000 -0.423300  1.481700
      Read-in Center 7847 is at   0.000000 -0.370400  1.481700
      Read-in Center 7848 is at   0.000000 -0.317500  1.481700
      Read-in Center 7849 is at   0.000000 -0.264600  1.481700
      Read-in Center 7850 is at   0.000000 -0.211700  1.481700
      Read-in Center 7851 is at   0.000000 -0.158800  1.481700
      Read-in Center 7852 is at   0.000000 -0.105800  1.481700
      Read-in Center 7853 is at   0.000000 -0.052900  1.481700
      Read-in Center 7854 is at   0.000000  0.000000  1.481700
      Read-in Center 7855 is at   0.000000  0.052900  1.481700
      Read-in Center 7856 is at   0.000000  0.105800  1.481700
      Read-in Center 7857 is at   0.000000  0.158800  1.481700
      Read-in Center 7858 is at   0.000000  0.211700  1.481700
      Read-in Center 7859 is at   0.000000  0.264600  1.481700
      Read-in Center 7860 is at   0.000000  0.317500  1.481700
      Read-in Center 7861 is at   0.000000  0.370400  1.481700
      Read-in Center 7862 is at   0.000000  0.423300  1.481700
      Read-in Center 7863 is at   0.000000  0.476300  1.481700
      Read-in Center 7864 is at   0.000000  0.529200  1.481700
      Read-in Center 7865 is at   0.000000  0.582100  1.481700
      Read-in Center 7866 is at   0.000000  0.635000  1.481700
      Read-in Center 7867 is at   0.000000  0.687900  1.481700
      Read-in Center 7868 is at   0.000000  0.740800  1.481700
      Read-in Center 7869 is at   0.000000  0.793800  1.481700
      Read-in Center 7870 is at   0.000000  0.846700  1.481700
      Read-in Center 7871 is at   0.000000  0.899600  1.481700
      Read-in Center 7872 is at   0.000000  0.952500  1.481700
      Read-in Center 7873 is at   0.000000  1.005400  1.481700
      Read-in Center 7874 is at   0.000000  1.058400  1.481700
      Read-in Center 7875 is at   0.000000  1.111300  1.481700
      Read-in Center 7876 is at   0.000000  1.164200  1.481700
      Read-in Center 7877 is at   0.000000  1.217100  1.481700
      Read-in Center 7878 is at   0.000000  1.270000  1.481700
      Read-in Center 7879 is at   0.000000  1.322900  1.481700
      Read-in Center 7880 is at   0.000000  1.375900  1.481700
      Read-in Center 7881 is at   0.000000  1.428800  1.481700
      Read-in Center 7882 is at   0.000000  1.481700  1.481700
      Read-in Center 7883 is at   0.000000  1.534600  1.481700
      Read-in Center 7884 is at   0.000000  1.587500  1.481700
      Read-in Center 7885 is at   0.000000  1.640400  1.481700
      Read-in Center 7886 is at   0.000000  1.693400  1.481700
      Read-in Center 7887 is at   0.000000  1.746300  1.481700
      Read-in Center 7888 is at   0.000000  1.799200  1.481700
      Read-in Center 7889 is at   0.000000  1.852100  1.481700
      Read-in Center 7890 is at   0.000000  1.905000  1.481700
      Read-in Center 7891 is at   0.000000  1.958000  1.481700
      Read-in Center 7892 is at   0.000000  2.010900  1.481700
      Read-in Center 7893 is at   0.000000  2.063800  1.481700
      Read-in Center 7894 is at   0.000000  2.116700  1.481700
      Read-in Center 7895 is at   0.000000  2.169600  1.481700
      Read-in Center 7896 is at   0.000000  2.222500  1.481700
      Read-in Center 7897 is at   0.000000  2.275500  1.481700
      Read-in Center 7898 is at   0.000000  2.328400  1.481700
      Read-in Center 7899 is at   0.000000  2.381300  1.481700
      Read-in Center 7900 is at   0.000000  2.434200  1.481700
      Read-in Center 7901 is at   0.000000  2.487100  1.481700
      Read-in Center 7902 is at   0.000000  2.540000  1.481700
      Read-in Center 7903 is at   0.000000  2.593000  1.481700
      Read-in Center 7904 is at   0.000000 -2.645900  1.534600
      Read-in Center 7905 is at   0.000000 -2.593000  1.534600
      Read-in Center 7906 is at   0.000000 -2.540000  1.534600
      Read-in Center 7907 is at   0.000000 -2.487100  1.534600
      Read-in Center 7908 is at   0.000000 -2.434200  1.534600
      Read-in Center 7909 is at   0.000000 -2.381300  1.534600
      Read-in Center 7910 is at   0.000000 -2.328400  1.534600
      Read-in Center 7911 is at   0.000000 -2.275500  1.534600
      Read-in Center 7912 is at   0.000000 -2.222500  1.534600
      Read-in Center 7913 is at   0.000000 -2.169600  1.534600
      Read-in Center 7914 is at   0.000000 -2.116700  1.534600
      Read-in Center 7915 is at   0.000000 -2.063800  1.534600
      Read-in Center 7916 is at   0.000000 -2.010900  1.534600
      Read-in Center 7917 is at   0.000000 -1.958000  1.534600
      Read-in Center 7918 is at   0.000000 -1.905000  1.534600
      Read-in Center 7919 is at   0.000000 -1.852100  1.534600
      Read-in Center 7920 is at   0.000000 -1.799200  1.534600
      Read-in Center 7921 is at   0.000000 -1.746300  1.534600
      Read-in Center 7922 is at   0.000000 -1.693400  1.534600
      Read-in Center 7923 is at   0.000000 -1.640400  1.534600
      Read-in Center 7924 is at   0.000000 -1.587500  1.534600
      Read-in Center 7925 is at   0.000000 -1.534600  1.534600
      Read-in Center 7926 is at   0.000000 -1.481700  1.534600
      Read-in Center 7927 is at   0.000000 -1.428800  1.534600
      Read-in Center 7928 is at   0.000000 -1.375900  1.534600
      Read-in Center 7929 is at   0.000000 -1.322900  1.534600
      Read-in Center 7930 is at   0.000000 -1.270000  1.534600
      Read-in Center 7931 is at   0.000000 -1.217100  1.534600
      Read-in Center 7932 is at   0.000000 -1.164200  1.534600
      Read-in Center 7933 is at   0.000000 -1.111300  1.534600
      Read-in Center 7934 is at   0.000000 -1.058400  1.534600
      Read-in Center 7935 is at   0.000000 -1.005400  1.534600
      Read-in Center 7936 is at   0.000000 -0.952500  1.534600
      Read-in Center 7937 is at   0.000000 -0.899600  1.534600
      Read-in Center 7938 is at   0.000000 -0.846700  1.534600
      Read-in Center 7939 is at   0.000000 -0.793800  1.534600
      Read-in Center 7940 is at   0.000000 -0.740800  1.534600
      Read-in Center 7941 is at   0.000000 -0.687900  1.534600
      Read-in Center 7942 is at   0.000000 -0.635000  1.534600
      Read-in Center 7943 is at   0.000000 -0.582100  1.534600
      Read-in Center 7944 is at   0.000000 -0.529200  1.534600
      Read-in Center 7945 is at   0.000000 -0.476300  1.534600
      Read-in Center 7946 is at   0.000000 -0.423300  1.534600
      Read-in Center 7947 is at   0.000000 -0.370400  1.534600
      Read-in Center 7948 is at   0.000000 -0.317500  1.534600
      Read-in Center 7949 is at   0.000000 -0.264600  1.534600
      Read-in Center 7950 is at   0.000000 -0.211700  1.534600
      Read-in Center 7951 is at   0.000000 -0.158800  1.534600
      Read-in Center 7952 is at   0.000000 -0.105800  1.534600
      Read-in Center 7953 is at   0.000000 -0.052900  1.534600
      Read-in Center 7954 is at   0.000000  0.000000  1.534600
      Read-in Center 7955 is at   0.000000  0.052900  1.534600
      Read-in Center 7956 is at   0.000000  0.105800  1.534600
      Read-in Center 7957 is at   0.000000  0.158800  1.534600
      Read-in Center 7958 is at   0.000000  0.211700  1.534600
      Read-in Center 7959 is at   0.000000  0.264600  1.534600
      Read-in Center 7960 is at   0.000000  0.317500  1.534600
      Read-in Center 7961 is at   0.000000  0.370400  1.534600
      Read-in Center 7962 is at   0.000000  0.423300  1.534600
      Read-in Center 7963 is at   0.000000  0.476300  1.534600
      Read-in Center 7964 is at   0.000000  0.529200  1.534600
      Read-in Center 7965 is at   0.000000  0.582100  1.534600
      Read-in Center 7966 is at   0.000000  0.635000  1.534600
      Read-in Center 7967 is at   0.000000  0.687900  1.534600
      Read-in Center 7968 is at   0.000000  0.740800  1.534600
      Read-in Center 7969 is at   0.000000  0.793800  1.534600
      Read-in Center 7970 is at   0.000000  0.846700  1.534600
      Read-in Center 7971 is at   0.000000  0.899600  1.534600
      Read-in Center 7972 is at   0.000000  0.952500  1.534600
      Read-in Center 7973 is at   0.000000  1.005400  1.534600
      Read-in Center 7974 is at   0.000000  1.058400  1.534600
      Read-in Center 7975 is at   0.000000  1.111300  1.534600
      Read-in Center 7976 is at   0.000000  1.164200  1.534600
      Read-in Center 7977 is at   0.000000  1.217100  1.534600
      Read-in Center 7978 is at   0.000000  1.270000  1.534600
      Read-in Center 7979 is at   0.000000  1.322900  1.534600
      Read-in Center 7980 is at   0.000000  1.375900  1.534600
      Read-in Center 7981 is at   0.000000  1.428800  1.534600
      Read-in Center 7982 is at   0.000000  1.481700  1.534600
      Read-in Center 7983 is at   0.000000  1.534600  1.534600
      Read-in Center 7984 is at   0.000000  1.587500  1.534600
      Read-in Center 7985 is at   0.000000  1.640400  1.534600
      Read-in Center 7986 is at   0.000000  1.693400  1.534600
      Read-in Center 7987 is at   0.000000  1.746300  1.534600
      Read-in Center 7988 is at   0.000000  1.799200  1.534600
      Read-in Center 7989 is at   0.000000  1.852100  1.534600
      Read-in Center 7990 is at   0.000000  1.905000  1.534600
      Read-in Center 7991 is at   0.000000  1.958000  1.534600
      Read-in Center 7992 is at   0.000000  2.010900  1.534600
      Read-in Center 7993 is at   0.000000  2.063800  1.534600
      Read-in Center 7994 is at   0.000000  2.116700  1.534600
      Read-in Center 7995 is at   0.000000  2.169600  1.534600
      Read-in Center 7996 is at   0.000000  2.222500  1.534600
      Read-in Center 7997 is at   0.000000  2.275500  1.534600
      Read-in Center 7998 is at   0.000000  2.328400  1.534600
      Read-in Center 7999 is at   0.000000  2.381300  1.534600
      Read-in Center 8000 is at   0.000000  2.434200  1.534600
      Read-in Center 8001 is at   0.000000  2.487100  1.534600
      Read-in Center 8002 is at   0.000000  2.540000  1.534600
      Read-in Center 8003 is at   0.000000  2.593000  1.534600
      Read-in Center 8004 is at   0.000000 -2.645900  1.587500
      Read-in Center 8005 is at   0.000000 -2.593000  1.587500
      Read-in Center 8006 is at   0.000000 -2.540000  1.587500
      Read-in Center 8007 is at   0.000000 -2.487100  1.587500
      Read-in Center 8008 is at   0.000000 -2.434200  1.587500
      Read-in Center 8009 is at   0.000000 -2.381300  1.587500
      Read-in Center 8010 is at   0.000000 -2.328400  1.587500
      Read-in Center 8011 is at   0.000000 -2.275500  1.587500
      Read-in Center 8012 is at   0.000000 -2.222500  1.587500
      Read-in Center 8013 is at   0.000000 -2.169600  1.587500
      Read-in Center 8014 is at   0.000000 -2.116700  1.587500
      Read-in Center 8015 is at   0.000000 -2.063800  1.587500
      Read-in Center 8016 is at   0.000000 -2.010900  1.587500
      Read-in Center 8017 is at   0.000000 -1.958000  1.587500
      Read-in Center 8018 is at   0.000000 -1.905000  1.587500
      Read-in Center 8019 is at   0.000000 -1.852100  1.587500
      Read-in Center 8020 is at   0.000000 -1.799200  1.587500
      Read-in Center 8021 is at   0.000000 -1.746300  1.587500
      Read-in Center 8022 is at   0.000000 -1.693400  1.587500
      Read-in Center 8023 is at   0.000000 -1.640400  1.587500
      Read-in Center 8024 is at   0.000000 -1.587500  1.587500
      Read-in Center 8025 is at   0.000000 -1.534600  1.587500
      Read-in Center 8026 is at   0.000000 -1.481700  1.587500
      Read-in Center 8027 is at   0.000000 -1.428800  1.587500
      Read-in Center 8028 is at   0.000000 -1.375900  1.587500
      Read-in Center 8029 is at   0.000000 -1.322900  1.587500
      Read-in Center 8030 is at   0.000000 -1.270000  1.587500
      Read-in Center 8031 is at   0.000000 -1.217100  1.587500
      Read-in Center 8032 is at   0.000000 -1.164200  1.587500
      Read-in Center 8033 is at   0.000000 -1.111300  1.587500
      Read-in Center 8034 is at   0.000000 -1.058400  1.587500
      Read-in Center 8035 is at   0.000000 -1.005400  1.587500
      Read-in Center 8036 is at   0.000000 -0.952500  1.587500
      Read-in Center 8037 is at   0.000000 -0.899600  1.587500
      Read-in Center 8038 is at   0.000000 -0.846700  1.587500
      Read-in Center 8039 is at   0.000000 -0.793800  1.587500
      Read-in Center 8040 is at   0.000000 -0.740800  1.587500
      Read-in Center 8041 is at   0.000000 -0.687900  1.587500
      Read-in Center 8042 is at   0.000000 -0.635000  1.587500
      Read-in Center 8043 is at   0.000000 -0.582100  1.587500
      Read-in Center 8044 is at   0.000000 -0.529200  1.587500
      Read-in Center 8045 is at   0.000000 -0.476300  1.587500
      Read-in Center 8046 is at   0.000000 -0.423300  1.587500
      Read-in Center 8047 is at   0.000000 -0.370400  1.587500
      Read-in Center 8048 is at   0.000000 -0.317500  1.587500
      Read-in Center 8049 is at   0.000000 -0.264600  1.587500
      Read-in Center 8050 is at   0.000000 -0.211700  1.587500
      Read-in Center 8051 is at   0.000000 -0.158800  1.587500
      Read-in Center 8052 is at   0.000000 -0.105800  1.587500
      Read-in Center 8053 is at   0.000000 -0.052900  1.587500
      Read-in Center 8054 is at   0.000000  0.000000  1.587500
      Read-in Center 8055 is at   0.000000  0.052900  1.587500
      Read-in Center 8056 is at   0.000000  0.105800  1.587500
      Read-in Center 8057 is at   0.000000  0.158800  1.587500
      Read-in Center 8058 is at   0.000000  0.211700  1.587500
      Read-in Center 8059 is at   0.000000  0.264600  1.587500
      Read-in Center 8060 is at   0.000000  0.317500  1.587500
      Read-in Center 8061 is at   0.000000  0.370400  1.587500
      Read-in Center 8062 is at   0.000000  0.423300  1.587500
      Read-in Center 8063 is at   0.000000  0.476300  1.587500
      Read-in Center 8064 is at   0.000000  0.529200  1.587500
      Read-in Center 8065 is at   0.000000  0.582100  1.587500
      Read-in Center 8066 is at   0.000000  0.635000  1.587500
      Read-in Center 8067 is at   0.000000  0.687900  1.587500
      Read-in Center 8068 is at   0.000000  0.740800  1.587500
      Read-in Center 8069 is at   0.000000  0.793800  1.587500
      Read-in Center 8070 is at   0.000000  0.846700  1.587500
      Read-in Center 8071 is at   0.000000  0.899600  1.587500
      Read-in Center 8072 is at   0.000000  0.952500  1.587500
      Read-in Center 8073 is at   0.000000  1.005400  1.587500
      Read-in Center 8074 is at   0.000000  1.058400  1.587500
      Read-in Center 8075 is at   0.000000  1.111300  1.587500
      Read-in Center 8076 is at   0.000000  1.164200  1.587500
      Read-in Center 8077 is at   0.000000  1.217100  1.587500
      Read-in Center 8078 is at   0.000000  1.270000  1.587500
      Read-in Center 8079 is at   0.000000  1.322900  1.587500
      Read-in Center 8080 is at   0.000000  1.375900  1.587500
      Read-in Center 8081 is at   0.000000  1.428800  1.587500
      Read-in Center 8082 is at   0.000000  1.481700  1.587500
      Read-in Center 8083 is at   0.000000  1.534600  1.587500
      Read-in Center 8084 is at   0.000000  1.587500  1.587500
      Read-in Center 8085 is at   0.000000  1.640400  1.587500
      Read-in Center 8086 is at   0.000000  1.693400  1.587500
      Read-in Center 8087 is at   0.000000  1.746300  1.587500
      Read-in Center 8088 is at   0.000000  1.799200  1.587500
      Read-in Center 8089 is at   0.000000  1.852100  1.587500
      Read-in Center 8090 is at   0.000000  1.905000  1.587500
      Read-in Center 8091 is at   0.000000  1.958000  1.587500
      Read-in Center 8092 is at   0.000000  2.010900  1.587500
      Read-in Center 8093 is at   0.000000  2.063800  1.587500
      Read-in Center 8094 is at   0.000000  2.116700  1.587500
      Read-in Center 8095 is at   0.000000  2.169600  1.587500
      Read-in Center 8096 is at   0.000000  2.222500  1.587500
      Read-in Center 8097 is at   0.000000  2.275500  1.587500
      Read-in Center 8098 is at   0.000000  2.328400  1.587500
      Read-in Center 8099 is at   0.000000  2.381300  1.587500
      Read-in Center 8100 is at   0.000000  2.434200  1.587500
      Read-in Center 8101 is at   0.000000  2.487100  1.587500
      Read-in Center 8102 is at   0.000000  2.540000  1.587500
      Read-in Center 8103 is at   0.000000  2.593000  1.587500
      Read-in Center 8104 is at   0.000000 -2.645900  1.640400
      Read-in Center 8105 is at   0.000000 -2.593000  1.640400
      Read-in Center 8106 is at   0.000000 -2.540000  1.640400
      Read-in Center 8107 is at   0.000000 -2.487100  1.640400
      Read-in Center 8108 is at   0.000000 -2.434200  1.640400
      Read-in Center 8109 is at   0.000000 -2.381300  1.640400
      Read-in Center 8110 is at   0.000000 -2.328400  1.640400
      Read-in Center 8111 is at   0.000000 -2.275500  1.640400
      Read-in Center 8112 is at   0.000000 -2.222500  1.640400
      Read-in Center 8113 is at   0.000000 -2.169600  1.640400
      Read-in Center 8114 is at   0.000000 -2.116700  1.640400
      Read-in Center 8115 is at   0.000000 -2.063800  1.640400
      Read-in Center 8116 is at   0.000000 -2.010900  1.640400
      Read-in Center 8117 is at   0.000000 -1.958000  1.640400
      Read-in Center 8118 is at   0.000000 -1.905000  1.640400
      Read-in Center 8119 is at   0.000000 -1.852100  1.640400
      Read-in Center 8120 is at   0.000000 -1.799200  1.640400
      Read-in Center 8121 is at   0.000000 -1.746300  1.640400
      Read-in Center 8122 is at   0.000000 -1.693400  1.640400
      Read-in Center 8123 is at   0.000000 -1.640400  1.640400
      Read-in Center 8124 is at   0.000000 -1.587500  1.640400
      Read-in Center 8125 is at   0.000000 -1.534600  1.640400
      Read-in Center 8126 is at   0.000000 -1.481700  1.640400
      Read-in Center 8127 is at   0.000000 -1.428800  1.640400
      Read-in Center 8128 is at   0.000000 -1.375900  1.640400
      Read-in Center 8129 is at   0.000000 -1.322900  1.640400
      Read-in Center 8130 is at   0.000000 -1.270000  1.640400
      Read-in Center 8131 is at   0.000000 -1.217100  1.640400
      Read-in Center 8132 is at   0.000000 -1.164200  1.640400
      Read-in Center 8133 is at   0.000000 -1.111300  1.640400
      Read-in Center 8134 is at   0.000000 -1.058400  1.640400
      Read-in Center 8135 is at   0.000000 -1.005400  1.640400
      Read-in Center 8136 is at   0.000000 -0.952500  1.640400
      Read-in Center 8137 is at   0.000000 -0.899600  1.640400
      Read-in Center 8138 is at   0.000000 -0.846700  1.640400
      Read-in Center 8139 is at   0.000000 -0.793800  1.640400
      Read-in Center 8140 is at   0.000000 -0.740800  1.640400
      Read-in Center 8141 is at   0.000000 -0.687900  1.640400
      Read-in Center 8142 is at   0.000000 -0.635000  1.640400
      Read-in Center 8143 is at   0.000000 -0.582100  1.640400
      Read-in Center 8144 is at   0.000000 -0.529200  1.640400
      Read-in Center 8145 is at   0.000000 -0.476300  1.640400
      Read-in Center 8146 is at   0.000000 -0.423300  1.640400
      Read-in Center 8147 is at   0.000000 -0.370400  1.640400
      Read-in Center 8148 is at   0.000000 -0.317500  1.640400
      Read-in Center 8149 is at   0.000000 -0.264600  1.640400
      Read-in Center 8150 is at   0.000000 -0.211700  1.640400
      Read-in Center 8151 is at   0.000000 -0.158800  1.640400
      Read-in Center 8152 is at   0.000000 -0.105800  1.640400
      Read-in Center 8153 is at   0.000000 -0.052900  1.640400
      Read-in Center 8154 is at   0.000000  0.000000  1.640400
      Read-in Center 8155 is at   0.000000  0.052900  1.640400
      Read-in Center 8156 is at   0.000000  0.105800  1.640400
      Read-in Center 8157 is at   0.000000  0.158800  1.640400
      Read-in Center 8158 is at   0.000000  0.211700  1.640400
      Read-in Center 8159 is at   0.000000  0.264600  1.640400
      Read-in Center 8160 is at   0.000000  0.317500  1.640400
      Read-in Center 8161 is at   0.000000  0.370400  1.640400
      Read-in Center 8162 is at   0.000000  0.423300  1.640400
      Read-in Center 8163 is at   0.000000  0.476300  1.640400
      Read-in Center 8164 is at   0.000000  0.529200  1.640400
      Read-in Center 8165 is at   0.000000  0.582100  1.640400
      Read-in Center 8166 is at   0.000000  0.635000  1.640400
      Read-in Center 8167 is at   0.000000  0.687900  1.640400
      Read-in Center 8168 is at   0.000000  0.740800  1.640400
      Read-in Center 8169 is at   0.000000  0.793800  1.640400
      Read-in Center 8170 is at   0.000000  0.846700  1.640400
      Read-in Center 8171 is at   0.000000  0.899600  1.640400
      Read-in Center 8172 is at   0.000000  0.952500  1.640400
      Read-in Center 8173 is at   0.000000  1.005400  1.640400
      Read-in Center 8174 is at   0.000000  1.058400  1.640400
      Read-in Center 8175 is at   0.000000  1.111300  1.640400
      Read-in Center 8176 is at   0.000000  1.164200  1.640400
      Read-in Center 8177 is at   0.000000  1.217100  1.640400
      Read-in Center 8178 is at   0.000000  1.270000  1.640400
      Read-in Center 8179 is at   0.000000  1.322900  1.640400
      Read-in Center 8180 is at   0.000000  1.375900  1.640400
      Read-in Center 8181 is at   0.000000  1.428800  1.640400
      Read-in Center 8182 is at   0.000000  1.481700  1.640400
      Read-in Center 8183 is at   0.000000  1.534600  1.640400
      Read-in Center 8184 is at   0.000000  1.587500  1.640400
      Read-in Center 8185 is at   0.000000  1.640400  1.640400
      Read-in Center 8186 is at   0.000000  1.693400  1.640400
      Read-in Center 8187 is at   0.000000  1.746300  1.640400
      Read-in Center 8188 is at   0.000000  1.799200  1.640400
      Read-in Center 8189 is at   0.000000  1.852100  1.640400
      Read-in Center 8190 is at   0.000000  1.905000  1.640400
      Read-in Center 8191 is at   0.000000  1.958000  1.640400
      Read-in Center 8192 is at   0.000000  2.010900  1.640400
      Read-in Center 8193 is at   0.000000  2.063800  1.640400
      Read-in Center 8194 is at   0.000000  2.116700  1.640400
      Read-in Center 8195 is at   0.000000  2.169600  1.640400
      Read-in Center 8196 is at   0.000000  2.222500  1.640400
      Read-in Center 8197 is at   0.000000  2.275500  1.640400
      Read-in Center 8198 is at   0.000000  2.328400  1.640400
      Read-in Center 8199 is at   0.000000  2.381300  1.640400
      Read-in Center 8200 is at   0.000000  2.434200  1.640400
      Read-in Center 8201 is at   0.000000  2.487100  1.640400
      Read-in Center 8202 is at   0.000000  2.540000  1.640400
      Read-in Center 8203 is at   0.000000  2.593000  1.640400
      Read-in Center 8204 is at   0.000000 -2.645900  1.693400
      Read-in Center 8205 is at   0.000000 -2.593000  1.693400
      Read-in Center 8206 is at   0.000000 -2.540000  1.693400
      Read-in Center 8207 is at   0.000000 -2.487100  1.693400
      Read-in Center 8208 is at   0.000000 -2.434200  1.693400
      Read-in Center 8209 is at   0.000000 -2.381300  1.693400
      Read-in Center 8210 is at   0.000000 -2.328400  1.693400
      Read-in Center 8211 is at   0.000000 -2.275500  1.693400
      Read-in Center 8212 is at   0.000000 -2.222500  1.693400
      Read-in Center 8213 is at   0.000000 -2.169600  1.693400
      Read-in Center 8214 is at   0.000000 -2.116700  1.693400
      Read-in Center 8215 is at   0.000000 -2.063800  1.693400
      Read-in Center 8216 is at   0.000000 -2.010900  1.693400
      Read-in Center 8217 is at   0.000000 -1.958000  1.693400
      Read-in Center 8218 is at   0.000000 -1.905000  1.693400
      Read-in Center 8219 is at   0.000000 -1.852100  1.693400
      Read-in Center 8220 is at   0.000000 -1.799200  1.693400
      Read-in Center 8221 is at   0.000000 -1.746300  1.693400
      Read-in Center 8222 is at   0.000000 -1.693400  1.693400
      Read-in Center 8223 is at   0.000000 -1.640400  1.693400
      Read-in Center 8224 is at   0.000000 -1.587500  1.693400
      Read-in Center 8225 is at   0.000000 -1.534600  1.693400
      Read-in Center 8226 is at   0.000000 -1.481700  1.693400
      Read-in Center 8227 is at   0.000000 -1.428800  1.693400
      Read-in Center 8228 is at   0.000000 -1.375900  1.693400
      Read-in Center 8229 is at   0.000000 -1.322900  1.693400
      Read-in Center 8230 is at   0.000000 -1.270000  1.693400
      Read-in Center 8231 is at   0.000000 -1.217100  1.693400
      Read-in Center 8232 is at   0.000000 -1.164200  1.693400
      Read-in Center 8233 is at   0.000000 -1.111300  1.693400
      Read-in Center 8234 is at   0.000000 -1.058400  1.693400
      Read-in Center 8235 is at   0.000000 -1.005400  1.693400
      Read-in Center 8236 is at   0.000000 -0.952500  1.693400
      Read-in Center 8237 is at   0.000000 -0.899600  1.693400
      Read-in Center 8238 is at   0.000000 -0.846700  1.693400
      Read-in Center 8239 is at   0.000000 -0.793800  1.693400
      Read-in Center 8240 is at   0.000000 -0.740800  1.693400
      Read-in Center 8241 is at   0.000000 -0.687900  1.693400
      Read-in Center 8242 is at   0.000000 -0.635000  1.693400
      Read-in Center 8243 is at   0.000000 -0.582100  1.693400
      Read-in Center 8244 is at   0.000000 -0.529200  1.693400
      Read-in Center 8245 is at   0.000000 -0.476300  1.693400
      Read-in Center 8246 is at   0.000000 -0.423300  1.693400
      Read-in Center 8247 is at   0.000000 -0.370400  1.693400
      Read-in Center 8248 is at   0.000000 -0.317500  1.693400
      Read-in Center 8249 is at   0.000000 -0.264600  1.693400
      Read-in Center 8250 is at   0.000000 -0.211700  1.693400
      Read-in Center 8251 is at   0.000000 -0.158800  1.693400
      Read-in Center 8252 is at   0.000000 -0.105800  1.693400
      Read-in Center 8253 is at   0.000000 -0.052900  1.693400
      Read-in Center 8254 is at   0.000000  0.000000  1.693400
      Read-in Center 8255 is at   0.000000  0.052900  1.693400
      Read-in Center 8256 is at   0.000000  0.105800  1.693400
      Read-in Center 8257 is at   0.000000  0.158800  1.693400
      Read-in Center 8258 is at   0.000000  0.211700  1.693400
      Read-in Center 8259 is at   0.000000  0.264600  1.693400
      Read-in Center 8260 is at   0.000000  0.317500  1.693400
      Read-in Center 8261 is at   0.000000  0.370400  1.693400
      Read-in Center 8262 is at   0.000000  0.423300  1.693400
      Read-in Center 8263 is at   0.000000  0.476300  1.693400
      Read-in Center 8264 is at   0.000000  0.529200  1.693400
      Read-in Center 8265 is at   0.000000  0.582100  1.693400
      Read-in Center 8266 is at   0.000000  0.635000  1.693400
      Read-in Center 8267 is at   0.000000  0.687900  1.693400
      Read-in Center 8268 is at   0.000000  0.740800  1.693400
      Read-in Center 8269 is at   0.000000  0.793800  1.693400
      Read-in Center 8270 is at   0.000000  0.846700  1.693400
      Read-in Center 8271 is at   0.000000  0.899600  1.693400
      Read-in Center 8272 is at   0.000000  0.952500  1.693400
      Read-in Center 8273 is at   0.000000  1.005400  1.693400
      Read-in Center 8274 is at   0.000000  1.058400  1.693400
      Read-in Center 8275 is at   0.000000  1.111300  1.693400
      Read-in Center 8276 is at   0.000000  1.164200  1.693400
      Read-in Center 8277 is at   0.000000  1.217100  1.693400
      Read-in Center 8278 is at   0.000000  1.270000  1.693400
      Read-in Center 8279 is at   0.000000  1.322900  1.693400
      Read-in Center 8280 is at   0.000000  1.375900  1.693400
      Read-in Center 8281 is at   0.000000  1.428800  1.693400
      Read-in Center 8282 is at   0.000000  1.481700  1.693400
      Read-in Center 8283 is at   0.000000  1.534600  1.693400
      Read-in Center 8284 is at   0.000000  1.587500  1.693400
      Read-in Center 8285 is at   0.000000  1.640400  1.693400
      Read-in Center 8286 is at   0.000000  1.693400  1.693400
      Read-in Center 8287 is at   0.000000  1.746300  1.693400
      Read-in Center 8288 is at   0.000000  1.799200  1.693400
      Read-in Center 8289 is at   0.000000  1.852100  1.693400
      Read-in Center 8290 is at   0.000000  1.905000  1.693400
      Read-in Center 8291 is at   0.000000  1.958000  1.693400
      Read-in Center 8292 is at   0.000000  2.010900  1.693400
      Read-in Center 8293 is at   0.000000  2.063800  1.693400
      Read-in Center 8294 is at   0.000000  2.116700  1.693400
      Read-in Center 8295 is at   0.000000  2.169600  1.693400
      Read-in Center 8296 is at   0.000000  2.222500  1.693400
      Read-in Center 8297 is at   0.000000  2.275500  1.693400
      Read-in Center 8298 is at   0.000000  2.328400  1.693400
      Read-in Center 8299 is at   0.000000  2.381300  1.693400
      Read-in Center 8300 is at   0.000000  2.434200  1.693400
      Read-in Center 8301 is at   0.000000  2.487100  1.693400
      Read-in Center 8302 is at   0.000000  2.540000  1.693400
      Read-in Center 8303 is at   0.000000  2.593000  1.693400
      Read-in Center 8304 is at   0.000000 -2.645900  1.746300
      Read-in Center 8305 is at   0.000000 -2.593000  1.746300
      Read-in Center 8306 is at   0.000000 -2.540000  1.746300
      Read-in Center 8307 is at   0.000000 -2.487100  1.746300
      Read-in Center 8308 is at   0.000000 -2.434200  1.746300
      Read-in Center 8309 is at   0.000000 -2.381300  1.746300
      Read-in Center 8310 is at   0.000000 -2.328400  1.746300
      Read-in Center 8311 is at   0.000000 -2.275500  1.746300
      Read-in Center 8312 is at   0.000000 -2.222500  1.746300
      Read-in Center 8313 is at   0.000000 -2.169600  1.746300
      Read-in Center 8314 is at   0.000000 -2.116700  1.746300
      Read-in Center 8315 is at   0.000000 -2.063800  1.746300
      Read-in Center 8316 is at   0.000000 -2.010900  1.746300
      Read-in Center 8317 is at   0.000000 -1.958000  1.746300
      Read-in Center 8318 is at   0.000000 -1.905000  1.746300
      Read-in Center 8319 is at   0.000000 -1.852100  1.746300
      Read-in Center 8320 is at   0.000000 -1.799200  1.746300
      Read-in Center 8321 is at   0.000000 -1.746300  1.746300
      Read-in Center 8322 is at   0.000000 -1.693400  1.746300
      Read-in Center 8323 is at   0.000000 -1.640400  1.746300
      Read-in Center 8324 is at   0.000000 -1.587500  1.746300
      Read-in Center 8325 is at   0.000000 -1.534600  1.746300
      Read-in Center 8326 is at   0.000000 -1.481700  1.746300
      Read-in Center 8327 is at   0.000000 -1.428800  1.746300
      Read-in Center 8328 is at   0.000000 -1.375900  1.746300
      Read-in Center 8329 is at   0.000000 -1.322900  1.746300
      Read-in Center 8330 is at   0.000000 -1.270000  1.746300
      Read-in Center 8331 is at   0.000000 -1.217100  1.746300
      Read-in Center 8332 is at   0.000000 -1.164200  1.746300
      Read-in Center 8333 is at   0.000000 -1.111300  1.746300
      Read-in Center 8334 is at   0.000000 -1.058400  1.746300
      Read-in Center 8335 is at   0.000000 -1.005400  1.746300
      Read-in Center 8336 is at   0.000000 -0.952500  1.746300
      Read-in Center 8337 is at   0.000000 -0.899600  1.746300
      Read-in Center 8338 is at   0.000000 -0.846700  1.746300
      Read-in Center 8339 is at   0.000000 -0.793800  1.746300
      Read-in Center 8340 is at   0.000000 -0.740800  1.746300
      Read-in Center 8341 is at   0.000000 -0.687900  1.746300
      Read-in Center 8342 is at   0.000000 -0.635000  1.746300
      Read-in Center 8343 is at   0.000000 -0.582100  1.746300
      Read-in Center 8344 is at   0.000000 -0.529200  1.746300
      Read-in Center 8345 is at   0.000000 -0.476300  1.746300
      Read-in Center 8346 is at   0.000000 -0.423300  1.746300
      Read-in Center 8347 is at   0.000000 -0.370400  1.746300
      Read-in Center 8348 is at   0.000000 -0.317500  1.746300
      Read-in Center 8349 is at   0.000000 -0.264600  1.746300
      Read-in Center 8350 is at   0.000000 -0.211700  1.746300
      Read-in Center 8351 is at   0.000000 -0.158800  1.746300
      Read-in Center 8352 is at   0.000000 -0.105800  1.746300
      Read-in Center 8353 is at   0.000000 -0.052900  1.746300
      Read-in Center 8354 is at   0.000000  0.000000  1.746300
      Read-in Center 8355 is at   0.000000  0.052900  1.746300
      Read-in Center 8356 is at   0.000000  0.105800  1.746300
      Read-in Center 8357 is at   0.000000  0.158800  1.746300
      Read-in Center 8358 is at   0.000000  0.211700  1.746300
      Read-in Center 8359 is at   0.000000  0.264600  1.746300
      Read-in Center 8360 is at   0.000000  0.317500  1.746300
      Read-in Center 8361 is at   0.000000  0.370400  1.746300
      Read-in Center 8362 is at   0.000000  0.423300  1.746300
      Read-in Center 8363 is at   0.000000  0.476300  1.746300
      Read-in Center 8364 is at   0.000000  0.529200  1.746300
      Read-in Center 8365 is at   0.000000  0.582100  1.746300
      Read-in Center 8366 is at   0.000000  0.635000  1.746300
      Read-in Center 8367 is at   0.000000  0.687900  1.746300
      Read-in Center 8368 is at   0.000000  0.740800  1.746300
      Read-in Center 8369 is at   0.000000  0.793800  1.746300
      Read-in Center 8370 is at   0.000000  0.846700  1.746300
      Read-in Center 8371 is at   0.000000  0.899600  1.746300
      Read-in Center 8372 is at   0.000000  0.952500  1.746300
      Read-in Center 8373 is at   0.000000  1.005400  1.746300
      Read-in Center 8374 is at   0.000000  1.058400  1.746300
      Read-in Center 8375 is at   0.000000  1.111300  1.746300
      Read-in Center 8376 is at   0.000000  1.164200  1.746300
      Read-in Center 8377 is at   0.000000  1.217100  1.746300
      Read-in Center 8378 is at   0.000000  1.270000  1.746300
      Read-in Center 8379 is at   0.000000  1.322900  1.746300
      Read-in Center 8380 is at   0.000000  1.375900  1.746300
      Read-in Center 8381 is at   0.000000  1.428800  1.746300
      Read-in Center 8382 is at   0.000000  1.481700  1.746300
      Read-in Center 8383 is at   0.000000  1.534600  1.746300
      Read-in Center 8384 is at   0.000000  1.587500  1.746300
      Read-in Center 8385 is at   0.000000  1.640400  1.746300
      Read-in Center 8386 is at   0.000000  1.693400  1.746300
      Read-in Center 8387 is at   0.000000  1.746300  1.746300
      Read-in Center 8388 is at   0.000000  1.799200  1.746300
      Read-in Center 8389 is at   0.000000  1.852100  1.746300
      Read-in Center 8390 is at   0.000000  1.905000  1.746300
      Read-in Center 8391 is at   0.000000  1.958000  1.746300
      Read-in Center 8392 is at   0.000000  2.010900  1.746300
      Read-in Center 8393 is at   0.000000  2.063800  1.746300
      Read-in Center 8394 is at   0.000000  2.116700  1.746300
      Read-in Center 8395 is at   0.000000  2.169600  1.746300
      Read-in Center 8396 is at   0.000000  2.222500  1.746300
      Read-in Center 8397 is at   0.000000  2.275500  1.746300
      Read-in Center 8398 is at   0.000000  2.328400  1.746300
      Read-in Center 8399 is at   0.000000  2.381300  1.746300
      Read-in Center 8400 is at   0.000000  2.434200  1.746300
      Read-in Center 8401 is at   0.000000  2.487100  1.746300
      Read-in Center 8402 is at   0.000000  2.540000  1.746300
      Read-in Center 8403 is at   0.000000  2.593000  1.746300
      Read-in Center 8404 is at   0.000000 -2.645900  1.799200
      Read-in Center 8405 is at   0.000000 -2.593000  1.799200
      Read-in Center 8406 is at   0.000000 -2.540000  1.799200
      Read-in Center 8407 is at   0.000000 -2.487100  1.799200
      Read-in Center 8408 is at   0.000000 -2.434200  1.799200
      Read-in Center 8409 is at   0.000000 -2.381300  1.799200
      Read-in Center 8410 is at   0.000000 -2.328400  1.799200
      Read-in Center 8411 is at   0.000000 -2.275500  1.799200
      Read-in Center 8412 is at   0.000000 -2.222500  1.799200
      Read-in Center 8413 is at   0.000000 -2.169600  1.799200
      Read-in Center 8414 is at   0.000000 -2.116700  1.799200
      Read-in Center 8415 is at   0.000000 -2.063800  1.799200
      Read-in Center 8416 is at   0.000000 -2.010900  1.799200
      Read-in Center 8417 is at   0.000000 -1.958000  1.799200
      Read-in Center 8418 is at   0.000000 -1.905000  1.799200
      Read-in Center 8419 is at   0.000000 -1.852100  1.799200
      Read-in Center 8420 is at   0.000000 -1.799200  1.799200
      Read-in Center 8421 is at   0.000000 -1.746300  1.799200
      Read-in Center 8422 is at   0.000000 -1.693400  1.799200
      Read-in Center 8423 is at   0.000000 -1.640400  1.799200
      Read-in Center 8424 is at   0.000000 -1.587500  1.799200
      Read-in Center 8425 is at   0.000000 -1.534600  1.799200
      Read-in Center 8426 is at   0.000000 -1.481700  1.799200
      Read-in Center 8427 is at   0.000000 -1.428800  1.799200
      Read-in Center 8428 is at   0.000000 -1.375900  1.799200
      Read-in Center 8429 is at   0.000000 -1.322900  1.799200
      Read-in Center 8430 is at   0.000000 -1.270000  1.799200
      Read-in Center 8431 is at   0.000000 -1.217100  1.799200
      Read-in Center 8432 is at   0.000000 -1.164200  1.799200
      Read-in Center 8433 is at   0.000000 -1.111300  1.799200
      Read-in Center 8434 is at   0.000000 -1.058400  1.799200
      Read-in Center 8435 is at   0.000000 -1.005400  1.799200
      Read-in Center 8436 is at   0.000000 -0.952500  1.799200
      Read-in Center 8437 is at   0.000000 -0.899600  1.799200
      Read-in Center 8438 is at   0.000000 -0.846700  1.799200
      Read-in Center 8439 is at   0.000000 -0.793800  1.799200
      Read-in Center 8440 is at   0.000000 -0.740800  1.799200
      Read-in Center 8441 is at   0.000000 -0.687900  1.799200
      Read-in Center 8442 is at   0.000000 -0.635000  1.799200
      Read-in Center 8443 is at   0.000000 -0.582100  1.799200
      Read-in Center 8444 is at   0.000000 -0.529200  1.799200
      Read-in Center 8445 is at   0.000000 -0.476300  1.799200
      Read-in Center 8446 is at   0.000000 -0.423300  1.799200
      Read-in Center 8447 is at   0.000000 -0.370400  1.799200
      Read-in Center 8448 is at   0.000000 -0.317500  1.799200
      Read-in Center 8449 is at   0.000000 -0.264600  1.799200
      Read-in Center 8450 is at   0.000000 -0.211700  1.799200
      Read-in Center 8451 is at   0.000000 -0.158800  1.799200
      Read-in Center 8452 is at   0.000000 -0.105800  1.799200
      Read-in Center 8453 is at   0.000000 -0.052900  1.799200
      Read-in Center 8454 is at   0.000000  0.000000  1.799200
      Read-in Center 8455 is at   0.000000  0.052900  1.799200
      Read-in Center 8456 is at   0.000000  0.105800  1.799200
      Read-in Center 8457 is at   0.000000  0.158800  1.799200
      Read-in Center 8458 is at   0.000000  0.211700  1.799200
      Read-in Center 8459 is at   0.000000  0.264600  1.799200
      Read-in Center 8460 is at   0.000000  0.317500  1.799200
      Read-in Center 8461 is at   0.000000  0.370400  1.799200
      Read-in Center 8462 is at   0.000000  0.423300  1.799200
      Read-in Center 8463 is at   0.000000  0.476300  1.799200
      Read-in Center 8464 is at   0.000000  0.529200  1.799200
      Read-in Center 8465 is at   0.000000  0.582100  1.799200
      Read-in Center 8466 is at   0.000000  0.635000  1.799200
      Read-in Center 8467 is at   0.000000  0.687900  1.799200
      Read-in Center 8468 is at   0.000000  0.740800  1.799200
      Read-in Center 8469 is at   0.000000  0.793800  1.799200
      Read-in Center 8470 is at   0.000000  0.846700  1.799200
      Read-in Center 8471 is at   0.000000  0.899600  1.799200
      Read-in Center 8472 is at   0.000000  0.952500  1.799200
      Read-in Center 8473 is at   0.000000  1.005400  1.799200
      Read-in Center 8474 is at   0.000000  1.058400  1.799200
      Read-in Center 8475 is at   0.000000  1.111300  1.799200
      Read-in Center 8476 is at   0.000000  1.164200  1.799200
      Read-in Center 8477 is at   0.000000  1.217100  1.799200
      Read-in Center 8478 is at   0.000000  1.270000  1.799200
      Read-in Center 8479 is at   0.000000  1.322900  1.799200
      Read-in Center 8480 is at   0.000000  1.375900  1.799200
      Read-in Center 8481 is at   0.000000  1.428800  1.799200
      Read-in Center 8482 is at   0.000000  1.481700  1.799200
      Read-in Center 8483 is at   0.000000  1.534600  1.799200
      Read-in Center 8484 is at   0.000000  1.587500  1.799200
      Read-in Center 8485 is at   0.000000  1.640400  1.799200
      Read-in Center 8486 is at   0.000000  1.693400  1.799200
      Read-in Center 8487 is at   0.000000  1.746300  1.799200
      Read-in Center 8488 is at   0.000000  1.799200  1.799200
      Read-in Center 8489 is at   0.000000  1.852100  1.799200
      Read-in Center 8490 is at   0.000000  1.905000  1.799200
      Read-in Center 8491 is at   0.000000  1.958000  1.799200
      Read-in Center 8492 is at   0.000000  2.010900  1.799200
      Read-in Center 8493 is at   0.000000  2.063800  1.799200
      Read-in Center 8494 is at   0.000000  2.116700  1.799200
      Read-in Center 8495 is at   0.000000  2.169600  1.799200
      Read-in Center 8496 is at   0.000000  2.222500  1.799200
      Read-in Center 8497 is at   0.000000  2.275500  1.799200
      Read-in Center 8498 is at   0.000000  2.328400  1.799200
      Read-in Center 8499 is at   0.000000  2.381300  1.799200
      Read-in Center 8500 is at   0.000000  2.434200  1.799200
      Read-in Center 8501 is at   0.000000  2.487100  1.799200
      Read-in Center 8502 is at   0.000000  2.540000  1.799200
      Read-in Center 8503 is at   0.000000  2.593000  1.799200
      Read-in Center 8504 is at   0.000000 -2.645900  1.852100
      Read-in Center 8505 is at   0.000000 -2.593000  1.852100
      Read-in Center 8506 is at   0.000000 -2.540000  1.852100
      Read-in Center 8507 is at   0.000000 -2.487100  1.852100
      Read-in Center 8508 is at   0.000000 -2.434200  1.852100
      Read-in Center 8509 is at   0.000000 -2.381300  1.852100
      Read-in Center 8510 is at   0.000000 -2.328400  1.852100
      Read-in Center 8511 is at   0.000000 -2.275500  1.852100
      Read-in Center 8512 is at   0.000000 -2.222500  1.852100
      Read-in Center 8513 is at   0.000000 -2.169600  1.852100
      Read-in Center 8514 is at   0.000000 -2.116700  1.852100
      Read-in Center 8515 is at   0.000000 -2.063800  1.852100
      Read-in Center 8516 is at   0.000000 -2.010900  1.852100
      Read-in Center 8517 is at   0.000000 -1.958000  1.852100
      Read-in Center 8518 is at   0.000000 -1.905000  1.852100
      Read-in Center 8519 is at   0.000000 -1.852100  1.852100
      Read-in Center 8520 is at   0.000000 -1.799200  1.852100
      Read-in Center 8521 is at   0.000000 -1.746300  1.852100
      Read-in Center 8522 is at   0.000000 -1.693400  1.852100
      Read-in Center 8523 is at   0.000000 -1.640400  1.852100
      Read-in Center 8524 is at   0.000000 -1.587500  1.852100
      Read-in Center 8525 is at   0.000000 -1.534600  1.852100
      Read-in Center 8526 is at   0.000000 -1.481700  1.852100
      Read-in Center 8527 is at   0.000000 -1.428800  1.852100
      Read-in Center 8528 is at   0.000000 -1.375900  1.852100
      Read-in Center 8529 is at   0.000000 -1.322900  1.852100
      Read-in Center 8530 is at   0.000000 -1.270000  1.852100
      Read-in Center 8531 is at   0.000000 -1.217100  1.852100
      Read-in Center 8532 is at   0.000000 -1.164200  1.852100
      Read-in Center 8533 is at   0.000000 -1.111300  1.852100
      Read-in Center 8534 is at   0.000000 -1.058400  1.852100
      Read-in Center 8535 is at   0.000000 -1.005400  1.852100
      Read-in Center 8536 is at   0.000000 -0.952500  1.852100
      Read-in Center 8537 is at   0.000000 -0.899600  1.852100
      Read-in Center 8538 is at   0.000000 -0.846700  1.852100
      Read-in Center 8539 is at   0.000000 -0.793800  1.852100
      Read-in Center 8540 is at   0.000000 -0.740800  1.852100
      Read-in Center 8541 is at   0.000000 -0.687900  1.852100
      Read-in Center 8542 is at   0.000000 -0.635000  1.852100
      Read-in Center 8543 is at   0.000000 -0.582100  1.852100
      Read-in Center 8544 is at   0.000000 -0.529200  1.852100
      Read-in Center 8545 is at   0.000000 -0.476300  1.852100
      Read-in Center 8546 is at   0.000000 -0.423300  1.852100
      Read-in Center 8547 is at   0.000000 -0.370400  1.852100
      Read-in Center 8548 is at   0.000000 -0.317500  1.852100
      Read-in Center 8549 is at   0.000000 -0.264600  1.852100
      Read-in Center 8550 is at   0.000000 -0.211700  1.852100
      Read-in Center 8551 is at   0.000000 -0.158800  1.852100
      Read-in Center 8552 is at   0.000000 -0.105800  1.852100
      Read-in Center 8553 is at   0.000000 -0.052900  1.852100
      Read-in Center 8554 is at   0.000000  0.000000  1.852100
      Read-in Center 8555 is at   0.000000  0.052900  1.852100
      Read-in Center 8556 is at   0.000000  0.105800  1.852100
      Read-in Center 8557 is at   0.000000  0.158800  1.852100
      Read-in Center 8558 is at   0.000000  0.211700  1.852100
      Read-in Center 8559 is at   0.000000  0.264600  1.852100
      Read-in Center 8560 is at   0.000000  0.317500  1.852100
      Read-in Center 8561 is at   0.000000  0.370400  1.852100
      Read-in Center 8562 is at   0.000000  0.423300  1.852100
      Read-in Center 8563 is at   0.000000  0.476300  1.852100
      Read-in Center 8564 is at   0.000000  0.529200  1.852100
      Read-in Center 8565 is at   0.000000  0.582100  1.852100
      Read-in Center 8566 is at   0.000000  0.635000  1.852100
      Read-in Center 8567 is at   0.000000  0.687900  1.852100
      Read-in Center 8568 is at   0.000000  0.740800  1.852100
      Read-in Center 8569 is at   0.000000  0.793800  1.852100
      Read-in Center 8570 is at   0.000000  0.846700  1.852100
      Read-in Center 8571 is at   0.000000  0.899600  1.852100
      Read-in Center 8572 is at   0.000000  0.952500  1.852100
      Read-in Center 8573 is at   0.000000  1.005400  1.852100
      Read-in Center 8574 is at   0.000000  1.058400  1.852100
      Read-in Center 8575 is at   0.000000  1.111300  1.852100
      Read-in Center 8576 is at   0.000000  1.164200  1.852100
      Read-in Center 8577 is at   0.000000  1.217100  1.852100
      Read-in Center 8578 is at   0.000000  1.270000  1.852100
      Read-in Center 8579 is at   0.000000  1.322900  1.852100
      Read-in Center 8580 is at   0.000000  1.375900  1.852100
      Read-in Center 8581 is at   0.000000  1.428800  1.852100
      Read-in Center 8582 is at   0.000000  1.481700  1.852100
      Read-in Center 8583 is at   0.000000  1.534600  1.852100
      Read-in Center 8584 is at   0.000000  1.587500  1.852100
      Read-in Center 8585 is at   0.000000  1.640400  1.852100
      Read-in Center 8586 is at   0.000000  1.693400  1.852100
      Read-in Center 8587 is at   0.000000  1.746300  1.852100
      Read-in Center 8588 is at   0.000000  1.799200  1.852100
      Read-in Center 8589 is at   0.000000  1.852100  1.852100
      Read-in Center 8590 is at   0.000000  1.905000  1.852100
      Read-in Center 8591 is at   0.000000  1.958000  1.852100
      Read-in Center 8592 is at   0.000000  2.010900  1.852100
      Read-in Center 8593 is at   0.000000  2.063800  1.852100
      Read-in Center 8594 is at   0.000000  2.116700  1.852100
      Read-in Center 8595 is at   0.000000  2.169600  1.852100
      Read-in Center 8596 is at   0.000000  2.222500  1.852100
      Read-in Center 8597 is at   0.000000  2.275500  1.852100
      Read-in Center 8598 is at   0.000000  2.328400  1.852100
      Read-in Center 8599 is at   0.000000  2.381300  1.852100
      Read-in Center 8600 is at   0.000000  2.434200  1.852100
      Read-in Center 8601 is at   0.000000  2.487100  1.852100
      Read-in Center 8602 is at   0.000000  2.540000  1.852100
      Read-in Center 8603 is at   0.000000  2.593000  1.852100
      Read-in Center 8604 is at   0.000000 -2.645900  1.905000
      Read-in Center 8605 is at   0.000000 -2.593000  1.905000
      Read-in Center 8606 is at   0.000000 -2.540000  1.905000
      Read-in Center 8607 is at   0.000000 -2.487100  1.905000
      Read-in Center 8608 is at   0.000000 -2.434200  1.905000
      Read-in Center 8609 is at   0.000000 -2.381300  1.905000
      Read-in Center 8610 is at   0.000000 -2.328400  1.905000
      Read-in Center 8611 is at   0.000000 -2.275500  1.905000
      Read-in Center 8612 is at   0.000000 -2.222500  1.905000
      Read-in Center 8613 is at   0.000000 -2.169600  1.905000
      Read-in Center 8614 is at   0.000000 -2.116700  1.905000
      Read-in Center 8615 is at   0.000000 -2.063800  1.905000
      Read-in Center 8616 is at   0.000000 -2.010900  1.905000
      Read-in Center 8617 is at   0.000000 -1.958000  1.905000
      Read-in Center 8618 is at   0.000000 -1.905000  1.905000
      Read-in Center 8619 is at   0.000000 -1.852100  1.905000
      Read-in Center 8620 is at   0.000000 -1.799200  1.905000
      Read-in Center 8621 is at   0.000000 -1.746300  1.905000
      Read-in Center 8622 is at   0.000000 -1.693400  1.905000
      Read-in Center 8623 is at   0.000000 -1.640400  1.905000
      Read-in Center 8624 is at   0.000000 -1.587500  1.905000
      Read-in Center 8625 is at   0.000000 -1.534600  1.905000
      Read-in Center 8626 is at   0.000000 -1.481700  1.905000
      Read-in Center 8627 is at   0.000000 -1.428800  1.905000
      Read-in Center 8628 is at   0.000000 -1.375900  1.905000
      Read-in Center 8629 is at   0.000000 -1.322900  1.905000
      Read-in Center 8630 is at   0.000000 -1.270000  1.905000
      Read-in Center 8631 is at   0.000000 -1.217100  1.905000
      Read-in Center 8632 is at   0.000000 -1.164200  1.905000
      Read-in Center 8633 is at   0.000000 -1.111300  1.905000
      Read-in Center 8634 is at   0.000000 -1.058400  1.905000
      Read-in Center 8635 is at   0.000000 -1.005400  1.905000
      Read-in Center 8636 is at   0.000000 -0.952500  1.905000
      Read-in Center 8637 is at   0.000000 -0.899600  1.905000
      Read-in Center 8638 is at   0.000000 -0.846700  1.905000
      Read-in Center 8639 is at   0.000000 -0.793800  1.905000
      Read-in Center 8640 is at   0.000000 -0.740800  1.905000
      Read-in Center 8641 is at   0.000000 -0.687900  1.905000
      Read-in Center 8642 is at   0.000000 -0.635000  1.905000
      Read-in Center 8643 is at   0.000000 -0.582100  1.905000
      Read-in Center 8644 is at   0.000000 -0.529200  1.905000
      Read-in Center 8645 is at   0.000000 -0.476300  1.905000
      Read-in Center 8646 is at   0.000000 -0.423300  1.905000
      Read-in Center 8647 is at   0.000000 -0.370400  1.905000
      Read-in Center 8648 is at   0.000000 -0.317500  1.905000
      Read-in Center 8649 is at   0.000000 -0.264600  1.905000
      Read-in Center 8650 is at   0.000000 -0.211700  1.905000
      Read-in Center 8651 is at   0.000000 -0.158800  1.905000
      Read-in Center 8652 is at   0.000000 -0.105800  1.905000
      Read-in Center 8653 is at   0.000000 -0.052900  1.905000
      Read-in Center 8654 is at   0.000000  0.000000  1.905000
      Read-in Center 8655 is at   0.000000  0.052900  1.905000
      Read-in Center 8656 is at   0.000000  0.105800  1.905000
      Read-in Center 8657 is at   0.000000  0.158800  1.905000
      Read-in Center 8658 is at   0.000000  0.211700  1.905000
      Read-in Center 8659 is at   0.000000  0.264600  1.905000
      Read-in Center 8660 is at   0.000000  0.317500  1.905000
      Read-in Center 8661 is at   0.000000  0.370400  1.905000
      Read-in Center 8662 is at   0.000000  0.423300  1.905000
      Read-in Center 8663 is at   0.000000  0.476300  1.905000
      Read-in Center 8664 is at   0.000000  0.529200  1.905000
      Read-in Center 8665 is at   0.000000  0.582100  1.905000
      Read-in Center 8666 is at   0.000000  0.635000  1.905000
      Read-in Center 8667 is at   0.000000  0.687900  1.905000
      Read-in Center 8668 is at   0.000000  0.740800  1.905000
      Read-in Center 8669 is at   0.000000  0.793800  1.905000
      Read-in Center 8670 is at   0.000000  0.846700  1.905000
      Read-in Center 8671 is at   0.000000  0.899600  1.905000
      Read-in Center 8672 is at   0.000000  0.952500  1.905000
      Read-in Center 8673 is at   0.000000  1.005400  1.905000
      Read-in Center 8674 is at   0.000000  1.058400  1.905000
      Read-in Center 8675 is at   0.000000  1.111300  1.905000
      Read-in Center 8676 is at   0.000000  1.164200  1.905000
      Read-in Center 8677 is at   0.000000  1.217100  1.905000
      Read-in Center 8678 is at   0.000000  1.270000  1.905000
      Read-in Center 8679 is at   0.000000  1.322900  1.905000
      Read-in Center 8680 is at   0.000000  1.375900  1.905000
      Read-in Center 8681 is at   0.000000  1.428800  1.905000
      Read-in Center 8682 is at   0.000000  1.481700  1.905000
      Read-in Center 8683 is at   0.000000  1.534600  1.905000
      Read-in Center 8684 is at   0.000000  1.587500  1.905000
      Read-in Center 8685 is at   0.000000  1.640400  1.905000
      Read-in Center 8686 is at   0.000000  1.693400  1.905000
      Read-in Center 8687 is at   0.000000  1.746300  1.905000
      Read-in Center 8688 is at   0.000000  1.799200  1.905000
      Read-in Center 8689 is at   0.000000  1.852100  1.905000
      Read-in Center 8690 is at   0.000000  1.905000  1.905000
      Read-in Center 8691 is at   0.000000  1.958000  1.905000
      Read-in Center 8692 is at   0.000000  2.010900  1.905000
      Read-in Center 8693 is at   0.000000  2.063800  1.905000
      Read-in Center 8694 is at   0.000000  2.116700  1.905000
      Read-in Center 8695 is at   0.000000  2.169600  1.905000
      Read-in Center 8696 is at   0.000000  2.222500  1.905000
      Read-in Center 8697 is at   0.000000  2.275500  1.905000
      Read-in Center 8698 is at   0.000000  2.328400  1.905000
      Read-in Center 8699 is at   0.000000  2.381300  1.905000
      Read-in Center 8700 is at   0.000000  2.434200  1.905000
      Read-in Center 8701 is at   0.000000  2.487100  1.905000
      Read-in Center 8702 is at   0.000000  2.540000  1.905000
      Read-in Center 8703 is at   0.000000  2.593000  1.905000
      Read-in Center 8704 is at   0.000000 -2.645900  1.958000
      Read-in Center 8705 is at   0.000000 -2.593000  1.958000
      Read-in Center 8706 is at   0.000000 -2.540000  1.958000
      Read-in Center 8707 is at   0.000000 -2.487100  1.958000
      Read-in Center 8708 is at   0.000000 -2.434200  1.958000
      Read-in Center 8709 is at   0.000000 -2.381300  1.958000
      Read-in Center 8710 is at   0.000000 -2.328400  1.958000
      Read-in Center 8711 is at   0.000000 -2.275500  1.958000
      Read-in Center 8712 is at   0.000000 -2.222500  1.958000
      Read-in Center 8713 is at   0.000000 -2.169600  1.958000
      Read-in Center 8714 is at   0.000000 -2.116700  1.958000
      Read-in Center 8715 is at   0.000000 -2.063800  1.958000
      Read-in Center 8716 is at   0.000000 -2.010900  1.958000
      Read-in Center 8717 is at   0.000000 -1.958000  1.958000
      Read-in Center 8718 is at   0.000000 -1.905000  1.958000
      Read-in Center 8719 is at   0.000000 -1.852100  1.958000
      Read-in Center 8720 is at   0.000000 -1.799200  1.958000
      Read-in Center 8721 is at   0.000000 -1.746300  1.958000
      Read-in Center 8722 is at   0.000000 -1.693400  1.958000
      Read-in Center 8723 is at   0.000000 -1.640400  1.958000
      Read-in Center 8724 is at   0.000000 -1.587500  1.958000
      Read-in Center 8725 is at   0.000000 -1.534600  1.958000
      Read-in Center 8726 is at   0.000000 -1.481700  1.958000
      Read-in Center 8727 is at   0.000000 -1.428800  1.958000
      Read-in Center 8728 is at   0.000000 -1.375900  1.958000
      Read-in Center 8729 is at   0.000000 -1.322900  1.958000
      Read-in Center 8730 is at   0.000000 -1.270000  1.958000
      Read-in Center 8731 is at   0.000000 -1.217100  1.958000
      Read-in Center 8732 is at   0.000000 -1.164200  1.958000
      Read-in Center 8733 is at   0.000000 -1.111300  1.958000
      Read-in Center 8734 is at   0.000000 -1.058400  1.958000
      Read-in Center 8735 is at   0.000000 -1.005400  1.958000
      Read-in Center 8736 is at   0.000000 -0.952500  1.958000
      Read-in Center 8737 is at   0.000000 -0.899600  1.958000
      Read-in Center 8738 is at   0.000000 -0.846700  1.958000
      Read-in Center 8739 is at   0.000000 -0.793800  1.958000
      Read-in Center 8740 is at   0.000000 -0.740800  1.958000
      Read-in Center 8741 is at   0.000000 -0.687900  1.958000
      Read-in Center 8742 is at   0.000000 -0.635000  1.958000
      Read-in Center 8743 is at   0.000000 -0.582100  1.958000
      Read-in Center 8744 is at   0.000000 -0.529200  1.958000
      Read-in Center 8745 is at   0.000000 -0.476300  1.958000
      Read-in Center 8746 is at   0.000000 -0.423300  1.958000
      Read-in Center 8747 is at   0.000000 -0.370400  1.958000
      Read-in Center 8748 is at   0.000000 -0.317500  1.958000
      Read-in Center 8749 is at   0.000000 -0.264600  1.958000
      Read-in Center 8750 is at   0.000000 -0.211700  1.958000
      Read-in Center 8751 is at   0.000000 -0.158800  1.958000
      Read-in Center 8752 is at   0.000000 -0.105800  1.958000
      Read-in Center 8753 is at   0.000000 -0.052900  1.958000
      Read-in Center 8754 is at   0.000000  0.000000  1.958000
      Read-in Center 8755 is at   0.000000  0.052900  1.958000
      Read-in Center 8756 is at   0.000000  0.105800  1.958000
      Read-in Center 8757 is at   0.000000  0.158800  1.958000
      Read-in Center 8758 is at   0.000000  0.211700  1.958000
      Read-in Center 8759 is at   0.000000  0.264600  1.958000
      Read-in Center 8760 is at   0.000000  0.317500  1.958000
      Read-in Center 8761 is at   0.000000  0.370400  1.958000
      Read-in Center 8762 is at   0.000000  0.423300  1.958000
      Read-in Center 8763 is at   0.000000  0.476300  1.958000
      Read-in Center 8764 is at   0.000000  0.529200  1.958000
      Read-in Center 8765 is at   0.000000  0.582100  1.958000
      Read-in Center 8766 is at   0.000000  0.635000  1.958000
      Read-in Center 8767 is at   0.000000  0.687900  1.958000
      Read-in Center 8768 is at   0.000000  0.740800  1.958000
      Read-in Center 8769 is at   0.000000  0.793800  1.958000
      Read-in Center 8770 is at   0.000000  0.846700  1.958000
      Read-in Center 8771 is at   0.000000  0.899600  1.958000
      Read-in Center 8772 is at   0.000000  0.952500  1.958000
      Read-in Center 8773 is at   0.000000  1.005400  1.958000
      Read-in Center 8774 is at   0.000000  1.058400  1.958000
      Read-in Center 8775 is at   0.000000  1.111300  1.958000
      Read-in Center 8776 is at   0.000000  1.164200  1.958000
      Read-in Center 8777 is at   0.000000  1.217100  1.958000
      Read-in Center 8778 is at   0.000000  1.270000  1.958000
      Read-in Center 8779 is at   0.000000  1.322900  1.958000
      Read-in Center 8780 is at   0.000000  1.375900  1.958000
      Read-in Center 8781 is at   0.000000  1.428800  1.958000
      Read-in Center 8782 is at   0.000000  1.481700  1.958000
      Read-in Center 8783 is at   0.000000  1.534600  1.958000
      Read-in Center 8784 is at   0.000000  1.587500  1.958000
      Read-in Center 8785 is at   0.000000  1.640400  1.958000
      Read-in Center 8786 is at   0.000000  1.693400  1.958000
      Read-in Center 8787 is at   0.000000  1.746300  1.958000
      Read-in Center 8788 is at   0.000000  1.799200  1.958000
      Read-in Center 8789 is at   0.000000  1.852100  1.958000
      Read-in Center 8790 is at   0.000000  1.905000  1.958000
      Read-in Center 8791 is at   0.000000  1.958000  1.958000
      Read-in Center 8792 is at   0.000000  2.010900  1.958000
      Read-in Center 8793 is at   0.000000  2.063800  1.958000
      Read-in Center 8794 is at   0.000000  2.116700  1.958000
      Read-in Center 8795 is at   0.000000  2.169600  1.958000
      Read-in Center 8796 is at   0.000000  2.222500  1.958000
      Read-in Center 8797 is at   0.000000  2.275500  1.958000
      Read-in Center 8798 is at   0.000000  2.328400  1.958000
      Read-in Center 8799 is at   0.000000  2.381300  1.958000
      Read-in Center 8800 is at   0.000000  2.434200  1.958000
      Read-in Center 8801 is at   0.000000  2.487100  1.958000
      Read-in Center 8802 is at   0.000000  2.540000  1.958000
      Read-in Center 8803 is at   0.000000  2.593000  1.958000
      Read-in Center 8804 is at   0.000000 -2.645900  2.010900
      Read-in Center 8805 is at   0.000000 -2.593000  2.010900
      Read-in Center 8806 is at   0.000000 -2.540000  2.010900
      Read-in Center 8807 is at   0.000000 -2.487100  2.010900
      Read-in Center 8808 is at   0.000000 -2.434200  2.010900
      Read-in Center 8809 is at   0.000000 -2.381300  2.010900
      Read-in Center 8810 is at   0.000000 -2.328400  2.010900
      Read-in Center 8811 is at   0.000000 -2.275500  2.010900
      Read-in Center 8812 is at   0.000000 -2.222500  2.010900
      Read-in Center 8813 is at   0.000000 -2.169600  2.010900
      Read-in Center 8814 is at   0.000000 -2.116700  2.010900
      Read-in Center 8815 is at   0.000000 -2.063800  2.010900
      Read-in Center 8816 is at   0.000000 -2.010900  2.010900
      Read-in Center 8817 is at   0.000000 -1.958000  2.010900
      Read-in Center 8818 is at   0.000000 -1.905000  2.010900
      Read-in Center 8819 is at   0.000000 -1.852100  2.010900
      Read-in Center 8820 is at   0.000000 -1.799200  2.010900
      Read-in Center 8821 is at   0.000000 -1.746300  2.010900
      Read-in Center 8822 is at   0.000000 -1.693400  2.010900
      Read-in Center 8823 is at   0.000000 -1.640400  2.010900
      Read-in Center 8824 is at   0.000000 -1.587500  2.010900
      Read-in Center 8825 is at   0.000000 -1.534600  2.010900
      Read-in Center 8826 is at   0.000000 -1.481700  2.010900
      Read-in Center 8827 is at   0.000000 -1.428800  2.010900
      Read-in Center 8828 is at   0.000000 -1.375900  2.010900
      Read-in Center 8829 is at   0.000000 -1.322900  2.010900
      Read-in Center 8830 is at   0.000000 -1.270000  2.010900
      Read-in Center 8831 is at   0.000000 -1.217100  2.010900
      Read-in Center 8832 is at   0.000000 -1.164200  2.010900
      Read-in Center 8833 is at   0.000000 -1.111300  2.010900
      Read-in Center 8834 is at   0.000000 -1.058400  2.010900
      Read-in Center 8835 is at   0.000000 -1.005400  2.010900
      Read-in Center 8836 is at   0.000000 -0.952500  2.010900
      Read-in Center 8837 is at   0.000000 -0.899600  2.010900
      Read-in Center 8838 is at   0.000000 -0.846700  2.010900
      Read-in Center 8839 is at   0.000000 -0.793800  2.010900
      Read-in Center 8840 is at   0.000000 -0.740800  2.010900
      Read-in Center 8841 is at   0.000000 -0.687900  2.010900
      Read-in Center 8842 is at   0.000000 -0.635000  2.010900
      Read-in Center 8843 is at   0.000000 -0.582100  2.010900
      Read-in Center 8844 is at   0.000000 -0.529200  2.010900
      Read-in Center 8845 is at   0.000000 -0.476300  2.010900
      Read-in Center 8846 is at   0.000000 -0.423300  2.010900
      Read-in Center 8847 is at   0.000000 -0.370400  2.010900
      Read-in Center 8848 is at   0.000000 -0.317500  2.010900
      Read-in Center 8849 is at   0.000000 -0.264600  2.010900
      Read-in Center 8850 is at   0.000000 -0.211700  2.010900
      Read-in Center 8851 is at   0.000000 -0.158800  2.010900
      Read-in Center 8852 is at   0.000000 -0.105800  2.010900
      Read-in Center 8853 is at   0.000000 -0.052900  2.010900
      Read-in Center 8854 is at   0.000000  0.000000  2.010900
      Read-in Center 8855 is at   0.000000  0.052900  2.010900
      Read-in Center 8856 is at   0.000000  0.105800  2.010900
      Read-in Center 8857 is at   0.000000  0.158800  2.010900
      Read-in Center 8858 is at   0.000000  0.211700  2.010900
      Read-in Center 8859 is at   0.000000  0.264600  2.010900
      Read-in Center 8860 is at   0.000000  0.317500  2.010900
      Read-in Center 8861 is at   0.000000  0.370400  2.010900
      Read-in Center 8862 is at   0.000000  0.423300  2.010900
      Read-in Center 8863 is at   0.000000  0.476300  2.010900
      Read-in Center 8864 is at   0.000000  0.529200  2.010900
      Read-in Center 8865 is at   0.000000  0.582100  2.010900
      Read-in Center 8866 is at   0.000000  0.635000  2.010900
      Read-in Center 8867 is at   0.000000  0.687900  2.010900
      Read-in Center 8868 is at   0.000000  0.740800  2.010900
      Read-in Center 8869 is at   0.000000  0.793800  2.010900
      Read-in Center 8870 is at   0.000000  0.846700  2.010900
      Read-in Center 8871 is at   0.000000  0.899600  2.010900
      Read-in Center 8872 is at   0.000000  0.952500  2.010900
      Read-in Center 8873 is at   0.000000  1.005400  2.010900
      Read-in Center 8874 is at   0.000000  1.058400  2.010900
      Read-in Center 8875 is at   0.000000  1.111300  2.010900
      Read-in Center 8876 is at   0.000000  1.164200  2.010900
      Read-in Center 8877 is at   0.000000  1.217100  2.010900
      Read-in Center 8878 is at   0.000000  1.270000  2.010900
      Read-in Center 8879 is at   0.000000  1.322900  2.010900
      Read-in Center 8880 is at   0.000000  1.375900  2.010900
      Read-in Center 8881 is at   0.000000  1.428800  2.010900
      Read-in Center 8882 is at   0.000000  1.481700  2.010900
      Read-in Center 8883 is at   0.000000  1.534600  2.010900
      Read-in Center 8884 is at   0.000000  1.587500  2.010900
      Read-in Center 8885 is at   0.000000  1.640400  2.010900
      Read-in Center 8886 is at   0.000000  1.693400  2.010900
      Read-in Center 8887 is at   0.000000  1.746300  2.010900
      Read-in Center 8888 is at   0.000000  1.799200  2.010900
      Read-in Center 8889 is at   0.000000  1.852100  2.010900
      Read-in Center 8890 is at   0.000000  1.905000  2.010900
      Read-in Center 8891 is at   0.000000  1.958000  2.010900
      Read-in Center 8892 is at   0.000000  2.010900  2.010900
      Read-in Center 8893 is at   0.000000  2.063800  2.010900
      Read-in Center 8894 is at   0.000000  2.116700  2.010900
      Read-in Center 8895 is at   0.000000  2.169600  2.010900
      Read-in Center 8896 is at   0.000000  2.222500  2.010900
      Read-in Center 8897 is at   0.000000  2.275500  2.010900
      Read-in Center 8898 is at   0.000000  2.328400  2.010900
      Read-in Center 8899 is at   0.000000  2.381300  2.010900
      Read-in Center 8900 is at   0.000000  2.434200  2.010900
      Read-in Center 8901 is at   0.000000  2.487100  2.010900
      Read-in Center 8902 is at   0.000000  2.540000  2.010900
      Read-in Center 8903 is at   0.000000  2.593000  2.010900
      Read-in Center 8904 is at   0.000000 -2.645900  2.063800
      Read-in Center 8905 is at   0.000000 -2.593000  2.063800
      Read-in Center 8906 is at   0.000000 -2.540000  2.063800
      Read-in Center 8907 is at   0.000000 -2.487100  2.063800
      Read-in Center 8908 is at   0.000000 -2.434200  2.063800
      Read-in Center 8909 is at   0.000000 -2.381300  2.063800
      Read-in Center 8910 is at   0.000000 -2.328400  2.063800
      Read-in Center 8911 is at   0.000000 -2.275500  2.063800
      Read-in Center 8912 is at   0.000000 -2.222500  2.063800
      Read-in Center 8913 is at   0.000000 -2.169600  2.063800
      Read-in Center 8914 is at   0.000000 -2.116700  2.063800
      Read-in Center 8915 is at   0.000000 -2.063800  2.063800
      Read-in Center 8916 is at   0.000000 -2.010900  2.063800
      Read-in Center 8917 is at   0.000000 -1.958000  2.063800
      Read-in Center 8918 is at   0.000000 -1.905000  2.063800
      Read-in Center 8919 is at   0.000000 -1.852100  2.063800
      Read-in Center 8920 is at   0.000000 -1.799200  2.063800
      Read-in Center 8921 is at   0.000000 -1.746300  2.063800
      Read-in Center 8922 is at   0.000000 -1.693400  2.063800
      Read-in Center 8923 is at   0.000000 -1.640400  2.063800
      Read-in Center 8924 is at   0.000000 -1.587500  2.063800
      Read-in Center 8925 is at   0.000000 -1.534600  2.063800
      Read-in Center 8926 is at   0.000000 -1.481700  2.063800
      Read-in Center 8927 is at   0.000000 -1.428800  2.063800
      Read-in Center 8928 is at   0.000000 -1.375900  2.063800
      Read-in Center 8929 is at   0.000000 -1.322900  2.063800
      Read-in Center 8930 is at   0.000000 -1.270000  2.063800
      Read-in Center 8931 is at   0.000000 -1.217100  2.063800
      Read-in Center 8932 is at   0.000000 -1.164200  2.063800
      Read-in Center 8933 is at   0.000000 -1.111300  2.063800
      Read-in Center 8934 is at   0.000000 -1.058400  2.063800
      Read-in Center 8935 is at   0.000000 -1.005400  2.063800
      Read-in Center 8936 is at   0.000000 -0.952500  2.063800
      Read-in Center 8937 is at   0.000000 -0.899600  2.063800
      Read-in Center 8938 is at   0.000000 -0.846700  2.063800
      Read-in Center 8939 is at   0.000000 -0.793800  2.063800
      Read-in Center 8940 is at   0.000000 -0.740800  2.063800
      Read-in Center 8941 is at   0.000000 -0.687900  2.063800
      Read-in Center 8942 is at   0.000000 -0.635000  2.063800
      Read-in Center 8943 is at   0.000000 -0.582100  2.063800
      Read-in Center 8944 is at   0.000000 -0.529200  2.063800
      Read-in Center 8945 is at   0.000000 -0.476300  2.063800
      Read-in Center 8946 is at   0.000000 -0.423300  2.063800
      Read-in Center 8947 is at   0.000000 -0.370400  2.063800
      Read-in Center 8948 is at   0.000000 -0.317500  2.063800
      Read-in Center 8949 is at   0.000000 -0.264600  2.063800
      Read-in Center 8950 is at   0.000000 -0.211700  2.063800
      Read-in Center 8951 is at   0.000000 -0.158800  2.063800
      Read-in Center 8952 is at   0.000000 -0.105800  2.063800
      Read-in Center 8953 is at   0.000000 -0.052900  2.063800
      Read-in Center 8954 is at   0.000000  0.000000  2.063800
      Read-in Center 8955 is at   0.000000  0.052900  2.063800
      Read-in Center 8956 is at   0.000000  0.105800  2.063800
      Read-in Center 8957 is at   0.000000  0.158800  2.063800
      Read-in Center 8958 is at   0.000000  0.211700  2.063800
      Read-in Center 8959 is at   0.000000  0.264600  2.063800
      Read-in Center 8960 is at   0.000000  0.317500  2.063800
      Read-in Center 8961 is at   0.000000  0.370400  2.063800
      Read-in Center 8962 is at   0.000000  0.423300  2.063800
      Read-in Center 8963 is at   0.000000  0.476300  2.063800
      Read-in Center 8964 is at   0.000000  0.529200  2.063800
      Read-in Center 8965 is at   0.000000  0.582100  2.063800
      Read-in Center 8966 is at   0.000000  0.635000  2.063800
      Read-in Center 8967 is at   0.000000  0.687900  2.063800
      Read-in Center 8968 is at   0.000000  0.740800  2.063800
      Read-in Center 8969 is at   0.000000  0.793800  2.063800
      Read-in Center 8970 is at   0.000000  0.846700  2.063800
      Read-in Center 8971 is at   0.000000  0.899600  2.063800
      Read-in Center 8972 is at   0.000000  0.952500  2.063800
      Read-in Center 8973 is at   0.000000  1.005400  2.063800
      Read-in Center 8974 is at   0.000000  1.058400  2.063800
      Read-in Center 8975 is at   0.000000  1.111300  2.063800
      Read-in Center 8976 is at   0.000000  1.164200  2.063800
      Read-in Center 8977 is at   0.000000  1.217100  2.063800
      Read-in Center 8978 is at   0.000000  1.270000  2.063800
      Read-in Center 8979 is at   0.000000  1.322900  2.063800
      Read-in Center 8980 is at   0.000000  1.375900  2.063800
      Read-in Center 8981 is at   0.000000  1.428800  2.063800
      Read-in Center 8982 is at   0.000000  1.481700  2.063800
      Read-in Center 8983 is at   0.000000  1.534600  2.063800
      Read-in Center 8984 is at   0.000000  1.587500  2.063800
      Read-in Center 8985 is at   0.000000  1.640400  2.063800
      Read-in Center 8986 is at   0.000000  1.693400  2.063800
      Read-in Center 8987 is at   0.000000  1.746300  2.063800
      Read-in Center 8988 is at   0.000000  1.799200  2.063800
      Read-in Center 8989 is at   0.000000  1.852100  2.063800
      Read-in Center 8990 is at   0.000000  1.905000  2.063800
      Read-in Center 8991 is at   0.000000  1.958000  2.063800
      Read-in Center 8992 is at   0.000000  2.010900  2.063800
      Read-in Center 8993 is at   0.000000  2.063800  2.063800
      Read-in Center 8994 is at   0.000000  2.116700  2.063800
      Read-in Center 8995 is at   0.000000  2.169600  2.063800
      Read-in Center 8996 is at   0.000000  2.222500  2.063800
      Read-in Center 8997 is at   0.000000  2.275500  2.063800
      Read-in Center 8998 is at   0.000000  2.328400  2.063800
      Read-in Center 8999 is at   0.000000  2.381300  2.063800
      Read-in Center 9000 is at   0.000000  2.434200  2.063800
      Read-in Center 9001 is at   0.000000  2.487100  2.063800
      Read-in Center 9002 is at   0.000000  2.540000  2.063800
      Read-in Center 9003 is at   0.000000  2.593000  2.063800
      Read-in Center 9004 is at   0.000000 -2.645900  2.116700
      Read-in Center 9005 is at   0.000000 -2.593000  2.116700
      Read-in Center 9006 is at   0.000000 -2.540000  2.116700
      Read-in Center 9007 is at   0.000000 -2.487100  2.116700
      Read-in Center 9008 is at   0.000000 -2.434200  2.116700
      Read-in Center 9009 is at   0.000000 -2.381300  2.116700
      Read-in Center 9010 is at   0.000000 -2.328400  2.116700
      Read-in Center 9011 is at   0.000000 -2.275500  2.116700
      Read-in Center 9012 is at   0.000000 -2.222500  2.116700
      Read-in Center 9013 is at   0.000000 -2.169600  2.116700
      Read-in Center 9014 is at   0.000000 -2.116700  2.116700
      Read-in Center 9015 is at   0.000000 -2.063800  2.116700
      Read-in Center 9016 is at   0.000000 -2.010900  2.116700
      Read-in Center 9017 is at   0.000000 -1.958000  2.116700
      Read-in Center 9018 is at   0.000000 -1.905000  2.116700
      Read-in Center 9019 is at   0.000000 -1.852100  2.116700
      Read-in Center 9020 is at   0.000000 -1.799200  2.116700
      Read-in Center 9021 is at   0.000000 -1.746300  2.116700
      Read-in Center 9022 is at   0.000000 -1.693400  2.116700
      Read-in Center 9023 is at   0.000000 -1.640400  2.116700
      Read-in Center 9024 is at   0.000000 -1.587500  2.116700
      Read-in Center 9025 is at   0.000000 -1.534600  2.116700
      Read-in Center 9026 is at   0.000000 -1.481700  2.116700
      Read-in Center 9027 is at   0.000000 -1.428800  2.116700
      Read-in Center 9028 is at   0.000000 -1.375900  2.116700
      Read-in Center 9029 is at   0.000000 -1.322900  2.116700
      Read-in Center 9030 is at   0.000000 -1.270000  2.116700
      Read-in Center 9031 is at   0.000000 -1.217100  2.116700
      Read-in Center 9032 is at   0.000000 -1.164200  2.116700
      Read-in Center 9033 is at   0.000000 -1.111300  2.116700
      Read-in Center 9034 is at   0.000000 -1.058400  2.116700
      Read-in Center 9035 is at   0.000000 -1.005400  2.116700
      Read-in Center 9036 is at   0.000000 -0.952500  2.116700
      Read-in Center 9037 is at   0.000000 -0.899600  2.116700
      Read-in Center 9038 is at   0.000000 -0.846700  2.116700
      Read-in Center 9039 is at   0.000000 -0.793800  2.116700
      Read-in Center 9040 is at   0.000000 -0.740800  2.116700
      Read-in Center 9041 is at   0.000000 -0.687900  2.116700
      Read-in Center 9042 is at   0.000000 -0.635000  2.116700
      Read-in Center 9043 is at   0.000000 -0.582100  2.116700
      Read-in Center 9044 is at   0.000000 -0.529200  2.116700
      Read-in Center 9045 is at   0.000000 -0.476300  2.116700
      Read-in Center 9046 is at   0.000000 -0.423300  2.116700
      Read-in Center 9047 is at   0.000000 -0.370400  2.116700
      Read-in Center 9048 is at   0.000000 -0.317500  2.116700
      Read-in Center 9049 is at   0.000000 -0.264600  2.116700
      Read-in Center 9050 is at   0.000000 -0.211700  2.116700
      Read-in Center 9051 is at   0.000000 -0.158800  2.116700
      Read-in Center 9052 is at   0.000000 -0.105800  2.116700
      Read-in Center 9053 is at   0.000000 -0.052900  2.116700
      Read-in Center 9054 is at   0.000000  0.000000  2.116700
      Read-in Center 9055 is at   0.000000  0.052900  2.116700
      Read-in Center 9056 is at   0.000000  0.105800  2.116700
      Read-in Center 9057 is at   0.000000  0.158800  2.116700
      Read-in Center 9058 is at   0.000000  0.211700  2.116700
      Read-in Center 9059 is at   0.000000  0.264600  2.116700
      Read-in Center 9060 is at   0.000000  0.317500  2.116700
      Read-in Center 9061 is at   0.000000  0.370400  2.116700
      Read-in Center 9062 is at   0.000000  0.423300  2.116700
      Read-in Center 9063 is at   0.000000  0.476300  2.116700
      Read-in Center 9064 is at   0.000000  0.529200  2.116700
      Read-in Center 9065 is at   0.000000  0.582100  2.116700
      Read-in Center 9066 is at   0.000000  0.635000  2.116700
      Read-in Center 9067 is at   0.000000  0.687900  2.116700
      Read-in Center 9068 is at   0.000000  0.740800  2.116700
      Read-in Center 9069 is at   0.000000  0.793800  2.116700
      Read-in Center 9070 is at   0.000000  0.846700  2.116700
      Read-in Center 9071 is at   0.000000  0.899600  2.116700
      Read-in Center 9072 is at   0.000000  0.952500  2.116700
      Read-in Center 9073 is at   0.000000  1.005400  2.116700
      Read-in Center 9074 is at   0.000000  1.058400  2.116700
      Read-in Center 9075 is at   0.000000  1.111300  2.116700
      Read-in Center 9076 is at   0.000000  1.164200  2.116700
      Read-in Center 9077 is at   0.000000  1.217100  2.116700
      Read-in Center 9078 is at   0.000000  1.270000  2.116700
      Read-in Center 9079 is at   0.000000  1.322900  2.116700
      Read-in Center 9080 is at   0.000000  1.375900  2.116700
      Read-in Center 9081 is at   0.000000  1.428800  2.116700
      Read-in Center 9082 is at   0.000000  1.481700  2.116700
      Read-in Center 9083 is at   0.000000  1.534600  2.116700
      Read-in Center 9084 is at   0.000000  1.587500  2.116700
      Read-in Center 9085 is at   0.000000  1.640400  2.116700
      Read-in Center 9086 is at   0.000000  1.693400  2.116700
      Read-in Center 9087 is at   0.000000  1.746300  2.116700
      Read-in Center 9088 is at   0.000000  1.799200  2.116700
      Read-in Center 9089 is at   0.000000  1.852100  2.116700
      Read-in Center 9090 is at   0.000000  1.905000  2.116700
      Read-in Center 9091 is at   0.000000  1.958000  2.116700
      Read-in Center 9092 is at   0.000000  2.010900  2.116700
      Read-in Center 9093 is at   0.000000  2.063800  2.116700
      Read-in Center 9094 is at   0.000000  2.116700  2.116700
      Read-in Center 9095 is at   0.000000  2.169600  2.116700
      Read-in Center 9096 is at   0.000000  2.222500  2.116700
      Read-in Center 9097 is at   0.000000  2.275500  2.116700
      Read-in Center 9098 is at   0.000000  2.328400  2.116700
      Read-in Center 9099 is at   0.000000  2.381300  2.116700
      Read-in Center 9100 is at   0.000000  2.434200  2.116700
      Read-in Center 9101 is at   0.000000  2.487100  2.116700
      Read-in Center 9102 is at   0.000000  2.540000  2.116700
      Read-in Center 9103 is at   0.000000  2.593000  2.116700
      Read-in Center 9104 is at   0.000000 -2.645900  2.169600
      Read-in Center 9105 is at   0.000000 -2.593000  2.169600
      Read-in Center 9106 is at   0.000000 -2.540000  2.169600
      Read-in Center 9107 is at   0.000000 -2.487100  2.169600
      Read-in Center 9108 is at   0.000000 -2.434200  2.169600
      Read-in Center 9109 is at   0.000000 -2.381300  2.169600
      Read-in Center 9110 is at   0.000000 -2.328400  2.169600
      Read-in Center 9111 is at   0.000000 -2.275500  2.169600
      Read-in Center 9112 is at   0.000000 -2.222500  2.169600
      Read-in Center 9113 is at   0.000000 -2.169600  2.169600
      Read-in Center 9114 is at   0.000000 -2.116700  2.169600
      Read-in Center 9115 is at   0.000000 -2.063800  2.169600
      Read-in Center 9116 is at   0.000000 -2.010900  2.169600
      Read-in Center 9117 is at   0.000000 -1.958000  2.169600
      Read-in Center 9118 is at   0.000000 -1.905000  2.169600
      Read-in Center 9119 is at   0.000000 -1.852100  2.169600
      Read-in Center 9120 is at   0.000000 -1.799200  2.169600
      Read-in Center 9121 is at   0.000000 -1.746300  2.169600
      Read-in Center 9122 is at   0.000000 -1.693400  2.169600
      Read-in Center 9123 is at   0.000000 -1.640400  2.169600
      Read-in Center 9124 is at   0.000000 -1.587500  2.169600
      Read-in Center 9125 is at   0.000000 -1.534600  2.169600
      Read-in Center 9126 is at   0.000000 -1.481700  2.169600
      Read-in Center 9127 is at   0.000000 -1.428800  2.169600
      Read-in Center 9128 is at   0.000000 -1.375900  2.169600
      Read-in Center 9129 is at   0.000000 -1.322900  2.169600
      Read-in Center 9130 is at   0.000000 -1.270000  2.169600
      Read-in Center 9131 is at   0.000000 -1.217100  2.169600
      Read-in Center 9132 is at   0.000000 -1.164200  2.169600
      Read-in Center 9133 is at   0.000000 -1.111300  2.169600
      Read-in Center 9134 is at   0.000000 -1.058400  2.169600
      Read-in Center 9135 is at   0.000000 -1.005400  2.169600
      Read-in Center 9136 is at   0.000000 -0.952500  2.169600
      Read-in Center 9137 is at   0.000000 -0.899600  2.169600
      Read-in Center 9138 is at   0.000000 -0.846700  2.169600
      Read-in Center 9139 is at   0.000000 -0.793800  2.169600
      Read-in Center 9140 is at   0.000000 -0.740800  2.169600
      Read-in Center 9141 is at   0.000000 -0.687900  2.169600
      Read-in Center 9142 is at   0.000000 -0.635000  2.169600
      Read-in Center 9143 is at   0.000000 -0.582100  2.169600
      Read-in Center 9144 is at   0.000000 -0.529200  2.169600
      Read-in Center 9145 is at   0.000000 -0.476300  2.169600
      Read-in Center 9146 is at   0.000000 -0.423300  2.169600
      Read-in Center 9147 is at   0.000000 -0.370400  2.169600
      Read-in Center 9148 is at   0.000000 -0.317500  2.169600
      Read-in Center 9149 is at   0.000000 -0.264600  2.169600
      Read-in Center 9150 is at   0.000000 -0.211700  2.169600
      Read-in Center 9151 is at   0.000000 -0.158800  2.169600
      Read-in Center 9152 is at   0.000000 -0.105800  2.169600
      Read-in Center 9153 is at   0.000000 -0.052900  2.169600
      Read-in Center 9154 is at   0.000000  0.000000  2.169600
      Read-in Center 9155 is at   0.000000  0.052900  2.169600
      Read-in Center 9156 is at   0.000000  0.105800  2.169600
      Read-in Center 9157 is at   0.000000  0.158800  2.169600
      Read-in Center 9158 is at   0.000000  0.211700  2.169600
      Read-in Center 9159 is at   0.000000  0.264600  2.169600
      Read-in Center 9160 is at   0.000000  0.317500  2.169600
      Read-in Center 9161 is at   0.000000  0.370400  2.169600
      Read-in Center 9162 is at   0.000000  0.423300  2.169600
      Read-in Center 9163 is at   0.000000  0.476300  2.169600
      Read-in Center 9164 is at   0.000000  0.529200  2.169600
      Read-in Center 9165 is at   0.000000  0.582100  2.169600
      Read-in Center 9166 is at   0.000000  0.635000  2.169600
      Read-in Center 9167 is at   0.000000  0.687900  2.169600
      Read-in Center 9168 is at   0.000000  0.740800  2.169600
      Read-in Center 9169 is at   0.000000  0.793800  2.169600
      Read-in Center 9170 is at   0.000000  0.846700  2.169600
      Read-in Center 9171 is at   0.000000  0.899600  2.169600
      Read-in Center 9172 is at   0.000000  0.952500  2.169600
      Read-in Center 9173 is at   0.000000  1.005400  2.169600
      Read-in Center 9174 is at   0.000000  1.058400  2.169600
      Read-in Center 9175 is at   0.000000  1.111300  2.169600
      Read-in Center 9176 is at   0.000000  1.164200  2.169600
      Read-in Center 9177 is at   0.000000  1.217100  2.169600
      Read-in Center 9178 is at   0.000000  1.270000  2.169600
      Read-in Center 9179 is at   0.000000  1.322900  2.169600
      Read-in Center 9180 is at   0.000000  1.375900  2.169600
      Read-in Center 9181 is at   0.000000  1.428800  2.169600
      Read-in Center 9182 is at   0.000000  1.481700  2.169600
      Read-in Center 9183 is at   0.000000  1.534600  2.169600
      Read-in Center 9184 is at   0.000000  1.587500  2.169600
      Read-in Center 9185 is at   0.000000  1.640400  2.169600
      Read-in Center 9186 is at   0.000000  1.693400  2.169600
      Read-in Center 9187 is at   0.000000  1.746300  2.169600
      Read-in Center 9188 is at   0.000000  1.799200  2.169600
      Read-in Center 9189 is at   0.000000  1.852100  2.169600
      Read-in Center 9190 is at   0.000000  1.905000  2.169600
      Read-in Center 9191 is at   0.000000  1.958000  2.169600
      Read-in Center 9192 is at   0.000000  2.010900  2.169600
      Read-in Center 9193 is at   0.000000  2.063800  2.169600
      Read-in Center 9194 is at   0.000000  2.116700  2.169600
      Read-in Center 9195 is at   0.000000  2.169600  2.169600
      Read-in Center 9196 is at   0.000000  2.222500  2.169600
      Read-in Center 9197 is at   0.000000  2.275500  2.169600
      Read-in Center 9198 is at   0.000000  2.328400  2.169600
      Read-in Center 9199 is at   0.000000  2.381300  2.169600
      Read-in Center 9200 is at   0.000000  2.434200  2.169600
      Read-in Center 9201 is at   0.000000  2.487100  2.169600
      Read-in Center 9202 is at   0.000000  2.540000  2.169600
      Read-in Center 9203 is at   0.000000  2.593000  2.169600
      Read-in Center 9204 is at   0.000000 -2.645900  2.222500
      Read-in Center 9205 is at   0.000000 -2.593000  2.222500
      Read-in Center 9206 is at   0.000000 -2.540000  2.222500
      Read-in Center 9207 is at   0.000000 -2.487100  2.222500
      Read-in Center 9208 is at   0.000000 -2.434200  2.222500
      Read-in Center 9209 is at   0.000000 -2.381300  2.222500
      Read-in Center 9210 is at   0.000000 -2.328400  2.222500
      Read-in Center 9211 is at   0.000000 -2.275500  2.222500
      Read-in Center 9212 is at   0.000000 -2.222500  2.222500
      Read-in Center 9213 is at   0.000000 -2.169600  2.222500
      Read-in Center 9214 is at   0.000000 -2.116700  2.222500
      Read-in Center 9215 is at   0.000000 -2.063800  2.222500
      Read-in Center 9216 is at   0.000000 -2.010900  2.222500
      Read-in Center 9217 is at   0.000000 -1.958000  2.222500
      Read-in Center 9218 is at   0.000000 -1.905000  2.222500
      Read-in Center 9219 is at   0.000000 -1.852100  2.222500
      Read-in Center 9220 is at   0.000000 -1.799200  2.222500
      Read-in Center 9221 is at   0.000000 -1.746300  2.222500
      Read-in Center 9222 is at   0.000000 -1.693400  2.222500
      Read-in Center 9223 is at   0.000000 -1.640400  2.222500
      Read-in Center 9224 is at   0.000000 -1.587500  2.222500
      Read-in Center 9225 is at   0.000000 -1.534600  2.222500
      Read-in Center 9226 is at   0.000000 -1.481700  2.222500
      Read-in Center 9227 is at   0.000000 -1.428800  2.222500
      Read-in Center 9228 is at   0.000000 -1.375900  2.222500
      Read-in Center 9229 is at   0.000000 -1.322900  2.222500
      Read-in Center 9230 is at   0.000000 -1.270000  2.222500
      Read-in Center 9231 is at   0.000000 -1.217100  2.222500
      Read-in Center 9232 is at   0.000000 -1.164200  2.222500
      Read-in Center 9233 is at   0.000000 -1.111300  2.222500
      Read-in Center 9234 is at   0.000000 -1.058400  2.222500
      Read-in Center 9235 is at   0.000000 -1.005400  2.222500
      Read-in Center 9236 is at   0.000000 -0.952500  2.222500
      Read-in Center 9237 is at   0.000000 -0.899600  2.222500
      Read-in Center 9238 is at   0.000000 -0.846700  2.222500
      Read-in Center 9239 is at   0.000000 -0.793800  2.222500
      Read-in Center 9240 is at   0.000000 -0.740800  2.222500
      Read-in Center 9241 is at   0.000000 -0.687900  2.222500
      Read-in Center 9242 is at   0.000000 -0.635000  2.222500
      Read-in Center 9243 is at   0.000000 -0.582100  2.222500
      Read-in Center 9244 is at   0.000000 -0.529200  2.222500
      Read-in Center 9245 is at   0.000000 -0.476300  2.222500
      Read-in Center 9246 is at   0.000000 -0.423300  2.222500
      Read-in Center 9247 is at   0.000000 -0.370400  2.222500
      Read-in Center 9248 is at   0.000000 -0.317500  2.222500
      Read-in Center 9249 is at   0.000000 -0.264600  2.222500
      Read-in Center 9250 is at   0.000000 -0.211700  2.222500
      Read-in Center 9251 is at   0.000000 -0.158800  2.222500
      Read-in Center 9252 is at   0.000000 -0.105800  2.222500
      Read-in Center 9253 is at   0.000000 -0.052900  2.222500
      Read-in Center 9254 is at   0.000000  0.000000  2.222500
      Read-in Center 9255 is at   0.000000  0.052900  2.222500
      Read-in Center 9256 is at   0.000000  0.105800  2.222500
      Read-in Center 9257 is at   0.000000  0.158800  2.222500
      Read-in Center 9258 is at   0.000000  0.211700  2.222500
      Read-in Center 9259 is at   0.000000  0.264600  2.222500
      Read-in Center 9260 is at   0.000000  0.317500  2.222500
      Read-in Center 9261 is at   0.000000  0.370400  2.222500
      Read-in Center 9262 is at   0.000000  0.423300  2.222500
      Read-in Center 9263 is at   0.000000  0.476300  2.222500
      Read-in Center 9264 is at   0.000000  0.529200  2.222500
      Read-in Center 9265 is at   0.000000  0.582100  2.222500
      Read-in Center 9266 is at   0.000000  0.635000  2.222500
      Read-in Center 9267 is at   0.000000  0.687900  2.222500
      Read-in Center 9268 is at   0.000000  0.740800  2.222500
      Read-in Center 9269 is at   0.000000  0.793800  2.222500
      Read-in Center 9270 is at   0.000000  0.846700  2.222500
      Read-in Center 9271 is at   0.000000  0.899600  2.222500
      Read-in Center 9272 is at   0.000000  0.952500  2.222500
      Read-in Center 9273 is at   0.000000  1.005400  2.222500
      Read-in Center 9274 is at   0.000000  1.058400  2.222500
      Read-in Center 9275 is at   0.000000  1.111300  2.222500
      Read-in Center 9276 is at   0.000000  1.164200  2.222500
      Read-in Center 9277 is at   0.000000  1.217100  2.222500
      Read-in Center 9278 is at   0.000000  1.270000  2.222500
      Read-in Center 9279 is at   0.000000  1.322900  2.222500
      Read-in Center 9280 is at   0.000000  1.375900  2.222500
      Read-in Center 9281 is at   0.000000  1.428800  2.222500
      Read-in Center 9282 is at   0.000000  1.481700  2.222500
      Read-in Center 9283 is at   0.000000  1.534600  2.222500
      Read-in Center 9284 is at   0.000000  1.587500  2.222500
      Read-in Center 9285 is at   0.000000  1.640400  2.222500
      Read-in Center 9286 is at   0.000000  1.693400  2.222500
      Read-in Center 9287 is at   0.000000  1.746300  2.222500
      Read-in Center 9288 is at   0.000000  1.799200  2.222500
      Read-in Center 9289 is at   0.000000  1.852100  2.222500
      Read-in Center 9290 is at   0.000000  1.905000  2.222500
      Read-in Center 9291 is at   0.000000  1.958000  2.222500
      Read-in Center 9292 is at   0.000000  2.010900  2.222500
      Read-in Center 9293 is at   0.000000  2.063800  2.222500
      Read-in Center 9294 is at   0.000000  2.116700  2.222500
      Read-in Center 9295 is at   0.000000  2.169600  2.222500
      Read-in Center 9296 is at   0.000000  2.222500  2.222500
      Read-in Center 9297 is at   0.000000  2.275500  2.222500
      Read-in Center 9298 is at   0.000000  2.328400  2.222500
      Read-in Center 9299 is at   0.000000  2.381300  2.222500
      Read-in Center 9300 is at   0.000000  2.434200  2.222500
      Read-in Center 9301 is at   0.000000  2.487100  2.222500
      Read-in Center 9302 is at   0.000000  2.540000  2.222500
      Read-in Center 9303 is at   0.000000  2.593000  2.222500
      Read-in Center 9304 is at   0.000000 -2.645900  2.275500
      Read-in Center 9305 is at   0.000000 -2.593000  2.275500
      Read-in Center 9306 is at   0.000000 -2.540000  2.275500
      Read-in Center 9307 is at   0.000000 -2.487100  2.275500
      Read-in Center 9308 is at   0.000000 -2.434200  2.275500
      Read-in Center 9309 is at   0.000000 -2.381300  2.275500
      Read-in Center 9310 is at   0.000000 -2.328400  2.275500
      Read-in Center 9311 is at   0.000000 -2.275500  2.275500
      Read-in Center 9312 is at   0.000000 -2.222500  2.275500
      Read-in Center 9313 is at   0.000000 -2.169600  2.275500
      Read-in Center 9314 is at   0.000000 -2.116700  2.275500
      Read-in Center 9315 is at   0.000000 -2.063800  2.275500
      Read-in Center 9316 is at   0.000000 -2.010900  2.275500
      Read-in Center 9317 is at   0.000000 -1.958000  2.275500
      Read-in Center 9318 is at   0.000000 -1.905000  2.275500
      Read-in Center 9319 is at   0.000000 -1.852100  2.275500
      Read-in Center 9320 is at   0.000000 -1.799200  2.275500
      Read-in Center 9321 is at   0.000000 -1.746300  2.275500
      Read-in Center 9322 is at   0.000000 -1.693400  2.275500
      Read-in Center 9323 is at   0.000000 -1.640400  2.275500
      Read-in Center 9324 is at   0.000000 -1.587500  2.275500
      Read-in Center 9325 is at   0.000000 -1.534600  2.275500
      Read-in Center 9326 is at   0.000000 -1.481700  2.275500
      Read-in Center 9327 is at   0.000000 -1.428800  2.275500
      Read-in Center 9328 is at   0.000000 -1.375900  2.275500
      Read-in Center 9329 is at   0.000000 -1.322900  2.275500
      Read-in Center 9330 is at   0.000000 -1.270000  2.275500
      Read-in Center 9331 is at   0.000000 -1.217100  2.275500
      Read-in Center 9332 is at   0.000000 -1.164200  2.275500
      Read-in Center 9333 is at   0.000000 -1.111300  2.275500
      Read-in Center 9334 is at   0.000000 -1.058400  2.275500
      Read-in Center 9335 is at   0.000000 -1.005400  2.275500
      Read-in Center 9336 is at   0.000000 -0.952500  2.275500
      Read-in Center 9337 is at   0.000000 -0.899600  2.275500
      Read-in Center 9338 is at   0.000000 -0.846700  2.275500
      Read-in Center 9339 is at   0.000000 -0.793800  2.275500
      Read-in Center 9340 is at   0.000000 -0.740800  2.275500
      Read-in Center 9341 is at   0.000000 -0.687900  2.275500
      Read-in Center 9342 is at   0.000000 -0.635000  2.275500
      Read-in Center 9343 is at   0.000000 -0.582100  2.275500
      Read-in Center 9344 is at   0.000000 -0.529200  2.275500
      Read-in Center 9345 is at   0.000000 -0.476300  2.275500
      Read-in Center 9346 is at   0.000000 -0.423300  2.275500
      Read-in Center 9347 is at   0.000000 -0.370400  2.275500
      Read-in Center 9348 is at   0.000000 -0.317500  2.275500
      Read-in Center 9349 is at   0.000000 -0.264600  2.275500
      Read-in Center 9350 is at   0.000000 -0.211700  2.275500
      Read-in Center 9351 is at   0.000000 -0.158800  2.275500
      Read-in Center 9352 is at   0.000000 -0.105800  2.275500
      Read-in Center 9353 is at   0.000000 -0.052900  2.275500
      Read-in Center 9354 is at   0.000000  0.000000  2.275500
      Read-in Center 9355 is at   0.000000  0.052900  2.275500
      Read-in Center 9356 is at   0.000000  0.105800  2.275500
      Read-in Center 9357 is at   0.000000  0.158800  2.275500
      Read-in Center 9358 is at   0.000000  0.211700  2.275500
      Read-in Center 9359 is at   0.000000  0.264600  2.275500
      Read-in Center 9360 is at   0.000000  0.317500  2.275500
      Read-in Center 9361 is at   0.000000  0.370400  2.275500
      Read-in Center 9362 is at   0.000000  0.423300  2.275500
      Read-in Center 9363 is at   0.000000  0.476300  2.275500
      Read-in Center 9364 is at   0.000000  0.529200  2.275500
      Read-in Center 9365 is at   0.000000  0.582100  2.275500
      Read-in Center 9366 is at   0.000000  0.635000  2.275500
      Read-in Center 9367 is at   0.000000  0.687900  2.275500
      Read-in Center 9368 is at   0.000000  0.740800  2.275500
      Read-in Center 9369 is at   0.000000  0.793800  2.275500
      Read-in Center 9370 is at   0.000000  0.846700  2.275500
      Read-in Center 9371 is at   0.000000  0.899600  2.275500
      Read-in Center 9372 is at   0.000000  0.952500  2.275500
      Read-in Center 9373 is at   0.000000  1.005400  2.275500
      Read-in Center 9374 is at   0.000000  1.058400  2.275500
      Read-in Center 9375 is at   0.000000  1.111300  2.275500
      Read-in Center 9376 is at   0.000000  1.164200  2.275500
      Read-in Center 9377 is at   0.000000  1.217100  2.275500
      Read-in Center 9378 is at   0.000000  1.270000  2.275500
      Read-in Center 9379 is at   0.000000  1.322900  2.275500
      Read-in Center 9380 is at   0.000000  1.375900  2.275500
      Read-in Center 9381 is at   0.000000  1.428800  2.275500
      Read-in Center 9382 is at   0.000000  1.481700  2.275500
      Read-in Center 9383 is at   0.000000  1.534600  2.275500
      Read-in Center 9384 is at   0.000000  1.587500  2.275500
      Read-in Center 9385 is at   0.000000  1.640400  2.275500
      Read-in Center 9386 is at   0.000000  1.693400  2.275500
      Read-in Center 9387 is at   0.000000  1.746300  2.275500
      Read-in Center 9388 is at   0.000000  1.799200  2.275500
      Read-in Center 9389 is at   0.000000  1.852100  2.275500
      Read-in Center 9390 is at   0.000000  1.905000  2.275500
      Read-in Center 9391 is at   0.000000  1.958000  2.275500
      Read-in Center 9392 is at   0.000000  2.010900  2.275500
      Read-in Center 9393 is at   0.000000  2.063800  2.275500
      Read-in Center 9394 is at   0.000000  2.116700  2.275500
      Read-in Center 9395 is at   0.000000  2.169600  2.275500
      Read-in Center 9396 is at   0.000000  2.222500  2.275500
      Read-in Center 9397 is at   0.000000  2.275500  2.275500
      Read-in Center 9398 is at   0.000000  2.328400  2.275500
      Read-in Center 9399 is at   0.000000  2.381300  2.275500
      Read-in Center 9400 is at   0.000000  2.434200  2.275500
      Read-in Center 9401 is at   0.000000  2.487100  2.275500
      Read-in Center 9402 is at   0.000000  2.540000  2.275500
      Read-in Center 9403 is at   0.000000  2.593000  2.275500
      Read-in Center 9404 is at   0.000000 -2.645900  2.328400
      Read-in Center 9405 is at   0.000000 -2.593000  2.328400
      Read-in Center 9406 is at   0.000000 -2.540000  2.328400
      Read-in Center 9407 is at   0.000000 -2.487100  2.328400
      Read-in Center 9408 is at   0.000000 -2.434200  2.328400
      Read-in Center 9409 is at   0.000000 -2.381300  2.328400
      Read-in Center 9410 is at   0.000000 -2.328400  2.328400
      Read-in Center 9411 is at   0.000000 -2.275500  2.328400
      Read-in Center 9412 is at   0.000000 -2.222500  2.328400
      Read-in Center 9413 is at   0.000000 -2.169600  2.328400
      Read-in Center 9414 is at   0.000000 -2.116700  2.328400
      Read-in Center 9415 is at   0.000000 -2.063800  2.328400
      Read-in Center 9416 is at   0.000000 -2.010900  2.328400
      Read-in Center 9417 is at   0.000000 -1.958000  2.328400
      Read-in Center 9418 is at   0.000000 -1.905000  2.328400
      Read-in Center 9419 is at   0.000000 -1.852100  2.328400
      Read-in Center 9420 is at   0.000000 -1.799200  2.328400
      Read-in Center 9421 is at   0.000000 -1.746300  2.328400
      Read-in Center 9422 is at   0.000000 -1.693400  2.328400
      Read-in Center 9423 is at   0.000000 -1.640400  2.328400
      Read-in Center 9424 is at   0.000000 -1.587500  2.328400
      Read-in Center 9425 is at   0.000000 -1.534600  2.328400
      Read-in Center 9426 is at   0.000000 -1.481700  2.328400
      Read-in Center 9427 is at   0.000000 -1.428800  2.328400
      Read-in Center 9428 is at   0.000000 -1.375900  2.328400
      Read-in Center 9429 is at   0.000000 -1.322900  2.328400
      Read-in Center 9430 is at   0.000000 -1.270000  2.328400
      Read-in Center 9431 is at   0.000000 -1.217100  2.328400
      Read-in Center 9432 is at   0.000000 -1.164200  2.328400
      Read-in Center 9433 is at   0.000000 -1.111300  2.328400
      Read-in Center 9434 is at   0.000000 -1.058400  2.328400
      Read-in Center 9435 is at   0.000000 -1.005400  2.328400
      Read-in Center 9436 is at   0.000000 -0.952500  2.328400
      Read-in Center 9437 is at   0.000000 -0.899600  2.328400
      Read-in Center 9438 is at   0.000000 -0.846700  2.328400
      Read-in Center 9439 is at   0.000000 -0.793800  2.328400
      Read-in Center 9440 is at   0.000000 -0.740800  2.328400
      Read-in Center 9441 is at   0.000000 -0.687900  2.328400
      Read-in Center 9442 is at   0.000000 -0.635000  2.328400
      Read-in Center 9443 is at   0.000000 -0.582100  2.328400
      Read-in Center 9444 is at   0.000000 -0.529200  2.328400
      Read-in Center 9445 is at   0.000000 -0.476300  2.328400
      Read-in Center 9446 is at   0.000000 -0.423300  2.328400
      Read-in Center 9447 is at   0.000000 -0.370400  2.328400
      Read-in Center 9448 is at   0.000000 -0.317500  2.328400
      Read-in Center 9449 is at   0.000000 -0.264600  2.328400
      Read-in Center 9450 is at   0.000000 -0.211700  2.328400
      Read-in Center 9451 is at   0.000000 -0.158800  2.328400
      Read-in Center 9452 is at   0.000000 -0.105800  2.328400
      Read-in Center 9453 is at   0.000000 -0.052900  2.328400
      Read-in Center 9454 is at   0.000000  0.000000  2.328400
      Read-in Center 9455 is at   0.000000  0.052900  2.328400
      Read-in Center 9456 is at   0.000000  0.105800  2.328400
      Read-in Center 9457 is at   0.000000  0.158800  2.328400
      Read-in Center 9458 is at   0.000000  0.211700  2.328400
      Read-in Center 9459 is at   0.000000  0.264600  2.328400
      Read-in Center 9460 is at   0.000000  0.317500  2.328400
      Read-in Center 9461 is at   0.000000  0.370400  2.328400
      Read-in Center 9462 is at   0.000000  0.423300  2.328400
      Read-in Center 9463 is at   0.000000  0.476300  2.328400
      Read-in Center 9464 is at   0.000000  0.529200  2.328400
      Read-in Center 9465 is at   0.000000  0.582100  2.328400
      Read-in Center 9466 is at   0.000000  0.635000  2.328400
      Read-in Center 9467 is at   0.000000  0.687900  2.328400
      Read-in Center 9468 is at   0.000000  0.740800  2.328400
      Read-in Center 9469 is at   0.000000  0.793800  2.328400
      Read-in Center 9470 is at   0.000000  0.846700  2.328400
      Read-in Center 9471 is at   0.000000  0.899600  2.328400
      Read-in Center 9472 is at   0.000000  0.952500  2.328400
      Read-in Center 9473 is at   0.000000  1.005400  2.328400
      Read-in Center 9474 is at   0.000000  1.058400  2.328400
      Read-in Center 9475 is at   0.000000  1.111300  2.328400
      Read-in Center 9476 is at   0.000000  1.164200  2.328400
      Read-in Center 9477 is at   0.000000  1.217100  2.328400
      Read-in Center 9478 is at   0.000000  1.270000  2.328400
      Read-in Center 9479 is at   0.000000  1.322900  2.328400
      Read-in Center 9480 is at   0.000000  1.375900  2.328400
      Read-in Center 9481 is at   0.000000  1.428800  2.328400
      Read-in Center 9482 is at   0.000000  1.481700  2.328400
      Read-in Center 9483 is at   0.000000  1.534600  2.328400
      Read-in Center 9484 is at   0.000000  1.587500  2.328400
      Read-in Center 9485 is at   0.000000  1.640400  2.328400
      Read-in Center 9486 is at   0.000000  1.693400  2.328400
      Read-in Center 9487 is at   0.000000  1.746300  2.328400
      Read-in Center 9488 is at   0.000000  1.799200  2.328400
      Read-in Center 9489 is at   0.000000  1.852100  2.328400
      Read-in Center 9490 is at   0.000000  1.905000  2.328400
      Read-in Center 9491 is at   0.000000  1.958000  2.328400
      Read-in Center 9492 is at   0.000000  2.010900  2.328400
      Read-in Center 9493 is at   0.000000  2.063800  2.328400
      Read-in Center 9494 is at   0.000000  2.116700  2.328400
      Read-in Center 9495 is at   0.000000  2.169600  2.328400
      Read-in Center 9496 is at   0.000000  2.222500  2.328400
      Read-in Center 9497 is at   0.000000  2.275500  2.328400
      Read-in Center 9498 is at   0.000000  2.328400  2.328400
      Read-in Center 9499 is at   0.000000  2.381300  2.328400
      Read-in Center 9500 is at   0.000000  2.434200  2.328400
      Read-in Center 9501 is at   0.000000  2.487100  2.328400
      Read-in Center 9502 is at   0.000000  2.540000  2.328400
      Read-in Center 9503 is at   0.000000  2.593000  2.328400
      Read-in Center 9504 is at   0.000000 -2.645900  2.381300
      Read-in Center 9505 is at   0.000000 -2.593000  2.381300
      Read-in Center 9506 is at   0.000000 -2.540000  2.381300
      Read-in Center 9507 is at   0.000000 -2.487100  2.381300
      Read-in Center 9508 is at   0.000000 -2.434200  2.381300
      Read-in Center 9509 is at   0.000000 -2.381300  2.381300
      Read-in Center 9510 is at   0.000000 -2.328400  2.381300
      Read-in Center 9511 is at   0.000000 -2.275500  2.381300
      Read-in Center 9512 is at   0.000000 -2.222500  2.381300
      Read-in Center 9513 is at   0.000000 -2.169600  2.381300
      Read-in Center 9514 is at   0.000000 -2.116700  2.381300
      Read-in Center 9515 is at   0.000000 -2.063800  2.381300
      Read-in Center 9516 is at   0.000000 -2.010900  2.381300
      Read-in Center 9517 is at   0.000000 -1.958000  2.381300
      Read-in Center 9518 is at   0.000000 -1.905000  2.381300
      Read-in Center 9519 is at   0.000000 -1.852100  2.381300
      Read-in Center 9520 is at   0.000000 -1.799200  2.381300
      Read-in Center 9521 is at   0.000000 -1.746300  2.381300
      Read-in Center 9522 is at   0.000000 -1.693400  2.381300
      Read-in Center 9523 is at   0.000000 -1.640400  2.381300
      Read-in Center 9524 is at   0.000000 -1.587500  2.381300
      Read-in Center 9525 is at   0.000000 -1.534600  2.381300
      Read-in Center 9526 is at   0.000000 -1.481700  2.381300
      Read-in Center 9527 is at   0.000000 -1.428800  2.381300
      Read-in Center 9528 is at   0.000000 -1.375900  2.381300
      Read-in Center 9529 is at   0.000000 -1.322900  2.381300
      Read-in Center 9530 is at   0.000000 -1.270000  2.381300
      Read-in Center 9531 is at   0.000000 -1.217100  2.381300
      Read-in Center 9532 is at   0.000000 -1.164200  2.381300
      Read-in Center 9533 is at   0.000000 -1.111300  2.381300
      Read-in Center 9534 is at   0.000000 -1.058400  2.381300
      Read-in Center 9535 is at   0.000000 -1.005400  2.381300
      Read-in Center 9536 is at   0.000000 -0.952500  2.381300
      Read-in Center 9537 is at   0.000000 -0.899600  2.381300
      Read-in Center 9538 is at   0.000000 -0.846700  2.381300
      Read-in Center 9539 is at   0.000000 -0.793800  2.381300
      Read-in Center 9540 is at   0.000000 -0.740800  2.381300
      Read-in Center 9541 is at   0.000000 -0.687900  2.381300
      Read-in Center 9542 is at   0.000000 -0.635000  2.381300
      Read-in Center 9543 is at   0.000000 -0.582100  2.381300
      Read-in Center 9544 is at   0.000000 -0.529200  2.381300
      Read-in Center 9545 is at   0.000000 -0.476300  2.381300
      Read-in Center 9546 is at   0.000000 -0.423300  2.381300
      Read-in Center 9547 is at   0.000000 -0.370400  2.381300
      Read-in Center 9548 is at   0.000000 -0.317500  2.381300
      Read-in Center 9549 is at   0.000000 -0.264600  2.381300
      Read-in Center 9550 is at   0.000000 -0.211700  2.381300
      Read-in Center 9551 is at   0.000000 -0.158800  2.381300
      Read-in Center 9552 is at   0.000000 -0.105800  2.381300
      Read-in Center 9553 is at   0.000000 -0.052900  2.381300
      Read-in Center 9554 is at   0.000000  0.000000  2.381300
      Read-in Center 9555 is at   0.000000  0.052900  2.381300
      Read-in Center 9556 is at   0.000000  0.105800  2.381300
      Read-in Center 9557 is at   0.000000  0.158800  2.381300
      Read-in Center 9558 is at   0.000000  0.211700  2.381300
      Read-in Center 9559 is at   0.000000  0.264600  2.381300
      Read-in Center 9560 is at   0.000000  0.317500  2.381300
      Read-in Center 9561 is at   0.000000  0.370400  2.381300
      Read-in Center 9562 is at   0.000000  0.423300  2.381300
      Read-in Center 9563 is at   0.000000  0.476300  2.381300
      Read-in Center 9564 is at   0.000000  0.529200  2.381300
      Read-in Center 9565 is at   0.000000  0.582100  2.381300
      Read-in Center 9566 is at   0.000000  0.635000  2.381300
      Read-in Center 9567 is at   0.000000  0.687900  2.381300
      Read-in Center 9568 is at   0.000000  0.740800  2.381300
      Read-in Center 9569 is at   0.000000  0.793800  2.381300
      Read-in Center 9570 is at   0.000000  0.846700  2.381300
      Read-in Center 9571 is at   0.000000  0.899600  2.381300
      Read-in Center 9572 is at   0.000000  0.952500  2.381300
      Read-in Center 9573 is at   0.000000  1.005400  2.381300
      Read-in Center 9574 is at   0.000000  1.058400  2.381300
      Read-in Center 9575 is at   0.000000  1.111300  2.381300
      Read-in Center 9576 is at   0.000000  1.164200  2.381300
      Read-in Center 9577 is at   0.000000  1.217100  2.381300
      Read-in Center 9578 is at   0.000000  1.270000  2.381300
      Read-in Center 9579 is at   0.000000  1.322900  2.381300
      Read-in Center 9580 is at   0.000000  1.375900  2.381300
      Read-in Center 9581 is at   0.000000  1.428800  2.381300
      Read-in Center 9582 is at   0.000000  1.481700  2.381300
      Read-in Center 9583 is at   0.000000  1.534600  2.381300
      Read-in Center 9584 is at   0.000000  1.587500  2.381300
      Read-in Center 9585 is at   0.000000  1.640400  2.381300
      Read-in Center 9586 is at   0.000000  1.693400  2.381300
      Read-in Center 9587 is at   0.000000  1.746300  2.381300
      Read-in Center 9588 is at   0.000000  1.799200  2.381300
      Read-in Center 9589 is at   0.000000  1.852100  2.381300
      Read-in Center 9590 is at   0.000000  1.905000  2.381300
      Read-in Center 9591 is at   0.000000  1.958000  2.381300
      Read-in Center 9592 is at   0.000000  2.010900  2.381300
      Read-in Center 9593 is at   0.000000  2.063800  2.381300
      Read-in Center 9594 is at   0.000000  2.116700  2.381300
      Read-in Center 9595 is at   0.000000  2.169600  2.381300
      Read-in Center 9596 is at   0.000000  2.222500  2.381300
      Read-in Center 9597 is at   0.000000  2.275500  2.381300
      Read-in Center 9598 is at   0.000000  2.328400  2.381300
      Read-in Center 9599 is at   0.000000  2.381300  2.381300
      Read-in Center 9600 is at   0.000000  2.434200  2.381300
      Read-in Center 9601 is at   0.000000  2.487100  2.381300
      Read-in Center 9602 is at   0.000000  2.540000  2.381300
      Read-in Center 9603 is at   0.000000  2.593000  2.381300
      Read-in Center 9604 is at   0.000000 -2.645900  2.434200
      Read-in Center 9605 is at   0.000000 -2.593000  2.434200
      Read-in Center 9606 is at   0.000000 -2.540000  2.434200
      Read-in Center 9607 is at   0.000000 -2.487100  2.434200
      Read-in Center 9608 is at   0.000000 -2.434200  2.434200
      Read-in Center 9609 is at   0.000000 -2.381300  2.434200
      Read-in Center 9610 is at   0.000000 -2.328400  2.434200
      Read-in Center 9611 is at   0.000000 -2.275500  2.434200
      Read-in Center 9612 is at   0.000000 -2.222500  2.434200
      Read-in Center 9613 is at   0.000000 -2.169600  2.434200
      Read-in Center 9614 is at   0.000000 -2.116700  2.434200
      Read-in Center 9615 is at   0.000000 -2.063800  2.434200
      Read-in Center 9616 is at   0.000000 -2.010900  2.434200
      Read-in Center 9617 is at   0.000000 -1.958000  2.434200
      Read-in Center 9618 is at   0.000000 -1.905000  2.434200
      Read-in Center 9619 is at   0.000000 -1.852100  2.434200
      Read-in Center 9620 is at   0.000000 -1.799200  2.434200
      Read-in Center 9621 is at   0.000000 -1.746300  2.434200
      Read-in Center 9622 is at   0.000000 -1.693400  2.434200
      Read-in Center 9623 is at   0.000000 -1.640400  2.434200
      Read-in Center 9624 is at   0.000000 -1.587500  2.434200
      Read-in Center 9625 is at   0.000000 -1.534600  2.434200
      Read-in Center 9626 is at   0.000000 -1.481700  2.434200
      Read-in Center 9627 is at   0.000000 -1.428800  2.434200
      Read-in Center 9628 is at   0.000000 -1.375900  2.434200
      Read-in Center 9629 is at   0.000000 -1.322900  2.434200
      Read-in Center 9630 is at   0.000000 -1.270000  2.434200
      Read-in Center 9631 is at   0.000000 -1.217100  2.434200
      Read-in Center 9632 is at   0.000000 -1.164200  2.434200
      Read-in Center 9633 is at   0.000000 -1.111300  2.434200
      Read-in Center 9634 is at   0.000000 -1.058400  2.434200
      Read-in Center 9635 is at   0.000000 -1.005400  2.434200
      Read-in Center 9636 is at   0.000000 -0.952500  2.434200
      Read-in Center 9637 is at   0.000000 -0.899600  2.434200
      Read-in Center 9638 is at   0.000000 -0.846700  2.434200
      Read-in Center 9639 is at   0.000000 -0.793800  2.434200
      Read-in Center 9640 is at   0.000000 -0.740800  2.434200
      Read-in Center 9641 is at   0.000000 -0.687900  2.434200
      Read-in Center 9642 is at   0.000000 -0.635000  2.434200
      Read-in Center 9643 is at   0.000000 -0.582100  2.434200
      Read-in Center 9644 is at   0.000000 -0.529200  2.434200
      Read-in Center 9645 is at   0.000000 -0.476300  2.434200
      Read-in Center 9646 is at   0.000000 -0.423300  2.434200
      Read-in Center 9647 is at   0.000000 -0.370400  2.434200
      Read-in Center 9648 is at   0.000000 -0.317500  2.434200
      Read-in Center 9649 is at   0.000000 -0.264600  2.434200
      Read-in Center 9650 is at   0.000000 -0.211700  2.434200
      Read-in Center 9651 is at   0.000000 -0.158800  2.434200
      Read-in Center 9652 is at   0.000000 -0.105800  2.434200
      Read-in Center 9653 is at   0.000000 -0.052900  2.434200
      Read-in Center 9654 is at   0.000000  0.000000  2.434200
      Read-in Center 9655 is at   0.000000  0.052900  2.434200
      Read-in Center 9656 is at   0.000000  0.105800  2.434200
      Read-in Center 9657 is at   0.000000  0.158800  2.434200
      Read-in Center 9658 is at   0.000000  0.211700  2.434200
      Read-in Center 9659 is at   0.000000  0.264600  2.434200
      Read-in Center 9660 is at   0.000000  0.317500  2.434200
      Read-in Center 9661 is at   0.000000  0.370400  2.434200
      Read-in Center 9662 is at   0.000000  0.423300  2.434200
      Read-in Center 9663 is at   0.000000  0.476300  2.434200
      Read-in Center 9664 is at   0.000000  0.529200  2.434200
      Read-in Center 9665 is at   0.000000  0.582100  2.434200
      Read-in Center 9666 is at   0.000000  0.635000  2.434200
      Read-in Center 9667 is at   0.000000  0.687900  2.434200
      Read-in Center 9668 is at   0.000000  0.740800  2.434200
      Read-in Center 9669 is at   0.000000  0.793800  2.434200
      Read-in Center 9670 is at   0.000000  0.846700  2.434200
      Read-in Center 9671 is at   0.000000  0.899600  2.434200
      Read-in Center 9672 is at   0.000000  0.952500  2.434200
      Read-in Center 9673 is at   0.000000  1.005400  2.434200
      Read-in Center 9674 is at   0.000000  1.058400  2.434200
      Read-in Center 9675 is at   0.000000  1.111300  2.434200
      Read-in Center 9676 is at   0.000000  1.164200  2.434200
      Read-in Center 9677 is at   0.000000  1.217100  2.434200
      Read-in Center 9678 is at   0.000000  1.270000  2.434200
      Read-in Center 9679 is at   0.000000  1.322900  2.434200
      Read-in Center 9680 is at   0.000000  1.375900  2.434200
      Read-in Center 9681 is at   0.000000  1.428800  2.434200
      Read-in Center 9682 is at   0.000000  1.481700  2.434200
      Read-in Center 9683 is at   0.000000  1.534600  2.434200
      Read-in Center 9684 is at   0.000000  1.587500  2.434200
      Read-in Center 9685 is at   0.000000  1.640400  2.434200
      Read-in Center 9686 is at   0.000000  1.693400  2.434200
      Read-in Center 9687 is at   0.000000  1.746300  2.434200
      Read-in Center 9688 is at   0.000000  1.799200  2.434200
      Read-in Center 9689 is at   0.000000  1.852100  2.434200
      Read-in Center 9690 is at   0.000000  1.905000  2.434200
      Read-in Center 9691 is at   0.000000  1.958000  2.434200
      Read-in Center 9692 is at   0.000000  2.010900  2.434200
      Read-in Center 9693 is at   0.000000  2.063800  2.434200
      Read-in Center 9694 is at   0.000000  2.116700  2.434200
      Read-in Center 9695 is at   0.000000  2.169600  2.434200
      Read-in Center 9696 is at   0.000000  2.222500  2.434200
      Read-in Center 9697 is at   0.000000  2.275500  2.434200
      Read-in Center 9698 is at   0.000000  2.328400  2.434200
      Read-in Center 9699 is at   0.000000  2.381300  2.434200
      Read-in Center 9700 is at   0.000000  2.434200  2.434200
      Read-in Center 9701 is at   0.000000  2.487100  2.434200
      Read-in Center 9702 is at   0.000000  2.540000  2.434200
      Read-in Center 9703 is at   0.000000  2.593000  2.434200
      Read-in Center 9704 is at   0.000000 -2.645900  2.487100
      Read-in Center 9705 is at   0.000000 -2.593000  2.487100
      Read-in Center 9706 is at   0.000000 -2.540000  2.487100
      Read-in Center 9707 is at   0.000000 -2.487100  2.487100
      Read-in Center 9708 is at   0.000000 -2.434200  2.487100
      Read-in Center 9709 is at   0.000000 -2.381300  2.487100
      Read-in Center 9710 is at   0.000000 -2.328400  2.487100
      Read-in Center 9711 is at   0.000000 -2.275500  2.487100
      Read-in Center 9712 is at   0.000000 -2.222500  2.487100
      Read-in Center 9713 is at   0.000000 -2.169600  2.487100
      Read-in Center 9714 is at   0.000000 -2.116700  2.487100
      Read-in Center 9715 is at   0.000000 -2.063800  2.487100
      Read-in Center 9716 is at   0.000000 -2.010900  2.487100
      Read-in Center 9717 is at   0.000000 -1.958000  2.487100
      Read-in Center 9718 is at   0.000000 -1.905000  2.487100
      Read-in Center 9719 is at   0.000000 -1.852100  2.487100
      Read-in Center 9720 is at   0.000000 -1.799200  2.487100
      Read-in Center 9721 is at   0.000000 -1.746300  2.487100
      Read-in Center 9722 is at   0.000000 -1.693400  2.487100
      Read-in Center 9723 is at   0.000000 -1.640400  2.487100
      Read-in Center 9724 is at   0.000000 -1.587500  2.487100
      Read-in Center 9725 is at   0.000000 -1.534600  2.487100
      Read-in Center 9726 is at   0.000000 -1.481700  2.487100
      Read-in Center 9727 is at   0.000000 -1.428800  2.487100
      Read-in Center 9728 is at   0.000000 -1.375900  2.487100
      Read-in Center 9729 is at   0.000000 -1.322900  2.487100
      Read-in Center 9730 is at   0.000000 -1.270000  2.487100
      Read-in Center 9731 is at   0.000000 -1.217100  2.487100
      Read-in Center 9732 is at   0.000000 -1.164200  2.487100
      Read-in Center 9733 is at   0.000000 -1.111300  2.487100
      Read-in Center 9734 is at   0.000000 -1.058400  2.487100
      Read-in Center 9735 is at   0.000000 -1.005400  2.487100
      Read-in Center 9736 is at   0.000000 -0.952500  2.487100
      Read-in Center 9737 is at   0.000000 -0.899600  2.487100
      Read-in Center 9738 is at   0.000000 -0.846700  2.487100
      Read-in Center 9739 is at   0.000000 -0.793800  2.487100
      Read-in Center 9740 is at   0.000000 -0.740800  2.487100
      Read-in Center 9741 is at   0.000000 -0.687900  2.487100
      Read-in Center 9742 is at   0.000000 -0.635000  2.487100
      Read-in Center 9743 is at   0.000000 -0.582100  2.487100
      Read-in Center 9744 is at   0.000000 -0.529200  2.487100
      Read-in Center 9745 is at   0.000000 -0.476300  2.487100
      Read-in Center 9746 is at   0.000000 -0.423300  2.487100
      Read-in Center 9747 is at   0.000000 -0.370400  2.487100
      Read-in Center 9748 is at   0.000000 -0.317500  2.487100
      Read-in Center 9749 is at   0.000000 -0.264600  2.487100
      Read-in Center 9750 is at   0.000000 -0.211700  2.487100
      Read-in Center 9751 is at   0.000000 -0.158800  2.487100
      Read-in Center 9752 is at   0.000000 -0.105800  2.487100
      Read-in Center 9753 is at   0.000000 -0.052900  2.487100
      Read-in Center 9754 is at   0.000000  0.000000  2.487100
      Read-in Center 9755 is at   0.000000  0.052900  2.487100
      Read-in Center 9756 is at   0.000000  0.105800  2.487100
      Read-in Center 9757 is at   0.000000  0.158800  2.487100
      Read-in Center 9758 is at   0.000000  0.211700  2.487100
      Read-in Center 9759 is at   0.000000  0.264600  2.487100
      Read-in Center 9760 is at   0.000000  0.317500  2.487100
      Read-in Center 9761 is at   0.000000  0.370400  2.487100
      Read-in Center 9762 is at   0.000000  0.423300  2.487100
      Read-in Center 9763 is at   0.000000  0.476300  2.487100
      Read-in Center 9764 is at   0.000000  0.529200  2.487100
      Read-in Center 9765 is at   0.000000  0.582100  2.487100
      Read-in Center 9766 is at   0.000000  0.635000  2.487100
      Read-in Center 9767 is at   0.000000  0.687900  2.487100
      Read-in Center 9768 is at   0.000000  0.740800  2.487100
      Read-in Center 9769 is at   0.000000  0.793800  2.487100
      Read-in Center 9770 is at   0.000000  0.846700  2.487100
      Read-in Center 9771 is at   0.000000  0.899600  2.487100
      Read-in Center 9772 is at   0.000000  0.952500  2.487100
      Read-in Center 9773 is at   0.000000  1.005400  2.487100
      Read-in Center 9774 is at   0.000000  1.058400  2.487100
      Read-in Center 9775 is at   0.000000  1.111300  2.487100
      Read-in Center 9776 is at   0.000000  1.164200  2.487100
      Read-in Center 9777 is at   0.000000  1.217100  2.487100
      Read-in Center 9778 is at   0.000000  1.270000  2.487100
      Read-in Center 9779 is at   0.000000  1.322900  2.487100
      Read-in Center 9780 is at   0.000000  1.375900  2.487100
      Read-in Center 9781 is at   0.000000  1.428800  2.487100
      Read-in Center 9782 is at   0.000000  1.481700  2.487100
      Read-in Center 9783 is at   0.000000  1.534600  2.487100
      Read-in Center 9784 is at   0.000000  1.587500  2.487100
      Read-in Center 9785 is at   0.000000  1.640400  2.487100
      Read-in Center 9786 is at   0.000000  1.693400  2.487100
      Read-in Center 9787 is at   0.000000  1.746300  2.487100
      Read-in Center 9788 is at   0.000000  1.799200  2.487100
      Read-in Center 9789 is at   0.000000  1.852100  2.487100
      Read-in Center 9790 is at   0.000000  1.905000  2.487100
      Read-in Center 9791 is at   0.000000  1.958000  2.487100
      Read-in Center 9792 is at   0.000000  2.010900  2.487100
      Read-in Center 9793 is at   0.000000  2.063800  2.487100
      Read-in Center 9794 is at   0.000000  2.116700  2.487100
      Read-in Center 9795 is at   0.000000  2.169600  2.487100
      Read-in Center 9796 is at   0.000000  2.222500  2.487100
      Read-in Center 9797 is at   0.000000  2.275500  2.487100
      Read-in Center 9798 is at   0.000000  2.328400  2.487100
      Read-in Center 9799 is at   0.000000  2.381300  2.487100
      Read-in Center 9800 is at   0.000000  2.434200  2.487100
      Read-in Center 9801 is at   0.000000  2.487100  2.487100
      Read-in Center 9802 is at   0.000000  2.540000  2.487100
      Read-in Center 9803 is at   0.000000  2.593000  2.487100
      Read-in Center 9804 is at   0.000000 -2.645900  2.540000
      Read-in Center 9805 is at   0.000000 -2.593000  2.540000
      Read-in Center 9806 is at   0.000000 -2.540000  2.540000
      Read-in Center 9807 is at   0.000000 -2.487100  2.540000
      Read-in Center 9808 is at   0.000000 -2.434200  2.540000
      Read-in Center 9809 is at   0.000000 -2.381300  2.540000
      Read-in Center 9810 is at   0.000000 -2.328400  2.540000
      Read-in Center 9811 is at   0.000000 -2.275500  2.540000
      Read-in Center 9812 is at   0.000000 -2.222500  2.540000
      Read-in Center 9813 is at   0.000000 -2.169600  2.540000
      Read-in Center 9814 is at   0.000000 -2.116700  2.540000
      Read-in Center 9815 is at   0.000000 -2.063800  2.540000
      Read-in Center 9816 is at   0.000000 -2.010900  2.540000
      Read-in Center 9817 is at   0.000000 -1.958000  2.540000
      Read-in Center 9818 is at   0.000000 -1.905000  2.540000
      Read-in Center 9819 is at   0.000000 -1.852100  2.540000
      Read-in Center 9820 is at   0.000000 -1.799200  2.540000
      Read-in Center 9821 is at   0.000000 -1.746300  2.540000
      Read-in Center 9822 is at   0.000000 -1.693400  2.540000
      Read-in Center 9823 is at   0.000000 -1.640400  2.540000
      Read-in Center 9824 is at   0.000000 -1.587500  2.540000
      Read-in Center 9825 is at   0.000000 -1.534600  2.540000
      Read-in Center 9826 is at   0.000000 -1.481700  2.540000
      Read-in Center 9827 is at   0.000000 -1.428800  2.540000
      Read-in Center 9828 is at   0.000000 -1.375900  2.540000
      Read-in Center 9829 is at   0.000000 -1.322900  2.540000
      Read-in Center 9830 is at   0.000000 -1.270000  2.540000
      Read-in Center 9831 is at   0.000000 -1.217100  2.540000
      Read-in Center 9832 is at   0.000000 -1.164200  2.540000
      Read-in Center 9833 is at   0.000000 -1.111300  2.540000
      Read-in Center 9834 is at   0.000000 -1.058400  2.540000
      Read-in Center 9835 is at   0.000000 -1.005400  2.540000
      Read-in Center 9836 is at   0.000000 -0.952500  2.540000
      Read-in Center 9837 is at   0.000000 -0.899600  2.540000
      Read-in Center 9838 is at   0.000000 -0.846700  2.540000
      Read-in Center 9839 is at   0.000000 -0.793800  2.540000
      Read-in Center 9840 is at   0.000000 -0.740800  2.540000
      Read-in Center 9841 is at   0.000000 -0.687900  2.540000
      Read-in Center 9842 is at   0.000000 -0.635000  2.540000
      Read-in Center 9843 is at   0.000000 -0.582100  2.540000
      Read-in Center 9844 is at   0.000000 -0.529200  2.540000
      Read-in Center 9845 is at   0.000000 -0.476300  2.540000
      Read-in Center 9846 is at   0.000000 -0.423300  2.540000
      Read-in Center 9847 is at   0.000000 -0.370400  2.540000
      Read-in Center 9848 is at   0.000000 -0.317500  2.540000
      Read-in Center 9849 is at   0.000000 -0.264600  2.540000
      Read-in Center 9850 is at   0.000000 -0.211700  2.540000
      Read-in Center 9851 is at   0.000000 -0.158800  2.540000
      Read-in Center 9852 is at   0.000000 -0.105800  2.540000
      Read-in Center 9853 is at   0.000000 -0.052900  2.540000
      Read-in Center 9854 is at   0.000000  0.000000  2.540000
      Read-in Center 9855 is at   0.000000  0.052900  2.540000
      Read-in Center 9856 is at   0.000000  0.105800  2.540000
      Read-in Center 9857 is at   0.000000  0.158800  2.540000
      Read-in Center 9858 is at   0.000000  0.211700  2.540000
      Read-in Center 9859 is at   0.000000  0.264600  2.540000
      Read-in Center 9860 is at   0.000000  0.317500  2.540000
      Read-in Center 9861 is at   0.000000  0.370400  2.540000
      Read-in Center 9862 is at   0.000000  0.423300  2.540000
      Read-in Center 9863 is at   0.000000  0.476300  2.540000
      Read-in Center 9864 is at   0.000000  0.529200  2.540000
      Read-in Center 9865 is at   0.000000  0.582100  2.540000
      Read-in Center 9866 is at   0.000000  0.635000  2.540000
      Read-in Center 9867 is at   0.000000  0.687900  2.540000
      Read-in Center 9868 is at   0.000000  0.740800  2.540000
      Read-in Center 9869 is at   0.000000  0.793800  2.540000
      Read-in Center 9870 is at   0.000000  0.846700  2.540000
      Read-in Center 9871 is at   0.000000  0.899600  2.540000
      Read-in Center 9872 is at   0.000000  0.952500  2.540000
      Read-in Center 9873 is at   0.000000  1.005400  2.540000
      Read-in Center 9874 is at   0.000000  1.058400  2.540000
      Read-in Center 9875 is at   0.000000  1.111300  2.540000
      Read-in Center 9876 is at   0.000000  1.164200  2.540000
      Read-in Center 9877 is at   0.000000  1.217100  2.540000
      Read-in Center 9878 is at   0.000000  1.270000  2.540000
      Read-in Center 9879 is at   0.000000  1.322900  2.540000
      Read-in Center 9880 is at   0.000000  1.375900  2.540000
      Read-in Center 9881 is at   0.000000  1.428800  2.540000
      Read-in Center 9882 is at   0.000000  1.481700  2.540000
      Read-in Center 9883 is at   0.000000  1.534600  2.540000
      Read-in Center 9884 is at   0.000000  1.587500  2.540000
      Read-in Center 9885 is at   0.000000  1.640400  2.540000
      Read-in Center 9886 is at   0.000000  1.693400  2.540000
      Read-in Center 9887 is at   0.000000  1.746300  2.540000
      Read-in Center 9888 is at   0.000000  1.799200  2.540000
      Read-in Center 9889 is at   0.000000  1.852100  2.540000
      Read-in Center 9890 is at   0.000000  1.905000  2.540000
      Read-in Center 9891 is at   0.000000  1.958000  2.540000
      Read-in Center 9892 is at   0.000000  2.010900  2.540000
      Read-in Center 9893 is at   0.000000  2.063800  2.540000
      Read-in Center 9894 is at   0.000000  2.116700  2.540000
      Read-in Center 9895 is at   0.000000  2.169600  2.540000
      Read-in Center 9896 is at   0.000000  2.222500  2.540000
      Read-in Center 9897 is at   0.000000  2.275500  2.540000
      Read-in Center 9898 is at   0.000000  2.328400  2.540000
      Read-in Center 9899 is at   0.000000  2.381300  2.540000
      Read-in Center 9900 is at   0.000000  2.434200  2.540000
      Read-in Center 9901 is at   0.000000  2.487100  2.540000
      Read-in Center 9902 is at   0.000000  2.540000  2.540000
      Read-in Center 9903 is at   0.000000  2.593000  2.540000
      Read-in Center 9904 is at   0.000000 -2.645900  2.593000
      Read-in Center 9905 is at   0.000000 -2.593000  2.593000
      Read-in Center 9906 is at   0.000000 -2.540000  2.593000
      Read-in Center 9907 is at   0.000000 -2.487100  2.593000
      Read-in Center 9908 is at   0.000000 -2.434200  2.593000
      Read-in Center 9909 is at   0.000000 -2.381300  2.593000
      Read-in Center 9910 is at   0.000000 -2.328400  2.593000
      Read-in Center 9911 is at   0.000000 -2.275500  2.593000
      Read-in Center 9912 is at   0.000000 -2.222500  2.593000
      Read-in Center 9913 is at   0.000000 -2.169600  2.593000
      Read-in Center 9914 is at   0.000000 -2.116700  2.593000
      Read-in Center 9915 is at   0.000000 -2.063800  2.593000
      Read-in Center 9916 is at   0.000000 -2.010900  2.593000
      Read-in Center 9917 is at   0.000000 -1.958000  2.593000
      Read-in Center 9918 is at   0.000000 -1.905000  2.593000
      Read-in Center 9919 is at   0.000000 -1.852100  2.593000
      Read-in Center 9920 is at   0.000000 -1.799200  2.593000
      Read-in Center 9921 is at   0.000000 -1.746300  2.593000
      Read-in Center 9922 is at   0.000000 -1.693400  2.593000
      Read-in Center 9923 is at   0.000000 -1.640400  2.593000
      Read-in Center 9924 is at   0.000000 -1.587500  2.593000
      Read-in Center 9925 is at   0.000000 -1.534600  2.593000
      Read-in Center 9926 is at   0.000000 -1.481700  2.593000
      Read-in Center 9927 is at   0.000000 -1.428800  2.593000
      Read-in Center 9928 is at   0.000000 -1.375900  2.593000
      Read-in Center 9929 is at   0.000000 -1.322900  2.593000
      Read-in Center 9930 is at   0.000000 -1.270000  2.593000
      Read-in Center 9931 is at   0.000000 -1.217100  2.593000
      Read-in Center 9932 is at   0.000000 -1.164200  2.593000
      Read-in Center 9933 is at   0.000000 -1.111300  2.593000
      Read-in Center 9934 is at   0.000000 -1.058400  2.593000
      Read-in Center 9935 is at   0.000000 -1.005400  2.593000
      Read-in Center 9936 is at   0.000000 -0.952500  2.593000
      Read-in Center 9937 is at   0.000000 -0.899600  2.593000
      Read-in Center 9938 is at   0.000000 -0.846700  2.593000
      Read-in Center 9939 is at   0.000000 -0.793800  2.593000
      Read-in Center 9940 is at   0.000000 -0.740800  2.593000
      Read-in Center 9941 is at   0.000000 -0.687900  2.593000
      Read-in Center 9942 is at   0.000000 -0.635000  2.593000
      Read-in Center 9943 is at   0.000000 -0.582100  2.593000
      Read-in Center 9944 is at   0.000000 -0.529200  2.593000
      Read-in Center 9945 is at   0.000000 -0.476300  2.593000
      Read-in Center 9946 is at   0.000000 -0.423300  2.593000
      Read-in Center 9947 is at   0.000000 -0.370400  2.593000
      Read-in Center 9948 is at   0.000000 -0.317500  2.593000
      Read-in Center 9949 is at   0.000000 -0.264600  2.593000
      Read-in Center 9950 is at   0.000000 -0.211700  2.593000
      Read-in Center 9951 is at   0.000000 -0.158800  2.593000
      Read-in Center 9952 is at   0.000000 -0.105800  2.593000
      Read-in Center 9953 is at   0.000000 -0.052900  2.593000
      Read-in Center 9954 is at   0.000000  0.000000  2.593000
      Read-in Center 9955 is at   0.000000  0.052900  2.593000
      Read-in Center 9956 is at   0.000000  0.105800  2.593000
      Read-in Center 9957 is at   0.000000  0.158800  2.593000
      Read-in Center 9958 is at   0.000000  0.211700  2.593000
      Read-in Center 9959 is at   0.000000  0.264600  2.593000
      Read-in Center 9960 is at   0.000000  0.317500  2.593000
      Read-in Center 9961 is at   0.000000  0.370400  2.593000
      Read-in Center 9962 is at   0.000000  0.423300  2.593000
      Read-in Center 9963 is at   0.000000  0.476300  2.593000
      Read-in Center 9964 is at   0.000000  0.529200  2.593000
      Read-in Center 9965 is at   0.000000  0.582100  2.593000
      Read-in Center 9966 is at   0.000000  0.635000  2.593000
      Read-in Center 9967 is at   0.000000  0.687900  2.593000
      Read-in Center 9968 is at   0.000000  0.740800  2.593000
      Read-in Center 9969 is at   0.000000  0.793800  2.593000
      Read-in Center 9970 is at   0.000000  0.846700  2.593000
      Read-in Center 9971 is at   0.000000  0.899600  2.593000
      Read-in Center 9972 is at   0.000000  0.952500  2.593000
      Read-in Center 9973 is at   0.000000  1.005400  2.593000
      Read-in Center 9974 is at   0.000000  1.058400  2.593000
      Read-in Center 9975 is at   0.000000  1.111300  2.593000
      Read-in Center 9976 is at   0.000000  1.164200  2.593000
      Read-in Center 9977 is at   0.000000  1.217100  2.593000
      Read-in Center 9978 is at   0.000000  1.270000  2.593000
      Read-in Center 9979 is at   0.000000  1.322900  2.593000
      Read-in Center 9980 is at   0.000000  1.375900  2.593000
      Read-in Center 9981 is at   0.000000  1.428800  2.593000
      Read-in Center 9982 is at   0.000000  1.481700  2.593000
      Read-in Center 9983 is at   0.000000  1.534600  2.593000
      Read-in Center 9984 is at   0.000000  1.587500  2.593000
      Read-in Center 9985 is at   0.000000  1.640400  2.593000
      Read-in Center 9986 is at   0.000000  1.693400  2.593000
      Read-in Center 9987 is at   0.000000  1.746300  2.593000
      Read-in Center 9988 is at   0.000000  1.799200  2.593000
      Read-in Center 9989 is at   0.000000  1.852100  2.593000
      Read-in Center 9990 is at   0.000000  1.905000  2.593000
      Read-in Center 9991 is at   0.000000  1.958000  2.593000
      Read-in Center 9992 is at   0.000000  2.010900  2.593000
      Read-in Center 9993 is at   0.000000  2.063800  2.593000
      Read-in Center 9994 is at   0.000000  2.116700  2.593000
      Read-in Center 9995 is at   0.000000  2.169600  2.593000
      Read-in Center 9996 is at   0.000000  2.222500  2.593000
      Read-in Center 9997 is at   0.000000  2.275500  2.593000
      Read-in Center 9998 is at   0.000000  2.328400  2.593000
      Read-in Center 9999 is at   0.000000  2.381300  2.593000
      Read-in Center **** is at   0.000000  2.434200  2.593000
      Read-in Center **** is at   0.000000  2.487100  2.593000
      Read-in Center **** is at   0.000000  2.540000  2.593000
      Read-in Center **** is at   0.000000  2.593000  2.593000
 -----------------------------------------------------------------

              Electrostatic Properties (Atomic Units)

 -----------------------------------------------------------------
    Center     Electric         -------- Electric Field --------
               Potential          X             Y             Z
 -----------------------------------------------------------------
    1 Atom    -22.332428
    2 Atom     -0.980704
    3 Atom     -0.980704
    4           0.015314
    5           0.015685
    6           0.016060
    7           0.016439
    8           0.016822
    9           0.017207
   10           0.017594
   11           0.017983
   12           0.018373
   13           0.018762
   14           0.019150
   15           0.019536
   16           0.019918
   17           0.020296
   18           0.020669
   19           0.021035
   20           0.021393
   21           0.021742
   22           0.022081
   23           0.022409
   24           0.022724
   25           0.023026
   26           0.023313
   27           0.023584
   28           0.023839
   29           0.024077
   30           0.024298
   31           0.024500
   32           0.024684
   33           0.024849
   34           0.024997
   35           0.025126
   36           0.025238
   37           0.025333
   38           0.025412
   39           0.025476
   40           0.025526
   41           0.025564
   42           0.025591
   43           0.025608
   44           0.025617
   45           0.025620
   46           0.025617
   47           0.025611
   48           0.025602
   49           0.025593
   50           0.025584
   51           0.025575
   52           0.025569
   53           0.025565
   54           0.025563
   55           0.025565
   56           0.025569
   57           0.025575
   58           0.025584
   59           0.025593
   60           0.025602
   61           0.025611
   62           0.025617
   63           0.025620
   64           0.025617
   65           0.025608
   66           0.025591
   67           0.025564
   68           0.025526
   69           0.025476
   70           0.025412
   71           0.025333
   72           0.025238
   73           0.025126
   74           0.024997
   75           0.024849
   76           0.024684
   77           0.024500
   78           0.024298
   79           0.024077
   80           0.023839
   81           0.023584
   82           0.023313
   83           0.023026
   84           0.022724
   85           0.022409
   86           0.022081
   87           0.021742
   88           0.021393
   89           0.021035
   90           0.020669
   91           0.020296
   92           0.019918
   93           0.019536
   94           0.019150
   95           0.018762
   96           0.018373
   97           0.017983
   98           0.017594
   99           0.017207
  100           0.016822
  101           0.016439
  102           0.016060
  103           0.015685
  104           0.015663
  105           0.016053
  106           0.016448
  107           0.016848
  108           0.017251
  109           0.017658
  110           0.018068
  111           0.018479
  112           0.018893
  113           0.019305
  114           0.019717
  115           0.020127
  116           0.020534
  117           0.020936
  118           0.021334
  119           0.021724
  120           0.022106
  121           0.022479
  122           0.022841
  123           0.023191
  124           0.023528
  125           0.023849
  126           0.024155
  127           0.024444
  128           0.024715
  129           0.024968
  130           0.025201
  131           0.025414
  132           0.025606
  133           0.025779
  134           0.025931
  135           0.026064
  136           0.026177
  137           0.026271
  138           0.026348
  139           0.026408
  140           0.026453
  141           0.026485
  142           0.026505
  143           0.026514
  144           0.026515
  145           0.026509
  146           0.026498
  147           0.026484
  148           0.026468
  149           0.026451
  150           0.026436
  151           0.026423
  152           0.026413
  153           0.026407
  154           0.026405
  155           0.026407
  156           0.026413
  157           0.026423
  158           0.026436
  159           0.026451
  160           0.026468
  161           0.026484
  162           0.026498
  163           0.026509
  164           0.026515
  165           0.026514
  166           0.026505
  167           0.026485
  168           0.026453
  169           0.026408
  170           0.026348
  171           0.026271
  172           0.026177
  173           0.026064
  174           0.025931
  175           0.025779
  176           0.025606
  177           0.025414
  178           0.025201
  179           0.024968
  180           0.024715
  181           0.024444
  182           0.024155
  183           0.023849
  184           0.023528
  185           0.023191
  186           0.022841
  187           0.022479
  188           0.022106
  189           0.021724
  190           0.021334
  191           0.020936
  192           0.020534
  193           0.020127
  194           0.019717
  195           0.019305
  196           0.018893
  197           0.018479
  198           0.018068
  199           0.017658
  200           0.017251
  201           0.016848
  202           0.016448
  203           0.016053
  204           0.016023
  205           0.016432
  206           0.016849
  207           0.017270
  208           0.017696
  209           0.018126
  210           0.018560
  211           0.018996
  212           0.019434
  213           0.019873
  214           0.020310
  215           0.020747
  216           0.021180
  217           0.021610
  218           0.022034
  219           0.022451
  220           0.022860
  221           0.023258
  222           0.023645
  223           0.024020
  224           0.024380
  225           0.024724
  226           0.025051
  227           0.025359
  228           0.025648
  229           0.025917
  230           0.026164
  231           0.026388
  232           0.026591
  233           0.026771
  234           0.026928
  235           0.027064
  236           0.027178
  237           0.027271
  238           0.027345
  239           0.027400
  240           0.027439
  241           0.027463
  242           0.027474
  243           0.027474
  244           0.027465
  245           0.027449
  246           0.027427
  247           0.027404
  248           0.027379
  249           0.027354
  250           0.027332
  251           0.027314
  252           0.027300
  253           0.027291
  254           0.027288
  255           0.027291
  256           0.027300
  257           0.027314
  258           0.027332
  259           0.027354
  260           0.027379
  261           0.027404
  262           0.027427
  263           0.027449
  264           0.027465
  265           0.027474
  266           0.027474
  267           0.027463
  268           0.027439
  269           0.027400
  270           0.027345
  271           0.027271
  272           0.027178
  273           0.027064
  274           0.026928
  275           0.026771
  276           0.026591
  277           0.026388
  278           0.026164
  279           0.025917
  280           0.025648
  281           0.025359
  282           0.025051
  283           0.024724
  284           0.024380
  285           0.024020
  286           0.023645
  287           0.023258
  288           0.022860
  289           0.022451
  290           0.022034
  291           0.021610
  292           0.021180
  293           0.020747
  294           0.020310
  295           0.019873
  296           0.019434
  297           0.018996
  298           0.018560
  299           0.018126
  300           0.017696
  301           0.017270
  302           0.016849
  303           0.016432
  304           0.016391
  305           0.016822
  306           0.017260
  307           0.017705
  308           0.018155
  309           0.018609
  310           0.019069
  311           0.019531
  312           0.019997
  313           0.020463
  314           0.020929
  315           0.021394
  316           0.021856
  317           0.022315
  318           0.022769
  319           0.023215
  320           0.023652
  321           0.024079
  322           0.024494
  323           0.024897
  324           0.025282
  325           0.025651
  326           0.026001
  327           0.026331
  328           0.026640
  329           0.026926
  330           0.027188
  331           0.027426
  332           0.027639
  333           0.027828
  334           0.027991
  335           0.028129
  336           0.028244
  337           0.028335
  338           0.028405
  339           0.028454
  340           0.028485
  341           0.028500
  342           0.028500
  343           0.028489
  344           0.028468
  345           0.028439
  346           0.028406
  347           0.028371
  348           0.028335
  349           0.028301
  350           0.028271
  351           0.028246
  352           0.028227
  353           0.028215
  354           0.028211
  355           0.028215
  356           0.028227
  357           0.028246
  358           0.028271
  359           0.028301
  360           0.028335
  361           0.028371
  362           0.028406
  363           0.028439
  364           0.028468
  365           0.028489
  366           0.028500
  367           0.028500
  368           0.028485
  369           0.028454
  370           0.028405
  371           0.028335
  372           0.028244
  373           0.028129
  374           0.027991
  375           0.027828
  376           0.027639
  377           0.027426
  378           0.027188
  379           0.026926
  380           0.026640
  381           0.026331
  382           0.026001
  383           0.025651
  384           0.025282
  385           0.024897
  386           0.024494
  387           0.024079
  388           0.023652
  389           0.023215
  390           0.022769
  391           0.022315
  392           0.021856
  393           0.021394
  394           0.020929
  395           0.020463
  396           0.019997
  397           0.019531
  398           0.019069
  399           0.018609
  400           0.018155
  401           0.017705
  402           0.017260
  403           0.016822
  404           0.016768
  405           0.017221
  406           0.017684
  407           0.018153
  408           0.018628
  409           0.019110
  410           0.019596
  411           0.020087
  412           0.020582
  413           0.021078
  414           0.021575
  415           0.022071
  416           0.022565
  417           0.023055
  418           0.023542
  419           0.024020
  420           0.024489
  421           0.024948
  422           0.025394
  423           0.025826
  424           0.026240
  425           0.026637
  426           0.027013
  427           0.027367
  428           0.027697
  429           0.028004
  430           0.028283
  431           0.028536
  432           0.028761
  433           0.028958
  434           0.029127
  435           0.029269
  436           0.029383
  437           0.029472
  438           0.029537
  439           0.029579
  440           0.029600
  441           0.029604
  442           0.029591
  443           0.029566
  444           0.029531
  445           0.029488
  446           0.029441
  447           0.029391
  448           0.029343
  449           0.029298
  450           0.029258
  451           0.029224
  452           0.029200
  453           0.029184
  454           0.029179
  455           0.029184
  456           0.029200
  457           0.029224
  458           0.029258
  459           0.029298
  460           0.029343
  461           0.029391
  462           0.029441
  463           0.029488
  464           0.029531
  465           0.029566
  466           0.029591
  467           0.029604
  468           0.029600
  469           0.029579
  470           0.029537
  471           0.029472
  472           0.029383
  473           0.029269
  474           0.029127
  475           0.028958
  476           0.028761
  477           0.028536
  478           0.028283
  479           0.028004
  480           0.027697
  481           0.027367
  482           0.027013
  483           0.026637
  484           0.026240
  485           0.025826
  486           0.025394
  487           0.024948
  488           0.024489
  489           0.024020
  490           0.023542
  491           0.023055
  492           0.022565
  493           0.022071
  494           0.021575
  495           0.021078
  496           0.020582
  497           0.020087
  498           0.019596
  499           0.019110
  500           0.018628
  501           0.018153
  502           0.017684
  503           0.017221
  504           0.017154
  505           0.017632
  506           0.018119
  507           0.018614
  508           0.019117
  509           0.019627
  510           0.020143
  511           0.020664
  512           0.021191
  513           0.021719
  514           0.022249
  515           0.022779
  516           0.023308
  517           0.023833
  518           0.024355
  519           0.024869
  520           0.025373
  521           0.025867
  522           0.026347
  523           0.026812
  524           0.027259
  525           0.027686
  526           0.028092
  527           0.028473
  528           0.028828
  529           0.029157
  530           0.029455
  531           0.029724
  532           0.029962
  533           0.030169
  534           0.030345
  535           0.030490
  536           0.030604
  537           0.030690
  538           0.030748
  539           0.030781
  540           0.030791
  541           0.030781
  542           0.030754
  543           0.030712
  544           0.030660
  545           0.030599
  546           0.030535
  547           0.030469
  548           0.030406
  549           0.030347
  550           0.030295
  551           0.030252
  552           0.030221
  553           0.030201
  554           0.030195
  555           0.030201
  556           0.030221
  557           0.030252
  558           0.030295
  559           0.030347
  560           0.030406
  561           0.030469
  562           0.030535
  563           0.030599
  564           0.030660
  565           0.030712
  566           0.030754
  567           0.030781
  568           0.030791
  569           0.030781
  570           0.030748
  571           0.030690
  572           0.030604
  573           0.030490
  574           0.030345
  575           0.030169
  576           0.029962
  577           0.029724
  578           0.029455
  579           0.029157
  580           0.028828
  581           0.028473
  582           0.028092
  583           0.027686
  584           0.027259
  585           0.026812
  586           0.026347
  587           0.025867
  588           0.025373
  589           0.024869
  590           0.024355
  591           0.023833
  592           0.023308
  593           0.022779
  594           0.022249
  595           0.021719
  596           0.021191
  597           0.020664
  598           0.020143
  599           0.019627
  600           0.019117
  601           0.018614
  602           0.018119
  603           0.017632
  604           0.017549
  605           0.018052
  606           0.018566
  607           0.019089
  608           0.019621
  609           0.020162
  610           0.020709
  611           0.021264
  612           0.021824
  613           0.022388
  614           0.022954
  615           0.023521
  616           0.024087
  617           0.024651
  618           0.025212
  619           0.025764
  620           0.026308
  621           0.026840
  622           0.027358
  623           0.027861
  624           0.028344
  625           0.028806
  626           0.029244
  627           0.029655
  628           0.030039
  629           0.030392
  630           0.030713
  631           0.031001
  632           0.031254
  633           0.031472
  634           0.031654
  635           0.031803
  636           0.031917
  637           0.031998
  638           0.032049
  639           0.032071
  640           0.032067
  641           0.032040
  642           0.031995
  643           0.031934
  644           0.031861
  645           0.031780
  646           0.031695
  647           0.031610
  648           0.031528
  649           0.031453
  650           0.031387
  651           0.031334
  652           0.031294
  653           0.031269
  654           0.031261
  655           0.031269
  656           0.031294
  657           0.031334
  658           0.031387
  659           0.031453
  660           0.031528
  661           0.031610
  662           0.031695
  663           0.031780
  664           0.031861
  665           0.031934
  666           0.031995
  667           0.032040
  668           0.032067
  669           0.032071
  670           0.032049
  671           0.031998
  672           0.031917
  673           0.031803
  674           0.031654
  675           0.031472
  676           0.031254
  677           0.031001
  678           0.030713
  679           0.030392
  680           0.030039
  681           0.029655
  682           0.029244
  683           0.028806
  684           0.028344
  685           0.027861
  686           0.027358
  687           0.026840
  688           0.026308
  689           0.025764
  690           0.025212
  691           0.024651
  692           0.024087
  693           0.023521
  694           0.022954
  695           0.022388
  696           0.021824
  697           0.021264
  698           0.020709
  699           0.020162
  700           0.019621
  701           0.019089
  702           0.018566
  703           0.018052
  704           0.017953
  705           0.018483
  706           0.019025
  707           0.019578
  708           0.020141
  709           0.020714
  710           0.021296
  711           0.021885
  712           0.022483
  713           0.023084
  714           0.023690
  715           0.024297
  716           0.024905
  717           0.025511
  718           0.026114
  719           0.026710
  720           0.027297
  721           0.027872
  722           0.028433
  723           0.028978
  724           0.029501
  725           0.030002
  726           0.030477
  727           0.030923
  728           0.031338
  729           0.031720
  730           0.032066
  731           0.032375
  732           0.032644
  733           0.032875
  734           0.033066
  735           0.033218
  736           0.033331
  737           0.033407
  738           0.033448
  739           0.033457
  740           0.033437
  741           0.033391
  742           0.033324
  743           0.033239
  744           0.033142
  745           0.033037
  746           0.032927
  747           0.032819
  748           0.032716
  749           0.032621
  750           0.032539
  751           0.032472
  752           0.032422
  753           0.032392
  754           0.032381
  755           0.032392
  756           0.032422
  757           0.032472
  758           0.032539
  759           0.032621
  760           0.032716
  761           0.032819
  762           0.032927
  763           0.033037
  764           0.033142
  765           0.033239
  766           0.033324
  767           0.033391
  768           0.033437
  769           0.033457
  770           0.033448
  771           0.033407
  772           0.033331
  773           0.033218
  774           0.033066
  775           0.032875
  776           0.032644
  777           0.032375
  778           0.032066
  779           0.031720
  780           0.031338
  781           0.030923
  782           0.030477
  783           0.030002
  784           0.029501
  785           0.028978
  786           0.028433
  787           0.027872
  788           0.027297
  789           0.026710
  790           0.026114
  791           0.025511
  792           0.024905
  793           0.024297
  794           0.023690
  795           0.023084
  796           0.022483
  797           0.021885
  798           0.021296
  799           0.020714
  800           0.020141
  801           0.019578
  802           0.019025
  803           0.018483
  804           0.018366
  805           0.018924
  806           0.019497
  807           0.020081
  808           0.020678
  809           0.021286
  810           0.021904
  811           0.022532
  812           0.023169
  813           0.023812
  814           0.024460
  815           0.025112
  816           0.025765
  817           0.026417
  818           0.027068
  819           0.027711
  820           0.028346
  821           0.028970
  822           0.029578
  823           0.030170
  824           0.030739
  825           0.031284
  826           0.031801
  827           0.032286
  828           0.032738
  829           0.033153
  830           0.033527
  831           0.033860
  832           0.034150
  833           0.034395
  834           0.034595
  835           0.034751
  836           0.034863
  837           0.034933
  838           0.034963
  839           0.034956
  840           0.034916
  841           0.034847
  842           0.034754
  843           0.034641
  844           0.034515
  845           0.034380
  846           0.034242
  847           0.034106
  848           0.033977
  849           0.033860
  850           0.033758
  851           0.033675
  852           0.033613
  853           0.033576
  854           0.033563
  855           0.033576
  856           0.033613
  857           0.033675
  858           0.033758
  859           0.033860
  860           0.033977
  861           0.034106
  862           0.034242
  863           0.034380
  864           0.034515
  865           0.034641
  866           0.034754
  867           0.034847
  868           0.034916
  869           0.034956
  870           0.034963
  871           0.034933
  872           0.034863
  873           0.034751
  874           0.034595
  875           0.034395
  876           0.034150
  877           0.033860
  878           0.033527
  879           0.033153
  880           0.032738
  881           0.032286
  882           0.031801
  883           0.031284
  884           0.030739
  885           0.030170
  886           0.029578
  887           0.028970
  888           0.028346
  889           0.027711
  890           0.027068
  891           0.026417
  892           0.025765
  893           0.025112
  894           0.024460
  895           0.023812
  896           0.023169
  897           0.022532
  898           0.021904
  899           0.021286
  900           0.020678
  901           0.020081
  902           0.019497
  903           0.018924
  904           0.018786
  905           0.019375
  906           0.019979
  907           0.020597
  908           0.021229
  909           0.021874
  910           0.022532
  911           0.023200
  912           0.023881
  913           0.024568
  914           0.025263
  915           0.025963
  916           0.026666
  917           0.027369
  918           0.028072
  919           0.028768
  920           0.029456
  921           0.030133
  922           0.030795
  923           0.031440
  924           0.032061
  925           0.032656
  926           0.033220
  927           0.033751
  928           0.034244
  929           0.034697
  930           0.035104
  931           0.035465
  932           0.035778
  933           0.036040
  934           0.036251
  935           0.036411
  936           0.036522
  937           0.036584
  938           0.036601
  939           0.036575
  940           0.036512
  941           0.036415
  942           0.036292
  943           0.036146
  944           0.035986
  945           0.035816
  946           0.035643
  947           0.035474
  948           0.035315
  949           0.035170
  950           0.035044
  951           0.034942
  952           0.034867
  953           0.034821
  954           0.034805
  955           0.034821
  956           0.034867
  957           0.034942
  958           0.035044
  959           0.035170
  960           0.035315
  961           0.035474
  962           0.035643
  963           0.035816
  964           0.035986
  965           0.036146
  966           0.036292
  967           0.036415
  968           0.036512
  969           0.036575
  970           0.036601
  971           0.036584
  972           0.036522
  973           0.036411
  974           0.036251
  975           0.036040
  976           0.035778
  977           0.035465
  978           0.035104
  979           0.034697
  980           0.034244
  981           0.033751
  982           0.033220
  983           0.032656
  984           0.032061
  985           0.031440
  986           0.030795
  987           0.030133
  988           0.029456
  989           0.028768
  990           0.028072
  991           0.027369
  992           0.026666
  993           0.025963
  994           0.025263
  995           0.024568
  996           0.023881
  997           0.023200
  998           0.022532
  999           0.021874
 1000           0.021229
 1001           0.020597
 1002           0.019979
 1003           0.019375
 1004           0.019214
 1005           0.019834
 1006           0.020472
 1007           0.021126
 1008           0.021796
 1009           0.022481
 1010           0.023180
 1011           0.023893
 1012           0.024620
 1013           0.025356
 1014           0.026101
 1015           0.026854
 1016           0.027611
 1017           0.028370
 1018           0.029131
 1019           0.029886
 1020           0.030634
 1021           0.031371
 1022           0.032093
 1023           0.032797
 1024           0.033477
 1025           0.034129
 1026           0.034748
 1027           0.035330
 1028           0.035871
 1029           0.036368
 1030           0.036814
 1031           0.037208
 1032           0.037547
 1033           0.037829
 1034           0.038053
 1035           0.038219
 1036           0.038328
 1037           0.038381
 1038           0.038382
 1039           0.038335
 1040           0.038244
 1041           0.038115
 1042           0.037955
 1043           0.037770
 1044           0.037569
 1045           0.037357
 1046           0.037143
 1047           0.036935
 1048           0.036739
 1049           0.036562
 1050           0.036408
 1051           0.036284
 1052           0.036192
 1053           0.036136
 1054           0.036117
 1055           0.036136
 1056           0.036192
 1057           0.036284
 1058           0.036408
 1059           0.036562
 1060           0.036739
 1061           0.036935
 1062           0.037143
 1063           0.037357
 1064           0.037569
 1065           0.037770
 1066           0.037955
 1067           0.038115
 1068           0.038244
 1069           0.038335
 1070           0.038382
 1071           0.038381
 1072           0.038328
 1073           0.038219
 1074           0.038053
 1075           0.037829
 1076           0.037547
 1077           0.037208
 1078           0.036814
 1079           0.036368
 1080           0.035871
 1081           0.035330
 1082           0.034748
 1083           0.034129
 1084           0.033477
 1085           0.032797
 1086           0.032093
 1087           0.031371
 1088           0.030634
 1089           0.029886
 1090           0.029131
 1091           0.028370
 1092           0.027611
 1093           0.026854
 1094           0.026101
 1095           0.025356
 1096           0.024620
 1097           0.023893
 1098           0.023180
 1099           0.022481
 1100           0.021796
 1101           0.021126
 1102           0.020472
 1103           0.019834
 1104           0.019649
 1105           0.020302
 1106           0.020976
 1107           0.021668
 1108           0.022377
 1109           0.023105
 1110           0.023850
 1111           0.024610
 1112           0.025387
 1113           0.026176
 1114           0.026976
 1115           0.027786
 1116           0.028603
 1117           0.029424
 1118           0.030248
 1119           0.031069
 1120           0.031883
 1121           0.032688
 1122           0.033478
 1123           0.034250
 1124           0.034996
 1125           0.035713
 1126           0.036395
 1127           0.037037
 1128           0.037634
 1129           0.038181
 1130           0.038673
 1131           0.039105
 1132           0.039476
 1133           0.039782
 1134           0.040022
 1135           0.040195
 1136           0.040303
 1137           0.040346
 1138           0.040329
 1139           0.040256
 1140           0.040132
 1141           0.039965
 1142           0.039761
 1143           0.039530
 1144           0.039279
 1145           0.039018
 1146           0.038756
 1147           0.038501
 1148           0.038261
 1149           0.038045
 1150           0.037858
 1151           0.037706
 1152           0.037595
 1153           0.037527
 1154           0.037504
 1155           0.037527
 1156           0.037595
 1157           0.037706
 1158           0.037858
 1159           0.038045
 1160           0.038261
 1161           0.038501
 1162           0.038756
 1163           0.039018
 1164           0.039279
 1165           0.039530
 1166           0.039761
 1167           0.039965
 1168           0.040132
 1169           0.040256
 1170           0.040329
 1171           0.040346
 1172           0.040303
 1173           0.040195
 1174           0.040022
 1175           0.039782
 1176           0.039476
 1177           0.039105
 1178           0.038673
 1179           0.038181
 1180           0.037634
 1181           0.037037
 1182           0.036395
 1183           0.035713
 1184           0.034996
 1185           0.034250
 1186           0.033478
 1187           0.032688
 1188           0.031883
 1189           0.031069
 1190           0.030248
 1191           0.029424
 1192           0.028603
 1193           0.027786
 1194           0.026976
 1195           0.026176
 1196           0.025387
 1197           0.024610
 1198           0.023850
 1199           0.023105
 1200           0.022377
 1201           0.021668
 1202           0.020976
 1203           0.020302
 1204           0.020089
 1205           0.020778
 1206           0.021490
 1207           0.022221
 1208           0.022974
 1209           0.023747
 1210           0.024540
 1211           0.025351
 1212           0.026182
 1213           0.027028
 1214           0.027888
 1215           0.028761
 1216           0.029643
 1217           0.030533
 1218           0.031428
 1219           0.032321
 1220           0.033210
 1221           0.034090
 1222           0.034957
 1223           0.035806
 1224           0.036629
 1225           0.037420
 1226           0.038175
 1227           0.038886
 1228           0.039548
 1229           0.040155
 1230           0.040700
 1231           0.041179
 1232           0.041588
 1233           0.041923
 1234           0.042182
 1235           0.042365
 1236           0.042471
 1237           0.042503
 1238           0.042465
 1239           0.042362
 1240           0.042200
 1241           0.041987
 1242           0.041732
 1243           0.041445
 1244           0.041136
 1245           0.040816
 1246           0.040495
 1247           0.040184
 1248           0.039892
 1249           0.039629
 1250           0.039403
 1251           0.039219
 1252           0.039084
 1253           0.039001
 1254           0.038974
 1255           0.039001
 1256           0.039084
 1257           0.039219
 1258           0.039403
 1259           0.039629
 1260           0.039892
 1261           0.040184
 1262           0.040495
 1263           0.040816
 1264           0.041136
 1265           0.041445
 1266           0.041732
 1267           0.041987
 1268           0.042200
 1269           0.042362
 1270           0.042465
 1271           0.042503
 1272           0.042471
 1273           0.042365
 1274           0.042182
 1275           0.041923
 1276           0.041588
 1277           0.041179
 1278           0.040700
 1279           0.040155
 1280           0.039548
 1281           0.038886
 1282           0.038175
 1283           0.037420
 1284           0.036629
 1285           0.035806
 1286           0.034957
 1287           0.034090
 1288           0.033210
 1289           0.032321
 1290           0.031428
 1291           0.030533
 1292           0.029643
 1293           0.028761
 1294           0.027888
 1295           0.027028
 1296           0.026182
 1297           0.025351
 1298           0.024540
 1299           0.023747
 1300           0.022974
 1301           0.022221
 1302           0.021490
 1303           0.020778
 1304           0.020535
 1305           0.021261
 1306           0.022012
 1307           0.022786
 1308           0.023584
 1309           0.024406
 1310           0.025250
 1311           0.026117
 1312           0.027006
 1313           0.027914
 1314           0.028839
 1315           0.029780
 1316           0.030735
 1317           0.031700
 1318           0.032673
 1319           0.033647
 1320           0.034620
 1321           0.035586
 1322           0.036539
 1323           0.037476
 1324           0.038386
 1325           0.039263
 1326           0.040102
 1327           0.040894
 1328           0.041633
 1329           0.042311
 1330           0.042919
 1331           0.043454
 1332           0.043908
 1333           0.044279
 1334           0.044562
 1335           0.044758
 1336           0.044864
 1337           0.044884
 1338           0.044822
 1339           0.044684
 1340           0.044476
 1341           0.044209
 1342           0.043894
 1343           0.043540
 1344           0.043162
 1345           0.042771
 1346           0.042379
 1347           0.042001
 1348           0.041648
 1349           0.041329
 1350           0.041055
 1351           0.040833
 1352           0.040669
 1353           0.040569
 1354           0.040536
 1355           0.040569
 1356           0.040669
 1357           0.040833
 1358           0.041055
 1359           0.041329
 1360           0.041648
 1361           0.042001
 1362           0.042379
 1363           0.042771
 1364           0.043162
 1365           0.043540
 1366           0.043894
 1367           0.044209
 1368           0.044476
 1369           0.044684
 1370           0.044822
 1371           0.044884
 1372           0.044864
 1373           0.044758
 1374           0.044562
 1375           0.044279
 1376           0.043908
 1377           0.043454
 1378           0.042919
 1379           0.042311
 1380           0.041633
 1381           0.040894
 1382           0.040102
 1383           0.039263
 1384           0.038386
 1385           0.037476
 1386           0.036539
 1387           0.035586
 1388           0.034620
 1389           0.033647
 1390           0.032673
 1391           0.031700
 1392           0.030735
 1393           0.029780
 1394           0.028839
 1395           0.027914
 1396           0.027006
 1397           0.026117
 1398           0.025250
 1399           0.024406
 1400           0.023584
 1401           0.022786
 1402           0.022012
 1403           0.021261
 1404           0.020986
 1405           0.021751
 1406           0.022544
 1407           0.023363
 1408           0.024209
 1409           0.025082
 1410           0.025982
 1411           0.026908
 1412           0.027860
 1413           0.028834
 1414           0.029831
 1415           0.030847
 1416           0.031881
 1417           0.032929
 1418           0.033990
 1419           0.035055
 1420           0.036121
 1421           0.037183
 1422           0.038236
 1423           0.039273
 1424           0.040283
 1425           0.041261
 1426           0.042198
 1427           0.043085
 1428           0.043914
 1429           0.044676
 1430           0.045362
 1431           0.045963
 1432           0.046474
 1433           0.046889
 1434           0.047203
 1435           0.047415
 1436           0.047524
 1437           0.047531
 1438           0.047442
 1439           0.047263
 1440           0.047002
 1441           0.046671
 1442           0.046282
 1443           0.045849
 1444           0.045388
 1445           0.044912
 1446           0.044436
 1447           0.043978
 1448           0.043549
 1449           0.043163
 1450           0.042832
 1451           0.042563
 1452           0.042366
 1453           0.042246
 1454           0.042205
 1455           0.042246
 1456           0.042366
 1457           0.042563
 1458           0.042832
 1459           0.043163
 1460           0.043549
 1461           0.043978
 1462           0.044436
 1463           0.044912
 1464           0.045388
 1465           0.045849
 1466           0.046282
 1467           0.046671
 1468           0.047002
 1469           0.047263
 1470           0.047442
 1471           0.047531
 1472           0.047524
 1473           0.047415
 1474           0.047203
 1475           0.046889
 1476           0.046474
 1477           0.045963
 1478           0.045362
 1479           0.044676
 1480           0.043914
 1481           0.043085
 1482           0.042198
 1483           0.041261
 1484           0.040283
 1485           0.039273
 1486           0.038236
 1487           0.037183
 1488           0.036121
 1489           0.035055
 1490           0.033990
 1491           0.032929
 1492           0.031881
 1493           0.030847
 1494           0.029831
 1495           0.028834
 1496           0.027860
 1497           0.026908
 1498           0.025982
 1499           0.025082
 1500           0.024209
 1501           0.023363
 1502           0.022544
 1503           0.021751
 1504           0.021440
 1505           0.022245
 1506           0.023082
 1507           0.023948
 1508           0.024845
 1509           0.025772
 1510           0.026731
 1511           0.027720
 1512           0.028740
 1513           0.029787
 1514           0.030860
 1515           0.031959
 1516           0.033079
 1517           0.034219
 1518           0.035376
 1519           0.036542
 1520           0.037714
 1521           0.038886
 1522           0.040050
 1523           0.041201
 1524           0.042328
 1525           0.043421
 1526           0.044473
 1527           0.045472
 1528           0.046408
 1529           0.047271
 1530           0.048049
 1531           0.048732
 1532           0.049313
 1533           0.049783
 1534           0.050137
 1535           0.050371
 1536           0.050485
 1537           0.050480
 1538           0.050360
 1539           0.050134
 1540           0.049810
 1541           0.049404
 1542           0.048928
 1543           0.048401
 1544           0.047839
 1545           0.047261
 1546           0.046685
 1547           0.046129
 1548           0.045611
 1549           0.045145
 1550           0.044745
 1551           0.044421
 1552           0.044183
 1553           0.044038
 1554           0.043989
 1555           0.044038
 1556           0.044183
 1557           0.044421
 1558           0.044745
 1559           0.045145
 1560           0.045611
 1561           0.046129
 1562           0.046685
 1563           0.047261
 1564           0.047839
 1565           0.048401
 1566           0.048928
 1567           0.049404
 1568           0.049810
 1569           0.050134
 1570           0.050360
 1571           0.050480
 1572           0.050485
 1573           0.050371
 1574           0.050137
 1575           0.049783
 1576           0.049313
 1577           0.048732
 1578           0.048049
 1579           0.047271
 1580           0.046408
 1581           0.045472
 1582           0.044473
 1583           0.043421
 1584           0.042328
 1585           0.041201
 1586           0.040050
 1587           0.038886
 1588           0.037714
 1589           0.036542
 1590           0.035376
 1591           0.034219
 1592           0.033079
 1593           0.031959
 1594           0.030860
 1595           0.029787
 1596           0.028740
 1597           0.027720
 1598           0.026731
 1599           0.025772
 1600           0.024845
 1601           0.023948
 1602           0.023082
 1603           0.022245
 1604           0.021894
 1605           0.022742
 1606           0.023624
 1607           0.024540
 1608           0.025491
 1609           0.026477
 1610           0.027498
 1611           0.028554
 1612           0.029648
 1613           0.030772
 1614           0.031930
 1615           0.033117
 1616           0.034333
 1617           0.035574
 1618           0.036839
 1619           0.038118
 1620           0.039408
 1621           0.040703
 1622           0.041995
 1623           0.043278
 1624           0.044538
 1625           0.045767
 1626           0.046953
 1627           0.048084
 1628           0.049148
 1629           0.050133
 1630           0.051023
 1631           0.051807
 1632           0.052473
 1633           0.053013
 1634           0.053419
 1635           0.053685
 1636           0.053808
 1637           0.053791
 1638           0.053638
 1639           0.053357
 1640           0.052960
 1641           0.052464
 1642           0.051885
 1643           0.051244
 1644           0.050562
 1645           0.049861
 1646           0.049163
 1647           0.048491
 1648           0.047865
 1649           0.047302
 1650           0.046818
 1651           0.046428
 1652           0.046141
 1653           0.045966
 1654           0.045907
 1655           0.045966
 1656           0.046141
 1657           0.046428
 1658           0.046818
 1659           0.047302
 1660           0.047865
 1661           0.048491
 1662           0.049163
 1663           0.049861
 1664           0.050562
 1665           0.051244
 1666           0.051885
 1667           0.052464
 1668           0.052960
 1669           0.053357
 1670           0.053638
 1671           0.053791
 1672           0.053808
 1673           0.053685
 1674           0.053419
 1675           0.053013
 1676           0.052473
 1677           0.051807
 1678           0.051023
 1679           0.050133
 1680           0.049148
 1681           0.048084
 1682           0.046953
 1683           0.045767
 1684           0.044538
 1685           0.043278
 1686           0.041995
 1687           0.040703
 1688           0.039408
 1689           0.038118
 1690           0.036839
 1691           0.035574
 1692           0.034333
 1693           0.033117
 1694           0.031930
 1695           0.030772
 1696           0.029648
 1697           0.028554
 1698           0.027498
 1699           0.026477
 1700           0.025491
 1701           0.024540
 1702           0.023624
 1703           0.022742
 1704           0.022349
 1705           0.023240
 1706           0.024171
 1707           0.025139
 1708           0.026146
 1709           0.027193
 1710           0.028281
 1711           0.029410
 1712           0.030581
 1713           0.031790
 1714           0.033038
 1715           0.034323
 1716           0.035644
 1717           0.036996
 1718           0.038380
 1719           0.039785
 1720           0.041209
 1721           0.042643
 1722           0.044081
 1723           0.045516
 1724           0.046931
 1725           0.048318
 1726           0.049663
 1727           0.050951
 1728           0.052168
 1729           0.053300
 1730           0.054327
 1731           0.055235
 1732           0.056010
 1733           0.056639
 1734           0.057111
 1735           0.057420
 1736           0.057561
 1737           0.057534
 1738           0.057344
 1739           0.057001
 1740           0.056519
 1741           0.055915
 1742           0.055213
 1743           0.054435
 1744           0.053608
 1745           0.052759
 1746           0.051913
 1747           0.051100
 1748           0.050343
 1749           0.049663
 1750           0.049080
 1751           0.048609
 1752           0.048263
 1753           0.048052
 1754           0.047982
 1755           0.048052
 1756           0.048263
 1757           0.048609
 1758           0.049080
 1759           0.049663
 1760           0.050343
 1761           0.051100
 1762           0.051913
 1763           0.052759
 1764           0.053608
 1765           0.054435
 1766           0.055213
 1767           0.055915
 1768           0.056519
 1769           0.057001
 1770           0.057344
 1771           0.057534
 1772           0.057561
 1773           0.057420
 1774           0.057111
 1775           0.056639
 1776           0.056010
 1777           0.055235
 1778           0.054327
 1779           0.053300
 1780           0.052168
 1781           0.050951
 1782           0.049663
 1783           0.048318
 1784           0.046931
 1785           0.045516
 1786           0.044081
 1787           0.042643
 1788           0.041209
 1789           0.039785
 1790           0.038380
 1791           0.036996
 1792           0.035644
 1793           0.034323
 1794           0.033038
 1795           0.031790
 1796           0.030581
 1797           0.029410
 1798           0.028281
 1799           0.027193
 1800           0.026146
 1801           0.025139
 1802           0.024171
 1803           0.023240
 1804           0.022801
 1805           0.023738
 1806           0.024719
 1807           0.025741
 1808           0.026808
 1809           0.027920
 1810           0.029079
 1811           0.030284
 1812           0.031539
 1813           0.032839
 1814           0.034185
 1815           0.035576
 1816           0.037011
 1817           0.038487
 1818           0.040003
 1819           0.041549
 1820           0.043122
 1821           0.044716
 1822           0.046321
 1823           0.047930
 1824           0.049525
 1825           0.051097
 1826           0.052629
 1827           0.054105
 1828           0.055507
 1829           0.056818
 1830           0.058013
 1831           0.059076
 1832           0.059987
 1833           0.060730
 1834           0.061291
 1835           0.061659
 1836           0.061826
 1837           0.061794
 1838           0.061566
 1839           0.061152
 1840           0.060568
 1841           0.059838
 1842           0.058987
 1843           0.058044
 1844           0.057041
 1845           0.056012
 1846           0.054988
 1847           0.054004
 1848           0.053087
 1849           0.052266
 1850           0.051561
 1851           0.050993
 1852           0.050576
 1853           0.050323
 1854           0.050237
 1855           0.050323
 1856           0.050576
 1857           0.050993
 1858           0.051561
 1859           0.052266
 1860           0.053087
 1861           0.054004
 1862           0.054988
 1863           0.056012
 1864           0.057041
 1865           0.058044
 1866           0.058987
 1867           0.059838
 1868           0.060568
 1869           0.061152
 1870           0.061566
 1871           0.061794
 1872           0.061826
 1873           0.061659
 1874           0.061291
 1875           0.060730
 1876           0.059987
 1877           0.059076
 1878           0.058013
 1879           0.056818
 1880           0.055507
 1881           0.054105
 1882           0.052629
 1883           0.051097
 1884           0.049525
 1885           0.047930
 1886           0.046321
 1887           0.044716
 1888           0.043122
 1889           0.041549
 1890           0.040003
 1891           0.038487
 1892           0.037011
 1893           0.035576
 1894           0.034185
 1895           0.032839
 1896           0.031539
 1897           0.030284
 1898           0.029079
 1899           0.027920
 1900           0.026808
 1901           0.025741
 1902           0.024719
 1903           0.023738
 1904           0.023251
 1905           0.024235
 1906           0.025268
 1907           0.026347
 1908           0.027476
 1909           0.028656
 1910           0.029890
 1911           0.031177
 1912           0.032521
 1913           0.033918
 1914           0.035370
 1915           0.036877
 1916           0.038438
 1917           0.040049
 1918           0.041713
 1919           0.043417
 1920           0.045159
 1921           0.046933
 1922           0.048729
 1923           0.050540
 1924           0.052346
 1925           0.054135
 1926           0.055890
 1927           0.057591
 1928           0.059217
 1929           0.060747
 1930           0.062151
 1931           0.063407
 1932           0.064491
 1933           0.065382
 1934           0.066059
 1935           0.066508
 1936           0.066719
 1937           0.066688
 1938           0.066420
 1939           0.065926
 1940           0.065224
 1941           0.064342
 1942           0.063311
 1943           0.062167
 1944           0.060950
 1945           0.059701
 1946           0.058458
 1947           0.057265
 1948           0.056155
 1949           0.055161
 1950           0.054309
 1951           0.053624
 1952           0.053120
 1953           0.052814
 1954           0.052711
 1955           0.052814
 1956           0.053120
 1957           0.053624
 1958           0.054309
 1959           0.055161
 1960           0.056155
 1961           0.057265
 1962           0.058458
 1963           0.059701
 1964           0.060950
 1965           0.062167
 1966           0.063311
 1967           0.064342
 1968           0.065224
 1969           0.065926
 1970           0.066420
 1971           0.066688
 1972           0.066719
 1973           0.066508
 1974           0.066059
 1975           0.065382
 1976           0.064491
 1977           0.063407
 1978           0.062151
 1979           0.060747
 1980           0.059217
 1981           0.057591
 1982           0.055890
 1983           0.054135
 1984           0.052346
 1985           0.050540
 1986           0.048729
 1987           0.046933
 1988           0.045159
 1989           0.043417
 1990           0.041713
 1991           0.040049
 1992           0.038438
 1993           0.036877
 1994           0.035370
 1995           0.033918
 1996           0.032521
 1997           0.031177
 1998           0.029890
 1999           0.028656
 2000           0.027476
 2001           0.026347
 2002           0.025268
 2003           0.024235
 2004           0.023693
 2005           0.024726
 2006           0.025812
 2007           0.026950
 2008           0.028144
 2009           0.029396
 2010           0.030708
 2011           0.032081
 2012           0.033521
 2013           0.035022
 2014           0.036589
 2015           0.038221
 2016           0.039918
 2017           0.041679
 2018           0.043504
 2019           0.045384
 2020           0.047316
 2021           0.049294
 2022           0.051309
 2023           0.053353
 2024           0.055403
 2025           0.057448
 2026           0.059468
 2027           0.061439
 2028           0.063336
 2029           0.065134
 2030           0.066796
 2031           0.068295
 2032           0.069600
 2033           0.070681
 2034           0.071512
 2035           0.072074
 2036           0.072348
 2037           0.072331
 2038           0.072023
 2039           0.071439
 2040           0.070598
 2041           0.069534
 2042           0.068286
 2043           0.066899
 2044           0.065420
 2045           0.063903
 2046           0.062393
 2047           0.060944
 2048           0.059598
 2049           0.058394
 2050           0.057364
 2051           0.056535
 2052           0.055928
 2053           0.055558
 2054           0.055434
 2055           0.055558
 2056           0.055928
 2057           0.056535
 2058           0.057364
 2059           0.058394
 2060           0.059598
 2061           0.060944
 2062           0.062393
 2063           0.063903
 2064           0.065420
 2065           0.066899
 2066           0.068286
 2067           0.069534
 2068           0.070598
 2069           0.071439
 2070           0.072023
 2071           0.072331
 2072           0.072348
 2073           0.072074
 2074           0.071512
 2075           0.070681
 2076           0.069600
 2077           0.068295
 2078           0.066796
 2079           0.065134
 2080           0.063336
 2081           0.061439
 2082           0.059468
 2083           0.057448
 2084           0.055403
 2085           0.053353
 2086           0.051309
 2087           0.049294
 2088           0.047316
 2089           0.045384
 2090           0.043504
 2091           0.041679
 2092           0.039918
 2093           0.038221
 2094           0.036589
 2095           0.035022
 2096           0.033521
 2097           0.032081
 2098           0.030708
 2099           0.029396
 2100           0.028144
 2101           0.026950
 2102           0.025812
 2103           0.024726
 2104           0.024126
 2105           0.025209
 2106           0.026350
 2107           0.027550
 2108           0.028811
 2109           0.030137
 2110           0.031532
 2111           0.032996
 2112           0.034537
 2113           0.036149
 2114           0.037838
 2115           0.039606
 2116           0.041452
 2117           0.043376
 2118           0.045381
 2119           0.047457
 2120           0.049603
 2121           0.051812
 2122           0.054077
 2123           0.056389
 2124           0.058725
 2125           0.061071
 2126           0.063405
 2127           0.065701
 2128           0.067928
 2129           0.070056
 2130           0.072040
 2131           0.073846
 2132           0.075433
 2133           0.076762
 2134           0.077799
 2135           0.078514
 2136           0.078884
 2137           0.078897
 2138           0.078555
 2139           0.077869
 2140           0.076865
 2141           0.075583
 2142           0.074071
 2143           0.072384
 2144           0.070584
 2145           0.068735
 2146           0.066897
 2147           0.065135
 2148           0.063500
 2149           0.062038
 2150           0.060791
 2151           0.059788
 2152           0.059054
 2153           0.058608
 2154           0.058459
 2155           0.058608
 2156           0.059054
 2157           0.059788
 2158           0.060791
 2159           0.062038
 2160           0.063500
 2161           0.065135
 2162           0.066897
 2163           0.068735
 2164           0.070584
 2165           0.072384
 2166           0.074071
 2167           0.075583
 2168           0.076865
 2169           0.077869
 2170           0.078555
 2171           0.078897
 2172           0.078884
 2173           0.078514
 2174           0.077799
 2175           0.076762
 2176           0.075433
 2177           0.073846
 2178           0.072040
 2179           0.070056
 2180           0.067928
 2181           0.065701
 2182           0.063405
 2183           0.061071
 2184           0.058725
 2185           0.056389
 2186           0.054077
 2187           0.051812
 2188           0.049603
 2189           0.047457
 2190           0.045381
 2191           0.043376
 2192           0.041452
 2193           0.039606
 2194           0.037838
 2195           0.036149
 2196           0.034537
 2197           0.032996
 2198           0.031532
 2199           0.030137
 2200           0.028811
 2201           0.027550
 2202           0.026350
 2203           0.025209
 2204           0.024548
 2205           0.025681
 2206           0.026879
 2207           0.028141
 2208           0.029472
 2209           0.030876
 2210           0.032357
 2211           0.033917
 2212           0.035564
 2213           0.037295
 2214           0.039115
 2215           0.041028
 2216           0.043036
 2217           0.045139
 2218           0.047342
 2219           0.049636
 2220           0.052021
 2221           0.054493
 2222           0.057043
 2223           0.059665
 2224           0.062334
 2225           0.065036
 2226           0.067745
 2227           0.070432
 2228           0.073061
 2229           0.075597
 2230           0.077983
 2231           0.080177
 2232           0.082126
 2233           0.083780
 2234           0.085090
 2235           0.086018
 2236           0.086526
 2237           0.086597
 2238           0.086228
 2239           0.085431
 2240           0.084235
 2241           0.082689
 2242           0.080852
 2243           0.078796
 2244           0.076598
 2245           0.074340
 2246           0.072095
 2247           0.069946
 2248           0.067955
 2249           0.066179
 2250           0.064665
 2251           0.063452
 2252           0.062564
 2253           0.062025
 2254           0.061845
 2255           0.062025
 2256           0.062564
 2257           0.063452
 2258           0.064665
 2259           0.066179
 2260           0.067955
 2261           0.069946
 2262           0.072095
 2263           0.074340
 2264           0.076598
 2265           0.078796
 2266           0.080852
 2267           0.082689
 2268           0.084235
 2269           0.085431
 2270           0.086228
 2271           0.086597
 2272           0.086526
 2273           0.086018
 2274           0.085090
 2275           0.083780
 2276           0.082126
 2277           0.080177
 2278           0.077983
 2279           0.075597
 2280           0.073061
 2281           0.070432
 2282           0.067745
 2283           0.065036
 2284           0.062334
 2285           0.059665
 2286           0.057043
 2287           0.054493
 2288           0.052021
 2289           0.049636
 2290           0.047342
 2291           0.045139
 2292           0.043036
 2293           0.041028
 2294           0.039115
 2295           0.037295
 2296           0.035564
 2297           0.033917
 2298           0.032357
 2299           0.030876
 2300           0.029472
 2301           0.028141
 2302           0.026879
 2303           0.025681
 2304           0.024955
 2305           0.026140
 2306           0.027395
 2307           0.028722
 2308           0.030124
 2309           0.031608
 2310           0.033178
 2311           0.034839
 2312           0.036598
 2313           0.038454
 2314           0.040414
 2315           0.042484
 2316           0.044666
 2317           0.046964
 2318           0.049385
 2319           0.051920
 2320           0.054574
 2321           0.057341
 2322           0.060217
 2323           0.063197
 2324           0.066254
 2325           0.069374
 2326           0.072530
 2327           0.075689
 2328           0.078809
 2329           0.081849
 2330           0.084740
 2331           0.087427
 2332           0.089844
 2333           0.091923
 2334           0.093601
 2335           0.094819
 2336           0.095527
 2337           0.095697
 2338           0.095317
 2339           0.094399
 2340           0.092975
 2341           0.091107
 2342           0.088870
 2343           0.086355
 2344           0.083661
 2345           0.080891
 2346           0.078140
 2347           0.075511
 2348           0.073081
 2349           0.070918
 2350           0.069080
 2351           0.067608
 2352           0.066535
 2353           0.065884
 2354           0.065665
 2355           0.065884
 2356           0.066535
 2357           0.067608
 2358           0.069080
 2359           0.070918
 2360           0.073081
 2361           0.075511
 2362           0.078140
 2363           0.080891
 2364           0.083661
 2365           0.086355
 2366           0.088870
 2367           0.091107
 2368           0.092975
 2369           0.094399
 2370           0.095317
 2371           0.095697
 2372           0.095527
 2373           0.094819
 2374           0.093601
 2375           0.091923
 2376           0.089844
 2377           0.087427
 2378           0.084740
 2379           0.081849
 2380           0.078809
 2381           0.075689
 2382           0.072530
 2383           0.069374
 2384           0.066254
 2385           0.063197
 2386           0.060217
 2387           0.057341
 2388           0.054574
 2389           0.051920
 2390           0.049385
 2391           0.046964
 2392           0.044666
 2393           0.042484
 2394           0.040414
 2395           0.038454
 2396           0.036598
 2397           0.034839
 2398           0.033178
 2399           0.031608
 2400           0.030124
 2401           0.028722
 2402           0.027395
 2403           0.026140
 2404           0.025344
 2405           0.026581
 2406           0.027895
 2407           0.029286
 2408           0.030762
 2409           0.032329
 2410           0.033991
 2411           0.035756
 2412           0.037632
 2413           0.039620
 2414           0.041729
 2415           0.043965
 2416           0.046335
 2417           0.048845
 2418           0.051504
 2419           0.054306
 2420           0.057258
 2421           0.060360
 2422           0.063606
 2423           0.066997
 2424           0.070505
 2425           0.074118
 2426           0.077806
 2427           0.081534
 2428           0.085255
 2429           0.088919
 2430           0.092444
 2431           0.095760
 2432           0.098783
 2433           0.101423
 2434           0.103592
 2435           0.105212
 2436           0.106208
 2437           0.106535
 2438           0.106171
 2439           0.105123
 2440           0.103430
 2441           0.101167
 2442           0.098431
 2443           0.095339
 2444           0.092021
 2445           0.088609
 2446           0.085224
 2447           0.081996
 2448           0.079021
 2449           0.076382
 2450           0.074146
 2451           0.072361
 2452           0.071062
 2453           0.070275
 2454           0.070012
 2455           0.070275
 2456           0.071062
 2457           0.072361
 2458           0.074146
 2459           0.076382
 2460           0.079021
 2461           0.081996
 2462           0.085224
 2463           0.088609
 2464           0.092021
 2465           0.095339
 2466           0.098431
 2467           0.101167
 2468           0.103430
 2469           0.105123
 2470           0.106171
 2471           0.106535
 2472           0.106208
 2473           0.105212
 2474           0.103592
 2475           0.101423
 2476           0.098783
 2477           0.095760
 2478           0.092444
 2479           0.088919
 2480           0.085255
 2481           0.081534
 2482           0.077806
 2483           0.074118
 2484           0.070505
 2485           0.066997
 2486           0.063606
 2487           0.060360
 2488           0.057258
 2489           0.054306
 2490           0.051504
 2491           0.048845
 2492           0.046335
 2493           0.043965
 2494           0.041729
 2495           0.039620
 2496           0.037632
 2497           0.035756
 2498           0.033991
 2499           0.032329
 2500           0.030762
 2501           0.029286
 2502           0.027895
 2503           0.026581
 2504           0.025712
 2505           0.027001
 2506           0.028374
 2507           0.029831
 2508           0.031382
 2509           0.033032
 2510           0.034790
 2511           0.036662
 2512           0.038661
 2513           0.040787
 2514           0.043053
 2515           0.045467
 2516           0.048039
 2517           0.050777
 2518           0.053697
 2519           0.056793
 2520           0.060077
 2521           0.063553
 2522           0.067221
 2523           0.071084
 2524           0.075117
 2525           0.079308
 2526           0.083632
 2527           0.088047
 2528           0.092504
 2529           0.096943
 2530           0.101267
 2531           0.105389
 2532           0.109200
 2533           0.112583
 2534           0.115418
 2535           0.117594
 2536           0.119003
 2537           0.119575
 2538           0.119269
 2539           0.118088
 2540           0.116072
 2541           0.113319
 2542           0.109952
 2543           0.106127
 2544           0.102014
 2545           0.097785
 2546           0.093599
 2547           0.089620
 2548           0.085966
 2549           0.082737
 2550           0.080012
 2551           0.077846
 2552           0.076273
 2553           0.075323
 2554           0.075005
 2555           0.075323
 2556           0.076273
 2557           0.077846
 2558           0.080012
 2559           0.082737
 2560           0.085966
 2561           0.089620
 2562           0.093599
 2563           0.097785
 2564           0.102014
 2565           0.106127
 2566           0.109952
 2567           0.113319
 2568           0.116072
 2569           0.118088
 2570           0.119269
 2571           0.119575
 2572           0.119003
 2573           0.117594
 2574           0.115418
 2575           0.112583
 2576           0.109200
 2577           0.105389
 2578           0.101267
 2579           0.096943
 2580           0.092504
 2581           0.088047
 2582           0.083632
 2583           0.079308
 2584           0.075117
 2585           0.071084
 2586           0.067221
 2587           0.063553
 2588           0.060077
 2589           0.056793
 2590           0.053697
 2591           0.050777
 2592           0.048039
 2593           0.045467
 2594           0.043053
 2595           0.040787
 2596           0.038661
 2597           0.036662
 2598           0.034790
 2599           0.033032
 2600           0.031382
 2601           0.029831
 2602           0.028374
 2603           0.027001
 2604           0.026054
 2605           0.027395
 2606           0.028826
 2607           0.030350
 2608           0.031975
 2609           0.033711
 2610           0.035566
 2611           0.037548
 2612           0.039673
 2613           0.041942
 2614           0.044372
 2615           0.046973
 2616           0.049759
 2617           0.052743
 2618           0.055943
 2619           0.059360
 2620           0.063010
 2621           0.066904
 2622           0.071046
 2623           0.075447
 2624           0.080085
 2625           0.084954
 2626           0.090028
 2627           0.095270
 2628           0.100623
 2629           0.106022
 2630           0.111352
 2631           0.116504
 2632           0.121341
 2633           0.125708
 2634           0.129443
 2635           0.132389
 2636           0.134389
 2637           0.135334
 2638           0.135156
 2639           0.133842
 2640           0.131440
 2641           0.128071
 2642           0.123900
 2643           0.119135
 2644           0.114001
 2645           0.108730
 2646           0.103527
 2647           0.098602
 2648           0.094101
 2649           0.090146
 2650           0.086824
 2651           0.084194
 2652           0.082293
 2653           0.081148
 2654           0.080765
 2655           0.081148
 2656           0.082293
 2657           0.084194
 2658           0.086824
 2659           0.090146
 2660           0.094101
 2661           0.098602
 2662           0.103527
 2663           0.108730
 2664           0.114001
 2665           0.119135
 2666           0.123900
 2667           0.128071
 2668           0.131440
 2669           0.133842
 2670           0.135156
 2671           0.135334
 2672           0.134389
 2673           0.132389
 2674           0.129443
 2675           0.125708
 2676           0.121341
 2677           0.116504
 2678           0.111352
 2679           0.106022
 2680           0.100623
 2681           0.095270
 2682           0.090028
 2683           0.084954
 2684           0.080085
 2685           0.075447
 2686           0.071046
 2687           0.066904
 2688           0.063010
 2689           0.059360
 2690           0.055943
 2691           0.052743
 2692           0.049759
 2693           0.046973
 2694           0.044372
 2695           0.041942
 2696           0.039673
 2697           0.037548
 2698           0.035566
 2699           0.033711
 2700           0.031975
 2701           0.030350
 2702           0.028826
 2703           0.027395
 2704           0.026368
 2705           0.027759
 2706           0.029248
 2707           0.030837
 2708           0.032538
 2709           0.034359
 2710           0.036312
 2711           0.038406
 2712           0.040660
 2713           0.043077
 2714           0.045677
 2715           0.048474
 2716           0.051487
 2717           0.054731
 2718           0.058233
 2719           0.061998
 2720           0.066050
 2721           0.070406
 2722           0.075080
 2723           0.080092
 2724           0.085425
 2725           0.091082
 2726           0.097045
 2727           0.103277
 2728           0.109723
 2729           0.116312
 2730           0.122909
 2731           0.129384
 2732           0.135562
 2733           0.141242
 2734           0.146202
 2735           0.150220
 2736           0.153068
 2737           0.154577
 2738           0.154633
 2739           0.153200
 2740           0.150327
 2741           0.146174
 2742           0.140960
 2743           0.134969
 2744           0.128510
 2745           0.121891
 2746           0.115385
 2747           0.109260
 2748           0.103699
 2749           0.098845
 2750           0.094794
 2751           0.091606
 2752           0.089312
 2753           0.087935
 2754           0.087476
 2755           0.087935
 2756           0.089312
 2757           0.091606
 2758           0.094794
 2759           0.098845
 2760           0.103699
 2761           0.109260
 2762           0.115385
 2763           0.121891
 2764           0.128510
 2765           0.134969
 2766           0.140960
 2767           0.146174
 2768           0.150327
 2769           0.153200
 2770           0.154633
 2771           0.154577
 2772           0.153068
 2773           0.150220
 2774           0.146202
 2775           0.141242
 2776           0.135562
 2777           0.129384
 2778           0.122909
 2779           0.116312
 2780           0.109723
 2781           0.103277
 2782           0.097045
 2783           0.091082
 2784           0.085425
 2785           0.080092
 2786           0.075080
 2787           0.070406
 2788           0.066050
 2789           0.061998
 2790           0.058233
 2791           0.054731
 2792           0.051487
 2793           0.048474
 2794           0.045677
 2795           0.043077
 2796           0.040660
 2797           0.038406
 2798           0.036312
 2799           0.034359
 2800           0.032538
 2801           0.030837
 2802           0.029248
 2803           0.027759
 2804           0.026649
 2805           0.028089
 2806           0.029635
 2807           0.031288
 2808           0.033063
 2809           0.034969
 2810           0.037020
 2811           0.039227
 2812           0.041612
 2813           0.044180
 2814           0.046955
 2815           0.049955
 2816           0.053204
 2817           0.056724
 2818           0.060548
 2819           0.064687
 2820           0.069176
 2821           0.074041
 2822           0.079307
 2823           0.085006
 2824           0.091133
 2825           0.097704
 2826           0.104711
 2827           0.112127
 2828           0.119899
 2829           0.127960
 2830           0.136154
 2831           0.144328
 2832           0.152266
 2833           0.159704
 2834           0.166339
 2835           0.171859
 2836           0.175929
 2837           0.178287
 2838           0.178744
 2839           0.177222
 2840           0.173770
 2841           0.168597
 2842           0.162008
 2843           0.154395
 2844           0.146188
 2845           0.137806
 2846           0.129617
 2847           0.121965
 2848           0.115074
 2849           0.109109
 2850           0.104172
 2851           0.100316
 2852           0.097557
 2853           0.095909
 2854           0.095360
 2855           0.095909
 2856           0.097557
 2857           0.100316
 2858           0.104172
 2859           0.109109
 2860           0.115074
 2861           0.121965
 2862           0.129617
 2863           0.137806
 2864           0.146188
 2865           0.154395
 2866           0.162008
 2867           0.168597
 2868           0.173770
 2869           0.177222
 2870           0.178744
 2871           0.178287
 2872           0.175929
 2873           0.171859
 2874           0.166339
 2875           0.159704
 2876           0.152266
 2877           0.144328
 2878           0.136154
 2879           0.127960
 2880           0.119899
 2881           0.112127
 2882           0.104711
 2883           0.097704
 2884           0.091133
 2885           0.085006
 2886           0.079307
 2887           0.074041
 2888           0.069176
 2889           0.064687
 2890           0.060548
 2891           0.056724
 2892           0.053204
 2893           0.049955
 2894           0.046955
 2895           0.044180
 2896           0.041612
 2897           0.039227
 2898           0.037020
 2899           0.034969
 2900           0.033063
 2901           0.031288
 2902           0.029635
 2903           0.028089
 2904           0.026893
 2905           0.028380
 2906           0.029980
 2907           0.031697
 2908           0.033544
 2909           0.035534
 2910           0.037681
 2911           0.040001
 2912           0.042516
 2913           0.045237
 2914           0.048190
 2915           0.051400
 2916           0.054893
 2917           0.058701
 2918           0.062864
 2919           0.067402
 2920           0.072361
 2921           0.077780
 2922           0.083698
 2923           0.090166
 2924           0.097191
 2925           0.104812
 2926           0.113037
 2927           0.121858
 2928           0.131236
 2929           0.141110
 2930           0.151315
 2931           0.161674
 2932           0.171924
 2933           0.181725
 2934           0.190667
 2935           0.198307
 2936           0.204148
 2937           0.207785
 2938           0.208903
 2939           0.207352
 2940           0.203174
 2941           0.196648
 2942           0.188207
 2943           0.178410
 2944           0.167864
 2945           0.157153
 2946           0.146772
 2947           0.137164
 2948           0.128603
 2949           0.121271
 2950           0.115265
 2951           0.110618
 2952           0.107318
 2953           0.105357
 2954           0.104706
 2955           0.105357
 2956           0.107318
 2957           0.110618
 2958           0.115265
 2959           0.121271
 2960           0.128603
 2961           0.137164
 2962           0.146772
 2963           0.157153
 2964           0.167864
 2965           0.178410
 2966           0.188207
 2967           0.196648
 2968           0.203174
 2969           0.207352
 2970           0.208903
 2971           0.207785
 2972           0.204148
 2973           0.198307
 2974           0.190667
 2975           0.181725
 2976           0.171924
 2977           0.161674
 2978           0.151315
 2979           0.141110
 2980           0.131236
 2981           0.121858
 2982           0.113037
 2983           0.104812
 2984           0.097191
 2985           0.090166
 2986           0.083698
 2987           0.077780
 2988           0.072361
 2989           0.067402
 2990           0.062864
 2991           0.058701
 2992           0.054893
 2993           0.051400
 2994           0.048190
 2995           0.045237
 2996           0.042516
 2997           0.040001
 2998           0.037681
 2999           0.035534
 3000           0.033544
 3001           0.031697
 3002           0.029980
 3003           0.028380
 3004           0.027096
 3005           0.028628
 3006           0.030280
 3007           0.032056
 3008           0.033973
 3009           0.036044
 3010           0.038286
 3011           0.040716
 3012           0.043361
 3013           0.046234
 3014           0.049367
 3015           0.052788
 3016           0.056532
 3017           0.060637
 3018           0.065153
 3019           0.070111
 3020           0.075570
 3021           0.081585
 3022           0.088214
 3023           0.095530
 3024           0.103563
 3025           0.112377
 3026           0.122013
 3027           0.132489
 3028           0.143795
 3029           0.155895
 3030           0.168621
 3031           0.181788
 3032           0.195081
 3033           0.208072
 3034           0.220210
 3035           0.230863
 3036           0.239294
 3037           0.244865
 3038           0.247055
 3039           0.245584
 3040           0.240474
 3041           0.232106
 3042           0.221115
 3043           0.208319
 3044           0.194596
 3045           0.180771
 3046           0.167513
 3047           0.155395
 3048           0.144739
 3049           0.135736
 3050           0.128456
 3051           0.122888
 3052           0.118973
 3053           0.116661
 3054           0.115897
 3055           0.116661
 3056           0.118973
 3057           0.122888
 3058           0.128456
 3059           0.135736
 3060           0.144739
 3061           0.155395
 3062           0.167513
 3063           0.180771
 3064           0.194596
 3065           0.208319
 3066           0.221115
 3067           0.232106
 3068           0.240474
 3069           0.245584
 3070           0.247055
 3071           0.244865
 3072           0.239294
 3073           0.230863
 3074           0.220210
 3075           0.208072
 3076           0.195081
 3077           0.181788
 3078           0.168621
 3079           0.155895
 3080           0.143795
 3081           0.132489
 3082           0.122013
 3083           0.112377
 3084           0.103563
 3085           0.095530
 3086           0.088214
 3087           0.081585
 3088           0.075570
 3089           0.070111
 3090           0.065153
 3091           0.060637
 3092           0.056532
 3093           0.052788
 3094           0.049367
 3095           0.046234
 3096           0.043361
 3097           0.040716
 3098           0.038286
 3099           0.036044
 3100           0.033973
 3101           0.032056
 3102           0.030280
 3103           0.028628
 3104           0.027254
 3105           0.028827
 3106           0.030528
 3107           0.032361
 3108           0.034343
 3109           0.036492
 3110           0.038825
 3111           0.041363
 3112           0.044135
 3113           0.047158
 3114           0.050468
 3115           0.054101
 3116           0.058098
 3117           0.062505
 3118           0.067385
 3119           0.072780
 3120           0.078765
 3121           0.085414
 3122           0.092808
 3123           0.101050
 3124           0.110198
 3125           0.120357
 3126           0.131609
 3127           0.144021
 3128           0.157628
 3129           0.172445
 3130           0.188325
 3131           0.205094
 3132           0.222404
 3133           0.239724
 3134           0.256323
 3135           0.271308
 3136           0.283575
 3137           0.292106
 3138           0.296028
 3139           0.294829
 3140           0.288479
 3141           0.277522
 3142           0.262922
 3143           0.245914
 3144           0.227798
 3145           0.209751
 3146           0.192687
 3147           0.177337
 3148           0.164062
 3149           0.153032
 3150           0.144257
 3151           0.137643
 3152           0.133048
 3153           0.130360
 3154           0.129475
 3155           0.130360
 3156           0.133048
 3157           0.137643
 3158           0.144257
 3159           0.153032
 3160           0.164062
 3161           0.177337
 3162           0.192687
 3163           0.209751
 3164           0.227798
 3165           0.245914
 3166           0.262922
 3167           0.277522
 3168           0.288479
 3169           0.294829
 3170           0.296028
 3171           0.292106
 3172           0.283575
 3173           0.271308
 3174           0.256323
 3175           0.239724
 3176           0.222404
 3177           0.205094
 3178           0.188325
 3179           0.172445
 3180           0.157628
 3181           0.144021
 3182           0.131609
 3183           0.120357
 3184           0.110198
 3185           0.101050
 3186           0.092808
 3187           0.085414
 3188           0.078765
 3189           0.072780
 3190           0.067385
 3191           0.062505
 3192           0.058098
 3193           0.054101
 3194           0.050468
 3195           0.047158
 3196           0.044135
 3197           0.041363
 3198           0.038825
 3199           0.036492
 3200           0.034343
 3201           0.032361
 3202           0.030528
 3203           0.028827
 3204           0.027362
 3205           0.028974
 3206           0.030718
 3207           0.032603
 3208           0.034647
 3209           0.036868
 3210           0.039287
 3211           0.041926
 3212           0.044820
 3213           0.047987
 3214           0.051470
 3215           0.055311
 3216           0.059558
 3217           0.064267
 3218           0.069514
 3219           0.075353
 3220           0.081879
 3221           0.089189
 3222           0.097391
 3223           0.106624
 3224           0.116984
 3225           0.128629
 3226           0.141700
 3227           0.156334
 3228           0.172643
 3229           0.190731
 3230           0.210510
 3231           0.231863
 3232           0.254439
 3233           0.277621
 3234           0.300464
 3235           0.321713
 3236           0.339708
 3237           0.352811
 3238           0.359537
 3239           0.358941
 3240           0.350883
 3241           0.336176
 3242           0.316339
 3243           0.293302
 3244           0.269031
 3245           0.245226
 3246           0.223126
 3247           0.203638
 3248           0.187129
 3249           0.173694
 3250           0.163218
 3251           0.155467
 3252           0.150167
 3253           0.147102
 3254           0.146099
 3255           0.147102
 3256           0.150167
 3257           0.155467
 3258           0.163218
 3259           0.173694
 3260           0.187129
 3261           0.203638
 3262           0.223126
 3263           0.245226
 3264           0.269031
 3265           0.293302
 3266           0.316339
 3267           0.336176
 3268           0.350883
 3269           0.358941
 3270           0.359537
 3271           0.352811
 3272           0.339708
 3273           0.321713
 3274           0.300464
 3275           0.277621
 3276           0.254439
 3277           0.231863
 3278           0.210510
 3279           0.190731
 3280           0.172643
 3281           0.156334
 3282           0.141700
 3283           0.128629
 3284           0.116984
 3285           0.106624
 3286           0.097391
 3287           0.089189
 3288           0.081879
 3289           0.075353
 3290           0.069514
 3291           0.064267
 3292           0.059558
 3293           0.055311
 3294           0.051470
 3295           0.047987
 3296           0.044820
 3297           0.041926
 3298           0.039287
 3299           0.036868
 3300           0.034647
 3301           0.032603
 3302           0.030718
 3303           0.028974
 3304           0.027418
 3305           0.029063
 3306           0.030846
 3307           0.032778
 3308           0.034878
 3309           0.037165
 3310           0.039662
 3311           0.042396
 3312           0.045403
 3313           0.048707
 3314           0.052355
 3315           0.056395
 3316           0.060884
 3317           0.065889
 3318           0.071499
 3319           0.077783
 3320           0.084857
 3321           0.092842
 3322           0.101881
 3323           0.112156
 3324           0.123811
 3325           0.137069
 3326           0.152154
 3327           0.169300
 3328           0.188739
 3329           0.210714
 3330           0.235264
 3331           0.262404
 3332           0.291863
 3333           0.322990
 3334           0.354623
 3335           0.385037
 3336           0.411737
 3337           0.432048
 3338           0.443381
 3339           0.443970
 3340           0.433433
 3341           0.413059
 3342           0.385361
 3343           0.353470
 3344           0.320423
 3345           0.288685
 3346           0.259905
 3347           0.235147
 3348           0.214702
 3349           0.198485
 3350           0.186154
 3351           0.177243
 3352           0.171275
 3353           0.167876
 3354           0.166773
 3355           0.167876
 3356           0.171275
 3357           0.177243
 3358           0.186154
 3359           0.198485
 3360           0.214702
 3361           0.235147
 3362           0.259905
 3363           0.288685
 3364           0.320423
 3365           0.353470
 3366           0.385361
 3367           0.413059
 3368           0.433433
 3369           0.443970
 3370           0.443381
 3371           0.432048
 3372           0.411737
 3373           0.385037
 3374           0.354623
 3375           0.322990
 3376           0.291863
 3377           0.262404
 3378           0.235264
 3379           0.210714
 3380           0.188739
 3381           0.169300
 3382           0.152154
 3383           0.137069
 3384           0.123811
 3385           0.112156
 3386           0.101881
 3387           0.092842
 3388           0.084857
 3389           0.077783
 3390           0.071499
 3391           0.065889
 3392           0.060884
 3393           0.056395
 3394           0.052355
 3395           0.048707
 3396           0.045403
 3397           0.042396
 3398           0.039662
 3399           0.037165
 3400           0.034878
 3401           0.032778
 3402           0.030846
 3403           0.029063
 3404           0.027418
 3405           0.029090
 3406           0.030907
 3407           0.032879
 3408           0.035028
 3409           0.037373
 3410           0.039941
 3411           0.042760
 3412           0.045871
 3413           0.049300
 3414           0.053101
 3415           0.057328
 3416           0.062047
 3417           0.067335
 3418           0.073295
 3419           0.080012
 3420           0.087626
 3421           0.096287
 3422           0.106173
 3423           0.117516
 3424           0.130520
 3425           0.145491
 3426           0.162755
 3427           0.182681
 3428           0.205667
 3429           0.232171
 3430           0.262458
 3431           0.296804
 3432           0.335166
 3433           0.377015
 3434           0.421053
 3435           0.465015
 3436           0.505186
 3437           0.537141
 3438           0.556238
 3439           0.559064
 3440           0.544701
 3441           0.515343
 3442           0.475408
 3443           0.430171
 3444           0.384418
 3445           0.341696
 3446           0.304093
 3447           0.272720
 3448           0.247605
 3449           0.228305
 3450           0.214089
 3451           0.204132
 3452           0.197647
 3453           0.194034
 3454           0.192876
 3455           0.194034
 3456           0.197647
 3457           0.204132
 3458           0.214089
 3459           0.228305
 3460           0.247605
 3461           0.272720
 3462           0.304093
 3463           0.341696
 3464           0.384418
 3465           0.430171
 3466           0.475408
 3467           0.515343
 3468           0.544701
 3469           0.559064
 3470           0.556238
 3471           0.537141
 3472           0.505186
 3473           0.465015
 3474           0.421053
 3475           0.377015
 3476           0.335166
 3477           0.296804
 3478           0.262458
 3479           0.232171
 3480           0.205667
 3481           0.182681
 3482           0.162755
 3483           0.145491
 3484           0.130520
 3485           0.117516
 3486           0.106173
 3487           0.096287
 3488           0.087626
 3489           0.080012
 3490           0.073295
 3491           0.067335
 3492           0.062047
 3493           0.057328
 3494           0.053101
 3495           0.049300
 3496           0.045871
 3497           0.042760
 3498           0.039941
 3499           0.037373
 3500           0.035028
 3501           0.032879
 3502           0.030907
 3503           0.029090
 3504           0.027357
 3505           0.029052
 3506           0.030897
 3507           0.032902
 3508           0.035091
 3509           0.037486
 3510           0.040114
 3511           0.043007
 3512           0.046208
 3513           0.049749
 3514           0.053687
 3515           0.058084
 3516           0.063014
 3517           0.068563
 3518           0.074851
 3519           0.081979
 3520           0.090110
 3521           0.099425
 3522           0.110143
 3523           0.122552
 3524           0.136922
 3525           0.153656
 3526           0.173210
 3527           0.196123
 3528           0.223020
 3529           0.254666
 3530           0.291681
 3531           0.334806
 3532           0.384494
 3533           0.440656
 3534           0.502167
 3535           0.566337
 3536           0.627796
 3537           0.679178
 3538           0.711868
 3539           0.718924
 3540           0.698193
 3541           0.653770
 3542           0.593966
 3543           0.528071
 3544           0.463697
 3545           0.405770
 3546           0.356656
 3547           0.317181
 3548           0.286754
 3549           0.264275
 3550           0.248388
 3551           0.237726
 3552           0.231063
 3553           0.227476
 3554           0.226349
 3555           0.227476
 3556           0.231063
 3557           0.237726
 3558           0.248388
 3559           0.264275
 3560           0.286754
 3561           0.317181
 3562           0.356656
 3563           0.405770
 3564           0.463697
 3565           0.528071
 3566           0.593966
 3567           0.653770
 3568           0.698193
 3569           0.718924
 3570           0.711868
 3571           0.679178
 3572           0.627796
 3573           0.566337
 3574           0.502167
 3575           0.440656
 3576           0.384494
 3577           0.334806
 3578           0.291681
 3579           0.254666
 3580           0.223020
 3581           0.196123
 3582           0.173210
 3583           0.153656
 3584           0.136922
 3585           0.122552
 3586           0.110143
 3587           0.099425
 3588           0.090110
 3589           0.081979
 3590           0.074851
 3591           0.068563
 3592           0.063014
 3593           0.058084
 3594           0.053687
 3595           0.049749
 3596           0.046208
 3597           0.043007
 3598           0.040114
 3599           0.037486
 3600           0.035091
 3601           0.032902
 3602           0.030897
 3603           0.029052
 3604           0.027235
 3605           0.028945
 3606           0.030810
 3607           0.032840
 3608           0.035061
 3609           0.037495
 3610           0.040173
 3611           0.043127
 3612           0.046404
 3613           0.050039
 3614           0.054095
 3615           0.058640
 3616           0.063754
 3617           0.069537
 3618           0.076121
 3619           0.083624
 3620           0.092232
 3621           0.102159
 3622           0.113665
 3623           0.127098
 3624           0.142804
 3625           0.161292
 3626           0.183170
 3627           0.209182
 3628           0.240242
 3629           0.277527
 3630           0.322183
 3631           0.375688
 3632           0.439421
 3633           0.514348
 3634           0.600282
 3635           0.694789
 3636           0.790718
 3637           0.875940
 3638           0.933816
 3639           0.949215
 3640           0.916841
 3641           0.845294
 3642           0.751599
 3643           0.652828
 3644           0.560932
 3645           0.482131
 3646           0.418337
 3647           0.369326
 3648           0.333255
 3649           0.307906
 3650           0.290977
 3651           0.280324
 3652           0.274113
 3653           0.270976
 3654           0.270029
 3655           0.270976
 3656           0.274113
 3657           0.280324
 3658           0.290977
 3659           0.307906
 3660           0.333255
 3661           0.369326
 3662           0.418337
 3663           0.482131
 3664           0.560932
 3665           0.652828
 3666           0.751599
 3667           0.845294
 3668           0.916841
 3669           0.949215
 3670           0.933816
 3671           0.875940
 3672           0.790718
 3673           0.694789
 3674           0.600282
 3675           0.514348
 3676           0.439421
 3677           0.375688
 3678           0.322183
 3679           0.277527
 3680           0.240242
 3681           0.209182
 3682           0.183170
 3683           0.161292
 3684           0.142804
 3685           0.127098
 3686           0.113665
 3687           0.102159
 3688           0.092232
 3689           0.083624
 3690           0.076121
 3691           0.069537
 3692           0.063754
 3693           0.058640
 3694           0.054095
 3695           0.050039
 3696           0.046404
 3697           0.043127
 3698           0.040173
 3699           0.037495
 3700           0.035061
 3701           0.032840
 3702           0.030810
 3703           0.028945
 3704           0.027047
 3705           0.028766
 3706           0.030644
 3707           0.032691
 3708           0.034933
 3709           0.037395
 3710           0.040109
 3711           0.043108
 3712           0.046444
 3713           0.050154
 3714           0.054304
 3715           0.058969
 3716           0.064237
 3717           0.070216
 3718           0.077051
 3719           0.084877
 3720           0.093903
 3721           0.104373
 3722           0.116588
 3723           0.130956
 3724           0.147900
 3725           0.168045
 3726           0.192159
 3727           0.221222
 3728           0.256486
 3729           0.299639
 3730           0.352530
 3731           0.417705
 3732           0.498059
 3733           0.596625
 3734           0.715740
 3735           0.855324
 3736           1.007893
 3737           1.154748
 3738           1.262777
 3739           1.295879
 3740           1.239723
 3741           1.115274
 3742           0.961006
 3743           0.809126
 3744           0.677003
 3745           0.570428
 3746           0.488821
 3747           0.429398
 3748           0.388069
 3749           0.360890
 3750           0.344214
 3751           0.334846
 3752           0.330146
 3753           0.328153
 3754           0.327628
 3755           0.328153
 3756           0.330146
 3757           0.334846
 3758           0.344214
 3759           0.360890
 3760           0.388069
 3761           0.429398
 3762           0.488821
 3763           0.570428
 3764           0.677003
 3765           0.809126
 3766           0.961006
 3767           1.115274
 3768           1.239723
 3769           1.295879
 3770           1.262777
 3771           1.154748
 3772           1.007893
 3773           0.855324
 3774           0.715740
 3775           0.596625
 3776           0.498059
 3777           0.417705
 3778           0.352530
 3779           0.299639
 3780           0.256486
 3781           0.221222
 3782           0.192159
 3783           0.168045
 3784           0.147900
 3785           0.130956
 3786           0.116588
 3787           0.104373
 3788           0.093903
 3789           0.084877
 3790           0.077051
 3791           0.070216
 3792           0.064237
 3793           0.058969
 3794           0.054304
 3795           0.050154
 3796           0.046444
 3797           0.043108
 3798           0.040109
 3799           0.037395
 3800           0.034933
 3801           0.032691
 3802           0.030644
 3803           0.028766
 3804           0.026793
 3805           0.028514
 3806           0.030396
 3807           0.032450
 3808           0.034703
 3809           0.037181
 3810           0.039916
 3811           0.042945
 3812           0.046320
 3813           0.050081
 3814           0.054299
 3815           0.059053
 3816           0.064437
 3817           0.070567
 3818           0.077600
 3819           0.085685
 3820           0.095052
 3821           0.105971
 3822           0.118784
 3823           0.133952
 3824           0.151973
 3825           0.173586
 3826           0.199720
 3827           0.231600
 3828           0.270844
 3829           0.319711
 3830           0.380907
 3831           0.458361
 3832           0.557155
 3833           0.683778
 3834           0.845827
 3835           1.050505
 3836           1.296635
 3837           1.561918
 3838           1.781017
 3839           1.857014
 3840           1.744177
 3841           1.504999
 3842           1.236962
 3843           0.999162
 3844           0.810109
 3845           0.668595
 3846           0.567044
 3847           0.497610
 3848           0.452652
 3849           0.425802
 3850           0.411654
 3851           0.405667
 3852           0.404164
 3853           0.404405
 3854           0.404678
 3855           0.404405
 3856           0.404164
 3857           0.405667
 3858           0.411654
 3859           0.425802
 3860           0.452652
 3861           0.497610
 3862           0.567044
 3863           0.668595
 3864           0.810109
 3865           0.999162
 3866           1.236962
 3867           1.504999
 3868           1.744177
 3869           1.857014
 3870           1.781017
 3871           1.561918
 3872           1.296635
 3873           1.050505
 3874           0.845827
 3875           0.683778
 3876           0.557155
 3877           0.458361
 3878           0.380907
 3879           0.319711
 3880           0.270844
 3881           0.231600
 3882           0.199720
 3883           0.173586
 3884           0.151973
 3885           0.133952
 3886           0.118784
 3887           0.105971
 3888           0.095052
 3889           0.085685
 3890           0.077600
 3891           0.070567
 3892           0.064437
 3893           0.059053
 3894           0.054299
 3895           0.050081
 3896           0.046320
 3897           0.042945
 3898           0.039916
 3899           0.037181
 3900           0.034703
 3901           0.032450
 3902           0.030396
 3903           0.028514
 3904           0.026471
 3905           0.028186
 3906           0.030064
 3907           0.032115
 3908           0.034367
 3909           0.036848
 3910           0.039589
 3911           0.042630
 3912           0.046023
 3913           0.049811
 3914           0.054068
 3915           0.058875
 3916           0.064332
 3917           0.070563
 3918           0.077732
 3919           0.086000
 3920           0.095613
 3921           0.106865
 3922           0.120129
 3923           0.135916
 3924           0.154787
 3925           0.177580
 3926           0.205374
 3927           0.239618
 3928           0.282282
 3929           0.336200
 3930           0.404977
 3931           0.494102
 3932           0.611351
 3933           0.768040
 3934           0.980574
 3935           1.272232
 3936           1.667286
 3937           2.168933
 3938           2.671783
 3939           2.877545
 3940           2.595346
 3941           2.072939
 3942           1.585717
 3943           1.215074
 3944           0.952045
 3945           0.771261
 3946           0.650534
 3947           0.573877
 3948           0.528853
 3949           0.506115
 3950           0.498159
 3951           0.498815
 3952           0.503099
 3953           0.507237
 3954           0.508885
 3955           0.507237
 3956           0.503099
 3957           0.498815
 3958           0.498159
 3959           0.506115
 3960           0.528853
 3961           0.573877
 3962           0.650534
 3963           0.771261
 3964           0.952045
 3965           1.215074
 3966           1.585717
 3967           2.072939
 3968           2.595346
 3969           2.877545
 3970           2.671783
 3971           2.168933
 3972           1.667286
 3973           1.272232
 3974           0.980574
 3975           0.768040
 3976           0.611351
 3977           0.494102
 3978           0.404977
 3979           0.336200
 3980           0.282282
 3981           0.239618
 3982           0.205374
 3983           0.177580
 3984           0.154787
 3985           0.135916
 3986           0.120129
 3987           0.106865
 3988           0.095613
 3989           0.086000
 3990           0.077732
 3991           0.070563
 3992           0.064332
 3993           0.058875
 3994           0.054068
 3995           0.049811
 3996           0.046023
 3997           0.042630
 3998           0.039589
 3999           0.036848
 4000           0.034367
 4001           0.032115
 4002           0.030064
 4003           0.028186
 4004           0.026080
 4005           0.027783
 4006           0.029647
 4007           0.031685
 4008           0.033925
 4009           0.036395
 4010           0.039126
 4011           0.042160
 4012           0.045549
 4013           0.049337
 4014           0.053601
 4015           0.058423
 4016           0.063908
 4017           0.070182
 4018           0.077417
 4019           0.085782
 4020           0.095533
 4021           0.106983
 4022           0.120526
 4023           0.136709
 4024           0.156141
 4025           0.179738
 4026           0.208691
 4027           0.244629
 4028           0.289809
 4029           0.347542
 4030           0.422224
 4031           0.520773
 4032           0.653629
 4033           0.837381
 4034           1.099585
 4035           1.489046
 4036           2.090238
 4037           3.044972
 4038           4.408898
 4039           5.200598
 4040           4.179732
 4041           2.859236
 4042           1.975561
 4043           1.428401
 4044           1.085755
 4045           0.869704
 4046           0.735899
 4047           0.658336
 4048           0.619423
 4049           0.606737
 4050           0.610528
 4051           0.622730
 4052           0.636705
 4053           0.647264
 4054           0.651160
 4055           0.647264
 4056           0.636705
 4057           0.622730
 4058           0.610528
 4059           0.606737
 4060           0.619423
 4061           0.658336
 4062           0.735899
 4063           0.869704
 4064           1.085755
 4065           1.428401
 4066           1.975561
 4067           2.859236
 4068           4.179732
 4069           5.200598
 4070           4.408898
 4071           3.044972
 4072           2.090238
 4073           1.489046
 4074           1.099585
 4075           0.837381
 4076           0.653629
 4077           0.520773
 4078           0.422224
 4079           0.347542
 4080           0.289809
 4081           0.244629
 4082           0.208691
 4083           0.179738
 4084           0.156141
 4085           0.136709
 4086           0.120526
 4087           0.106983
 4088           0.095533
 4089           0.085782
 4090           0.077417
 4091           0.070182
 4092           0.063908
 4093           0.058423
 4094           0.053601
 4095           0.049337
 4096           0.045549
 4097           0.042160
 4098           0.039126
 4099           0.036395
 4100           0.033925
 4101           0.031685
 4102           0.029647
 4103           0.027783
 4104           0.025621
 4105           0.027303
 4106           0.029144
 4107           0.031160
 4108           0.033376
 4109           0.035820
 4110           0.038526
 4111           0.041532
 4112           0.044894
 4113           0.048656
 4114           0.052893
 4115           0.057692
 4116           0.063155
 4117           0.069413
 4118           0.076640
 4119           0.085008
 4120           0.094781
 4121           0.106278
 4122           0.119908
 4123           0.136235
 4124           0.155895
 4125           0.179848
 4126           0.209351
 4127           0.246135
 4128           0.292630
 4129           0.352436
 4130           0.430439
 4131           0.534483
 4132           0.676802
 4133           0.877799
 4134           1.173955
 4135           1.638138
 4136           2.430874
 4137           4.002895
 4138           8.024409
 4139          15.003308
 4140           7.064365
 4141           3.668856
 4142           2.287180
 4143           1.585079
 4144           1.187431
 4145           0.953935
 4146           0.820127
 4147           0.752144
 4148           0.728591
 4149           0.734584
 4150           0.758216
 4151           0.789219
 4152           0.818755
 4153           0.839579
 4154           0.847067
 4155           0.839579
 4156           0.818755
 4157           0.789219
 4158           0.758216
 4159           0.734584
 4160           0.728591
 4161           0.752144
 4162           0.820127
 4163           0.953935
 4164           1.187431
 4165           1.585079
 4166           2.287180
 4167           3.668856
 4168           7.064365
 4169          15.003308
 4170           8.024409
 4171           4.002895
 4172           2.430874
 4173           1.638138
 4174           1.173955
 4175           0.877799
 4176           0.676802
 4177           0.534483
 4178           0.430439
 4179           0.352436
 4180           0.292630
 4181           0.246135
 4182           0.209351
 4183           0.179848
 4184           0.155895
 4185           0.136235
 4186           0.119908
 4187           0.106278
 4188           0.094781
 4189           0.085008
 4190           0.076640
 4191           0.069413
 4192           0.063155
 4193           0.057692
 4194           0.052893
 4195           0.048656
 4196           0.044894
 4197           0.041532
 4198           0.038526
 4199           0.035820
 4200           0.033376
 4201           0.031160
 4202           0.029144
 4203           0.027303
 4204           0.025094
 4205           0.026746
 4206           0.028557
 4207           0.030539
 4208           0.032718
 4209           0.035123
 4210           0.037786
 4211           0.040747
 4212           0.044059
 4213           0.047766
 4214           0.051944
 4215           0.056677
 4216           0.062070
 4217           0.068250
 4218           0.075392
 4219           0.083668
 4220           0.093340
 4221           0.104729
 4222           0.118242
 4223           0.134446
 4224           0.153978
 4225           0.177802
 4226           0.207182
 4227           0.243864
 4228           0.290296
 4229           0.350117
 4230           0.428276
 4231           0.532734
 4232           0.675957
 4233           0.878838
 4234           1.179029
 4235           1.652737
 4236           2.472359
 4237           4.151765
 4238           9.050182
 4239          24.764236
 4240           7.775009
 4241           3.800718
 4242           2.343925
 4243           1.630016
 4244           1.237624
 4245           1.017467
 4246           0.902742
 4247           0.858609
 4248           0.862915
 4249           0.899544
 4250           0.954645
 4251           1.015248
 4252           1.069345
 4253           1.106447
 4254           1.119650
 4255           1.106447
 4256           1.069345
 4257           1.015248
 4258           0.954645
 4259           0.899544
 4260           0.862915
 4261           0.858609
 4262           0.902742
 4263           1.017467
 4264           1.237624
 4265           1.630016
 4266           2.343925
 4267           3.800718
 4268           7.775009
 4269          24.764236
 4270           9.050182
 4271           4.151765
 4272           2.472359
 4273           1.652737
 4274           1.179029
 4275           0.878838
 4276           0.675957
 4277           0.532734
 4278           0.428276
 4279           0.350117
 4280           0.290296
 4281           0.243864
 4282           0.207182
 4283           0.177802
 4284           0.153978
 4285           0.134446
 4286           0.118242
 4287           0.104729
 4288           0.093340
 4289           0.083668
 4290           0.075392
 4291           0.068250
 4292           0.062070
 4293           0.056677
 4294           0.051944
 4295           0.047766
 4296           0.044059
 4297           0.040747
 4298           0.037786
 4299           0.035123
 4300           0.032718
 4301           0.030539
 4302           0.028557
 4303           0.026746
 4304           0.024502
 4305           0.026117
 4306           0.027888
 4307           0.029826
 4308           0.031958
 4309           0.034310
 4310           0.036914
 4311           0.039810
 4312           0.043049
 4313           0.046674
 4314           0.050760
 4315           0.055389
 4316           0.060662
 4317           0.066705
 4318           0.073687
 4319           0.081776
 4320           0.091229
 4321           0.102355
 4322           0.115552
 4323           0.131367
 4324           0.150417
 4325           0.173630
 4326           0.202218
 4327           0.237846
 4328           0.282830
 4329           0.340586
 4330           0.415685
 4331           0.515364
 4332           0.650663
 4333           0.839393
 4334           1.111824
 4335           1.523606
 4336           2.178650
 4337           3.282279
 4338           5.062836
 4339           6.270652
 4340           4.763190
 4341           3.101999
 4342           2.117880
 4343           1.556435
 4344           1.235115
 4345           1.061435
 4346           0.986713
 4347           0.982883
 4348           1.030495
 4349           1.113630
 4350           1.216632
 4351           1.322900
 4352           1.415496
 4353           1.478454
 4354           1.500800
 4355           1.478454
 4356           1.415496
 4357           1.322900
 4358           1.216632
 4359           1.113630
 4360           1.030495
 4361           0.982883
 4362           0.986713
 4363           1.061435
 4364           1.235115
 4365           1.556435
 4366           2.117880
 4367           3.101999
 4368           4.763190
 4369           6.270652
 4370           5.062836
 4371           3.282279
 4372           2.178650
 4373           1.523606
 4374           1.111824
 4375           0.839393
 4376           0.650663
 4377           0.515364
 4378           0.415685
 4379           0.340586
 4380           0.282830
 4381           0.237846
 4382           0.202218
 4383           0.173630
 4384           0.150417
 4385           0.131367
 4386           0.115552
 4387           0.102355
 4388           0.091229
 4389           0.081776
 4390           0.073687
 4391           0.066705
 4392           0.060662
 4393           0.055389
 4394           0.050760
 4395           0.046674
 4396           0.043049
 4397           0.039810
 4398           0.036914
 4399           0.034310
 4400           0.031958
 4401           0.029826
 4402           0.027888
 4403           0.026117
 4404           0.023846
 4405           0.025417
 4406           0.027140
 4407           0.029024
 4408           0.031097
 4409           0.033383
 4410           0.035913
 4411           0.038726
 4412           0.041870
 4413           0.045388
 4414           0.049351
 4415           0.053837
 4416           0.058945
 4417           0.064793
 4418           0.071545
 4419           0.079359
 4420           0.088479
 4421           0.099199
 4422           0.111892
 4423           0.127074
 4424           0.145316
 4425           0.167477
 4426           0.194667
 4427           0.228388
 4428           0.270704
 4429           0.324598
 4430           0.393929
 4431           0.484619
 4432           0.605219
 4433           0.768487
 4434           0.993646
 4435           1.309768
 4436           1.752612
 4437           2.344357
 4438           2.981081
 4439           3.266263
 4440           2.916768
 4441           2.299662
 4442           1.773395
 4443           1.414400
 4444           1.199202
 4445           1.095476
 4446           1.079042
 4447           1.132653
 4448           1.242081
 4449           1.392861
 4450           1.567514
 4451           1.744214
 4452           1.897888
 4453           2.002760
 4454           2.040105
 4455           2.002760
 4456           1.897888
 4457           1.744214
 4458           1.567514
 4459           1.392861
 4460           1.242081
 4461           1.132653
 4462           1.079042
 4463           1.095476
 4464           1.199202
 4465           1.414400
 4466           1.773395
 4467           2.299662
 4468           2.916768
 4469           3.266263
 4470           2.981081
 4471           2.344357
 4472           1.752612
 4473           1.309768
 4474           0.993646
 4475           0.768487
 4476           0.605219
 4477           0.484619
 4478           0.393929
 4479           0.324598
 4480           0.270704
 4481           0.228388
 4482           0.194667
 4483           0.167477
 4484           0.145316
 4485           0.127074
 4486           0.111892
 4487           0.099199
 4488           0.088479
 4489           0.079359
 4490           0.071545
 4491           0.064793
 4492           0.058945
 4493           0.053837
 4494           0.049351
 4495           0.045388
 4496           0.041870
 4497           0.038726
 4498           0.035913
 4499           0.033383
 4500           0.031097
 4501           0.029024
 4502           0.027140
 4503           0.025417
 4504           0.023129
 4505           0.024650
 4506           0.026315
 4507           0.028137
 4508           0.030140
 4509           0.032348
 4510           0.034790
 4511           0.037502
 4512           0.040532
 4513           0.043918
 4514           0.047729
 4515           0.052039
 4516           0.056939
 4517           0.062541
 4518           0.068998
 4519           0.076457
 4520           0.085145
 4521           0.095332
 4522           0.107360
 4523           0.121697
 4524           0.138855
 4525           0.159598
 4526           0.184895
 4527           0.216039
 4528           0.254759
 4529           0.303493
 4530           0.365231
 4531           0.444370
 4532           0.546757
 4533           0.680165
 4534           0.854324
 4535           1.079878
 4536           1.360100
 4537           1.675314
 4538           1.951152
 4539           2.064126
 4540           1.954104
 4541           1.708332
 4542           1.454386
 4543           1.263268
 4544           1.155859
 4545           1.132322
 4546           1.187674
 4547           1.315630
 4548           1.508708
 4549           1.755313
 4550           2.036382
 4551           2.322662
 4552           2.575405
 4553           2.750472
 4554           2.813375
 4555           2.750472
 4556           2.575405
 4557           2.322662
 4558           2.036382
 4559           1.755313
 4560           1.508708
 4561           1.315630
 4562           1.187674
 4563           1.132322
 4564           1.155859
 4565           1.263268
 4566           1.454386
 4567           1.708332
 4568           1.954104
 4569           2.064126
 4570           1.951152
 4571           1.675314
 4572           1.360100
 4573           1.079878
 4574           0.854324
 4575           0.680165
 4576           0.546757
 4577           0.444370
 4578           0.365231
 4579           0.303493
 4580           0.254759
 4581           0.216039
 4582           0.184895
 4583           0.159598
 4584           0.138855
 4585           0.121697
 4586           0.107360
 4587           0.095332
 4588           0.085145
 4589           0.076457
 4590           0.068998
 4591           0.062541
 4592           0.056939
 4593           0.052039
 4594           0.047729
 4595           0.043918
 4596           0.040532
 4597           0.037502
 4598           0.034790
 4599           0.032348
 4600           0.030140
 4601           0.028137
 4602           0.026315
 4603           0.024650
 4604           0.022355
 4605           0.023818
 4606           0.025420
 4607           0.027171
 4608           0.029094
 4609           0.031212
 4610           0.033553
 4611           0.036149
 4612           0.039046
 4613           0.042280
 4614           0.045913
 4615           0.050015
 4616           0.054671
 4617           0.059982
 4618           0.066090
 4619           0.073127
 4620           0.081298
 4621           0.090848
 4622           0.102079
 4623           0.115405
 4624           0.131268
 4625           0.150321
 4626           0.173378
 4627           0.201496
 4628           0.236047
 4629           0.278905
 4630           0.332211
 4631           0.398951
 4632           0.482689
 4633           0.587462
 4634           0.717031
 4635           0.873144
 4636           1.049850
 4637           1.228237
 4638           1.371004
 4639           1.435335
 4640           1.405824
 4641           1.312929
 4642           1.209099
 4643           1.137713
 4644           1.124267
 4645           1.181902
 4646           1.318186
 4647           1.536911
 4648           1.839150
 4649           2.218984
 4650           2.657782
 4651           3.117192
 4652           3.536021
 4653           3.834433
 4654           3.943446
 4655           3.834433
 4656           3.536021
 4657           3.117192
 4658           2.657782
 4659           2.218984
 4660           1.839150
 4661           1.536911
 4662           1.318186
 4663           1.181902
 4664           1.124267
 4665           1.137713
 4666           1.209099
 4667           1.312929
 4668           1.405824
 4669           1.435335
 4670           1.371004
 4671           1.228237
 4672           1.049850
 4673           0.873144
 4674           0.717031
 4675           0.587462
 4676           0.482689
 4677           0.398951
 4678           0.332211
 4679           0.278905
 4680           0.236047
 4681           0.201496
 4682           0.173378
 4683           0.150321
 4684           0.131268
 4685           0.115405
 4686           0.102079
 4687           0.090848
 4688           0.081298
 4689           0.073127
 4690           0.066090
 4691           0.059982
 4692           0.054671
 4693           0.050015
 4694           0.045913
 4695           0.042280
 4696           0.039046
 4697           0.036149
 4698           0.033553
 4699           0.031212
 4700           0.029094
 4701           0.027171
 4702           0.025420
 4703           0.023818
 4704           0.021527
 4705           0.022928
 4706           0.024459
 4707           0.026131
 4708           0.027966
 4709           0.029984
 4710           0.032211
 4711           0.034679
 4712           0.037427
 4713           0.040490
 4714           0.043924
 4715           0.047793
 4716           0.052172
 4717           0.057156
 4718           0.062870
 4719           0.069431
 4720           0.077022
 4721           0.085856
 4722           0.096195
 4723           0.108394
 4724           0.122820
 4725           0.140015
 4726           0.160635
 4727           0.185505
 4728           0.215664
 4729           0.252473
 4730           0.297356
 4731           0.352188
 4732           0.418912
 4733           0.499269
 4734           0.594054
 4735           0.701949
 4736           0.816620
 4737           0.926092
 4738           1.012976
 4739           1.061719
 4740           1.069585
 4741           1.051815
 4742           1.035247
 4743           1.047665
 4744           1.112250
 4745           1.247654
 4746           1.470622
 4747           1.795434
 4748           2.235345
 4749           2.796424
 4750           3.467889
 4751           4.205703
 4752           4.915484
 4753           5.446123
 4754           5.645678
 4755           5.446123
 4756           4.915484
 4757           4.205703
 4758           3.467889
 4759           2.796424
 4760           2.235345
 4761           1.795434
 4762           1.470622
 4763           1.247654
 4764           1.112250
 4765           1.047665
 4766           1.035247
 4767           1.051815
 4768           1.069585
 4769           1.061719
 4770           1.012976
 4771           0.926092
 4772           0.816620
 4773           0.701949
 4774           0.594054
 4775           0.499269
 4776           0.418912
 4777           0.352188
 4778           0.297356
 4779           0.252473
 4780           0.215664
 4781           0.185505
 4782           0.160635
 4783           0.140015
 4784           0.122820
 4785           0.108394
 4786           0.096195
 4787           0.085856
 4788           0.077022
 4789           0.069431
 4790           0.062870
 4791           0.057156
 4792           0.052172
 4793           0.047793
 4794           0.043924
 4795           0.040490
 4796           0.037427
 4797           0.034679
 4798           0.032211
 4799           0.029984
 4800           0.027966
 4801           0.026131
 4802           0.024459
 4803           0.022928
 4804           0.020649
 4805           0.021981
 4806           0.023436
 4807           0.025023
 4808           0.026761
 4809           0.028670
 4810           0.030774
 4811           0.033101
 4812           0.035687
 4813           0.038563
 4814           0.041780
 4815           0.045394
 4816           0.049474
 4817           0.054101
 4818           0.059387
 4819           0.065433
 4820           0.072397
 4821           0.080461
 4822           0.089847
 4823           0.100850
 4824           0.113766
 4825           0.129029
 4826           0.147149
 4827           0.168746
 4828           0.194571
 4829           0.225568
 4830           0.262618
 4831           0.306821
 4832           0.359121
 4833           0.420069
 4834           0.489333
 4835           0.565140
 4836           0.642991
 4837           0.716419
 4838           0.777896
 4839           0.822506
 4840           0.852210
 4841           0.877578
 4842           0.916251
 4843           0.989132
 4844           1.118087
 4845           1.325795
 4846           1.637697
 4847           2.080731
 4848           2.686599
 4849           3.485675
 4850           4.494822
 4851           5.686909
 4852           6.936064
 4853           7.950620
 4854           8.353035
 4855           7.950620
 4856           6.936064
 4857           5.686909
 4858           4.494822
 4859           3.485675
 4860           2.686599
 4861           2.080731
 4862           1.637697
 4863           1.325795
 4864           1.118087
 4865           0.989132
 4866           0.916251
 4867           0.877578
 4868           0.852210
 4869           0.822506
 4870           0.777896
 4871           0.716419
 4872           0.642991
 4873           0.565140
 4874           0.489333
 4875           0.420069
 4876           0.359121
 4877           0.306821
 4878           0.262618
 4879           0.225568
 4880           0.194571
 4881           0.168746
 4882           0.147149
 4883           0.129029
 4884           0.113766
 4885           0.100850
 4886           0.089847
 4887           0.080461
 4888           0.072397
 4889           0.065433
 4890           0.059387
 4891           0.054101
 4892           0.049474
 4893           0.045394
 4894           0.041780
 4895           0.038563
 4896           0.035687
 4897           0.033101
 4898           0.030774
 4899           0.028670
 4900           0.026761
 4901           0.025023
 4902           0.023436
 4903           0.021981
 4904           0.019729
 4905           0.020988
 4906           0.022361
 4907           0.023857
 4908           0.025493
 4909           0.027287
 4910           0.029259
 4911           0.031436
 4912           0.033850
 4913           0.036527
 4914           0.039514
 4915           0.042859
 4916           0.046621
 4917           0.050874
 4918           0.055711
 4919           0.061219
 4920           0.067532
 4921           0.074801
 4922           0.083210
 4923           0.092998
 4924           0.104398
 4925           0.117747
 4926           0.133428
 4927           0.151895
 4928           0.173668
 4929           0.199381
 4930           0.229546
 4931           0.264777
 4932           0.305483
 4933           0.351722
 4934           0.402942
 4935           0.457825
 4936           0.513748
 4937           0.567770
 4938           0.617331
 4939           0.662203
 4940           0.706508
 4941           0.759413
 4942           0.835329
 4943           0.952722
 4944           1.133856
 4945           1.405675
 4946           1.802845
 4947           2.367674
 4948           3.158052
 4949           4.247309
 4950           5.720608
 4951           7.643185
 4952           9.944893
 4953          12.113895
 4954          13.071432
 4955          12.113895
 4956           9.944893
 4957           7.643185
 4958           5.720608
 4959           4.247309
 4960           3.158052
 4961           2.367674
 4962           1.802845
 4963           1.405675
 4964           1.133856
 4965           0.952722
 4966           0.835329
 4967           0.759413
 4968           0.706508
 4969           0.662203
 4970           0.617331
 4971           0.567770
 4972           0.513748
 4973           0.457825
 4974           0.402942
 4975           0.351722
 4976           0.305483
 4977           0.264777
 4978           0.229546
 4979           0.199381
 4980           0.173668
 4981           0.151895
 4982           0.133428
 4983           0.117747
 4984           0.104398
 4985           0.092998
 4986           0.083210
 4987           0.074801
 4988           0.067532
 4989           0.061219
 4990           0.055711
 4991           0.050874
 4992           0.046621
 4993           0.042859
 4994           0.039514
 4995           0.036527
 4996           0.033850
 4997           0.031436
 4998           0.029259
 4999           0.027287
 5000           0.025493
 5001           0.023857
 5002           0.022361
 5003           0.020988
 5004           0.018771
 5005           0.019953
 5006           0.021240
 5007           0.022640
 5008           0.024168
 5009           0.025841
 5010           0.027676
 5011           0.029696
 5012           0.031930
 5013           0.034400
 5014           0.037147
 5015           0.040213
 5016           0.043648
 5017           0.047514
 5018           0.051892
 5019           0.056852
 5020           0.062506
 5021           0.068977
 5022           0.076413
 5023           0.085005
 5024           0.094929
 5025           0.106442
 5026           0.119826
 5027           0.135401
 5028           0.153520
 5029           0.174601
 5030           0.198926
 5031           0.226835
 5032           0.258496
 5033           0.293841
 5034           0.332473
 5035           0.373697
 5036           0.416327
 5037           0.459559
 5038           0.503443
 5039           0.550036
 5040           0.604607
 5041           0.675989
 5042           0.777572
 5043           0.927295
 5044           1.148673
 5045           1.472892
 5046           1.943661
 5047           2.619733
 5048           3.590511
 5049           4.990836
 5050           7.032526
 5051          10.045525
 5052          14.442471
 5053          19.911935
 5054          22.992619
 5055          19.911935
 5056          14.442471
 5057          10.045525
 5058           7.032526
 5059           4.990836
 5060           3.590511
 5061           2.619733
 5062           1.943661
 5063           1.472892
 5064           1.148673
 5065           0.927295
 5066           0.777572
 5067           0.675989
 5068           0.604607
 5069           0.550036
 5070           0.503443
 5071           0.459559
 5072           0.416327
 5073           0.373697
 5074           0.332473
 5075           0.293841
 5076           0.258496
 5077           0.226835
 5078           0.198926
 5079           0.174601
 5080           0.153520
 5081           0.135401
 5082           0.119826
 5083           0.106442
 5084           0.094929
 5085           0.085005
 5086           0.076413
 5087           0.068977
 5088           0.062506
 5089           0.056852
 5090           0.051892
 5091           0.047514
 5092           0.043648
 5093           0.040213
 5094           0.037147
 5095           0.034400
 5096           0.031930
 5097           0.029696
 5098           0.027676
 5099           0.025841
 5100           0.024168
 5101           0.022640
 5102           0.021240
 5103           0.019953
 5104           0.017779
 5105           0.018881
 5106           0.020080
 5107           0.021380
 5108           0.022797
 5109           0.024344
 5110           0.026037
 5111           0.027895
 5112           0.029944
 5113           0.032203
 5114           0.034705
 5115           0.037486
 5116           0.040590
 5117           0.044066
 5118           0.047984
 5119           0.052398
 5120           0.057400
 5121           0.063089
 5122           0.069581
 5123           0.077025
 5124           0.085551
 5125           0.095352
 5126           0.106631
 5127           0.119610
 5128           0.134527
 5129           0.151657
 5130           0.171156
 5131           0.193229
 5132           0.217966
 5133           0.245342
 5134           0.275213
 5135           0.307436
 5136           0.341806
 5137           0.378788
 5138           0.419837
 5139           0.468170
 5140           0.529660
 5141           0.613067
 5142           0.731439
 5143           0.902687
 5144           1.151415
 5145           1.511973
 5146           2.035002
 5147           2.792747
 5148           3.902855
 5149           5.563020
 5150           8.143219
 5151          12.435316
 5152          20.421391
 5153          36.686834
 5154          54.236481
 5155          36.686834
 5156          20.421391
 5157          12.435316
 5158           8.143219
 5159           5.563020
 5160           3.902855
 5161           2.792747
 5162           2.035002
 5163           1.511973
 5164           1.151415
 5165           0.902687
 5166           0.731439
 5167           0.613067
 5168           0.529660
 5169           0.468170
 5170           0.419837
 5171           0.378788
 5172           0.341806
 5173           0.307436
 5174           0.275213
 5175           0.245342
 5176           0.217966
 5177           0.193229
 5178           0.171156
 5179           0.151657
 5180           0.134527
 5181           0.119610
 5182           0.106631
 5183           0.095352
 5184           0.085551
 5185           0.077025
 5186           0.069581
 5187           0.063089
 5188           0.057400
 5189           0.052398
 5190           0.047984
 5191           0.044066
 5192           0.040590
 5193           0.037486
 5194           0.034705
 5195           0.032203
 5196           0.029944
 5197           0.027895
 5198           0.026037
 5199           0.024344
 5200           0.022797
 5201           0.021380
 5202           0.020080
 5203           0.018881
 5204           0.016759
 5205           0.017780
 5206           0.018887
 5207           0.020086
 5208           0.021389
 5209           0.022807
 5210           0.024356
 5211           0.026049
 5212           0.027911
 5213           0.029955
 5214           0.032211
 5215           0.034708
 5216           0.037481
 5217           0.040572
 5218           0.044036
 5219           0.047917
 5220           0.052287
 5221           0.057224
 5222           0.062818
 5223           0.069184
 5224           0.076415
 5225           0.084654
 5226           0.094045
 5227           0.104744
 5228           0.116912
 5229           0.130738
 5230           0.146316
 5231           0.163798
 5232           0.183278
 5233           0.204831
 5234           0.228551
 5235           0.254706
 5236           0.283722
 5237           0.316791
 5238           0.356100
 5239           0.405416
 5240           0.470833
 5241           0.560947
 5242           0.688442
 5243           0.870925
 5244           1.133119
 5245           1.510398
 5246           2.056125
 5247           2.848351
 5248           4.017753
 5249           5.793389
 5250           8.633998
 5251          13.650768
 5252          24.473610
 5253          60.477250
 5254         814.217721
 5255          60.477250
 5256          24.473610
 5257          13.650768
 5258           8.633998
 5259           5.793389
 5260           4.017753
 5261           2.848351
 5262           2.056125
 5263           1.510398
 5264           1.133119
 5265           0.870925
 5266           0.688442
 5267           0.560947
 5268           0.470833
 5269           0.405416
 5270           0.356100
 5271           0.316791
 5272           0.283722
 5273           0.254706
 5274           0.228551
 5275           0.204831
 5276           0.183278
 5277           0.163798
 5278           0.146316
 5279           0.130738
 5280           0.116912
 5281           0.104744
 5282           0.094045
 5283           0.084654
 5284           0.076415
 5285           0.069184
 5286           0.062818
 5287           0.057224
 5288           0.052287
 5289           0.047917
 5290           0.044036
 5291           0.040572
 5292           0.037481
 5293           0.034708
 5294           0.032211
 5295           0.029955
 5296           0.027911
 5297           0.026049
 5298           0.024356
 5299           0.022807
 5300           0.021389
 5301           0.020086
 5302           0.018887
 5303           0.017780
 5304           0.015715
 5305           0.016653
 5306           0.017667
 5307           0.018763
 5308           0.019950
 5309           0.021239
 5310           0.022641
 5311           0.024169
 5312           0.025843
 5313           0.027674
 5314           0.029684
 5315           0.031900
 5316           0.034348
 5317           0.037061
 5318           0.040085
 5319           0.043452
 5320           0.047219
 5321           0.051445
 5322           0.056198
 5323           0.061567
 5324           0.067616
 5325           0.074451
 5326           0.082175
 5327           0.090898
 5328           0.100735
 5329           0.111825
 5330           0.124241
 5331           0.138122
 5332           0.153599
 5333           0.170848
 5334           0.190149
 5335           0.212040
 5336           0.237327
 5337           0.267578
 5338           0.305312
 5339           0.354457
 5340           0.421015
 5341           0.513181
 5342           0.642946
 5343           0.826958
 5344           1.088705
 5345           1.461901
 5346           1.997389
 5347           2.768898
 5348           3.898218
 5349           5.593163
 5350           8.251089
 5351          12.755050
 5352          21.498841
 5353          41.495701
 5354          68.694135
 5355          41.495701
 5356          21.498841
 5357          12.755050
 5358           8.251089
 5359           5.593163
 5360           3.898218
 5361           2.768898
 5362           1.997389
 5363           1.461901
 5364           1.088705
 5365           0.826958
 5366           0.642946
 5367           0.513181
 5368           0.421015
 5369           0.354457
 5370           0.305312
 5371           0.267578
 5372           0.237327
 5373           0.212040
 5374           0.190149
 5375           0.170848
 5376           0.153599
 5377           0.138122
 5378           0.124241
 5379           0.111825
 5380           0.100735
 5381           0.090898
 5382           0.082175
 5383           0.074451
 5384           0.067616
 5385           0.061567
 5386           0.056198
 5387           0.051445
 5388           0.047219
 5389           0.043452
 5390           0.040085
 5391           0.037061
 5392           0.034348
 5393           0.031900
 5394           0.029684
 5395           0.027674
 5396           0.025843
 5397           0.024169
 5398           0.022641
 5399           0.021239
 5400           0.019950
 5401           0.018763
 5402           0.017667
 5403           0.016653
 5404           0.014657
 5405           0.015511
 5406           0.016432
 5407           0.017423
 5408           0.018495
 5409           0.019654
 5410           0.020911
 5411           0.022276
 5412           0.023764
 5413           0.025384
 5414           0.027156
 5415           0.029097
 5416           0.031231
 5417           0.033582
 5418           0.036185
 5419           0.039065
 5420           0.042265
 5421           0.045831
 5422           0.049811
 5423           0.054273
 5424           0.059262
 5425           0.064856
 5426           0.071130
 5427           0.078166
 5428           0.086051
 5429           0.094897
 5430           0.104777
 5431           0.115835
 5432           0.128243
 5433           0.142257
 5434           0.158278
 5435           0.177005
 5436           0.199435
 5437           0.227286
 5438           0.263142
 5439           0.310830
 5440           0.376002
 5441           0.466184
 5442           0.592281
 5443           0.769284
 5444           1.018146
 5445           1.368527
 5446           1.864230
 5447           2.566275
 5448           3.570363
 5449           5.025349
 5450           7.175223
 5451          10.429448
 5452          15.400932
 5453          22.068520
 5454          26.124719
 5455          22.068520
 5456          15.400932
 5457          10.429448
 5458           7.175223
 5459           5.025349
 5460           3.570363
 5461           2.566275
 5462           1.864230
 5463           1.368527
 5464           1.018146
 5465           0.769284
 5466           0.592281
 5467           0.466184
 5468           0.376002
 5469           0.310830
 5470           0.263142
 5471           0.227286
 5472           0.199435
 5473           0.177005
 5474           0.158278
 5475           0.142257
 5476           0.128243
 5477           0.115835
 5478           0.104777
 5479           0.094897
 5480           0.086051
 5481           0.078166
 5482           0.071130
 5483           0.064856
 5484           0.059262
 5485           0.054273
 5486           0.049811
 5487           0.045831
 5488           0.042265
 5489           0.039065
 5490           0.036185
 5491           0.033582
 5492           0.031231
 5493           0.029097
 5494           0.027156
 5495           0.025384
 5496           0.023764
 5497           0.022276
 5498           0.020911
 5499           0.019654
 5500           0.018495
 5501           0.017423
 5502           0.016432
 5503           0.015511
 5504           0.013588
 5505           0.014357
 5506           0.015185
 5507           0.016074
 5508           0.017030
 5509           0.018062
 5510           0.019175
 5511           0.020379
 5512           0.021685
 5513           0.023101
 5514           0.024640
 5515           0.026317
 5516           0.028149
 5517           0.030154
 5518           0.032361
 5519           0.034784
 5520           0.037458
 5521           0.040414
 5522           0.043691
 5523           0.047336
 5524           0.051382
 5525           0.055888
 5526           0.060909
 5527           0.066510
 5528           0.072762
 5529           0.079764
 5530           0.087597
 5531           0.096414
 5532           0.106419
 5533           0.117918
 5534           0.131382
 5535           0.147578
 5536           0.167568
 5537           0.193067
 5538           0.226550
 5539           0.271562
 5540           0.333208
 5541           0.418131
 5542           0.535822
 5543           0.699097
 5544           0.925493
 5545           1.239129
 5546           1.674293
 5547           2.275537
 5548           3.107052
 5549           4.254404
 5550           5.825104
 5551           7.921786
 5552          10.518408
 5553          13.068672
 5554          14.231234
 5555          13.068672
 5556          10.518408
 5557           7.921786
 5558           5.825104
 5559           4.254404
 5560           3.107052
 5561           2.275537
 5562           1.674293
 5563           1.239129
 5564           0.925493
 5565           0.699097
 5566           0.535822
 5567           0.418131
 5568           0.333208
 5569           0.271562
 5570           0.226550
 5571           0.193067
 5572           0.167568
 5573           0.147578
 5574           0.131382
 5575           0.117918
 5576           0.106419
 5577           0.096414
 5578           0.087597
 5579           0.079764
 5580           0.072762
 5581           0.066510
 5582           0.060909
 5583           0.055888
 5584           0.051382
 5585           0.047336
 5586           0.043691
 5587           0.040414
 5588           0.037458
 5589           0.034784
 5590           0.032361
 5591           0.030154
 5592           0.028149
 5593           0.026317
 5594           0.024640
 5595           0.023101
 5596           0.021685
 5597           0.020379
 5598           0.019175
 5599           0.018062
 5600           0.017030
 5601           0.016074
 5602           0.015185
 5603           0.014357
 5604           0.012514
 5605           0.013199
 5606           0.013935
 5607           0.014721
 5608           0.015565
 5609           0.016470
 5610           0.017443
 5611           0.018490
 5612           0.019620
 5613           0.020837
 5614           0.022153
 5615           0.023578
 5616           0.025123
 5617           0.026804
 5618           0.028638
 5619           0.030638
 5620           0.032827
 5621           0.035228
 5622           0.037869
 5623           0.040785
 5624           0.044000
 5625           0.047556
 5626           0.051500
 5627           0.055882
 5628           0.060768
 5629           0.066249
 5630           0.072413
 5631           0.079421
 5632           0.087496
 5633           0.096966
 5634           0.108328
 5635           0.122350
 5636           0.140074
 5637           0.163111
 5638           0.193716
 5639           0.235035
 5640           0.291508
 5641           0.368763
 5642           0.474708
 5643           0.619746
 5644           0.817681
 5645           1.086726
 5646           1.451443
 5647           1.940756
 5648           2.591806
 5649           3.443826
 5650           4.526132
 5651           5.824797
 5652           7.216072
 5653           8.372619
 5654           8.838576
 5655           8.372619
 5656           7.216072
 5657           5.824797
 5658           4.526132
 5659           3.443826
 5660           2.591806
 5661           1.940756
 5662           1.451443
 5663           1.086726
 5664           0.817681
 5665           0.619746
 5666           0.474708
 5667           0.368763
 5668           0.291508
 5669           0.235035
 5670           0.193716
 5671           0.163111
 5672           0.140074
 5673           0.122350
 5674           0.108328
 5675           0.096966
 5676           0.087496
 5677           0.079421
 5678           0.072413
 5679           0.066249
 5680           0.060768
 5681           0.055882
 5682           0.051500
 5683           0.047556
 5684           0.044000
 5685           0.040785
 5686           0.037869
 5687           0.035228
 5688           0.032827
 5689           0.030638
 5690           0.028638
 5691           0.026804
 5692           0.025123
 5693           0.023578
 5694           0.022153
 5695           0.020837
 5696           0.019620
 5697           0.018490
 5698           0.017443
 5699           0.016470
 5700           0.015565
 5701           0.014721
 5702           0.013935
 5703           0.013199
 5704           0.011438
 5705           0.012042
 5706           0.012686
 5707           0.013373
 5708           0.014106
 5709           0.014888
 5710           0.015725
 5711           0.016620
 5712           0.017581
 5713           0.018608
 5714           0.019711
 5715           0.020897
 5716           0.022173
 5717           0.023549
 5718           0.025038
 5719           0.026648
 5720           0.028394
 5721           0.030294
 5722           0.032365
 5723           0.034635
 5724           0.037119
 5725           0.039852
 5726           0.042870
 5727           0.046219
 5728           0.049958
 5729           0.054173
 5730           0.058959
 5731           0.064478
 5732           0.070957
 5733           0.078726
 5734           0.088271
 5735           0.100321
 5736           0.115840
 5737           0.136273
 5738           0.163593
 5739           0.200493
 5740           0.250690
 5741           0.318776
 5742           0.411061
 5743           0.535592
 5744           0.702627
 5745           0.925011
 5746           1.218965
 5747           1.601216
 5748           2.090186
 5749           2.698613
 5750           3.422914
 5751           4.223878
 5752           5.003996
 5753           5.594867
 5754           5.818938
 5755           5.594867
 5756           5.003996
 5757           4.223878
 5758           3.422914
 5759           2.698613
 5760           2.090186
 5761           1.601216
 5762           1.218965
 5763           0.925011
 5764           0.702627
 5765           0.535592
 5766           0.411061
 5767           0.318776
 5768           0.250690
 5769           0.200493
 5770           0.163593
 5771           0.136273
 5772           0.115840
 5773           0.100321
 5774           0.088271
 5775           0.078726
 5776           0.070957
 5777           0.064478
 5778           0.058959
 5779           0.054173
 5780           0.049958
 5781           0.046219
 5782           0.042870
 5783           0.039852
 5784           0.037119
 5785           0.034635
 5786           0.032365
 5787           0.030294
 5788           0.028394
 5789           0.026648
 5790           0.025038
 5791           0.023549
 5792           0.022173
 5793           0.020897
 5794           0.019711
 5795           0.018608
 5796           0.017581
 5797           0.016620
 5798           0.015725
 5799           0.014888
 5800           0.014106
 5801           0.013373
 5802           0.012686
 5803           0.012042
 5804           0.010367
 5805           0.010890
 5806           0.011446
 5807           0.012035
 5808           0.012661
 5809           0.013325
 5810           0.014030
 5811           0.014779
 5812           0.015578
 5813           0.016425
 5814           0.017327
 5815           0.018288
 5816           0.019312
 5817           0.020406
 5818           0.021578
 5819           0.022831
 5820           0.024176
 5821           0.025625
 5822           0.027189
 5823           0.028887
 5824           0.030732
 5825           0.032752
 5826           0.034976
 5827           0.037446
 5828           0.040216
 5829           0.043368
 5830           0.046999
 5831           0.051266
 5832           0.056390
 5833           0.062685
 5834           0.070603
 5835           0.080804
 5836           0.094145
 5837           0.111873
 5838           0.135658
 5839           0.167724
 5840           0.211083
 5841           0.269341
 5842           0.347343
 5843           0.451035
 5844           0.587650
 5845           0.765698
 5846           0.995101
 5847           1.284353
 5848           1.640819
 5849           2.064863
 5850           2.543414
 5851           3.041381
 5852           3.496530
 5853           3.822564
 5854           3.942131
 5855           3.822564
 5856           3.496530
 5857           3.041381
 5858           2.543414
 5859           2.064863
 5860           1.640819
 5861           1.284353
 5862           0.995101
 5863           0.765698
 5864           0.587650
 5865           0.451035
 5866           0.347343
 5867           0.269341
 5868           0.211083
 5869           0.167724
 5870           0.135658
 5871           0.111873
 5872           0.094145
 5873           0.080804
 5874           0.070603
 5875           0.062685
 5876           0.056390
 5877           0.051266
 5878           0.046999
 5879           0.043368
 5880           0.040216
 5881           0.037446
 5882           0.034976
 5883           0.032752
 5884           0.030732
 5885           0.028887
 5886           0.027189
 5887           0.025625
 5888           0.024176
 5889           0.022831
 5890           0.021578
 5891           0.020406
 5892           0.019312
 5893           0.018288
 5894           0.017327
 5895           0.016425
 5896           0.015578
 5897           0.014779
 5898           0.014030
 5899           0.013325
 5900           0.012661
 5901           0.012035
 5902           0.011446
 5903           0.010890
 5904           0.009303
 5905           0.009748
 5906           0.010217
 5907           0.010712
 5908           0.011233
 5909           0.011783
 5910           0.012363
 5911           0.012973
 5912           0.013617
 5913           0.014294
 5914           0.015007
 5915           0.015758
 5916           0.016548
 5917           0.017381
 5918           0.018262
 5919           0.019191
 5920           0.020174
 5921           0.021218
 5922           0.022332
 5923           0.023527
 5924           0.024815
 5925           0.026216
 5926           0.027757
 5927           0.029475
 5928           0.031420
 5929           0.033669
 5930           0.036316
 5931           0.039509
 5932           0.043453
 5933           0.048437
 5934           0.054863
 5935           0.063309
 5936           0.074508
 5937           0.089505
 5938           0.109665
 5939           0.136772
 5940           0.173188
 5941           0.221649
 5942           0.285738
 5943           0.369669
 5944           0.478304
 5945           0.616947
 5946           0.791216
 5947           1.004659
 5948           1.258952
 5949           1.549944
 5950           1.864452
 5951           2.177107
 5952           2.450557
 5953           2.639512
 5954           2.707401
 5955           2.639512
 5956           2.450557
 5957           2.177107
 5958           1.864452
 5959           1.549944
 5960           1.258952
 5961           1.004659
 5962           0.791216
 5963           0.616947
 5964           0.478304
 5965           0.369669
 5966           0.285738
 5967           0.221649
 5968           0.173188
 5969           0.136772
 5970           0.109665
 5971           0.089505
 5972           0.074508
 5973           0.063309
 5974           0.054863
 5975           0.048437
 5976           0.043453
 5977           0.039509
 5978           0.036316
 5979           0.033669
 5980           0.031420
 5981           0.029475
 5982           0.027757
 5983           0.026216
 5984           0.024815
 5985           0.023527
 5986           0.022332
 5987           0.021218
 5988           0.020174
 5989           0.019191
 5990           0.018262
 5991           0.017381
 5992           0.016548
 5993           0.015758
 5994           0.015007
 5995           0.014294
 5996           0.013617
 5997           0.012973
 5998           0.012363
 5999           0.011783
 6000           0.011233
 6001           0.010712
 6002           0.010217
 6003           0.009748
 6004           0.008254
 6005           0.008623
 6006           0.009009
 6007           0.009413
 6008           0.009836
 6009           0.010277
 6010           0.010736
 6011           0.011215
 6012           0.011715
 6013           0.012232
 6014           0.012769
 6015           0.013325
 6016           0.013901
 6017           0.014496
 6018           0.015112
 6019           0.015748
 6020           0.016408
 6021           0.017093
 6022           0.017810
 6023           0.018565
 6024           0.019368
 6025           0.020236
 6026           0.021190
 6027           0.022265
 6028           0.023507
 6029           0.024984
 6030           0.026787
 6031           0.029049
 6032           0.031955
 6033           0.035760
 6034           0.040813
 6035           0.047603
 6036           0.056740
 6037           0.069074
 6038           0.085693
 6039           0.107986
 6040           0.137753
 6041           0.177008
 6042           0.228315
 6043           0.294554
 6044           0.378853
 6045           0.484334
 6046           0.613915
 6047           0.768507
 6048           0.947302
 6049           1.145318
 6050           1.352035
 6051           1.550515
 6052           1.718686
 6053           1.832062
 6054           1.872246
 6055           1.832062
 6056           1.718686
 6057           1.550515
 6058           1.352035
 6059           1.145318
 6060           0.947302
 6061           0.768507
 6062           0.613915
 6063           0.484334
 6064           0.378853
 6065           0.294554
 6066           0.228315
 6067           0.177008
 6068           0.137753
 6069           0.107986
 6070           0.085693
 6071           0.069074
 6072           0.056740
 6073           0.047603
 6074           0.040813
 6075           0.035760
 6076           0.031955
 6077           0.029049
 6078           0.026787
 6079           0.024984
 6080           0.023507
 6081           0.022265
 6082           0.021190
 6083           0.020236
 6084           0.019368
 6085           0.018565
 6086           0.017810
 6087           0.017093
 6088           0.016408
 6089           0.015748
 6090           0.015112
 6091           0.014496
 6092           0.013901
 6093           0.013325
 6094           0.012769
 6095           0.012232
 6096           0.011715
 6097           0.011215
 6098           0.010736
 6099           0.010277
 6100           0.009836
 6101           0.009413
 6102           0.009009
 6103           0.008623
 6104           0.007222
 6105           0.007517
 6106           0.007825
 6107           0.008142
 6108           0.008470
 6109           0.008807
 6110           0.009154
 6111           0.009510
 6112           0.009873
 6113           0.010242
 6114           0.010616
 6115           0.010993
 6116           0.011371
 6117           0.011749
 6118           0.012127
 6119           0.012500
 6120           0.012871
 6121           0.013239
 6122           0.013606
 6123           0.013977
 6124           0.014359
 6125           0.014765
 6126           0.015215
 6127           0.015738
 6128           0.016375
 6129           0.017188
 6130           0.018258
 6131           0.019705
 6132           0.021688
 6133           0.024427
 6134           0.028217
 6135           0.033460
 6136           0.040652
 6137           0.050465
 6138           0.063746
 6139           0.081549
 6140           0.105210
 6141           0.136167
 6142           0.176203
 6143           0.227221
 6144           0.291151
 6145           0.369717
 6146           0.464266
 6147           0.574480
 6148           0.698740
 6149           0.832657
 6150           0.968606
 6151           1.095670
 6152           1.200813
 6153           1.270446
 6154           1.294890
 6155           1.270446
 6156           1.200813
 6157           1.095670
 6158           0.968606
 6159           0.832657
 6160           0.698740
 6161           0.574480
 6162           0.464266
 6163           0.369717
 6164           0.291151
 6165           0.227221
 6166           0.176203
 6167           0.136167
 6168           0.105210
 6169           0.081549
 6170           0.063746
 6171           0.050465
 6172           0.040652
 6173           0.033460
 6174           0.028217
 6175           0.024427
 6176           0.021688
 6177           0.019705
 6178           0.018258
 6179           0.017188
 6180           0.016375
 6181           0.015738
 6182           0.015215
 6183           0.014765
 6184           0.014359
 6185           0.013977
 6186           0.013606
 6187           0.013239
 6188           0.012871
 6189           0.012500
 6190           0.012127
 6191           0.011749
 6192           0.011371
 6193           0.010993
 6194           0.010616
 6195           0.010242
 6196           0.009873
 6197           0.009510
 6198           0.009154
 6199           0.008807
 6200           0.008470
 6201           0.008142
 6202           0.007825
 6203           0.007517
 6204           0.006210
 6205           0.006436
 6206           0.006667
 6207           0.006902
 6208           0.007140
 6209           0.007380
 6210           0.007621
 6211           0.007861
 6212           0.008098
 6213           0.008329
 6214           0.008553
 6215           0.008765
 6216           0.008964
 6217           0.009146
 6218           0.009309
 6219           0.009448
 6220           0.009562
 6221           0.009650
 6222           0.009710
 6223           0.009747
 6224           0.009764
 6225           0.009772
 6226           0.009788
 6227           0.009835
 6228           0.009952
 6229           0.010191
 6230           0.010624
 6231           0.011354
 6232           0.012517
 6233           0.014296
 6234           0.016934
 6235           0.020756
 6236           0.026155
 6237           0.033652
 6238           0.043890
 6239           0.057649
 6240           0.075899
 6241           0.099640
 6242           0.130080
 6243           0.168436
 6244           0.215851
 6245           0.273209
 6246           0.341007
 6247           0.418479
 6248           0.503966
 6249           0.594050
 6250           0.683473
 6251           0.765310
 6252           0.831817
 6253           0.875281
 6254           0.890431
 6255           0.875281
 6256           0.831817
 6257           0.765310
 6258           0.683473
 6259           0.594050
 6260           0.503966
 6261           0.418479
 6262           0.341007
 6263           0.273209
 6264           0.215851
 6265           0.168436
 6266           0.130080
 6267           0.099640
 6268           0.075899
 6269           0.057649
 6270           0.043890
 6271           0.033652
 6272           0.026155
 6273           0.020756
 6274           0.016934
 6275           0.014296
 6276           0.012517
 6277           0.011354
 6278           0.010624
 6279           0.010191
 6280           0.009952
 6281           0.009835
 6282           0.009788
 6283           0.009772
 6284           0.009764
 6285           0.009747
 6286           0.009710
 6287           0.009650
 6288           0.009562
 6289           0.009448
 6290           0.009309
 6291           0.009146
 6292           0.008964
 6293           0.008765
 6294           0.008553
 6295           0.008329
 6296           0.008098
 6297           0.007861
 6298           0.007621
 6299           0.007380
 6300           0.007140
 6301           0.006902
 6302           0.006667
 6303           0.006436
 6304           0.005221
 6305           0.005381
 6306           0.005540
 6307           0.005697
 6308           0.005851
 6309           0.006000
 6310           0.006142
 6311           0.006274
 6312           0.006394
 6313           0.006499
 6314           0.006584
 6315           0.006647
 6316           0.006684
 6317           0.006689
 6318           0.006659
 6319           0.006590
 6320           0.006478
 6321           0.006319
 6322           0.006113
 6323           0.005859
 6324           0.005561
 6325           0.005226
 6326           0.004869
 6327           0.004508
 6328           0.004177
 6329           0.003920
 6330           0.003799
 6331           0.003900
 6332           0.004337
 6333           0.005260
 6334           0.006864
 6335           0.009405
 6336           0.013190
 6337           0.018615
 6338           0.026156
 6339           0.036384
 6340           0.049982
 6341           0.067635
 6342           0.090137
 6343           0.118249
 6344           0.152617
 6345           0.193646
 6346           0.241418
 6347           0.295107
 6348           0.353309
 6349           0.413533
 6350           0.472254
 6351           0.525114
 6352           0.567480
 6353           0.594889
 6354           0.604392
 6355           0.594889
 6356           0.567480
 6357           0.525114
 6358           0.472254
 6359           0.413533
 6360           0.353309
 6361           0.295107
 6362           0.241418
 6363           0.193646
 6364           0.152617
 6365           0.118249
 6366           0.090137
 6367           0.067635
 6368           0.049982
 6369           0.036384
 6370           0.026156
 6371           0.018615
 6372           0.013190
 6373           0.009405
 6374           0.006864
 6375           0.005260
 6376           0.004337
 6377           0.003900
 6378           0.003799
 6379           0.003920
 6380           0.004177
 6381           0.004508
 6382           0.004869
 6383           0.005226
 6384           0.005561
 6385           0.005859
 6386           0.006113
 6387           0.006319
 6388           0.006478
 6389           0.006590
 6390           0.006659
 6391           0.006689
 6392           0.006684
 6393           0.006647
 6394           0.006584
 6395           0.006499
 6396           0.006394
 6397           0.006274
 6398           0.006142
 6399           0.006000
 6400           0.005851
 6401           0.005697
 6402           0.005540
 6403           0.005381
 6404           0.004259
 6405           0.004356
 6406           0.004447
 6407           0.004531
 6408           0.004606
 6409           0.004670
 6410           0.004720
 6411           0.004752
 6412           0.004765
 6413           0.004753
 6414           0.004713
 6415           0.004641
 6416           0.004531
 6417           0.004378
 6418           0.004178
 6419           0.003924
 6420           0.003613
 6421           0.003240
 6422           0.002802
 6423           0.002296
 6424           0.001728
 6425           0.001099
 6426           0.000422
 6427          -0.000287
 6428          -0.001002
 6429          -0.001686
 6430          -0.002286
 6431          -0.002732
 6432          -0.002931
 6433          -0.002763
 6434          -0.002073
 6435          -0.000664
 6436           0.001702
 6437           0.005319
 6438           0.010535
 6439           0.017757
 6440           0.027463
 6441           0.040107
 6442           0.056207
 6443           0.076221
 6444           0.100499
 6445           0.129194
 6446           0.162210
 6447           0.198825
 6448           0.237959
 6449           0.277871
 6450           0.316244
 6451           0.350347
 6452           0.377390
 6453           0.394752
 6454           0.400747
 6455           0.394752
 6456           0.377390
 6457           0.350347
 6458           0.316244
 6459           0.277871
 6460           0.237959
 6461           0.198825
 6462           0.162210
 6463           0.129194
 6464           0.100499
 6465           0.076221
 6466           0.056207
 6467           0.040107
 6468           0.027463
 6469           0.017757
 6470           0.010535
 6471           0.005319
 6472           0.001702
 6473          -0.000664
 6474          -0.002073
 6475          -0.002763
 6476          -0.002931
 6477          -0.002732
 6478          -0.002286
 6479          -0.001686
 6480          -0.001002
 6481          -0.000287
 6482           0.000422
 6483           0.001099
 6484           0.001728
 6485           0.002296
 6486           0.002802
 6487           0.003240
 6488           0.003613
 6489           0.003924
 6490           0.004178
 6491           0.004378
 6492           0.004531
 6493           0.004641
 6494           0.004713
 6495           0.004753
 6496           0.004765
 6497           0.004752
 6498           0.004720
 6499           0.004670
 6500           0.004606
 6501           0.004531
 6502           0.004447
 6503           0.004356
 6504           0.003324
 6505           0.003361
 6506           0.003389
 6507           0.003405
 6508           0.003406
 6509           0.003390
 6510           0.003355
 6511           0.003296
 6512           0.003210
 6513           0.003092
 6514           0.002939
 6515           0.002744
 6516           0.002503
 6517           0.002210
 6518           0.001857
 6519           0.001442
 6520           0.000957
 6521           0.000397
 6522          -0.000241
 6523          -0.000962
 6524          -0.001763
 6525          -0.002642
 6526          -0.003591
 6527          -0.004597
 6528          -0.005639
 6529          -0.006688
 6530          -0.007699
 6531          -0.008617
 6532          -0.009367
 6533          -0.009855
 6534          -0.009959
 6535          -0.009529
 6536          -0.008385
 6537          -0.006312
 6538          -0.003057
 6539           0.001665
 6540           0.008181
 6541           0.016794
 6542           0.027829
 6543           0.041560
 6544           0.058167
 6545           0.077677
 6546           0.099943
 6547           0.124399
 6548           0.150261
 6549           0.176349
 6550           0.201164
 6551           0.223005
 6552           0.240185
 6553           0.251151
 6554           0.254927
 6555           0.251151
 6556           0.240185
 6557           0.223005
 6558           0.201164
 6559           0.176349
 6560           0.150261
 6561           0.124399
 6562           0.099943
 6563           0.077677
 6564           0.058167
 6565           0.041560
 6566           0.027829
 6567           0.016794
 6568           0.008181
 6569           0.001665
 6570          -0.003057
 6571          -0.006312
 6572          -0.008385
 6573          -0.009529
 6574          -0.009959
 6575          -0.009855
 6576          -0.009367
 6577          -0.008617
 6578          -0.007699
 6579          -0.006688
 6580          -0.005639
 6581          -0.004597
 6582          -0.003591
 6583          -0.002642
 6584          -0.001763
 6585          -0.000962
 6586          -0.000241
 6587           0.000397
 6588           0.000957
 6589           0.001442
 6590           0.001857
 6591           0.002210
 6592           0.002503
 6593           0.002744
 6594           0.002939
 6595           0.003092
 6596           0.003210
 6597           0.003296
 6598           0.003355
 6599           0.003390
 6600           0.003406
 6601           0.003405
 6602           0.003389
 6603           0.003361
 6604           0.002421
 6605           0.002403
 6606           0.002371
 6607           0.002323
 6608           0.002257
 6609           0.002168
 6610           0.002054
 6611           0.001912
 6612           0.001736
 6613           0.001523
 6614           0.001268
 6615           0.000964
 6616           0.000607
 6617           0.000189
 6618          -0.000295
 6619          -0.000852
 6620          -0.001487
 6621          -0.002207
 6622          -0.003016
 6623          -0.003921
 6624          -0.004917
 6625          -0.006007
 6626          -0.007184
 6627          -0.008439
 6628          -0.009756
 6629          -0.011113
 6630          -0.012471
 6631          -0.013788
 6632          -0.015005
 6633          -0.016049
 6634          -0.016825
 6635          -0.017223
 6636          -0.017105
 6637          -0.016314
 6638          -0.014670
 6639          -0.011974
 6640          -0.008004
 6641          -0.002561
 6642           0.004563
 6643           0.013524
 6644           0.024413
 6645           0.037207
 6646           0.051763
 6647           0.067666
 6648           0.084374
 6649           0.101103
 6650           0.116899
 6651           0.130707
 6652           0.141506
 6653           0.148371
 6654           0.150730
 6655           0.148371
 6656           0.141506
 6657           0.130707
 6658           0.116899
 6659           0.101103
 6660           0.084374
 6661           0.067666
 6662           0.051763
 6663           0.037207
 6664           0.024413
 6665           0.013524
 6666           0.004563
 6667          -0.002561
 6668          -0.008004
 6669          -0.011974
 6670          -0.014670
 6671          -0.016314
 6672          -0.017105
 6673          -0.017223
 6674          -0.016825
 6675          -0.016049
 6676          -0.015005
 6677          -0.013788
 6678          -0.012471
 6679          -0.011113
 6680          -0.009756
 6681          -0.008439
 6682          -0.007184
 6683          -0.006007
 6684          -0.004917
 6685          -0.003921
 6686          -0.003016
 6687          -0.002207
 6688          -0.001487
 6689          -0.000852
 6690          -0.000295
 6691           0.000189
 6692           0.000607
 6693           0.000964
 6694           0.001268
 6695           0.001523
 6696           0.001736
 6697           0.001912
 6698           0.002054
 6699           0.002168
 6700           0.002257
 6701           0.002323
 6702           0.002371
 6703           0.002403
 6704           0.001550
 6705           0.001481
 6706           0.001394
 6707           0.001287
 6708           0.001158
 6709           0.001002
 6710           0.000818
 6711           0.000599
 6712           0.000343
 6713           0.000044
 6714          -0.000303
 6715          -0.000703
 6716          -0.001163
 6717          -0.001689
 6718          -0.002289
 6719          -0.002968
 6720          -0.003732
 6721          -0.004588
 6722          -0.005542
 6723          -0.006601
 6724          -0.007761
 6725          -0.009027
 6726          -0.010394
 6727          -0.011855
 6728          -0.013399
 6729          -0.015011
 6730          -0.016657
 6731          -0.018305
 6732          -0.019912
 6733          -0.021418
 6734          -0.022753
 6735          -0.023835
 6736          -0.024558
 6737          -0.024811
 6738          -0.024464
 6739          -0.023378
 6740          -0.021404
 6741          -0.018407
 6742          -0.014259
 6743          -0.008864
 6744          -0.002182
 6745           0.005755
 6746           0.014830
 6747           0.024756
 6748           0.035171
 6749           0.045568
 6750           0.055349
 6751           0.063867
 6752           0.070506
 6753           0.074716
 6754           0.076160
 6755           0.074716
 6756           0.070506
 6757           0.063867
 6758           0.055349
 6759           0.045568
 6760           0.035171
 6761           0.024756
 6762           0.014830
 6763           0.005755
 6764          -0.002182
 6765          -0.008864
 6766          -0.014259
 6767          -0.018407
 6768          -0.021404
 6769          -0.023378
 6770          -0.024464
 6771          -0.024811
 6772          -0.024558
 6773          -0.023835
 6774          -0.022753
 6775          -0.021418
 6776          -0.019912
 6777          -0.018305
 6778          -0.016657
 6779          -0.015011
 6780          -0.013399
 6781          -0.011855
 6782          -0.010394
 6783          -0.009027
 6784          -0.007761
 6785          -0.006601
 6786          -0.005542
 6787          -0.004588
 6788          -0.003732
 6789          -0.002968
 6790          -0.002289
 6791          -0.001689
 6792          -0.001163
 6793          -0.000703
 6794          -0.000303
 6795           0.000044
 6796           0.000343
 6797           0.000599
 6798           0.000818
 6799           0.001002
 6800           0.001158
 6801           0.001287
 6802           0.001394
 6803           0.001481
 6804           0.000714
 6805           0.000596
 6806           0.000459
 6807           0.000298
 6808           0.000111
 6809          -0.000105
 6810          -0.000355
 6811          -0.000641
 6812          -0.000971
 6813          -0.001346
 6814          -0.001774
 6815          -0.002259
 6816          -0.002809
 6817          -0.003429
 6818          -0.004129
 6819          -0.004911
 6820          -0.005785
 6821          -0.006757
 6822          -0.007832
 6823          -0.009017
 6824          -0.010313
 6825          -0.011723
 6826          -0.013244
 6827          -0.014874
 6828          -0.016602
 6829          -0.018420
 6830          -0.020299
 6831          -0.022217
 6832          -0.024138
 6833          -0.026019
 6834          -0.027806
 6835          -0.029439
 6836          -0.030836
 6837          -0.031917
 6838          -0.032590
 6839          -0.032756
 6840          -0.032314
 6841          -0.031171
 6842          -0.029243
 6843          -0.026473
 6844          -0.022839
 6845          -0.018369
 6846          -0.013147
 6847          -0.007359
 6848          -0.001240
 6849           0.004892
 6850           0.010671
 6851           0.015704
 6852           0.019626
 6853           0.022112
 6854           0.022964
 6855           0.022112
 6856           0.019626
 6857           0.015704
 6858           0.010671
 6859           0.004892
 6860          -0.001240
 6861          -0.007359
 6862          -0.013147
 6863          -0.018369
 6864          -0.022839
 6865          -0.026473
 6866          -0.029243
 6867          -0.031171
 6868          -0.032314
 6869          -0.032756
 6870          -0.032590
 6871          -0.031917
 6872          -0.030836
 6873          -0.029439
 6874          -0.027806
 6875          -0.026019
 6876          -0.024138
 6877          -0.022217
 6878          -0.020299
 6879          -0.018420
 6880          -0.016602
 6881          -0.014874
 6882          -0.013244
 6883          -0.011723
 6884          -0.010313
 6885          -0.009017
 6886          -0.007832
 6887          -0.006757
 6888          -0.005785
 6889          -0.004911
 6890          -0.004129
 6891          -0.003429
 6892          -0.002809
 6893          -0.002259
 6894          -0.001774
 6895          -0.001346
 6896          -0.000971
 6897          -0.000641
 6898          -0.000355
 6899          -0.000105
 6900           0.000111
 6901           0.000298
 6902           0.000459
 6903           0.000596
 6904          -0.000088
 6905          -0.000249
 6906          -0.000434
 6907          -0.000644
 6908          -0.000883
 6909          -0.001155
 6910          -0.001462
 6911          -0.001810
 6912          -0.002204
 6913          -0.002647
 6914          -0.003146
 6915          -0.003706
 6916          -0.004333
 6917          -0.005035
 6918          -0.005819
 6919          -0.006690
 6920          -0.007655
 6921          -0.008722
 6922          -0.009896
 6923          -0.011186
 6924          -0.012591
 6925          -0.014115
 6926          -0.015760
 6927          -0.017522
 6928          -0.019396
 6929          -0.021376
 6930          -0.023438
 6931          -0.025567
 6932          -0.027736
 6933          -0.029912
 6934          -0.032053
 6935          -0.034115
 6936          -0.036036
 6937          -0.037757
 6938          -0.039211
 6939          -0.040327
 6940          -0.041036
 6941          -0.041267
 6942          -0.040966
 6943          -0.040093
 6944          -0.038632
 6945          -0.036602
 6946          -0.034054
 6947          -0.031102
 6948          -0.027892
 6949          -0.024612
 6950          -0.021485
 6951          -0.018738
 6952          -0.016588
 6953          -0.015221
 6954          -0.014752
 6955          -0.015221
 6956          -0.016588
 6957          -0.018738
 6958          -0.021485
 6959          -0.024612
 6960          -0.027892
 6961          -0.031102
 6962          -0.034054
 6963          -0.036602
 6964          -0.038632
 6965          -0.040093
 6966          -0.040966
 6967          -0.041267
 6968          -0.041036
 6969          -0.040327
 6970          -0.039211
 6971          -0.037757
 6972          -0.036036
 6973          -0.034115
 6974          -0.032053
 6975          -0.029912
 6976          -0.027736
 6977          -0.025567
 6978          -0.023438
 6979          -0.021376
 6980          -0.019396
 6981          -0.017522
 6982          -0.015760
 6983          -0.014115
 6984          -0.012591
 6985          -0.011186
 6986          -0.009896
 6987          -0.008722
 6988          -0.007655
 6989          -0.006690
 6990          -0.005819
 6991          -0.005035
 6992          -0.004333
 6993          -0.003706
 6994          -0.003146
 6995          -0.002647
 6996          -0.002204
 6997          -0.001810
 6998          -0.001462
 6999          -0.001155
 7000          -0.000883
 7001          -0.000644
 7002          -0.000434
 7003          -0.000249
 7004          -0.000856
 7005          -0.001057
 7006          -0.001285
 7007          -0.001540
 7008          -0.001826
 7009          -0.002147
 7010          -0.002507
 7011          -0.002909
 7012          -0.003360
 7013          -0.003863
 7014          -0.004424
 7015          -0.005048
 7016          -0.005742
 7017          -0.006513
 7018          -0.007368
 7019          -0.008312
 7020          -0.009353
 7021          -0.010497
 7022          -0.011752
 7023          -0.013124
 7024          -0.014615
 7025          -0.016229
 7026          -0.017968
 7027          -0.019831
 7028          -0.021815
 7029          -0.023917
 7030          -0.026118
 7031          -0.028406
 7032          -0.030762
 7033          -0.033159
 7034          -0.035567
 7035          -0.037953
 7036          -0.040268
 7037          -0.042468
 7038          -0.044502
 7039          -0.046320
 7040          -0.047870
 7041          -0.049098
 7042          -0.049964
 7043          -0.050437
 7044          -0.050501
 7045          -0.050161
 7046          -0.049446
 7047          -0.048415
 7048          -0.047153
 7049          -0.045768
 7050          -0.044384
 7051          -0.043131
 7052          -0.042131
 7053          -0.041488
 7054          -0.041266
 7055          -0.041488
 7056          -0.042131
 7057          -0.043131
 7058          -0.044384
 7059          -0.045768
 7060          -0.047153
 7061          -0.048415
 7062          -0.049446
 7063          -0.050161
 7064          -0.050501
 7065          -0.050437
 7066          -0.049964
 7067          -0.049098
 7068          -0.047870
 7069          -0.046320
 7070          -0.044502
 7071          -0.042468
 7072          -0.040268
 7073          -0.037953
 7074          -0.035567
 7075          -0.033159
 7076          -0.030762
 7077          -0.028406
 7078          -0.026118
 7079          -0.023917
 7080          -0.021815
 7081          -0.019831
 7082          -0.017968
 7083          -0.016229
 7084          -0.014615
 7085          -0.013124
 7086          -0.011752
 7087          -0.010497
 7088          -0.009353
 7089          -0.008312
 7090          -0.007368
 7091          -0.006513
 7092          -0.005742
 7093          -0.005048
 7094          -0.004424
 7095          -0.003863
 7096          -0.003360
 7097          -0.002909
 7098          -0.002507
 7099          -0.002147
 7100          -0.001826
 7101          -0.001540
 7102          -0.001285
 7103          -0.001057
 7104          -0.001586
 7105          -0.001824
 7106          -0.002090
 7107          -0.002386
 7108          -0.002714
 7109          -0.003080
 7110          -0.003485
 7111          -0.003936
 7112          -0.004437
 7113          -0.004991
 7114          -0.005605
 7115          -0.006284
 7116          -0.007034
 7117          -0.007863
 7118          -0.008777
 7119          -0.009780
 7120          -0.010882
 7121          -0.012087
 7122          -0.013404
 7123          -0.014840
 7124          -0.016394
 7125          -0.018074
 7126          -0.019882
 7127          -0.021818
 7128          -0.023880
 7129          -0.026068
 7130          -0.028366
 7131          -0.030766
 7132          -0.033253
 7133          -0.035807
 7134          -0.038404
 7135          -0.041022
 7136          -0.043616
 7137          -0.046157
 7138          -0.048606
 7139          -0.050920
 7140          -0.053064
 7141          -0.054990
 7142          -0.056667
 7143          -0.058066
 7144          -0.059170
 7145          -0.059973
 7146          -0.060484
 7147          -0.060728
 7148          -0.060747
 7149          -0.060596
 7150          -0.060340
 7151          -0.060048
 7152          -0.059784
 7153          -0.059603
 7154          -0.059538
 7155          -0.059603
 7156          -0.059784
 7157          -0.060048
 7158          -0.060340
 7159          -0.060596
 7160          -0.060747
 7161          -0.060728
 7162          -0.060484
 7163          -0.059973
 7164          -0.059170
 7165          -0.058066
 7166          -0.056667
 7167          -0.054990
 7168          -0.053064
 7169          -0.050920
 7170          -0.048606
 7171          -0.046157
 7172          -0.043616
 7173          -0.041022
 7174          -0.038404
 7175          -0.035807
 7176          -0.033253
 7177          -0.030766
 7178          -0.028366
 7179          -0.026068
 7180          -0.023880
 7181          -0.021818
 7182          -0.019882
 7183          -0.018074
 7184          -0.016394
 7185          -0.014840
 7186          -0.013404
 7187          -0.012087
 7188          -0.010882
 7189          -0.009780
 7190          -0.008777
 7191          -0.007863
 7192          -0.007034
 7193          -0.006284
 7194          -0.005605
 7195          -0.004991
 7196          -0.004437
 7197          -0.003936
 7198          -0.003485
 7199          -0.003080
 7200          -0.002714
 7201          -0.002386
 7202          -0.002090
 7203          -0.001824
 7204          -0.002279
 7205          -0.002551
 7206          -0.002852
 7207          -0.003183
 7208          -0.003550
 7209          -0.003954
 7210          -0.004400
 7211          -0.004893
 7212          -0.005437
 7213          -0.006035
 7214          -0.006694
 7215          -0.007419
 7216          -0.008217
 7217          -0.009092
 7218          -0.010054
 7219          -0.011104
 7220          -0.012253
 7221          -0.013505
 7222          -0.014868
 7223          -0.016350
 7224          -0.017950
 7225          -0.019675
 7226          -0.021529
 7227          -0.023512
 7228          -0.025624
 7229          -0.027867
 7230          -0.030225
 7231          -0.032696
 7232          -0.035266
 7233          -0.037921
 7234          -0.040641
 7235          -0.043410
 7236          -0.046191
 7237          -0.048960
 7238          -0.051686
 7239          -0.054335
 7240          -0.056879
 7241          -0.059273
 7242          -0.061492
 7243          -0.063509
 7244          -0.065302
 7245          -0.066856
 7246          -0.068168
 7247          -0.069238
 7248          -0.070081
 7249          -0.070719
 7250          -0.071181
 7251          -0.071497
 7252          -0.071697
 7253          -0.071806
 7254          -0.071840
 7255          -0.071806
 7256          -0.071697
 7257          -0.071497
 7258          -0.071181
 7259          -0.070719
 7260          -0.070081
 7261          -0.069238
 7262          -0.068168
 7263          -0.066856
 7264          -0.065302
 7265          -0.063509
 7266          -0.061492
 7267          -0.059273
 7268          -0.056879
 7269          -0.054335
 7270          -0.051686
 7271          -0.048960
 7272          -0.046191
 7273          -0.043410
 7274          -0.040641
 7275          -0.037921
 7276          -0.035266
 7277          -0.032696
 7278          -0.030225
 7279          -0.027867
 7280          -0.025624
 7281          -0.023512
 7282          -0.021529
 7283          -0.019675
 7284          -0.017950
 7285          -0.016350
 7286          -0.014868
 7287          -0.013505
 7288          -0.012253
 7289          -0.011104
 7290          -0.010054
 7291          -0.009092
 7292          -0.008217
 7293          -0.007419
 7294          -0.006694
 7295          -0.006035
 7296          -0.005437
 7297          -0.004893
 7298          -0.004400
 7299          -0.003954
 7300          -0.003550
 7301          -0.003183
 7302          -0.002852
 7303          -0.002551
 7304          -0.002936
 7305          -0.003238
 7306          -0.003569
 7307          -0.003933
 7308          -0.004333
 7309          -0.004772
 7310          -0.005253
 7311          -0.005782
 7312          -0.006362
 7313          -0.006998
 7314          -0.007695
 7315          -0.008458
 7316          -0.009293
 7317          -0.010206
 7318          -0.011205
 7319          -0.012292
 7320          -0.013476
 7321          -0.014762
 7322          -0.016157
 7323          -0.017669
 7324          -0.019298
 7325          -0.021050
 7326          -0.022930
 7327          -0.024938
 7328          -0.027075
 7329          -0.029344
 7330          -0.031732
 7331          -0.034237
 7332          -0.036849
 7333          -0.039557
 7334          -0.042344
 7335          -0.045198
 7336          -0.048088
 7337          -0.050993
 7338          -0.053887
 7339          -0.056743
 7340          -0.059535
 7341          -0.062223
 7342          -0.064784
 7343          -0.067189
 7344          -0.069415
 7345          -0.071440
 7346          -0.073254
 7347          -0.074838
 7348          -0.076193
 7349          -0.077320
 7350          -0.078223
 7351          -0.078912
 7352          -0.079396
 7353          -0.079682
 7354          -0.079777
 7355          -0.079682
 7356          -0.079396
 7357          -0.078912
 7358          -0.078223
 7359          -0.077320
 7360          -0.076193
 7361          -0.074838
 7362          -0.073254
 7363          -0.071440
 7364          -0.069415
 7365          -0.067189
 7366          -0.064784
 7367          -0.062223
 7368          -0.059535
 7369          -0.056743
 7370          -0.053887
 7371          -0.050993
 7372          -0.048088
 7373          -0.045198
 7374          -0.042344
 7375          -0.039557
 7376          -0.036849
 7377          -0.034237
 7378          -0.031732
 7379          -0.029344
 7380          -0.027075
 7381          -0.024938
 7382          -0.022930
 7383          -0.021050
 7384          -0.019298
 7385          -0.017669
 7386          -0.016157
 7387          -0.014762
 7388          -0.013476
 7389          -0.012292
 7390          -0.011205
 7391          -0.010206
 7392          -0.009293
 7393          -0.008458
 7394          -0.007695
 7395          -0.006998
 7396          -0.006362
 7397          -0.005782
 7398          -0.005253
 7399          -0.004772
 7400          -0.004333
 7401          -0.003933
 7402          -0.003569
 7403          -0.003238
 7404          -0.003556
 7405          -0.003885
 7406          -0.004244
 7407          -0.004636
 7408          -0.005065
 7409          -0.005533
 7410          -0.006045
 7411          -0.006604
 7412          -0.007216
 7413          -0.007883
 7414          -0.008611
 7415          -0.009405
 7416          -0.010270
 7417          -0.011211
 7418          -0.012238
 7419          -0.013351
 7420          -0.014559
 7421          -0.015868
 7422          -0.017282
 7423          -0.018811
 7424          -0.020454
 7425          -0.022218
 7426          -0.024105
 7427          -0.026119
 7428          -0.028260
 7429          -0.030532
 7430          -0.032922
 7431          -0.035431
 7432          -0.038050
 7433          -0.040770
 7434          -0.043577
 7435          -0.046463
 7436          -0.049397
 7437          -0.052364
 7438          -0.055341
 7439          -0.058303
 7440          -0.061228
 7441          -0.064078
 7442          -0.066830
 7443          -0.069457
 7444          -0.071932
 7445          -0.074230
 7446          -0.076335
 7447          -0.078220
 7448          -0.079873
 7449          -0.081284
 7450          -0.082445
 7451          -0.083351
 7452          -0.084000
 7453          -0.084389
 7454          -0.084519
 7455          -0.084389
 7456          -0.084000
 7457          -0.083351
 7458          -0.082445
 7459          -0.081284
 7460          -0.079873
 7461          -0.078220
 7462          -0.076335
 7463          -0.074230
 7464          -0.071932
 7465          -0.069457
 7466          -0.066830
 7467          -0.064078
 7468          -0.061228
 7469          -0.058303
 7470          -0.055341
 7471          -0.052364
 7472          -0.049397
 7473          -0.046463
 7474          -0.043577
 7475          -0.040770
 7476          -0.038050
 7477          -0.035431
 7478          -0.032922
 7479          -0.030532
 7480          -0.028260
 7481          -0.026119
 7482          -0.024105
 7483          -0.022218
 7484          -0.020454
 7485          -0.018811
 7486          -0.017282
 7487          -0.015868
 7488          -0.014559
 7489          -0.013351
 7490          -0.012238
 7491          -0.011211
 7492          -0.010270
 7493          -0.009405
 7494          -0.008611
 7495          -0.007883
 7496          -0.007216
 7497          -0.006604
 7498          -0.006045
 7499          -0.005533
 7500          -0.005065
 7501          -0.004636
 7502          -0.004244
 7503          -0.003885
 7504          -0.004141
 7505          -0.004493
 7506          -0.004876
 7507          -0.005293
 7508          -0.005747
 7509          -0.006240
 7510          -0.006778
 7511          -0.007363
 7512          -0.008000
 7513          -0.008692
 7514          -0.009445
 7515          -0.010262
 7516          -0.011150
 7517          -0.012113
 7518          -0.013159
 7519          -0.014290
 7520          -0.015512
 7521          -0.016832
 7522          -0.018255
 7523          -0.019789
 7524          -0.021433
 7525          -0.023194
 7526          -0.025075
 7527          -0.027078
 7528          -0.029205
 7529          -0.031459
 7530          -0.033829
 7531          -0.036316
 7532          -0.038913
 7533          -0.041612
 7534          -0.044401
 7535          -0.047274
 7536          -0.050202
 7537          -0.053173
 7538          -0.056165
 7539          -0.059157
 7540          -0.062128
 7541          -0.065041
 7542          -0.067875
 7543          -0.070603
 7544          -0.073195
 7545          -0.075628
 7546          -0.077878
 7547          -0.079915
 7548          -0.081721
 7549          -0.083280
 7550          -0.084575
 7551          -0.085595
 7552          -0.086331
 7553          -0.086774
 7554          -0.086922
 7555          -0.086774
 7556          -0.086331
 7557          -0.085595
 7558          -0.084575
 7559          -0.083280
 7560          -0.081721
 7561          -0.079915
 7562          -0.077878
 7563          -0.075628
 7564          -0.073195
 7565          -0.070603
 7566          -0.067875
 7567          -0.065041
 7568          -0.062128
 7569          -0.059157
 7570          -0.056165
 7571          -0.053173
 7572          -0.050202
 7573          -0.047274
 7574          -0.044401
 7575          -0.041612
 7576          -0.038913
 7577          -0.036316
 7578          -0.033829
 7579          -0.031459
 7580          -0.029205
 7581          -0.027078
 7582          -0.025075
 7583          -0.023194
 7584          -0.021433
 7585          -0.019789
 7586          -0.018255
 7587          -0.016832
 7588          -0.015512
 7589          -0.014290
 7590          -0.013159
 7591          -0.012113
 7592          -0.011150
 7593          -0.010262
 7594          -0.009445
 7595          -0.008692
 7596          -0.008000
 7597          -0.007363
 7598          -0.006778
 7599          -0.006240
 7600          -0.005747
 7601          -0.005293
 7602          -0.004876
 7603          -0.004493
 7604          -0.004691
 7605          -0.005063
 7606          -0.005468
 7607          -0.005906
 7608          -0.006381
 7609          -0.006896
 7610          -0.007455
 7611          -0.008060
 7612          -0.008719
 7613          -0.009431
 7614          -0.010202
 7615          -0.011037
 7616          -0.011941
 7617          -0.012918
 7618          -0.013976
 7619          -0.015116
 7620          -0.016345
 7621          -0.017668
 7622          -0.019090
 7623          -0.020618
 7624          -0.022253
 7625          -0.023999
 7626          -0.025860
 7627          -0.027838
 7628          -0.029935
 7629          -0.032154
 7630          -0.034485
 7631          -0.036929
 7632          -0.039481
 7633          -0.042132
 7634          -0.044873
 7635          -0.047698
 7636          -0.050581
 7637          -0.053510
 7638          -0.056467
 7639          -0.059431
 7640          -0.062383
 7641          -0.065287
 7642          -0.068123
 7643          -0.070865
 7644          -0.073483
 7645          -0.075951
 7646          -0.078247
 7647          -0.080335
 7648          -0.082197
 7649          -0.083812
 7650          -0.085160
 7651          -0.086225
 7652          -0.086997
 7653          -0.087462
 7654          -0.087618
 7655          -0.087462
 7656          -0.086997
 7657          -0.086225
 7658          -0.085160
 7659          -0.083812
 7660          -0.082197
 7661          -0.080335
 7662          -0.078247
 7663          -0.075951
 7664          -0.073483
 7665          -0.070865
 7666          -0.068123
 7667          -0.065287
 7668          -0.062383
 7669          -0.059431
 7670          -0.056467
 7671          -0.053510
 7672          -0.050581
 7673          -0.047698
 7674          -0.044873
 7675          -0.042132
 7676          -0.039481
 7677          -0.036929
 7678          -0.034485
 7679          -0.032154
 7680          -0.029935
 7681          -0.027838
 7682          -0.025860
 7683          -0.023999
 7684          -0.022253
 7685          -0.020618
 7686          -0.019090
 7687          -0.017668
 7688          -0.016345
 7689          -0.015116
 7690          -0.013976
 7691          -0.012918
 7692          -0.011941
 7693          -0.011037
 7694          -0.010202
 7695          -0.009431
 7696          -0.008719
 7697          -0.008060
 7698          -0.007455
 7699          -0.006896
 7700          -0.006381
 7701          -0.005906
 7702          -0.005468
 7703          -0.005063
 7704          -0.005204
 7705          -0.005594
 7706          -0.006017
 7707          -0.006473
 7708          -0.006966
 7709          -0.007499
 7710          -0.008075
 7711          -0.008697
 7712          -0.009371
 7713          -0.010098
 7714          -0.010883
 7715          -0.011731
 7716          -0.012644
 7717          -0.013630
 7718          -0.014693
 7719          -0.015835
 7720          -0.017063
 7721          -0.018380
 7722          -0.019793
 7723          -0.021308
 7724          -0.022922
 7725          -0.024644
 7726          -0.026475
 7727          -0.028417
 7728          -0.030471
 7729          -0.032642
 7730          -0.034919
 7731          -0.037304
 7732          -0.039791
 7733          -0.042374
 7734          -0.045043
 7735          -0.047793
 7736          -0.050601
 7737          -0.053456
 7738          -0.056339
 7739          -0.059232
 7740          -0.062118
 7741          -0.064962
 7742          -0.067744
 7743          -0.070438
 7744          -0.073017
 7745          -0.075455
 7746          -0.077727
 7747          -0.079800
 7748          -0.081653
 7749          -0.083263
 7750          -0.084610
 7751          -0.085676
 7752          -0.086450
 7753          -0.086917
 7754          -0.087073
 7755          -0.086917
 7756          -0.086450
 7757          -0.085676
 7758          -0.084610
 7759          -0.083263
 7760          -0.081653
 7761          -0.079800
 7762          -0.077727
 7763          -0.075455
 7764          -0.073017
 7765          -0.070438
 7766          -0.067744
 7767          -0.064962
 7768          -0.062118
 7769          -0.059232
 7770          -0.056339
 7771          -0.053456
 7772          -0.050601
 7773          -0.047793
 7774          -0.045043
 7775          -0.042374
 7776          -0.039791
 7777          -0.037304
 7778          -0.034919
 7779          -0.032642
 7780          -0.030471
 7781          -0.028417
 7782          -0.026475
 7783          -0.024644
 7784          -0.022922
 7785          -0.021308
 7786          -0.019793
 7787          -0.018380
 7788          -0.017063
 7789          -0.015835
 7790          -0.014693
 7791          -0.013630
 7792          -0.012644
 7793          -0.011731
 7794          -0.010883
 7795          -0.010098
 7796          -0.009371
 7797          -0.008697
 7798          -0.008075
 7799          -0.007499
 7800          -0.006966
 7801          -0.006473
 7802          -0.006017
 7803          -0.005594
 7804          -0.005684
 7805          -0.006089
 7806          -0.006527
 7807          -0.006998
 7808          -0.007505
 7809          -0.008052
 7810          -0.008642
 7811          -0.009277
 7812          -0.009963
 7813          -0.010700
 7814          -0.011494
 7815          -0.012349
 7816          -0.013267
 7817          -0.014255
 7818          -0.015318
 7819          -0.016456
 7820          -0.017676
 7821          -0.018982
 7822          -0.020378
 7823          -0.021871
 7824          -0.023459
 7825          -0.025147
 7826          -0.026939
 7827          -0.028835
 7828          -0.030837
 7829          -0.032949
 7830          -0.035160
 7831          -0.037472
 7832          -0.039881
 7833          -0.042379
 7834          -0.044959
 7835          -0.047616
 7836          -0.050326
 7837          -0.053081
 7838          -0.055864
 7839          -0.058656
 7840          -0.061442
 7841          -0.064189
 7842          -0.066878
 7843          -0.069484
 7844          -0.071981
 7845          -0.074343
 7846          -0.076547
 7847          -0.078560
 7848          -0.080361
 7849          -0.081927
 7850          -0.083239
 7851          -0.084278
 7852          -0.085032
 7853          -0.085488
 7854          -0.085641
 7855          -0.085488
 7856          -0.085032
 7857          -0.084278
 7858          -0.083239
 7859          -0.081927
 7860          -0.080361
 7861          -0.078560
 7862          -0.076547
 7863          -0.074343
 7864          -0.071981
 7865          -0.069484
 7866          -0.066878
 7867          -0.064189
 7868          -0.061442
 7869          -0.058656
 7870          -0.055864
 7871          -0.053081
 7872          -0.050326
 7873          -0.047616
 7874          -0.044959
 7875          -0.042379
 7876          -0.039881
 7877          -0.037472
 7878          -0.035160
 7879          -0.032949
 7880          -0.030837
 7881          -0.028835
 7882          -0.026939
 7883          -0.025147
 7884          -0.023459
 7885          -0.021871
 7886          -0.020378
 7887          -0.018982
 7888          -0.017676
 7889          -0.016456
 7890          -0.015318
 7891          -0.014255
 7892          -0.013267
 7893          -0.012349
 7894          -0.011494
 7895          -0.010700
 7896          -0.009963
 7897          -0.009277
 7898          -0.008642
 7899          -0.008052
 7900          -0.007505
 7901          -0.006998
 7902          -0.006527
 7903          -0.006089
 7904          -0.006130
 7905          -0.006548
 7906          -0.006998
 7907          -0.007481
 7908          -0.008000
 7909          -0.008558
 7910          -0.009158
 7911          -0.009802
 7912          -0.010496
 7913          -0.011240
 7914          -0.012039
 7915          -0.012896
 7916          -0.013815
 7917          -0.014800
 7918          -0.015857
 7919          -0.016986
 7920          -0.018192
 7921          -0.019481
 7922          -0.020854
 7923          -0.022319
 7924          -0.023874
 7925          -0.025522
 7926          -0.027268
 7927          -0.029111
 7928          -0.031053
 7929          -0.033097
 7930          -0.035234
 7931          -0.037464
 7932          -0.039784
 7933          -0.042186
 7934          -0.044664
 7935          -0.047213
 7936          -0.049811
 7937          -0.052449
 7938          -0.055112
 7939          -0.057783
 7940          -0.060447
 7941          -0.063073
 7942          -0.065643
 7943          -0.068134
 7944          -0.070520
 7945          -0.072778
 7946          -0.074885
 7947          -0.076810
 7948          -0.078532
 7949          -0.080030
 7950          -0.081285
 7951          -0.082280
 7952          -0.083002
 7953          -0.083438
 7954          -0.083584
 7955          -0.083438
 7956          -0.083002
 7957          -0.082280
 7958          -0.081285
 7959          -0.080030
 7960          -0.078532
 7961          -0.076810
 7962          -0.074885
 7963          -0.072778
 7964          -0.070520
 7965          -0.068134
 7966          -0.065643
 7967          -0.063073
 7968          -0.060447
 7969          -0.057783
 7970          -0.055112
 7971          -0.052449
 7972          -0.049811
 7973          -0.047213
 7974          -0.044664
 7975          -0.042186
 7976          -0.039784
 7977          -0.037464
 7978          -0.035234
 7979          -0.033097
 7980          -0.031053
 7981          -0.029111
 7982          -0.027268
 7983          -0.025522
 7984          -0.023874
 7985          -0.022319
 7986          -0.020854
 7987          -0.019481
 7988          -0.018192
 7989          -0.016986
 7990          -0.015857
 7991          -0.014800
 7992          -0.013815
 7993          -0.012896
 7994          -0.012039
 7995          -0.011240
 7996          -0.010496
 7997          -0.009802
 7998          -0.009158
 7999          -0.008558
 8000          -0.008000
 8001          -0.007481
 8002          -0.006998
 8003          -0.006548
 8004          -0.006544
 8005          -0.006972
 8006          -0.007432
 8007          -0.007925
 8008          -0.008453
 8009          -0.009019
 8010          -0.009626
 8011          -0.010276
 8012          -0.010974
 8013          -0.011721
 8014          -0.012520
 8015          -0.013376
 8016          -0.014291
 8017          -0.015269
 8018          -0.016316
 8019          -0.017431
 8020          -0.018620
 8021          -0.019886
 8022          -0.021232
 8023          -0.022664
 8024          -0.024180
 8025          -0.025784
 8026          -0.027477
 8027          -0.029262
 8028          -0.031138
 8029          -0.033109
 8030          -0.035164
 8031          -0.037306
 8032          -0.039529
 8033          -0.041828
 8034          -0.044195
 8035          -0.046627
 8036          -0.049102
 8037          -0.051613
 8038          -0.054145
 8039          -0.056682
 8040          -0.059209
 8041          -0.061699
 8042          -0.064135
 8043          -0.066494
 8044          -0.068753
 8045          -0.070889
 8046          -0.072882
 8047          -0.074701
 8048          -0.076329
 8049          -0.077745
 8050          -0.078931
 8051          -0.079871
 8052          -0.080553
 8053          -0.080965
 8054          -0.081103
 8055          -0.080965
 8056          -0.080553
 8057          -0.079871
 8058          -0.078931
 8059          -0.077745
 8060          -0.076329
 8061          -0.074701
 8062          -0.072882
 8063          -0.070889
 8064          -0.068753
 8065          -0.066494
 8066          -0.064135
 8067          -0.061699
 8068          -0.059209
 8069          -0.056682
 8070          -0.054145
 8071          -0.051613
 8072          -0.049102
 8073          -0.046627
 8074          -0.044195
 8075          -0.041828
 8076          -0.039529
 8077          -0.037306
 8078          -0.035164
 8079          -0.033109
 8080          -0.031138
 8081          -0.029262
 8082          -0.027477
 8083          -0.025784
 8084          -0.024180
 8085          -0.022664
 8086          -0.021232
 8087          -0.019886
 8088          -0.018620
 8089          -0.017431
 8090          -0.016316
 8091          -0.015269
 8092          -0.014291
 8093          -0.013376
 8094          -0.012520
 8095          -0.011721
 8096          -0.010974
 8097          -0.010276
 8098          -0.009626
 8099          -0.009019
 8100          -0.008453
 8101          -0.007925
 8102          -0.007432
 8103          -0.006972
 8104          -0.006926
 8105          -0.007362
 8106          -0.007830
 8107          -0.008330
 8108          -0.008864
 8109          -0.009436
 8110          -0.010047
 8111          -0.010700
 8112          -0.011400
 8113          -0.012146
 8114          -0.012943
 8115          -0.013794
 8116          -0.014702
 8117          -0.015669
 8118          -0.016702
 8119          -0.017800
 8120          -0.018967
 8121          -0.020206
 8122          -0.021521
 8123          -0.022916
 8124          -0.024389
 8125          -0.025944
 8126          -0.027581
 8127          -0.029303
 8128          -0.031109
 8129          -0.033002
 8130          -0.034972
 8131          -0.037020
 8132          -0.039143
 8133          -0.041334
 8134          -0.043586
 8135          -0.045895
 8136          -0.048242
 8137          -0.050619
 8138          -0.053013
 8139          -0.055408
 8140          -0.057793
 8141          -0.060139
 8142          -0.062431
 8143          -0.064649
 8144          -0.066771
 8145          -0.068777
 8146          -0.070647
 8147          -0.072353
 8148          -0.073878
 8149          -0.075205
 8150          -0.076315
 8151          -0.077195
 8152          -0.077833
 8153          -0.078219
 8154          -0.078348
 8155          -0.078219
 8156          -0.077833
 8157          -0.077195
 8158          -0.076315
 8159          -0.075205
 8160          -0.073878
 8161          -0.072353
 8162          -0.070647
 8163          -0.068777
 8164          -0.066771
 8165          -0.064649
 8166          -0.062431
 8167          -0.060139
 8168          -0.057793
 8169          -0.055408
 8170          -0.053013
 8171          -0.050619
 8172          -0.048242
 8173          -0.045895
 8174          -0.043586
 8175          -0.041334
 8176          -0.039143
 8177          -0.037020
 8178          -0.034972
 8179          -0.033002
 8180          -0.031109
 8181          -0.029303
 8182          -0.027581
 8183          -0.025944
 8184          -0.024389
 8185          -0.022916
 8186          -0.021521
 8187          -0.020206
 8188          -0.018967
 8189          -0.017800
 8190          -0.016702
 8191          -0.015669
 8192          -0.014702
 8193          -0.013794
 8194          -0.012943
 8195          -0.012146
 8196          -0.011400
 8197          -0.010700
 8198          -0.010047
 8199          -0.009436
 8200          -0.008864
 8201          -0.008330
 8202          -0.007830
 8203          -0.007362
 8204          -0.007279
 8205          -0.007721
 8206          -0.008195
 8207          -0.008699
 8208          -0.009238
 8209          -0.009812
 8210          -0.010425
 8211          -0.011079
 8212          -0.011777
 8213          -0.012520
 8214          -0.013312
 8215          -0.014155
 8216          -0.015052
 8217          -0.016006
 8218          -0.017021
 8219          -0.018098
 8220          -0.019240
 8221          -0.020450
 8222          -0.021729
 8223          -0.023084
 8224          -0.024511
 8225          -0.026014
 8226          -0.027593
 8227          -0.029249
 8228          -0.030982
 8229          -0.032794
 8230          -0.034676
 8231          -0.036629
 8232          -0.038648
 8233          -0.040728
 8234          -0.042863
 8235          -0.045047
 8236          -0.047263
 8237          -0.049504
 8238          -0.051758
 8239          -0.054009
 8240          -0.056247
 8241          -0.058445
 8242          -0.060591
 8243          -0.062665
 8244          -0.064648
 8245          -0.066519
 8246          -0.068262
 8247          -0.069851
 8248          -0.071272
 8249          -0.072506
 8250          -0.073538
 8251          -0.074356
 8252          -0.074948
 8253          -0.075307
 8254          -0.075427
 8255          -0.075307
 8256          -0.074948
 8257          -0.074356
 8258          -0.073538
 8259          -0.072506
 8260          -0.071272
 8261          -0.069851
 8262          -0.068262
 8263          -0.066519
 8264          -0.064648
 8265          -0.062665
 8266          -0.060591
 8267          -0.058445
 8268          -0.056247
 8269          -0.054009
 8270          -0.051758
 8271          -0.049504
 8272          -0.047263
 8273          -0.045047
 8274          -0.042863
 8275          -0.040728
 8276          -0.038648
 8277          -0.036629
 8278          -0.034676
 8279          -0.032794
 8280          -0.030982
 8281          -0.029249
 8282          -0.027593
 8283          -0.026014
 8284          -0.024511
 8285          -0.023084
 8286          -0.021729
 8287          -0.020450
 8288          -0.019240
 8289          -0.018098
 8290          -0.017021
 8291          -0.016006
 8292          -0.015052
 8293          -0.014155
 8294          -0.013312
 8295          -0.012520
 8296          -0.011777
 8297          -0.011079
 8298          -0.010425
 8299          -0.009812
 8300          -0.009238
 8301          -0.008699
 8302          -0.008195
 8303          -0.007721
 8304          -0.007601
 8305          -0.008049
 8306          -0.008526
 8307          -0.009033
 8308          -0.009574
 8309          -0.010149
 8310          -0.010761
 8311          -0.011413
 8312          -0.012107
 8313          -0.012844
 8314          -0.013628
 8315          -0.014460
 8316          -0.015344
 8317          -0.016282
 8318          -0.017278
 8319          -0.018331
 8320          -0.019445
 8321          -0.020622
 8322          -0.021865
 8323          -0.023178
 8324          -0.024556
 8325          -0.026004
 8326          -0.027523
 8327          -0.029112
 8328          -0.030771
 8329          -0.032501
 8330          -0.034295
 8331          -0.036151
 8332          -0.038067
 8333          -0.040037
 8334          -0.042053
 8335          -0.044113
 8336          -0.046199
 8337          -0.048305
 8338          -0.050419
 8339          -0.052527
 8340          -0.054620
 8341          -0.056673
 8342          -0.058673
 8343          -0.060605
 8344          -0.062448
 8345          -0.064187
 8346          -0.065804
 8347          -0.067278
 8348          -0.068594
 8349          -0.069736
 8350          -0.070691
 8351          -0.071447
 8352          -0.071995
 8353          -0.072326
 8354          -0.072437
 8355          -0.072326
 8356          -0.071995
 8357          -0.071447
 8358          -0.070691
 8359          -0.069736
 8360          -0.068594
 8361          -0.067278
 8362          -0.065804
 8363          -0.064187
 8364          -0.062448
 8365          -0.060605
 8366          -0.058673
 8367          -0.056673
 8368          -0.054620
 8369          -0.052527
 8370          -0.050419
 8371          -0.048305
 8372          -0.046199
 8373          -0.044113
 8374          -0.042053
 8375          -0.040037
 8376          -0.038067
 8377          -0.036151
 8378          -0.034295
 8379          -0.032501
 8380          -0.030771
 8381          -0.029112
 8382          -0.027523
 8383          -0.026004
 8384          -0.024556
 8385          -0.023178
 8386          -0.021865
 8387          -0.020622
 8388          -0.019445
 8389          -0.018331
 8390          -0.017278
 8391          -0.016282
 8392          -0.015344
 8393          -0.014460
 8394          -0.013628
 8395          -0.012844
 8396          -0.012107
 8397          -0.011413
 8398          -0.010761
 8399          -0.010149
 8400          -0.009574
 8401          -0.009033
 8402          -0.008526
 8403          -0.008049
 8404          -0.007896
 8405          -0.008346
 8406          -0.008825
 8407          -0.009334
 8408          -0.009874
 8409          -0.010448
 8410          -0.011058
 8411          -0.011705
 8412          -0.012394
 8413          -0.013123
 8414          -0.013896
 8415          -0.014716
 8416          -0.015584
 8417          -0.016503
 8418          -0.017477
 8419          -0.018505
 8420          -0.019589
 8421          -0.020732
 8422          -0.021936
 8423          -0.023204
 8424          -0.024533
 8425          -0.025926
 8426          -0.027383
 8427          -0.028904
 8428          -0.030488
 8429          -0.032137
 8430          -0.033842
 8431          -0.035603
 8432          -0.037417
 8433          -0.039277
 8434          -0.041178
 8435          -0.043115
 8436          -0.045074
 8437          -0.047047
 8438          -0.049025
 8439          -0.050993
 8440          -0.052944
 8441          -0.054855
 8442          -0.056714
 8443          -0.058506
 8444          -0.060215
 8445          -0.061824
 8446          -0.063320
 8447          -0.064681
 8448          -0.065896
 8449          -0.066949
 8450          -0.067829
 8451          -0.068526
 8452          -0.069030
 8453          -0.069335
 8454          -0.069437
 8455          -0.069335
 8456          -0.069030
 8457          -0.068526
 8458          -0.067829
 8459          -0.066949
 8460          -0.065896
 8461          -0.064681
 8462          -0.063320
 8463          -0.061824
 8464          -0.060215
 8465          -0.058506
 8466          -0.056714
 8467          -0.054855
 8468          -0.052944
 8469          -0.050993
 8470          -0.049025
 8471          -0.047047
 8472          -0.045074
 8473          -0.043115
 8474          -0.041178
 8475          -0.039277
 8476          -0.037417
 8477          -0.035603
 8478          -0.033842
 8479          -0.032137
 8480          -0.030488
 8481          -0.028904
 8482          -0.027383
 8483          -0.025926
 8484          -0.024533
 8485          -0.023204
 8486          -0.021936
 8487          -0.020732
 8488          -0.019589
 8489          -0.018505
 8490          -0.017477
 8491          -0.016503
 8492          -0.015584
 8493          -0.014716
 8494          -0.013896
 8495          -0.013123
 8496          -0.012394
 8497          -0.011705
 8498          -0.011058
 8499          -0.010448
 8500          -0.009874
 8501          -0.009334
 8502          -0.008825
 8503          -0.008346
 8504          -0.008164
 8505          -0.008615
 8506          -0.009095
 8507          -0.009603
 8508          -0.010141
 8509          -0.010712
 8510          -0.011318
 8511          -0.011959
 8512          -0.012640
 8513          -0.013359
 8514          -0.014120
 8515          -0.014925
 8516          -0.015776
 8517          -0.016675
 8518          -0.017625
 8519          -0.018625
 8520          -0.019678
 8521          -0.020786
 8522          -0.021949
 8523          -0.023172
 8524          -0.024451
 8525          -0.025788
 8526          -0.027183
 8527          -0.028636
 8528          -0.030146
 8529          -0.031714
 8530          -0.033332
 8531          -0.034999
 8532          -0.036712
 8533          -0.038466
 8534          -0.040254
 8535          -0.042073
 8536          -0.043908
 8537          -0.045753
 8538          -0.047598
 8539          -0.049433
 8540          -0.051247
 8541          -0.053021
 8542          -0.054744
 8543          -0.056404
 8544          -0.057983
 8545          -0.059469
 8546          -0.060849
 8547          -0.062103
 8548          -0.063220
 8549          -0.064189
 8550          -0.064998
 8551          -0.065638
 8552          -0.066101
 8553          -0.066380
 8554          -0.066474
 8555          -0.066380
 8556          -0.066101
 8557          -0.065638
 8558          -0.064998
 8559          -0.064189
 8560          -0.063220
 8561          -0.062103
 8562          -0.060849
 8563          -0.059469
 8564          -0.057983
 8565          -0.056404
 8566          -0.054744
 8567          -0.053021
 8568          -0.051247
 8569          -0.049433
 8570          -0.047598
 8571          -0.045753
 8572          -0.043908
 8573          -0.042073
 8574          -0.040254
 8575          -0.038466
 8576          -0.036712
 8577          -0.034999
 8578          -0.033332
 8579          -0.031714
 8580          -0.030146
 8581          -0.028636
 8582          -0.027183
 8583          -0.025788
 8584          -0.024451
 8585          -0.023172
 8586          -0.021949
 8587          -0.020786
 8588          -0.019678
 8589          -0.018625
 8590          -0.017625
 8591          -0.016675
 8592          -0.015776
 8593          -0.014925
 8594          -0.014120
 8595          -0.013359
 8596          -0.012640
 8597          -0.011959
 8598          -0.011318
 8599          -0.010712
 8600          -0.010141
 8601          -0.009603
 8602          -0.009095
 8603          -0.008615
 8604          -0.008407
 8605          -0.008858
 8606          -0.009336
 8607          -0.009842
 8608          -0.010377
 8609          -0.010944
 8610          -0.011543
 8611          -0.012177
 8612          -0.012848
 8613          -0.013556
 8614          -0.014303
 8615          -0.015092
 8616          -0.015924
 8617          -0.016801
 8618          -0.017726
 8619          -0.018697
 8620          -0.019718
 8621          -0.020789
 8622          -0.021911
 8623          -0.023089
 8624          -0.024316
 8625          -0.025597
 8626          -0.026931
 8627          -0.028317
 8628          -0.029754
 8629          -0.031243
 8630          -0.032775
 8631          -0.034352
 8632          -0.035967
 8633          -0.037617
 8634          -0.039296
 8635          -0.041001
 8636          -0.042718
 8637          -0.044440
 8638          -0.046159
 8639          -0.047865
 8640          -0.049549
 8641          -0.051194
 8642          -0.052789
 8643          -0.054322
 8644          -0.055780
 8645          -0.057149
 8646          -0.058418
 8647          -0.059571
 8648          -0.060598
 8649          -0.061487
 8650          -0.062229
 8651          -0.062815
 8652          -0.063239
 8653          -0.063496
 8654          -0.063581
 8655          -0.063496
 8656          -0.063239
 8657          -0.062815
 8658          -0.062229
 8659          -0.061487
 8660          -0.060598
 8661          -0.059571
 8662          -0.058418
 8663          -0.057149
 8664          -0.055780
 8665          -0.054322
 8666          -0.052789
 8667          -0.051194
 8668          -0.049549
 8669          -0.047865
 8670          -0.046159
 8671          -0.044440
 8672          -0.042718
 8673          -0.041001
 8674          -0.039296
 8675          -0.037617
 8676          -0.035967
 8677          -0.034352
 8678          -0.032775
 8679          -0.031243
 8680          -0.029754
 8681          -0.028317
 8682          -0.026931
 8683          -0.025597
 8684          -0.024316
 8685          -0.023089
 8686          -0.021911
 8687          -0.020789
 8688          -0.019718
 8689          -0.018697
 8690          -0.017726
 8691          -0.016801
 8692          -0.015924
 8693          -0.015092
 8694          -0.014303
 8695          -0.013556
 8696          -0.012848
 8697          -0.012177
 8698          -0.011543
 8699          -0.010944
 8700          -0.010377
 8701          -0.009842
 8702          -0.009336
 8703          -0.008858
 8704          -0.008625
 8705          -0.009075
 8706          -0.009551
 8707          -0.010053
 8708          -0.010584
 8709          -0.011145
 8710          -0.011736
 8711          -0.012361
 8712          -0.013021
 8713          -0.013716
 8714          -0.014449
 8715          -0.015220
 8716          -0.016032
 8717          -0.016886
 8718          -0.017784
 8719          -0.018726
 8720          -0.019713
 8721          -0.020747
 8722          -0.021828
 8723          -0.022959
 8724          -0.024136
 8725          -0.025362
 8726          -0.026635
 8727          -0.027955
 8728          -0.029320
 8729          -0.030732
 8730          -0.032181
 8731          -0.033669
 8732          -0.035191
 8733          -0.036742
 8734          -0.038316
 8735          -0.039912
 8736          -0.041515
 8737          -0.043121
 8738          -0.044720
 8739          -0.046305
 8740          -0.047866
 8741          -0.049388
 8742          -0.050862
 8743          -0.052276
 8744          -0.053619
 8745          -0.054879
 8746          -0.056046
 8747          -0.057104
 8748          -0.058045
 8749          -0.058860
 8750          -0.059539
 8751          -0.060076
 8752          -0.060464
 8753          -0.060698
 8754          -0.060777
 8755          -0.060698
 8756          -0.060464
 8757          -0.060076
 8758          -0.059539
 8759          -0.058860
 8760          -0.058045
 8761          -0.057104
 8762          -0.056046
 8763          -0.054879
 8764          -0.053619
 8765          -0.052276
 8766          -0.050862
 8767          -0.049388
 8768          -0.047866
 8769          -0.046305
 8770          -0.044720
 8771          -0.043121
 8772          -0.041515
 8773          -0.039912
 8774          -0.038316
 8775          -0.036742
 8776          -0.035191
 8777          -0.033669
 8778          -0.032181
 8779          -0.030732
 8780          -0.029320
 8781          -0.027955
 8782          -0.026635
 8783          -0.025362
 8784          -0.024136
 8785          -0.022959
 8786          -0.021828
 8787          -0.020747
 8788          -0.019713
 8789          -0.018726
 8790          -0.017784
 8791          -0.016886
 8792          -0.016032
 8793          -0.015220
 8794          -0.014449
 8795          -0.013716
 8796          -0.013021
 8797          -0.012361
 8798          -0.011736
 8799          -0.011145
 8800          -0.010584
 8801          -0.010053
 8802          -0.009551
 8803          -0.009075
 8804          -0.008820
 8805          -0.009267
 8806          -0.009740
 8807          -0.010238
 8808          -0.010762
 8809          -0.011316
 8810          -0.011899
 8811          -0.012514
 8812          -0.013162
 8813          -0.013843
 8814          -0.014559
 8815          -0.015312
 8816          -0.016102
 8817          -0.016932
 8818          -0.017804
 8819          -0.018715
 8820          -0.019669
 8821          -0.020665
 8822          -0.021705
 8823          -0.022790
 8824          -0.023918
 8825          -0.025089
 8826          -0.026302
 8827          -0.027558
 8828          -0.028854
 8829          -0.030191
 8830          -0.031561
 8831          -0.032964
 8832          -0.034396
 8833          -0.035852
 8834          -0.037328
 8835          -0.038820
 8836          -0.040316
 8837          -0.041811
 8838          -0.043299
 8839          -0.044769
 8840          -0.046215
 8841          -0.047622
 8842          -0.048983
 8843          -0.050287
 8844          -0.051524
 8845          -0.052682
 8846          -0.053754
 8847          -0.054725
 8848          -0.055587
 8849          -0.056333
 8850          -0.056955
 8851          -0.057446
 8852          -0.057801
 8853          -0.058015
 8854          -0.058086
 8855          -0.058015
 8856          -0.057801
 8857          -0.057446
 8858          -0.056955
 8859          -0.056333
 8860          -0.055587
 8861          -0.054725
 8862          -0.053754
 8863          -0.052682
 8864          -0.051524
 8865          -0.050287
 8866          -0.048983
 8867          -0.047622
 8868          -0.046215
 8869          -0.044769
 8870          -0.043299
 8871          -0.041811
 8872          -0.040316
 8873          -0.038820
 8874          -0.037328
 8875          -0.035852
 8876          -0.034396
 8877          -0.032964
 8878          -0.031561
 8879          -0.030191
 8880          -0.028854
 8881          -0.027558
 8882          -0.026302
 8883          -0.025089
 8884          -0.023918
 8885          -0.022790
 8886          -0.021705
 8887          -0.020665
 8888          -0.019669
 8889          -0.018715
 8890          -0.017804
 8891          -0.016932
 8892          -0.016102
 8893          -0.015312
 8894          -0.014559
 8895          -0.013843
 8896          -0.013162
 8897          -0.012514
 8898          -0.011899
 8899          -0.011316
 8900          -0.010762
 8901          -0.010238
 8902          -0.009740
 8903          -0.009267
 8904          -0.008993
 8905          -0.009436
 8906          -0.009905
 8907          -0.010397
 8908          -0.010915
 8909          -0.011460
 8910          -0.012034
 8911          -0.012637
 8912          -0.013272
 8913          -0.013938
 8914          -0.014637
 8915          -0.015371
 8916          -0.016140
 8917          -0.016945
 8918          -0.017789
 8919          -0.018670
 8920          -0.019589
 8921          -0.020548
 8922          -0.021547
 8923          -0.022588
 8924          -0.023666
 8925          -0.024783
 8926          -0.025939
 8927          -0.027132
 8928          -0.028361
 8929          -0.029626
 8930          -0.030920
 8931          -0.032242
 8932          -0.033588
 8933          -0.034955
 8934          -0.036336
 8935          -0.037731
 8936          -0.039126
 8937          -0.040518
 8938          -0.041900
 8939          -0.043264
 8940          -0.044603
 8941          -0.045904
 8942          -0.047160
 8943          -0.048361
 8944          -0.049499
 8945          -0.050564
 8946          -0.051547
 8947          -0.052438
 8948          -0.053228
 8949          -0.053911
 8950          -0.054479
 8951          -0.054928
 8952          -0.055252
 8953          -0.055448
 8954          -0.055513
 8955          -0.055448
 8956          -0.055252
 8957          -0.054928
 8958          -0.054479
 8959          -0.053911
 8960          -0.053228
 8961          -0.052438
 8962          -0.051547
 8963          -0.050564
 8964          -0.049499
 8965          -0.048361
 8966          -0.047160
 8967          -0.045904
 8968          -0.044603
 8969          -0.043264
 8970          -0.041900
 8971          -0.040518
 8972          -0.039126
 8973          -0.037731
 8974          -0.036336
 8975          -0.034955
 8976          -0.033588
 8977          -0.032242
 8978          -0.030920
 8979          -0.029626
 8980          -0.028361
 8981          -0.027132
 8982          -0.025939
 8983          -0.024783
 8984          -0.023666
 8985          -0.022588
 8986          -0.021547
 8987          -0.020548
 8988          -0.019589
 8989          -0.018670
 8990          -0.017789
 8991          -0.016945
 8992          -0.016140
 8993          -0.015371
 8994          -0.014637
 8995          -0.013938
 8996          -0.013272
 8997          -0.012637
 8998          -0.012034
 8999          -0.011460
 9000          -0.010915
 9001          -0.010397
 9002          -0.009905
 9003          -0.009436
 9004          -0.009145
 9005          -0.009584
 9006          -0.010047
 9007          -0.010532
 9008          -0.011043
 9009          -0.011579
 9010          -0.012142
 9011          -0.012733
 9012          -0.013355
 9013          -0.014005
 9014          -0.014686
 9015          -0.015400
 9016          -0.016146
 9017          -0.016926
 9018          -0.017743
 9019          -0.018593
 9020          -0.019479
 9021          -0.020401
 9022          -0.021359
 9023          -0.022355
 9024          -0.023386
 9025          -0.024451
 9026          -0.025551
 9027          -0.026684
 9028          -0.027848
 9029          -0.029044
 9030          -0.030265
 9031          -0.031510
 9032          -0.032775
 9033          -0.034057
 9034          -0.035350
 9035          -0.036652
 9036          -0.037953
 9037          -0.039249
 9038          -0.040532
 9039          -0.041797
 9040          -0.043036
 9041          -0.044238
 9042          -0.045397
 9043          -0.046504
 9044          -0.047551
 9045          -0.048529
 9046          -0.049432
 9047          -0.050248
 9048          -0.050972
 9049          -0.051597
 9050          -0.052118
 9051          -0.052528
 9052          -0.052824
 9053          -0.053003
 9054          -0.053062
 9055          -0.053003
 9056          -0.052824
 9057          -0.052528
 9058          -0.052118
 9059          -0.051597
 9060          -0.050972
 9061          -0.050248
 9062          -0.049432
 9063          -0.048529
 9064          -0.047551
 9065          -0.046504
 9066          -0.045397
 9067          -0.044238
 9068          -0.043036
 9069          -0.041797
 9070          -0.040532
 9071          -0.039249
 9072          -0.037953
 9073          -0.036652
 9074          -0.035350
 9075          -0.034057
 9076          -0.032775
 9077          -0.031510
 9078          -0.030265
 9079          -0.029044
 9080          -0.027848
 9081          -0.026684
 9082          -0.025551
 9083          -0.024451
 9084          -0.023386
 9085          -0.022355
 9086          -0.021359
 9087          -0.020401
 9088          -0.019479
 9089          -0.018593
 9090          -0.017743
 9091          -0.016926
 9092          -0.016146
 9093          -0.015400
 9094          -0.014686
 9095          -0.014005
 9096          -0.013355
 9097          -0.012733
 9098          -0.012142
 9099          -0.011579
 9100          -0.011043
 9101          -0.010532
 9102          -0.010047
 9103          -0.009584
 9104          -0.009277
 9105          -0.009711
 9106          -0.010168
 9107          -0.010646
 9108          -0.011148
 9109          -0.011674
 9110          -0.012226
 9111          -0.012805
 9112          -0.013411
 9113          -0.014045
 9114          -0.014709
 9115          -0.015402
 9116          -0.016125
 9117          -0.016881
 9118          -0.017669
 9119          -0.018489
 9120          -0.019341
 9121          -0.020227
 9122          -0.021145
 9123          -0.022098
 9124          -0.023082
 9125          -0.024097
 9126          -0.025142
 9127          -0.026218
 9128          -0.027320
 9129          -0.028451
 9130          -0.029602
 9131          -0.030773
 9132          -0.031962
 9133          -0.033163
 9134          -0.034373
 9135          -0.035589
 9136          -0.036802
 9137          -0.038008
 9138          -0.039200
 9139          -0.040372
 9140          -0.041519
 9141          -0.042630
 9142          -0.043699
 9143          -0.044720
 9144          -0.045683
 9145          -0.046582
 9146          -0.047411
 9147          -0.048159
 9148          -0.048823
 9149          -0.049395
 9150          -0.049871
 9151          -0.050246
 9152          -0.050517
 9153          -0.050680
 9154          -0.050735
 9155          -0.050680
 9156          -0.050517
 9157          -0.050246
 9158          -0.049871
 9159          -0.049395
 9160          -0.048823
 9161          -0.048159
 9162          -0.047411
 9163          -0.046582
 9164          -0.045683
 9165          -0.044720
 9166          -0.043699
 9167          -0.042630
 9168          -0.041519
 9169          -0.040372
 9170          -0.039200
 9171          -0.038008
 9172          -0.036802
 9173          -0.035589
 9174          -0.034373
 9175          -0.033163
 9176          -0.031962
 9177          -0.030773
 9178          -0.029602
 9179          -0.028451
 9180          -0.027320
 9181          -0.026218
 9182          -0.025142
 9183          -0.024097
 9184          -0.023082
 9185          -0.022098
 9186          -0.021145
 9187          -0.020227
 9188          -0.019341
 9189          -0.018489
 9190          -0.017669
 9191          -0.016881
 9192          -0.016125
 9193          -0.015402
 9194          -0.014709
 9195          -0.014045
 9196          -0.013411
 9197          -0.012805
 9198          -0.012226
 9199          -0.011674
 9200          -0.011148
 9201          -0.010646
 9202          -0.010168
 9203          -0.009711
 9204          -0.009391
 9205          -0.009819
 9206          -0.010269
 9207          -0.010739
 9208          -0.011232
 9209          -0.011748
 9210          -0.012288
 9211          -0.012853
 9212          -0.013445
 9213          -0.014062
 9214          -0.014707
 9215          -0.015379
 9216          -0.016080
 9217          -0.016810
 9218          -0.017571
 9219          -0.018360
 9220          -0.019180
 9221          -0.020029
 9222          -0.020909
 9223          -0.021820
 9224          -0.022758
 9225          -0.023724
 9226          -0.024718
 9227          -0.025737
 9228          -0.026781
 9229          -0.027849
 9230          -0.028934
 9231          -0.030037
 9232          -0.031153
 9233          -0.032279
 9234          -0.033411
 9235          -0.034547
 9236          -0.035677
 9237          -0.036798
 9238          -0.037906
 9239          -0.038992
 9240          -0.040055
 9241          -0.041082
 9242          -0.042069
 9243          -0.043009
 9244          -0.043896
 9245          -0.044722
 9246          -0.045484
 9247          -0.046170
 9248          -0.046778
 9249          -0.047303
 9250          -0.047738
 9251          -0.048081
 9252          -0.048329
 9253          -0.048478
 9254          -0.048528
 9255          -0.048478
 9256          -0.048329
 9257          -0.048081
 9258          -0.047738
 9259          -0.047303
 9260          -0.046778
 9261          -0.046170
 9262          -0.045484
 9263          -0.044722
 9264          -0.043896
 9265          -0.043009
 9266          -0.042069
 9267          -0.041082
 9268          -0.040055
 9269          -0.038992
 9270          -0.037906
 9271          -0.036798
 9272          -0.035677
 9273          -0.034547
 9274          -0.033411
 9275          -0.032279
 9276          -0.031153
 9277          -0.030037
 9278          -0.028934
 9279          -0.027849
 9280          -0.026781
 9281          -0.025737
 9282          -0.024718
 9283          -0.023724
 9284          -0.022758
 9285          -0.021820
 9286          -0.020909
 9287          -0.020029
 9288          -0.019180
 9289          -0.018360
 9290          -0.017571
 9291          -0.016810
 9292          -0.016080
 9293          -0.015379
 9294          -0.014707
 9295          -0.014062
 9296          -0.013445
 9297          -0.012853
 9298          -0.012288
 9299          -0.011748
 9300          -0.011232
 9301          -0.010739
 9302          -0.010269
 9303          -0.009819
 9304          -0.009488
 9305          -0.009910
 9306          -0.010352
 9307          -0.010813
 9308          -0.011296
 9309          -0.011801
 9310          -0.012329
 9311          -0.012880
 9312          -0.013457
 9313          -0.014057
 9314          -0.014682
 9315          -0.015334
 9316          -0.016012
 9317          -0.016717
 9318          -0.017451
 9319          -0.018210
 9320          -0.018997
 9321          -0.019811
 9322          -0.020653
 9323          -0.021523
 9324          -0.022417
 9325          -0.023337
 9326          -0.024280
 9327          -0.025247
 9328          -0.026234
 9329          -0.027242
 9330          -0.028265
 9331          -0.029302
 9332          -0.030350
 9333          -0.031405
 9334          -0.032464
 9335          -0.033525
 9336          -0.034578
 9337          -0.035622
 9338          -0.036650
 9339          -0.037658
 9340          -0.038642
 9341          -0.039591
 9342          -0.040503
 9343          -0.041370
 9344          -0.042186
 9345          -0.042946
 9346          -0.043646
 9347          -0.044276
 9348          -0.044834
 9349          -0.045315
 9350          -0.045714
 9351          -0.046028
 9352          -0.046254
 9353          -0.046391
 9354          -0.046437
 9355          -0.046391
 9356          -0.046254
 9357          -0.046028
 9358          -0.045714
 9359          -0.045315
 9360          -0.044834
 9361          -0.044276
 9362          -0.043646
 9363          -0.042946
 9364          -0.042186
 9365          -0.041370
 9366          -0.040503
 9367          -0.039591
 9368          -0.038642
 9369          -0.037658
 9370          -0.036650
 9371          -0.035622
 9372          -0.034578
 9373          -0.033525
 9374          -0.032464
 9375          -0.031405
 9376          -0.030350
 9377          -0.029302
 9378          -0.028265
 9379          -0.027242
 9380          -0.026234
 9381          -0.025247
 9382          -0.024280
 9383          -0.023337
 9384          -0.022417
 9385          -0.021523
 9386          -0.020653
 9387          -0.019811
 9388          -0.018997
 9389          -0.018210
 9390          -0.017451
 9391          -0.016717
 9392          -0.016012
 9393          -0.015334
 9394          -0.014682
 9395          -0.014057
 9396          -0.013457
 9397          -0.012880
 9398          -0.012329
 9399          -0.011801
 9400          -0.011296
 9401          -0.010813
 9402          -0.010352
 9403          -0.009910
 9404          -0.009568
 9405          -0.009983
 9406          -0.010417
 9407          -0.010870
 9408          -0.011342
 9409          -0.011836
 9410          -0.012351
 9411          -0.012888
 9412          -0.013449
 9413          -0.014032
 9414          -0.014638
 9415          -0.015269
 9416          -0.015925
 9417          -0.016605
 9418          -0.017311
 9419          -0.018041
 9420          -0.018797
 9421          -0.019577
 9422          -0.020382
 9423          -0.021212
 9424          -0.022064
 9425          -0.022939
 9426          -0.023835
 9427          -0.024750
 9428          -0.025684
 9429          -0.026636
 9430          -0.027600
 9431          -0.028575
 9432          -0.029559
 9433          -0.030548
 9434          -0.031539
 9435          -0.032530
 9436          -0.033512
 9437          -0.034484
 9438          -0.035439
 9439          -0.036375
 9440          -0.037286
 9441          -0.038165
 9442          -0.039007
 9443          -0.039807
 9444          -0.040559
 9445          -0.041259
 9446          -0.041903
 9447          -0.042482
 9448          -0.042994
 9449          -0.043435
 9450          -0.043801
 9451          -0.044088
 9452          -0.044296
 9453          -0.044421
 9454          -0.044463
 9455          -0.044421
 9456          -0.044296
 9457          -0.044088
 9458          -0.043801
 9459          -0.043435
 9460          -0.042994
 9461          -0.042482
 9462          -0.041903
 9463          -0.041259
 9464          -0.040559
 9465          -0.039807
 9466          -0.039007
 9467          -0.038165
 9468          -0.037286
 9469          -0.036375
 9470          -0.035439
 9471          -0.034484
 9472          -0.033512
 9473          -0.032530
 9474          -0.031539
 9475          -0.030548
 9476          -0.029559
 9477          -0.028575
 9478          -0.027600
 9479          -0.026636
 9480          -0.025684
 9481          -0.024750
 9482          -0.023835
 9483          -0.022939
 9484          -0.022064
 9485          -0.021212
 9486          -0.020382
 9487          -0.019577
 9488          -0.018797
 9489          -0.018041
 9490          -0.017311
 9491          -0.016605
 9492          -0.015925
 9493          -0.015269
 9494          -0.014638
 9495          -0.014032
 9496          -0.013449
 9497          -0.012888
 9498          -0.012351
 9499          -0.011836
 9500          -0.011342
 9501          -0.010870
 9502          -0.010417
 9503          -0.009983
 9504          -0.009633
 9505          -0.010040
 9506          -0.010466
 9507          -0.010909
 9508          -0.011372
 9509          -0.011853
 9510          -0.012356
 9511          -0.012878
 9512          -0.013423
 9513          -0.013989
 9514          -0.014577
 9515          -0.015187
 9516          -0.015819
 9517          -0.016475
 9518          -0.017155
 9519          -0.017857
 9520          -0.018581
 9521          -0.019328
 9522          -0.020097
 9523          -0.020890
 9524          -0.021701
 9525          -0.022533
 9526          -0.023383
 9527          -0.024250
 9528          -0.025133
 9529          -0.026032
 9530          -0.026940
 9531          -0.027858
 9532          -0.028782
 9533          -0.029709
 9534          -0.030636
 9535          -0.031562
 9536          -0.032478
 9537          -0.033383
 9538          -0.034272
 9539          -0.035140
 9540          -0.035985
 9541          -0.036798
 9542          -0.037577
 9543          -0.038316
 9544          -0.039010
 9545          -0.039655
 9546          -0.040247
 9547          -0.040780
 9548          -0.041250
 9549          -0.041655
 9550          -0.041991
 9551          -0.042255
 9552          -0.042445
 9553          -0.042560
 9554          -0.042598
 9555          -0.042560
 9556          -0.042445
 9557          -0.042255
 9558          -0.041991
 9559          -0.041655
 9560          -0.041250
 9561          -0.040780
 9562          -0.040247
 9563          -0.039655
 9564          -0.039010
 9565          -0.038316
 9566          -0.037577
 9567          -0.036798
 9568          -0.035985
 9569          -0.035140
 9570          -0.034272
 9571          -0.033383
 9572          -0.032478
 9573          -0.031562
 9574          -0.030636
 9575          -0.029709
 9576          -0.028782
 9577          -0.027858
 9578          -0.026940
 9579          -0.026032
 9580          -0.025133
 9581          -0.024250
 9582          -0.023383
 9583          -0.022533
 9584          -0.021701
 9585          -0.020890
 9586          -0.020097
 9587          -0.019328
 9588          -0.018581
 9589          -0.017857
 9590          -0.017155
 9591          -0.016475
 9592          -0.015819
 9593          -0.015187
 9594          -0.014577
 9595          -0.013989
 9596          -0.013423
 9597          -0.012878
 9598          -0.012356
 9599          -0.011853
 9600          -0.011372
 9601          -0.010909
 9602          -0.010466
 9603          -0.010040
 9604          -0.009684
 9605          -0.010083
 9606          -0.010500
 9607          -0.010934
 9608          -0.011385
 9609          -0.011855
 9610          -0.012344
 9611          -0.012852
 9612          -0.013381
 9613          -0.013930
 9614          -0.014499
 9615          -0.015088
 9616          -0.015699
 9617          -0.016330
 9618          -0.016984
 9619          -0.017658
 9620          -0.018352
 9621          -0.019067
 9622          -0.019802
 9623          -0.020558
 9624          -0.021330
 9625          -0.022121
 9626          -0.022927
 9627          -0.023749
 9628          -0.024584
 9629          -0.025432
 9630          -0.026288
 9631          -0.027152
 9632          -0.028019
 9633          -0.028889
 9634          -0.029757
 9635          -0.030622
 9636          -0.031477
 9637          -0.032320
 9638          -0.033147
 9639          -0.033953
 9640          -0.034737
 9641          -0.035491
 9642          -0.036211
 9643          -0.036894
 9644          -0.037535
 9645          -0.038130
 9646          -0.038675
 9647          -0.039166
 9648          -0.039599
 9649          -0.039971
 9650          -0.040279
 9651          -0.040522
 9652          -0.040697
 9653          -0.040802
 9654          -0.040837
 9655          -0.040802
 9656          -0.040697
 9657          -0.040522
 9658          -0.040279
 9659          -0.039971
 9660          -0.039599
 9661          -0.039166
 9662          -0.038675
 9663          -0.038130
 9664          -0.037535
 9665          -0.036894
 9666          -0.036211
 9667          -0.035491
 9668          -0.034737
 9669          -0.033953
 9670          -0.033147
 9671          -0.032320
 9672          -0.031477
 9673          -0.030622
 9674          -0.029757
 9675          -0.028889
 9676          -0.028019
 9677          -0.027152
 9678          -0.026288
 9679          -0.025432
 9680          -0.024584
 9681          -0.023749
 9682          -0.022927
 9683          -0.022121
 9684          -0.021330
 9685          -0.020558
 9686          -0.019802
 9687          -0.019067
 9688          -0.018352
 9689          -0.017658
 9690          -0.016984
 9691          -0.016330
 9692          -0.015699
 9693          -0.015088
 9694          -0.014499
 9695          -0.013930
 9696          -0.013381
 9697          -0.012852
 9698          -0.012344
 9699          -0.011855
 9700          -0.011385
 9701          -0.010934
 9702          -0.010500
 9703          -0.010083
 9704          -0.009721
 9705          -0.010112
 9706          -0.010520
 9707          -0.010944
 9708          -0.011384
 9709          -0.011842
 9710          -0.012318
 9711          -0.012812
 9712          -0.013325
 9713          -0.013856
 9714          -0.014406
 9715          -0.014975
 9716          -0.015564
 9717          -0.016172
 9718          -0.016801
 9719          -0.017447
 9720          -0.018112
 9721          -0.018796
 9722          -0.019498
 9723          -0.020218
 9724          -0.020954
 9725          -0.021705
 9726          -0.022470
 9727          -0.023249
 9728          -0.024038
 9729          -0.024839
 9730          -0.025646
 9731          -0.026458
 9732          -0.027273
 9733          -0.028089
 9734          -0.028902
 9735          -0.029710
 9736          -0.030508
 9737          -0.031294
 9738          -0.032064
 9739          -0.032814
 9740          -0.033542
 9741          -0.034240
 9742          -0.034908
 9743          -0.035539
 9744          -0.036132
 9745          -0.036681
 9746          -0.037184
 9747          -0.037636
 9748          -0.038035
 9749          -0.038377
 9750          -0.038661
 9751          -0.038884
 9752          -0.039045
 9753          -0.039141
 9754          -0.039174
 9755          -0.039141
 9756          -0.039045
 9757          -0.038884
 9758          -0.038661
 9759          -0.038377
 9760          -0.038035
 9761          -0.037636
 9762          -0.037184
 9763          -0.036681
 9764          -0.036132
 9765          -0.035539
 9766          -0.034908
 9767          -0.034240
 9768          -0.033542
 9769          -0.032814
 9770          -0.032064
 9771          -0.031294
 9772          -0.030508
 9773          -0.029710
 9774          -0.028902
 9775          -0.028089
 9776          -0.027273
 9777          -0.026458
 9778          -0.025646
 9779          -0.024839
 9780          -0.024038
 9781          -0.023249
 9782          -0.022470
 9783          -0.021705
 9784          -0.020954
 9785          -0.020218
 9786          -0.019498
 9787          -0.018796
 9788          -0.018112
 9789          -0.017447
 9790          -0.016801
 9791          -0.016172
 9792          -0.015564
 9793          -0.014975
 9794          -0.014406
 9795          -0.013856
 9796          -0.013325
 9797          -0.012812
 9798          -0.012318
 9799          -0.011842
 9800          -0.011384
 9801          -0.010944
 9802          -0.010520
 9803          -0.010112
 9804          -0.009746
 9805          -0.010129
 9806          -0.010528
 9807          -0.010941
 9808          -0.011370
 9809          -0.011816
 9810          -0.012278
 9811          -0.012758
 9812          -0.013255
 9813          -0.013769
 9814          -0.014301
 9815          -0.014850
 9816          -0.015417
 9817          -0.016002
 9818          -0.016606
 9819          -0.017226
 9820          -0.017863
 9821          -0.018517
 9822          -0.019187
 9823          -0.019874
 9824          -0.020574
 9825          -0.021288
 9826          -0.022014
 9827          -0.022751
 9828          -0.023497
 9829          -0.024253
 9830          -0.025014
 9831          -0.025778
 9832          -0.026544
 9833          -0.027309
 9834          -0.028071
 9835          -0.028828
 9836          -0.029573
 9837          -0.030306
 9838          -0.031023
 9839          -0.031721
 9840          -0.032397
 9841          -0.033046
 9842          -0.033664
 9843          -0.034249
 9844          -0.034797
 9845          -0.035304
 9846          -0.035769
 9847          -0.036186
 9848          -0.036553
 9849          -0.036869
 9850          -0.037130
 9851          -0.037336
 9852          -0.037484
 9853          -0.037573
 9854          -0.037602
 9855          -0.037573
 9856          -0.037484
 9857          -0.037336
 9858          -0.037130
 9859          -0.036869
 9860          -0.036553
 9861          -0.036186
 9862          -0.035769
 9863          -0.035304
 9864          -0.034797
 9865          -0.034249
 9866          -0.033664
 9867          -0.033046
 9868          -0.032397
 9869          -0.031721
 9870          -0.031023
 9871          -0.030306
 9872          -0.029573
 9873          -0.028828
 9874          -0.028071
 9875          -0.027309
 9876          -0.026544
 9877          -0.025778
 9878          -0.025014
 9879          -0.024253
 9880          -0.023497
 9881          -0.022751
 9882          -0.022014
 9883          -0.021288
 9884          -0.020574
 9885          -0.019874
 9886          -0.019187
 9887          -0.018517
 9888          -0.017863
 9889          -0.017226
 9890          -0.016606
 9891          -0.016002
 9892          -0.015417
 9893          -0.014850
 9894          -0.014301
 9895          -0.013769
 9896          -0.013255
 9897          -0.012758
 9898          -0.012278
 9899          -0.011816
 9900          -0.011370
 9901          -0.010941
 9902          -0.010528
 9903          -0.010129
 9904          -0.009760
 9905          -0.010134
 9906          -0.010523
 9907          -0.010926
 9908          -0.011344
 9909          -0.011778
 9910          -0.012227
 9911          -0.012692
 9912          -0.013173
 9913          -0.013670
 9914          -0.014184
 9915          -0.014714
 9916          -0.015260
 9917          -0.015822
 9918          -0.016402
 9919          -0.016997
 9920          -0.017607
 9921          -0.018232
 9922          -0.018871
 9923          -0.019525
 9924          -0.020191
 9925          -0.020869
 9926          -0.021558
 9927          -0.022256
 9928          -0.022962
 9929          -0.023676
 9930          -0.024393
 9931          -0.025112
 9932          -0.025832
 9933          -0.026550
 9934          -0.027264
 9935          -0.027972
 9936          -0.028669
 9937          -0.029353
 9938          -0.030021
 9939          -0.030671
 9940          -0.031299
 9941          -0.031902
 9942          -0.032475
 9943          -0.033017
 9944          -0.033525
 9945          -0.033994
 9946          -0.034423
 9947          -0.034809
 9948          -0.035148
 9949          -0.035439
 9950          -0.035680
 9951          -0.035869
 9952          -0.036006
 9953          -0.036088
 9954          -0.036115
 9955          -0.036088
 9956          -0.036006
 9957          -0.035869
 9958          -0.035680
 9959          -0.035439
 9960          -0.035148
 9961          -0.034809
 9962          -0.034423
 9963          -0.033994
 9964          -0.033525
 9965          -0.033017
 9966          -0.032475
 9967          -0.031902
 9968          -0.031299
 9969          -0.030671
 9970          -0.030021
 9971          -0.029353
 9972          -0.028669
 9973          -0.027972
 9974          -0.027264
 9975          -0.026550
 9976          -0.025832
 9977          -0.025112
 9978          -0.024393
 9979          -0.023676
 9980          -0.022962
 9981          -0.022256
 9982          -0.021558
 9983          -0.020869
 9984          -0.020191
 9985          -0.019525
 9986          -0.018871
 9987          -0.018232
 9988          -0.017607
 9989          -0.016997
 9990          -0.016402
 9991          -0.015822
 9992          -0.015260
 9993          -0.014714
 9994          -0.014184
 9995          -0.013670
 9996          -0.013173
 9997          -0.012692
 9998          -0.012227
 9999          -0.011778
 ****          -0.011344
 ****          -0.010926
 ****          -0.010523
 ****          -0.010134
 -----------------------------------------------------------------
 1\1\GINC-DESKTOP-2Q9IDMK\SP\RMP2-FC\6-311++G(d,p)\H2O1\BRUNO\30-May-20
 19\0\\#mp2/6-311++G(d,p) density=current prop=(read,potential)\\water 
 pruebas\\0,1\O,0,0.,0.110861,0.\H,0,0.797899,-0.443444,0.\H,0,-0.79789
 9,-0.443443,0.\\Version=AM64L-G09RevA.02\State=1-A1\HF=-76.0508475\MP2
 =-76.2733925\RMSD=3.345e-09\Dipole=-0.0000008,-0.8247576,0.\Quadrupole
 =1.5495436,-0.3432588,-1.2062848,-0.0000019,0.,0.\PG=C02V [C2(O1),SGV(
 H2)]\\@


 THOSE WITH THE GOLD MAKE THE RULES.

             -- PETER'S GOLDEN RULE
 Job cpu time:  0 days  0 hours  0 minutes  1.8 seconds.
 File lengths (MBytes):  RWF=      5 Int=      0 D2E=      0 Chk=      1 Scr=      1
 Normal termination of Gaussian 09 at Thu May 30 15:11:10 2019.
