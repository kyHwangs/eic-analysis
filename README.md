# eic-analysis
Scripts for DR-EIC analysis

## Analysis
1. Making calib.csv file
2. Copy analysis.cc to <path-to-dual-readout>/anaylsis and "make" at build

    cp analysis.cc <path-to-dual-readout>/anaylsis
    cd <path-to-dual-readout>/build
    make
  
3. Run analysis code to analysis certain root file.
  
    cd analysis
    ./analysis <filename> <integer : low> <integer : high>
    

## Electron simulation for EM resolution estimate
For EM energy resolution estimation, 5 GeV ~ 110 GeV electron beams are used. (5, 10, 20, 30, 50, 70, 90, 110 GeV)
Beam spot is alligned to the center of 70th tower and has angle of (theta, phi) = (1.5 deg, 1 deg).
