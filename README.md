# eic-analysis
Scripts for DR-EIC analysis

## Analysis
First, Making csv file with equlization constant.
Edit csv file name in analysis.cc.
Move analysis.cc <path-to-dualreadout>/build/analysis

     cp analysis.cc <path-to-dual-readout>/anaylsis
     cd <path-to-dual-readout>/build
     make
     cd analysis
     ./analysis <root-filename> <integer : low> <integer : high>

From <root-filename>_Ecs.png, get "Mean" and "sigma" from gaussian fit, caculating scale factor, correcting equalization constant.
