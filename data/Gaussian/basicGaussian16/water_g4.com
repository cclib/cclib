%Chk="water.chk"
%Rwf="water.rwf"
%NProcShared=4
%Mem=10GB
#p   SCF=(DIIS) Symmetry=(Tight) G4

water_Gaussian4

0, 1
O          -0.00000        -0.11916         0.00000
H          -0.79065         0.47664         0.00000
H           0.79065         0.47664        -0.00000

