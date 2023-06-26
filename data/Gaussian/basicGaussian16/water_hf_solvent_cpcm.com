%NProcShared=1
%Mem=500MB
#p SP HF/STO-3G SCF=(DIIS) Symmetry=(Tight) SCRF=(CPCM, Solvent=Toluene)

water_Single_Point_HF_CPCM_Pople_Basis_Sets_STO_3G

0, 1
O          -0.00000        -0.11916         0.00000
H          -0.79065         0.47664         0.00000
H           0.79065         0.47664        -0.00000

