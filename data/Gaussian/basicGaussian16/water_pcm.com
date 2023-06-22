%Chk="water_ccsd.chk"
%Rwf="water_ccsd.rwf"
%NProcShared=2
%Mem=10GB
#p SP HF/STO-3G SCF=(DIIS) Symmetry=(Tight) SCRF=(PCM, Solvent=Toluene) Population=(Regular) Density=(Current)

water_Single_Point_HF_Toluene_Pople_Basis_Sets_STO_3G

0, 1
O          -0.00000         0.00000         0.11916
H          -0.00000         0.79065        -0.47664
H          -0.00000        -0.79065        -0.47664

