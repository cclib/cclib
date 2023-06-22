%Chk="water.chk"
%Rwf="water.rwf"
%NProcShared=2
%Mem=4GB
#p SP PBE1PBE/STO-3G SCF=(DIIS) Symmetry=(Tight) SCRF=(SMD, Solvent=Toluene) EmpiricalDispersion=(GD3BJ) Integral=(Ultrafine) Population=(Regular) Density=(Current)

water_Single_Point_Singlet_PBE0__GD3BJ__Toluene_Pople_Basis_Sets_STO_3G

0, 1
O          -0.00000        -0.11916         0.00000
H          -0.79065         0.47664         0.00000
H           0.79065         0.47664        -0.00000

