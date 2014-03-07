%chk=mo4ocl4
%NProcShared=2
%Mem=500MB
#p gfinput iop(6/7=3) iop(3/33=1)
#B3LYP/Gen pseudo=read
#SP

Mo4OCl4 restricted single point

-2 1
Mo     0.000325    -0.000325     0.051810
 O     0.000206    -0.000206     1.741504
Cl     1.711727    -1.711727    -0.659846
Cl    -1.710992    -1.711738    -0.660079
Cl    -1.711004     1.711004    -0.660311
Cl     1.711738     1.710992    -0.660079

Mo Cl 0
LanL2MB
****
O 0
6-31G(d)
****

Mo Cl 0
LanL2

