%chk=water_cisd
#P CIS(D,50-50,NStates=5)/STO-3G

Water

0 1
O
H 1 R1
H 1 R1 2 A1

R1=0.99
A1=106.0
