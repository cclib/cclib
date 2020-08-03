***, divinylbenzene geom opt
gprint,basis,orbital
gthresh,thrgrad=0.0002
basis=sto-3g
symmetry,nosym
geometry={angstrom;
 Q1
 c   1 cxx2    
 c   1 cxx3       2   60.000  
 c   1 cxx4       3   60.000     2 180.000
 c   1 cxx2       4   60.000     3 180.000
 c   1 cxx3       5   60.000     4 180.000
 c   1 cxx4       6   60.000     5 180.000
 c   2 cc8        3 ccc8         4 180.000
 c   5 cc8        6 ccc8         7 180.000
 h   7 hc10       2 hcc10        1 180.000
 h   4 hc10       5 hcc10        1 180.000
 h   3 hc12       2 hcc12        1 180.000
 h   6 hc12       5 hcc12        1 180.000
 c   8 cc14       2 ccc14        3 180.000
 c   9 cc14       5 ccc14        6 180.000
 h   8 hc16       2 hcc16        3   0.000
 h   9 hc16       5 hcc16        6   0.000
 h  15 hc18       9 hcc18        5 180.000
 h  14 hc18       8 hcc18        2 180.000
 h  15 hc20       9 hcc20        5   0.000
 h  14 hc20       8 hcc20        2   0.000
}
cxx2=       1.430000
cxx3=       1.403000
cxx4=       1.409000
cc8=        1.450000
ccc8=       120.000
hc10=       1.097000
hcc10=      120.000
hc12=       1.098720
hcc12=      120.000
cc14=       1.342300
ccc14=      120.000
hc16=       1.100000
hcc16=      120.000
hc18=       1.096600
hcc18=      109.471
hc20=       1.095900
hcc20=      126.000

optg,procedure=rundft
rundft={ks,b3lyp}
ks,b3lyp;
orbprint,9999
pop
