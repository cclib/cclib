***,h2o                   !A title
r=1.85,theta=104          !set geometry parameters
geometry={nosym,
          O;              !z-matrix geometry input
          H1,O,r;
          H2,O,r,H1,theta}      
hf                        !closed-shell scf
MP4                       !fourth-order Moeller-Plesset perturbation theory
