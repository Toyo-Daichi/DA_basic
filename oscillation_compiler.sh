#!/bin/bash
# attach by oscillation.f90

rm -rf *.cnf *.mod maintools

if [ -e initial_Osc.cnf ]; then rm -f initial_Osc.cnf ; fi
cat >> initial_Osc.cnf << EOF
  &initial_osc
    x_t(0) = ,
    v_t(0) = ,
    x_s(0) = ,
    v_s(0) = ,
    
    Pf(1,1) = ,
    Pf(1,2) = ,
    Pf(2,1) = ,
    Pf(2,2) = , 
    R(1,1)  = , 
    Kg(1:2, 1:1) = ,
    H(1,1), H(1,2)
 /
EOF

if [ -e conv.cnf ]; then rm -f conv.cnf ; fi
cat >> conv.cnf << EOF
 &conv_param
    outputfile  ='${outputname}',
    opt_veach   = .${boolen}.,
    xgrid_file  ='${refdir}/xgrd.data',
    ygrid_file  ='${refdir}/ygrd.data',
    zgrid_file  ='${refdir}/zgrd.data',
 /
EOF

# Compiler
gfortran -fbounds-check -o maintools
./maintools

echo 'Normal END'
exit