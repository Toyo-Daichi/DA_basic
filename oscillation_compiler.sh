#!/bin/bash
# attach by oscillation.f90

rm -rf *.cnf *.mod maintools
 
# initial setting
DA_METHOD='KF   ' #or 'EnKF ' or 'Ajoint'

x_tinit=5.0d0; v_tinit=0.0d0
x_sinit=4.0d0; v_sinit=1.0d0

Pf_init=( 1.0d0 0.0d0 0.0d0 1.0d0 )
R_init=( 0.1d0 )
Kg_init=( 0.0d0 )
H_init=( 1.0d0 0.0d0 )

#----------------------------------------------------------------------
# +++ Set conf files
#----------------------------------------------------------------------


# Compiler
gfortran -fbounds-check -o maintools oscillation.f90
./maintools << EOF
  &initial_osc
    x_tinit = ${x_tinit},　v_tinit = ${v_tinit},
    x_sinit = ${x_sinit},  v_sinit = ${v_sinit},
  &initial_que
    Pf_init = ${Pf_init[0]}, ${Pf_init[1]}, ${Pf_init[2]}, ${Pf_init[3]}, 
    R_init  = ${R_init[0]}, Kg_init = ${Kg_init[0]}, H_init = ${H_init[0]}, ${H_init[1]}
  &conv_param
    outputfile  = '${outputname}',
    opt_veach   = .${boolen}.,
  &da_setting
    da_method = ${DA_METHOD}
  &END
EOF

rm -rf *.cnf *.mod maintools
echo 'Normal END'
exit