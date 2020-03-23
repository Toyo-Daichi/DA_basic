#!/bin/bash
# attach by oscillation.f90

rm -rf *.cnf *.mod maintools
 
# initial setting
x_tinit=5.0d0; v_tinit=0.0d0
x_sinit=4.0d0; v_sinit=1.0d0

Pf_init=( 1.0d0 0.0d0 0.0d0 1.0d0 )
R_init=( 0.1d0 )
Kg_init=( 0.0d0 )
H_init=( 1.0d0 0.0d0 )

#----------------------------------------------------------------------
# +++ Set conf files
#----------------------------------------------------------------------

if [ -e initial_Osc.cnf ]; then rm -f initial_Osc.cnf ; fi
 cat >> initial_Osc.cnf << EOF
   &initial_osc
     x_tinit = ${x_tinit},
     v_tinit = ${v_tinit},
     x_sinit = ${x_sinit},
     v_sinit = ${v_sinit},
     
     Pf_init = ${Pf_init[0]}, ${Pf_init[1]}, ${Pf_init[2]}, ${Pf_init[3]}, 
     R_init  = ${R_init[0]}, 
     Kg_init  = ${Kg_init[0]},
     H_init  = ${H_init[0]}, ${H_init[1]}
  /
EOF
 
if [ -e conv.cnf ]; then rm -f conv.cnf ; fi
 cat > conv.cnf << EOF
  &conv_param
     outputfile  = '${outputname}',
     opt_veach   = .${boolen}.,
  /
EOF

# Compiler
gfortran -fbounds-check -o maintools oscillation.f90
./maintools

rm -rf *.cnf *.mod maintools
echo 'Normal END'
exit