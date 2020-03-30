#!/bin/bash
# attach by oscillation.f90

# set -ex
prg=maintools
rm -rf *.cnf *.mod ${prg}

#----------------------------------------------------------------------
# +++ Set intial setting
#----------------------------------------------------------------------
nt_asm=400 
nt_prd=400
obs_interval=40
DA_METHOD='KF' #or 'EnKF' or 'Ajoint'
mem=500

# +++ initial value
x_tinit=5.0d0; v_tinit=0.0d0
x_sinit=4.0d0; v_sinit=1.0d0

# +++ initial matrix
# Forecast error covariance matrix
Pf_init=( 1.0d0 0.0d0 
          0.0d0 1.0d0 ) 
# Observation error covariance matrix
R_init=( 0.1d0 )
# Kalman gain matrix
K_init=( 0.0d0 0.0d0 )
# Observation　operator
H_init=( 1.0d0 0.0d0 )

# +++ output info
boolen='.true.' # write putput
outputname=oscillation_${DA_METHOD}.csv
outputfile=./output/${outputname}

#----------------------------------------------------------------------
# +++ Compiler & Run exp.
#----------------------------------------------------------------------
gfortran -fbounds-check -o ${prg} kinddef.f90 oscillation.f90
./${prg} << EOF
  &set_parm
    nt_asm = ${nt_asm}, nt_prd = ${nt_prd},
    obs_interval = ${obs_interval},
  &da_setting
    da_method= ${DA_METHOD},
  &EnKF
    mem = ${mem},
  &initial_osc
    x_tinit = ${x_tinit},　v_tinit = ${v_tinit},
    x_sinit = ${x_sinit},  v_sinit = ${v_sinit},
  &initial_que
    Pf_init = ${Pf_init[0]}, ${Pf_init[1]}, ${Pf_init[2]}, ${Pf_init[3]}, 
    R_init  = ${R_init[0]}, Kg_init = ${Kg_init[0]}, H_init = ${H_init[0]}, ${H_init[1]}
  &output
    outputfile  = '${outputfile}',
    opt_veach   = .${boolen}.
  &END
EOF

rm -rf *.mod ${prg}
echo 'Normal END'
exit