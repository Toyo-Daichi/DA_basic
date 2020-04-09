#!/bin/bash
# attached by lorenz63.f90

# set -ex
prg=lorenz_maintools
rm -rf *.cnf *.mod ${prg}

#----------------------------------------------------------------------
# +++ Set intial setting
#----------------------------------------------------------------------
nt_asm=400 
nt_prd=400
obs_interval=40
DA_METHOD='KF' #or 'EnKF' or 'Ajoint'
mem=5000

# +++ initial value
x_tinit=5.0d0; y_tinit=0.0d0; z_tinit=0.0d0
x_sinit=4.0d0; y_sinit=1.0d0; z_sinit=0.0d0

# +++ initial matrix
# Forecast error covariance matrix
# Background error covariance matrix (= same matrix in this case)
Pf_init=( 1.0d0 0.0d0 0.0d0
          0.0d0 1.0d0 0.0d0
          0.0d0 0.0d0 1.0d0 )

# Observation error covariance matrix
R_init=( 0.1d0 0.0d0 
         0.0d0 0.2d0 )

# Kalman gain matrix
K_init=( 0.0d0 0.0d0 
         0.0d0 0.0d0
         0.0d0 0.0d0 )

# Observation　operator
H_init=( 1.0d0 0.0d0 0.0d0 
         0.0d0 1.0d0 0.0d0 )

# +++ output info
boolen='true' # write putput
outputname='lorenz63_'${DA_METHOD}'.csv'
outputfile='./output/'${outputname}

#----------------------------------------------------------------------
# +++ Run exp.
#----------------------------------------------------------------------
gfortran -fbounds-check -o ${prg} kinddef.f90 lorenz63.f90

./${prg} > ./log.${prg} << EOF
  &set_parm
    nt_asm = ${nt_asm}, nt_prd = ${nt_prd},
    obs_interval = ${obs_interval}
  /
  &da_setting
    da_method = '${DA_METHOD}'
  /
  &ensemble_size
    mems = ${mem}
  /
  &initial_osc
    x_tinit = ${x_tinit},  v_tinit = ${v_tinit},
    x_sinit = ${x_sinit},  v_sinit = ${v_sinit},
  /
  &initial_matrix
    Pf_init = ${Pf_init[0]}, ${Pf_init[1]}, ${Pf_init[2]},
              ${Pf_init[3]}, ${Pf_init[4]}, ${Pf_init[5]},
              ${Pf_init[6]}, ${Pf_init[7]}, ${Pf_init[8]},

    B_init  = ${Pf_init[0]}, ${Pf_init[1]}, ${Pf_init[2]},
              ${Pf_init[3]}, ${Pf_init[4]}, ${Pf_init[5]},
              ${Pf_init[6]}, ${Pf_init[7]}, ${Pf_init[8]},
 
    R_init  = ${R_init[0]}, ${R_init[1]},
              ${R_init[2]}, ${R_init[3]},

    Kg_init = ${Kg_init[0]}, ${Kg_init[1]},
              ${Kg_init[2]}, ${Kg_init[3]},
              ${Kg_init[4]}, ${Kg_init[5]},

    H_init = ${H_init[0]}, ${H_init[1]},
             ${H_init[2]}, ${H_init[3]},
             ${H_init[4]}, ${H_init[5]},
  /
  &output
    output_file  = '${outputfile}',
    !opt_veach   = '.${boolen}.'
  /
EOF

rm -rf *.mod ${prg}
echo 'Normal END'
exit
