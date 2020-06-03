#!/bin/bash
# attached by lorenz63.f90

set -ex
CDIR=`pwd`
prg=lorenz_maintools
today=$(date "+%Y%m%d%H%M")
echo ${today}
rm -rf *.mod ${prg}

#----------------------------------------------------------------------
# +++ Set intial setting
#----------------------------------------------------------------------
nt_asm=2500
nt_prd=2500
obs_interval=1
da_method='EnKF' #'KF' or 'EnKF'
intg_method='Runge-Kutta' #'Euler' or 'Runge-Kutta'
mem=1
enkf_method='none'
if [ ${da_method} = 'EnKF' ]; then 
  mem=50
  enkf_method='SRF' # 'PO' or "SRF"
fi

# +++ initial value
x_tinit=1.50880d0; y_tinit=-1.531271d0; z_tinit=25.46091d0
x_sinit=7.71000d0; y_sinit=-7.721271d0; z_sinit=20.46091d0

# +++ initial matrix
# forecast var & obs size
nx=3; ny=3

# Adaptive inflation mode
alpha=0.0d0

# +++ outqput info
boolen='true' # output logical format

# for KF
outputname=${da_method}'.csv'
outputfile_score='./output/lorenz63/'${outputname}
outputfile_errcov='./output/lorenz63/'errcov_${outputname}
exp_log=./log/${today}_${prg}_${da_method}.log

# for EnKF
if [ ${da_method} = 'EnKF' ]; then 
  outputname=${da_method}_${mem}'m'_${enkf_method}'.csv'
  outputfile_score='./output/lorenz63/'${outputname}
  outputfile_errcov='./output/lorenz63/'errcov_${outputname}
  exp_log=./log/${today}_${prg}_${da_method}_${enkf_method}.log
fi

#----------------------------------------------------------------------
# +++ Run exp.
#----------------------------------------------------------------------
cp ${CDIR}/common/common.f90 common_mod.f90
cp ${CDIR}/common/common_mtx.f90 common_mtx_mod.f90
cp ${CDIR}/common/SFMT.f90 SFMT_mod.f90
cp ${CDIR}/common/netlib.f netlib_mod.f


# +++ compile
gfortran -fbounds-check \
  SFMT_mod.f90 common_mod.f90 netlib_mod.f common_mtx_mod.f90 lorenz63_prm.f90 lorenz63_cal.f90 lorenz63_main.f90 \
  -o ${prg} -I/usr/local/include -llapack -lblas \
  -w # error message Suppression

./${prg} > ${exp_log} << EOF
  $day_info
    today = '${today}'
  /
  &set_dim
    nx = ${nx},
    ny = ${ny}
  /
  &set_parm
    nt_asm = ${nt_asm}, nt_prd = ${nt_prd},
    obs_interval = ${obs_interval}
  /
  &da_setting
    da_method = '${da_method}',
    alpha = ${alpha}
  /
  &intg_setting
    intg_method = '${intg_method}'
  /
  &enkf_setting
    mems = ${mem},
    enkf_method = '${enkf_method}'
  /
  &initial_score
    x_tinit = ${x_tinit},  y_tinit = ${y_tinit}, z_tinit = ${z_tinit}
    x_sinit = ${x_sinit},  y_sinit = ${y_sinit}, z_sinit = ${z_sinit}
  /
  &output
    output_file = '${outputfile_score}',
    output_file_errcov = '${outputfile_errcov}'
    opt_veach = .${boolen}.
  /
EOF

rm -rf *.mod *_mod.f *_mod.f90 ${prg}
echo 'Normal END'

exit

#----------------------------------------------------------------------
# +++ matrix memo !
#  The following cases are for 
#  three predictor var and two observation var.
#----------------------------------------------------------------------
# 
# Forecast error covariance matrix
# Background error covariance matrix
# Pf = (nx, nx)
Pf_init=( 1.0d0 0.0d0 0.0d0
          0.0d0 1.0d0 0.0d0
          0.0d0 0.0d0 1.0d0 )

# Observation error covariance matrix
# R = (ny, ny)
R_init=(  1.0d0 0.0d0 
          0.0d0 1.0d0 )

# Kalman gain matrix
# K = (nx, ny)
K_init=(  0.0d0 0.0d0 
          0.0d0 0.0d0 
          0.0d0 0.0d0 )

# Observation　operator
# H = (ny, nx)
H_init=(  1.0d0 0.0d0 0.0d0 
          0.0d0 1.0d0 0.0d0 )
