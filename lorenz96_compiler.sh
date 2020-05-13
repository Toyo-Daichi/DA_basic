#!/bin/bash
# attached by lorenz96.f90

set -ex
CDIR=`pwd`
tool='normal' #spinup or normal
ts_check='sim' # if spinup output is 'true' or 'sim'.
prg=lorenz96_${tool}_maintools
today=$(date "+%Y%m%d%H%M")
rm -rf *.mod ${prg}

#----------------------------------------------------------------------
# +++ Set intial setting
#----------------------------------------------------------------------
nx=24
dt=0.05d0
force=3.85d0
oneday=0.2d0

# +++ integral period(day)
spinup_period=365
normal_period=40

da_method='KF'
intg_method='Runge-Kutta'
mem=40

# +++ making obs. info
obs_xintv=2
obs_tintv=10

# +++ output info
out_boolen='true' # write putput
da_boolen='true'
outputname='lorenz96'
initial_true_file='./output/'${outputname}/'spinup_true_initial.csv'
initial_sim_file='./output/'${outputname}/'spinup_sim_initial.csv'
output_true_file='./output/'${outputname}/'normal_true_score.csv'
output_NoDA_file='./output/'${outputname}/'normal_NoDA_score.csv'
output_DA_file='./output/'${outputname}/'normal_'${da_method}'_DA_score.csv'
if [ ${da_method} = 'EnKF' ]; then 
  output_DA_file='./output/'${outputname}/${tool}'_'${da_method}${mem}'m_DA_score.csv'
fi
output_obs_file='./output/'${outputname}/${tool}'_obs_score.csv'
#----------------------------------------------------------------------
# +++ Run exp.
#----------------------------------------------------------------------
cp ${CDIR}/common/SFMT.f90 SFMT_mod.f90
cp ${CDIR}/common/common.f90 common_mod.f90

# +++ compile
gfortran -fbounds-check -I/usr/local/include -lm -lblas -llapack \
SFMT_mod.f90 common_mod.f90 lorenz96_prm.f90 lorenz96_cal.f90 lorenz96_main.f90 -o ${prg} 

./${prg} > ./log/${today}_${prg}_${da_method}.log << EOF
  &set_parm
    nx = ${nx},
    dt = ${dt},
    force = ${force},
    oneday = ${oneday}
  /
  &set_exp
    tool = '${tool}',
    ts_check = '${ts_check}',
    intg_method = '${intg_method}'
  /
  &set_da_exp
    da_veach  = .${da_boolen}.,
    mems      = ${mem},
    da_method = '${da_method}'
  /
  &set_period
    spinup_period = ${spinup_period},
    normal_period = ${normal_period}
  /
  &set_mobs
    obs_xintv = ${obs_xintv},
    obs_tintv = ${obs_tintv}
  /
  &output
    initial_true_file = '${initial_true_file}',
    initial_sim_file  = '${initial_sim_file}',
    output_true_file  = '${output_true_file}',
    output_DA_file    = '${output_DA_file}',
    output_NoDA_file  = '${output_NoDA_file}',
    output_obs_file   = '${output_obs_file}', 
    opt_veach = .${out_boolen}.
  /
EOF

rm -rf *.mod *_mod.f90 ${prg}
echo 'Normal END'
exit
