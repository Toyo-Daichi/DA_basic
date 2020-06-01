#!/bin/bash
# attached by lorenz96.f90

set -ex
CDIR=`pwd`
tool='normal' #spinup or normal
ts_check='true' # if spinup output is 'true' or 'sim'.
prg=lorenz96_${tool}_maintools
today=$(date "+%Y%m%d%H%M")
rm -rf *.mod ${prg}

#----------------------------------------------------------------------
# +++ Set intial setting
#----------------------------------------------------------------------

# +++ model dimension
nx=40
dt=0.05d0
force=8.0d0
oneday=0.2d0

# +++ integral period(day)
spinup_period=365
normal_period=40

da_method='EnKF'
intg_method='Runge-Kutta'
mem=1
enkf_method='none'
if [ ${da_method} = 'EnKF' ]; then 
  mem=15
  enkf_method='PO' # 'PO' or "SRF"
fi

# +++ adaptive inflation
alpha=0.0d0

# +++ making obs. info
obs_xintv=1
obs_tintv=2

# +++ output info
out_boolen='true' # write putput
da_boolen='true'
outputname='lorenz96'
initial_true_file='./output/'${outputname}/'spinup_true_initial_'${nx}'n.csv'
initial_sim_file='./output/'${outputname}/'spinup_sim_initial_'${nx}'n.csv'
output_true_file='./output/'${outputname}/'normal_true_score_'${nx}'n.csv'
output_sim_file='./output/'${outputname}/'normal_sim_score_'${nx}'n.csv'
output_anl_file='./output/'${outputname}/'normal_'${da_method}'_anl_score_'${nx}'n.csv'
output_errcov_file='./output/'${outputname}/'Error_matrix_'${da_method}'_'${nx}'n.csv'

if [ ${da_method} = 'EnKF' ]; then 
  output_anl_file='./output/'${outputname}/${tool}'_'${da_method}${mem}'m_anl_score_'${nx}'n.csv'
fi
output_obs_file='./output/'${outputname}/${tool}'_obs_score_'${nx}'.csv'
#----------------------------------------------------------------------
# +++ Run exp.
#----------------------------------------------------------------------
cp ${CDIR}/common/common.f90 common_mod.f90
cp ${CDIR}/common/common_mtx.f90 common_mtx_mod.f90
cp ${CDIR}/common/SFMT.f90 SFMT_mod.f90
cp ${CDIR}/common/netlib.f netlib_mod.f

# +++ compile
gfortran -fbounds-check  \
  SFMT_mod.f90 common_mod.f90 netlib_mod.f common_mtx_mod.f90 lorenz96_prm.f90 lorenz96_cal.f90 lorenz96_main.f90 \
  -o ${prg} -I/usr/local/include -lm -lblas -llapack \
  -w # error message Suppression

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
    da_method = '${da_method}',
    alpha = ${alpha}
  /
  &enkf_setting
    mems = ${mem},
    enkf_method = '${enkf_method}'
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
    initial_true_file  = '${initial_true_file}',
    initial_sim_file   = '${initial_sim_file}',
    output_true_file   = '${output_true_file}',
    output_anl_file    = '${output_anl_file}',
    output_sim_file    = '${output_sim_file}',
    output_obs_file    = '${output_obs_file}', 
    output_errcov_file = '${output_errcov_file}', 
    opt_veach = .${out_boolen}.
  /
EOF

rm -rf *.mod *_mod.f90 ${prg}
echo 'Normal END'

exit
