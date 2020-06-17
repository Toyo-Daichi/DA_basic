#!/bin/bash
# attached by lorenz96.f90

set -ex
CDIR=`pwd`
tool='normal' #spinup or normal
ts_check='sim' # if spinup output is 'true' or 'sim'.
prg=lorenz96_${tool}_maintools
today=$(date "+%Y%m%d%H%M")
rm -rf *.mod ${prg}

echo ${today}

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

# +++ exp. info
da_method='EnKF'
intg_method='Runge-Kutta'
mem=1
enkf_method='none'
if [ ${da_method} = 'EnKF' ]; then 
  mem=10
  enkf_method='SRF' # 'PO' or "SRF" or "ETKF"
fi

# +++ adaptive inflation & localization
alpha=0.02d0
localization_mode=1
shchur_length_scale=3

# +++ making obs. info
obs_tintv=1
# >> For OBSERVATION OPERATER(H)
obs_xintv=1
#  OBS x coordinate set is full veriosn -> obs_set=0
#  OBS x coordinate set is interval lack version -> obs_set=1
#  OBS x coordinate set is bias kacl(wind) version -> obs_set=2
obs_set=0
obs_wnd_point=0
if [ ${obs_xintv} -ge 2  ]; then obs_set=1 ;fi
if [ ${obs_xintv} -eq 99 ]; then 
  obs_set=2
  obs_wnd_point=26
fi

# +++ output info
out_boolen='true' # write putput
da_boolen='true'
output_dir='./output/lorenz96/'

# >> (1) spinup file set
initial_true_file=${output_dir}/'spinup_true_initial_'${nx}'ndim.csv'
initial_sim_file=${output_dir}/'spinup_sim_initial_'${nx}'ndim.csv'

# >> (2) experiment file set
# *** The three file names below do not change in any way.
output_true_file=${output_dir}/'normal_true_score_'${nx}'ndim.csv'
output_sim_file=${output_dir}/'normal_sim_score_'${nx}'ndim.csv'
output_obs_file=${output_dir}/'normal_obs_score_'${nx}'ndim.csv'
#
output_anl_file=${output_dir}/'normal_'${da_method}'_anl_score_'${nx}'ndim.csv'
output_anlinc_file=${output_dir}/'normal_'${da_method}'_anlinc_'${nx}'ndim.csv'
output_errcov_file=${output_dir}/'normal_'${da_method}'_errcov_'${nx}'ndim.csv'
exp_log=./log/${today}_${prg}_${da_method}.log

if [ ${da_method} = 'EnKF' ]; then 
  output_anl_file=${output_dir}/'normal_'${da_method}${mem}'m_anl_score_'${nx}'ndim.csv'
  output_anlinc_file=${output_dir}/'normal_'${da_method}${mem}'m_anlinc_'${nx}'ndim.csv'
  output_errcov_file=${output_dir}/'normal_'${da_method}${mem}'m_errcov_'${nx}'ndim.csv'
  exp_log=./log/${today}_${prg}_${da_method}_${enkf_method}.log
  if [ ${localization_mode} == 1 ]; then
  #output_anl_file=${output_dir}/'normal_'${da_method}${mem}'m_anl_score_'${nx}'ndim_loc_'${shchur_length_scale}'_alpha_'${alpha}'.csv'
  output_anl_file=${output_dir}/'normal_'${da_method}${mem}'m_anl_score_'${nx}'ndim_loc_.csv'
  output_anlinc_file=${output_dir}/'normal_'${da_method}${mem}'m_anlinc_'${nx}'ndim_loc.csv'
  output_errcov_file=${output_dir}/'normal_'${da_method}${mem}'m_errcov_'${nx}'ndim_loc.csv'
  exp_log=./log/${today}_${prg}_${da_method}_${enkf_method}_loc.log
  fi
fi

# >> *** wind effect experiment
input_wnd_errcov_file=''
if [ ${obs_xintv} -eq 99  ]; then
  input_wnd_errcov_file=${output_dir}'/input_KF_errcov_40ndim.csv'
fi
#----------------------------------------------------------------------
# +++ Run exp.
#----------------------------------------------------------------------
cp ${CDIR}/common/common.f90 common_mod.f90
cp ${CDIR}/common/common_mtx.f90 common_mtx_mod.f90
cp ${CDIR}/common/common_enkf.f90 common_enkf_mod.f90
cp ${CDIR}/common/SFMT.f90 SFMT_mod.f90
cp ${CDIR}/common/netlib.f netlib_mod.f

# +++ compile
gfortran -fbounds-check  \
  SFMT_mod.f90 common_mod.f90 netlib_mod.f common_mtx_mod.f90 common_enkf_mod.f90 lorenz96_prm.f90 lorenz96_cal.f90 lorenz96_main.f90 \
  -o ${prg} -I/usr/local/include -lm -lblas -llapack \
  -w # error message Suppression

for alpha in 0.01d0 0.02d0 0.03d0 0.04d0  0.05d0  0.06d0  0.07d0  0.08d0  0.09d0  0.10d0
do
for shchur_length_scale in 1 2 3 4 5 6 7 8 9 10
do
output_anl_file=${output_dir}/'normal_'${da_method}${mem}'m_anl_score_'${nx}'ndim_loc_'${shchur_length_scale}'_alpha_'${alpha}'.csv'


./${prg} > ${exp_log} << EOF
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
    da_method = '${da_method}',
    alpha = ${alpha}
  /
  &enkf_setting
    mems = ${mem},
    enkf_method = '${enkf_method}'
    localization_mode = ${localization_mode},
    shchur_length_scale = ${shchur_length_scale}
  /
  &set_period
    spinup_period = ${spinup_period},
    normal_period = ${normal_period}
  /
  &set_obs
    obs_set = ${obs_set},
    obs_xintv = ${obs_xintv},
    obs_tintv = ${obs_tintv},
    obs_wnd_point = ${obs_wnd_point}
  /
  &inoutput_file
    initial_true_file  = '${initial_true_file}',
    initial_sim_file   = '${initial_sim_file}',
    input_wnd_errcov_file = '${input_wnd_errcov_file}'
  /
  &exp_outputfile
    output_true_file   = '${output_true_file}',
    output_anl_file    = '${output_anl_file}',
    output_sim_file    = '${output_sim_file}',
    output_obs_file    = '${output_obs_file}', 
    output_errcov_file = '${output_errcov_file}', 
    output_anlinc_file = '${output_anlinc_file}',
    opt_veach = .${out_boolen}.
  /

EOF

#rm -rf *.mod *_mod.f* ${prg}
echo 'Normal END'
done
done
exit
