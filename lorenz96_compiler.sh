#!/bin/bash
# attached by lorenz96.f90

set -ex
CDIR=`pwd`
tool='normal' #spinup or normal
prg=lorenz96_${tool}_maintools
today=$(date "+%Y%m%d%H%M")
rm -rf *.mod ${prg}

#----------------------------------------------------------------------
# +++ Set intial setting
#----------------------------------------------------------------------
nx=24
dt=0.005d0
force=3.85d0
oneday=0.2d0

# +++ integral period
spinup_period=365
normal_period=40

da_method=''
intg_method='Runge-Kutta'

# +++ making obs. info
obs_xintv=2
obs_tintv=20

# +++ output info
out_boolen='false' # write putput
da_boolen='false'
outputname='lorenz96'
initial_true_file='./output/'${outputname}_spinup_initial'.csv'
initial_sim_file='./output/'${outputname}_spinup_initial'.csv'
outputfile='./output/'${outputname}'.csv'

#----------------------------------------------------------------------
# +++ Run exp.
#----------------------------------------------------------------------
cp ${CDIR}/common/SFMT.f90 SFMT_mod.f90
cp ${CDIR}/common/common.f90 common_mod.f90

# +++ compile
gfortran -fbounds-check \
SFMT_mod.f90 common_mod.f90 lorenz96_prm.f90 lorenz96_cal.f90 lorenz96_main.f90 -o ${prg} 

./${prg} > ./log/${today}_${prg}_${DA_METHOD}.log << EOF
  &set_parm
    nx = ${nx},
    dt = ${dt},
    force = ${force},
    oneday = ${oneday}
  /
  &set_exp
    tool = '${tool}'
    intg_method = '${intg_method}'
  /
  &set_da_exp
    da_veach = .${da_boolen}.,
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
    initial_true_file = '${initial_true_file}'
    initial_sim_file = '${initial_sim_file}'
    output_file  = '${outputfile}',
    opt_veach    = .${out_boolen}.
  /
EOF

rm -rf *.mod *_mod.f90 ${prg}
echo 'Normal END'
exit
