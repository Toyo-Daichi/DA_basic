#!/bin/bash
# attached by lorenz96.f90

set -ex
prg=lorenz96_maintools
today=$(date "+%Y%m%d%H%M")
rm -rf *.mod ${prg}

#----------------------------------------------------------------------
# +++ Set intial setting
#----------------------------------------------------------------------

# +++ output info
boolen='false' # write putput
outputname='lorenz96'
outputfile='./output/'${outputname}'.csv'

#----------------------------------------------------------------------
# +++ Run exp.
#----------------------------------------------------------------------
gfortran -fbounds-check SFMT.f90 kinddef.f90 lorenz96_prm.f90 lorenz96_main.f90 lorenz96_spinup.f90 -o ${prg} 

./${prg} > ./log/${today}_${prg}_${DA_METHOD}.log << EOF
  &output
    output_file  = '${outputfile}',
    opt_veach    = .${boolen}.
  /
EOF

rm -rf *.mod ${prg}
echo 'Normal END'
exit
