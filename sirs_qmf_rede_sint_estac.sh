#!/bin/bash

#label1=${1}

#label2=${2}

modulos='mod_types.f90 mod_rndgen_multiple.f90 geraRede.f90 mod_tools_redes.f90' # mod_epidemic.f90'

principal='sirs_qmf_rede_sint_estac.f90'

exe='qmf12g28a05_1'

#flags='-traceback -check all -heap-arrays -O3 -g -fp-stack-check -check uninit -check bounds -ftrapuv'

flags=''

rm ${exe} &> erroRM_exe.log 

ifort ${flags} ${modulos} ${principal} -o ${exe} #&> erroCompilacao.log

chmod +x ${exe}

rm -r *.mod &> erroRM.log

time ./${exe} #${label1} ${label2}
