#$ -S /bin/sh

modulos='mod_types.f90 mod_rndgen_multiple.f90 geraRede.f90 mod_tools_redes.f90 mod_epidemic.f90'

principal='sirs_qmf_rede_sint_estac.f90'

executavel='qmf_g23'

job=${executavel}

#flags='-traceback -check all -heap-arrays -O3 -g -fp-stack-check'
flags=''

rm ${exe} &> erroRM_exe.log 

ifort ${flags} ${modulos} ${principal} -o ${executavel} #&> erroCompilacao.log

chmod +x ${executavel}

rm -r *.mod &> erroRM.log


qsub -N ${job} -cwd executa.sh ${executavel}

