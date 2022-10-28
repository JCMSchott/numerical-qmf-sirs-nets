#!/bin/bash

job=${1}
bashi=${2}
exec=${3}

qsub -N ${job} -cwd executa.sh ${exec}

