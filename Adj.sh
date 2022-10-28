#!/bin/bash

modulos='mod_types.f90 mod_rndgen_multiple.f90 geraRede.f90 mod_tools_redes.f90'
principal='Adj.f90'
executavel='ex_Adj_g35'

flags= #'-traceback -check all -heap-arrays -O3 -g -fp-stack-check'

ifort ${modulos} ${principal} ${flags} -o ${executavel}

rm -r *.mod

chmod +x ${executavel}
