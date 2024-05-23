#!/bin/bash   
#SBATCH -J LCO.bulk.um_opt
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -o /home/a.boev//LCO//LCO.bulk.um_opt/sbatch.out
#SBATCH -e /home/a.boev//LCO//LCO.bulk.um_opt/sbatch.err
#SBATCH -p AMG-medium
cd /home/a.boev//LCO//LCO.bulk.um_opt/
module load Compiler/Intel/17u8; module load Q-Ch/VASP/5.4.4_OPT; module load ScriptLang/python/3.6i_2018u3
 ulimit -s unlimited

export PATH=$PATH:/home/a.boev/tools/
touch RUNNING
#Basic run:
cp 1.POSCAR POSCAR
mpirun vasp_std >LCO.bulk.um_opt.1.log
sleep 20
mv OUTCAR 1.OUTCAR
mv CONTCAR 1.CONTCAR
mv CHGCAR 1.CHGCAR

#Footer section: 
rm XDATCAR EIGENVAL PROCAR LOCPOT vasprun.xml OSZICAR WAVEDER AECCAR0 WAVECAR ELFCAR DOSCAR PARCHG AECCAR2 CHG 
rm RUNNING
