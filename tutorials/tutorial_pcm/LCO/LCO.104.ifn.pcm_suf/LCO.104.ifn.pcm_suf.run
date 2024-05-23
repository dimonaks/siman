#!/bin/bash   
#SBATCH -J LCO.104.ifn.pcm_suf
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -o /home/a.boev//LCO//LCO.104.ifn.pcm_suf/sbatch.out
#SBATCH -e /home/a.boev//LCO//LCO.104.ifn.pcm_suf/sbatch.err
#SBATCH -p AMG-medium
cd /home/a.boev//LCO//LCO.104.ifn.pcm_suf/
module load Compiler/Intel/17u8; module load Q-Ch/VASP/5.4.4_OPT; module load ScriptLang/python/3.6i_2018u3
 ulimit -s unlimited

export PATH=$PATH:/home/a.boev/tools/
touch RUNNING
#Basic run:
cp 1.POSCAR POSCAR
mpirun vasp_std >LCO.104.ifn.pcm_suf.1.log
sleep 20
mv OUTCAR 1.OUTCAR
mv CONTCAR 1.CONTCAR
mv CHG 1.CHG
mv CHGCAR 1.CHGCAR
mv LOCPOT 1.LOCPOT
mv WAVECAR 1.WAVECAR

#Footer section: 
rm XDATCAR EIGENVAL PROCAR LOCPOT vasprun.xml OSZICAR WAVEDER AECCAR0 WAVECAR ELFCAR DOSCAR PARCHG AECCAR2 CHG 
rm RUNNING
