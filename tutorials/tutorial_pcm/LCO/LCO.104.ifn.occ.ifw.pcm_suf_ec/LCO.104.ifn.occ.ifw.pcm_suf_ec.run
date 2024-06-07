#!/bin/bash   
#SBATCH -J LCO.104.ifn.occ.ifw.pcm_suf_ec
#SBATCH -t 72:00:00
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -o /home/a.boev//LCO//LCO.104.ifn.occ.ifw.pcm_suf_ec/sbatch.out
#SBATCH -e /home/a.boev//LCO//LCO.104.ifn.occ.ifw.pcm_suf_ec/sbatch.err
#SBATCH -p AMG-medium
cd /home/a.boev//LCO//LCO.104.ifn.occ.ifw.pcm_suf_ec/
module load Compiler/Intel/17u8; module load Q-Ch/VASP/5.4.4_SOL; module load ScriptLang/python/3.6i_2018u3
 ulimit -s unlimited

export PATH=$PATH:/home/a.boev/tools/
touch RUNNING
#Basic run:
cp 1.POSCAR POSCAR
mpirun vasp_std >LCO.104.ifn.occ.ifw.pcm_suf_ec.1.log
sleep 20
mv OUTCAR 1.OUTCAR
mv CONTCAR 1.CONTCAR
mv CHG 1.CHG
mv CHGCAR 1.CHGCAR
mv LOCPOT 1.LOCPOT

#Footer section: 
rm vasprun.xml XDATCAR PROCAR AECCAR0 EIGENVAL ELFCAR OSZICAR PARCHG WAVEDER DOSCAR AECCAR2 WAVECAR 
rm RUNNING
