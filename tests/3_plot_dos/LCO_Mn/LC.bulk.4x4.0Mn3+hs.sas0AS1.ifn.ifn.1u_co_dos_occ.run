#!/bin/bash   
#SBATCH -J LC.bulk.4x4.0Mn3+hs.sas0AS1.ifn.ifn.1u_co_dos_occ
#SBATCH -t 72:00:00
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -o /home/a.boev/vasp/surseg_tem//seg_paper/sol/as//LC.bulk.4x4.0Mn3+hs.sas0AS1.ifn.ifn.1u_co_dos_occ/sbatch.out
#SBATCH -e /home/a.boev/vasp/surseg_tem//seg_paper/sol/as//LC.bulk.4x4.0Mn3+hs.sas0AS1.ifn.ifn.1u_co_dos_occ/sbatch.err
#SBATCH -p AMG-medium
cd /home/a.boev/vasp/surseg_tem//seg_paper/sol/as//LC.bulk.4x4.0Mn3+hs.sas0AS1.ifn.ifn.1u_co_dos_occ/
module load Compiler/Intel/17u8; module load Q-Ch/VASP/5.4.4_OMC; module load ScriptLang/python/3.6i_2018u3
 ulimit -s unlimited

export PATH=$PATH:/home/a.boev/tools/
touch RUNNING
#Basic run:
cp 1.POSCAR POSCAR
mpirun vasp_std >LC.bulk.4x4.0Mn3+hs.sas0AS1.ifn.ifn.1u_co_dos_occ.1.log
sleep 20
mv OUTCAR 1.OUTCAR
mv CONTCAR 1.CONTCAR
mv CHGCAR 1.CHGCAR
mv DOSCAR 1.DOSCAR
mv vasprun.xml 1.vasprun.xml

#Footer section: 
rm CHG AECCAR0 AECCAR2 EIGENVAL PROCAR WAVECAR OSZICAR XDATCAR LOCPOT WAVEDER PARCHG 
rm RUNNING
