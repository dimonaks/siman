current
    - commit qe1:
        - 


    - Massive changes by AK
    Now LOCPOT is moved to 1.LOCPOT if 'l' in savefile is provided, otherwise removed
    showing atom numbers from zero in determine_symmetry_positions()
    fontsize argument was removed from plot_dos(). Use plot_param dictionary, which is transfered to fit_and_plot
    velocities and predictor are read from CONTCAR, allows to continiue MD

    - commit gaus:
        - determine_file_format moved to inout.py
        - 



0.9.6

    new regime for add() 
        calc_method = 'polaron'

    new module small_classes to where some of simple classes would be transferred

    new control parameter added to set
        periodic - defines if calculation is using periodic or open boundary conditions

    
    sleep changed from 5 to 1 for PBS

    to show neb xyz files use show = 'neb_geo'

    procmemgb key added for PBS schedule 
    -neb analysis on cluster is switched off

    update_incar() moved to functions.py


    add_neb()

        - corenum arg removed , use add_loop_dic dictionary
        - now if *end_pos_types_z* is used only provided elements are chosen as possible final positions 
        - new naming convention for neb calculations when *end_pos_types_z* parameter is used. Use *old_behaviour* = '261018' for old naming

    Structure ()

        - new method get_total_number_electrons(), works only with cl.init, when zval exist 

    CalculationVasp()

        - new fields after reading OUTCAR:

            - memory (float) -  total amount of memory  per job (Gb)
            - memory_max (float) -  Maximum memory used per job (Gb)
            - nelect (int) - total number of electrons
            - list_e_conv


    IMPORTANT NOTES:

        - To optimize the size of database file from time to time run siman using ::

         header.reorganize = 1

        - now select works correctly  with inheritance 





0.9.5
new regime for calc_method in add_loop:
    'Monte-Carlo' - Monter Carlo simulation above DFT

sshpass parameter in clusters now can take 'proxy' value. this allows to connect cluster with sshpass using proxy. 
TODO: now proxy server and path to password is hardcoded in run_on_server(), push_to_server() and download()

new attibute added to calculation objects: self.cluster - dictionary with all parameters for cluster


read_poscar() moved from classes.py to inout.py as a separate function


new type 'void' with Z = 300 added to system. It is not written in POSCAR. Be carefull, bugs are expected, 
since cl.add_potcar() and cl.calculate_nbands() are affected

BUGS:
bug related to the wrong reading of forces with selective dynamics was fixed (st.select was used from init, while in end the order of atoms can change)
now in the end structure is always created as new one. So, no i

but with incorrect writing of magmom in the case of atom reordering for poscar was fixed;
warning! if you provide magmom in set, than you should take care that the provided order is consistent with that in generated POSCAR


0.8.3
cl.fix_atoms() method added
bug fixed in cl.run()
calc_oxidation_states() improved, use silent = 0 to show charges
get_surface_area() added to Structure()
get_dipole() added to Structure()
cl.get_bader_ACF() improved
creating_sets tutorial updated
'jmol' option of show in res() added
calculation name in output changed to db['name']
st.jmol() now uses vasp format
new property st.outfile - path to output file added
mcif argument added to write_xyz() enabling visualization of magnetic structures with jmol=1 
cl.add_new_name() method added
cluster_tools/fit_tool.py new paths added for numpy support on cee cluster (local change for skoltech)
now modules are written for SLURM from project_conf.py, please add them as in example below:
    CLUSTERS['cee']['modules'] = 'module add prun/1.0; module add intel/16.0.2.181; module add impi/5.1.3.181\n'
}

New tutorial for calc_barriers subroutine


