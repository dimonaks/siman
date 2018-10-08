#!/bin/python
from functions import *
"""  
Programm Monte-Carlo

"""

init = read_lammps_poscar("conf.dump_init")

for i, tat in enumerate(init.typat):
    if tat == 1:
        init.typat1.append( (i, tat )  )
    elif tat == 2:
        init.typat2.append( (i, tat )  )

#print init.typat1
#print init.typat2
#print "Energy of init structure", calc_energy()

"""1. Generate two vacancies  """
add_vac_auto = False
add_Al = True
if add_vac_auto:
    ir1 = randint(513, len(init.typat)-1 ) 


    init.xred10.append(  init.xred[ir1]   )
    del init.typat[ir1]
    del init.xred[ir1]

    ir2 = randint(513, len(init.typat)-1 ) 

    init.xred10.append(  init.xred[ir2]   )
    del init.typat[ir2]
    del init.xred[ir2]
    
    ir3 = randint(513, len(init.typat)-1 ) 

    init.xred10.append(  init.xred[ir3]   )
    del init.typat[ir3]
    del init.xred[ir3]

    ir4 = randint(513, len(init.typat)-1 ) 

    init.xred10.append(  init.xred[ir4]   )
    del init.typat[ir4]
    del init.xred[ir4]

    ir5 = randint(513, len(init.typat)-1 ) 

    init.xred10.append(  init.xred[ir5]   )
    del init.typat[ir5]
    del init.xred[ir5]

    ir6 = randint(513, len(init.typat)-1 ) 

    init.xred10.append(  init.xred[ir6]   )
    del init.typat[ir6]
    del init.xred[ir6]

    ir7 = randint(513, len(init.typat)-1 ) 

    init.xred10.append(  init.xred[ir7]   )
    del init.typat[ir7]
    del init.xred[ir7]
    
    init.natom = len(init.typat)
    #print ir1, ir2
    #init.typat[ ir1 ] = 10
    #init.typat[ ir2 ] = 10
elif add_Al:


    for j in range(10):
        alnum, fenum = [], []
        for i, t in enumerate(init.typat):
            if t == 1: alnum.append(i)
            if t == 2: fenum.append(i)
        print 'list of iron numbers', fenum

        ir1 = randint(0, len(fenum)-1 ) #random Fe
        print 'random number in list of iron atoms ', ir1
        ir1 = fenum[ir1]
        print 'random number of iron atom in complete list of atoms', ir1
        init.typat[ir1] = 1


    antisites = copy.deepcopy(init)


"""2. Run initial lammps  """
write_lammps_poscar(init, "conf.dump")
write_lammps_poscar(init, "conf.dump_init_two_vacancies")
init.etot = calc_energy()
print "Energy after making two vacancies",init.etot

#cur = Structure()
nmcstep = 10000

#for i_mcstep in range(nmcstep):
#    write_lammps_poscar(init, "conf.dump"+str(i_mcstep) )

for i_mcstep in range(nmcstep):
    

    """3. Exchange one any vacancy and one any atom """
    
    cur = copy.deepcopy(init)
    #print i_mcstep%2
    #if i_mcstep%2 == 0:

    cur = exchange_atom_atom(cur, 12)
    #else:
    #    cur = exchange_atom_vac(cur)
     #   pass


    """4. Run lammps  """
    write_lammps_poscar(cur, "conf.dump" )

    cur.etot = calc_energy()
    
    if i_mcstep%2 == 0:
        print "Energy after exhanging atom and vacancy", cur.etot
    else:
        print "Energy after changing atom and atom", cur.etot
        pass

    """5. Check if to accept new structure  """
    if metropolis(init.etot, cur.etot):
        #print init.etot, cur.etot, "accept"
        init = copy.deepcopy(cur)


write_lammps_poscar(init, "conf.dump_last")
print "\n\nFinal energy is", init.etot


write_xyz(cur, path = 'conf.dump_last', repeat = 1)
write_xyz(antisites, path = 'before_calc.dump_init', repeat = 1)


#emfile = initposcar.read()
#initposcar_words = memfile.split()
#xred = read_lammps_xred("zs", 128, initposcar_words)

#print xred



#print runBash("ls")