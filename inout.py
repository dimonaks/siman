#Copyright Aksyonov D.A
from __future__ import division, unicode_literals, absolute_import 
import os
from header import printlog
from functions import element_name_inv

def makedir(path):
    """
    *path* - path to some file 
    Make dirname(path) directory if it does not exist
    """
    dirname = os.path.dirname(path)

    if dirname and not os.path.exists(dirname):
        os.makedirs(dirname)
        print_and_log("Directory", dirname, " was created", imp = 'y')
    return



def write_xyz(st, path = None, repeat = 1, shift = 1.0,  gbpos2 = None, gbwidth = 1 , 
    imp_positions = [], specialcommand = None, analysis = None, show_around = None, replications = None, nnumber = 6, topview = True,
    filename = None, file_name = None, full_cell = False, orientation = None, boundbox = 2, withgb = False,
    include_boundary = 2, rotate = None, imp_sub_positions = None, jmol = None, show_around_x = None,
    ):
    """Writes st structure in xyz format in the folder xyz/path

    if repeat == 2: produces jmol script
    shift - in rprimd[1][1] - shift of the second view
    gbpos2 - position of grain boundary in A
    gbwidth - atoms aroung gbpos2 will be colored differently

    imp_positions - (x1,x2,x3, element, label)- xcart and element name coordinates additionally to be added to structure; to visulaze all impurity positions: for jmol, additional key 's', 'i' can be added after element
    imp_sub_positions - list of atom numbers; the typat of these atoms is changed: not used now

    specialcommand - any command at the end of script

    analysis - additional processing, allows to show only specifice atoms, 
        'imp_surrounding' - shows Ti atoms only around impurity
        nnumber - number of neighbours to show
        show_around - choose atom number around which to show
        show_around_x - show atoms around point, has higher priority

    replications - list of replications, (2,2,2) 

    full_cell - returns atoms to cell and replicate boundary atoms

    jmol - 1,0 - allow to use jmol
    """
    if replications:
        st = replic(st, mul = replications, inv = 1 )
  
    def update_var():
        return st.rprimd, st.xcart, st.xred, st.typat, st.znucl, len(st.xred)

    rprimd, xcart, xred, typat, znucl, natom = update_var()


    if file_name:
        name = file_name
    elif filename:
        name = filename
    else:
        name = st.name



    if natom != len(xred) != len(xcart) != len(typat) or len(znucl) != max(typat): 
        printlog( "Error! write_xyz: check your arrays.\n\n"    )
    # print st.natom, len(st.xred), len(st.xcart), len(st.typat), len(st.znucl), max(st.typat)
    
    printlog("write_xyz(): Name is", name, important = 'n')
    if name == '': name = 'noname'


    if xcart == [] or len(xcart) != len(xred):
        printlog( "Warining! write_xyz: len(xcart) != len(xred) making xcart from xred.\n")
        xcart = xred2xcart(xred, rprimd)
        #print xcart[1]
    
    if path:
        basepath = path
    else:
        basepath = 'xyz/'



    xyzfile = os.path.join(basepath, name+".xyz")
    makedir(xyzfile)


    """Processing section"""


    if analysis == 'imp_surrounding':
        lxcart = []
        ltypat = []
        i=0
        for t, x in zip(typat, xcart):
            
            condition = False
            # print show_around, 'show'
            if show_around:
                # print i, condition
                condition = (i + 1 == show_around)
                # print i, condition

            else:
                condition = (t > 1) # compat with prev behav
            
            # print 'se', condition

            if condition: 
                # lxcart.append(x)
                # ltypat.append(t)
                # print x, ' x'
                x_t = local_surrounding(x, st, nnumber, control = 'atoms', periodic = True)
                # print x_t[1]
                lxcart+=x_t[0]
                ltypat+=x_t[1]
            i+=1
        
        if show_around_x:
            x = show_around_x
            x_t = local_surrounding(x, st, nnumber, control = 'atoms', periodic = True)
            # print x_t[1]
            lxcart+=x_t[0]
            ltypat+=x_t[1]            


        xcart = lxcart
        typat = ltypat
        natom = len(typat)
        # print natom, 'nat'



    """Include atoms on the edge of cell"""
    if full_cell:
        # print xred
        # print natom
        # st = return_atoms_to_cell(st)
        # print xred
        st = replic(st, mul = (1,1,2), inv = 0, cut_one_cell = 1, include_boundary = include_boundary)
        # print natom, st.natom

        # print st.xred

        rprimd, xcart, xred, typat, znucl, natom = update_var()
        
    # asdegf

    """Writing section"""   
    # printlog("Writing xyz: "+xyzfile, imp = 'y')

    #analyze imp_positions
    if imp_sub_positions == None:
        imp_sub_positions = []
    nsub = 0
    for pos in imp_positions:
        # if len(pos) > 4:
        #     if 's' not in pos[4]: continue # skip interstitial positions
        
        xs = np.asarray([pos[0],pos[1],pos[2]])
        nsub+=1
        # print xs
        for i, x in enumerate(xcart):
            # print np.linalg.norm( x-xs)
            if np.linalg.norm( x-xs) < 1:
                imp_sub_positions.append(i)

    if imp_sub_positions : 
        printlog( imp_sub_positions, ': numbers of found atoms to be changed ' )


    # for i in sorted(indices, reverse=True):
    #     del somelist[i]




    with open(xyzfile,'w') as f:
        for i in range(repeat):
            f.write(str(natom + len(imp_positions)-nsub + 3)+"\n") #+3 vectors
            f.write(name+"\n")

            if imp_positions: 
                for i, el in enumerate(imp_positions):
                    # if len(el) != 4: continue
                    f.write( "%s %.5f %.5f %.5f \n"%( el[3], el[0], el[1], el[2] ) )
                    # print 'composite -pointsize 60 label:{0:d} -geometry +{1:d}+{2:d} 1.png 2.png'.format(i, el[0], el[1])


            for i in range(natom):
                typ = typat[i] - 1
                
                z = int ( znucl[ typ ] )

                if i in imp_sub_positions: 
                    # f.write( "Be " )
                    continue
                else:
                    el = element_name_inv(z)
                    f.write( el+" " )


                f.write( "%.5f %.5f %.5f \n"%( xcart[i][0], xcart[i][1], xcart[i][2] ) )

            for r in st.rprimd:
                f.write('Tv {:.10f} {:.10f} {:.10f}\n'.format(*r)  )


    # os._exit(1)
    printlog('File', xyzfile, 'was written', imp = 'y')


    if jmol:
        """
        script mode for jmol. Create script file as well for elobarate visualization
        """
        
        """Choose gb atoms to change their color"""
        printlog( 'position of boundary 2', gbpos2)
        atomselection = ''

        #create consistent xcart_new list like it will be in Jmol
        xcart_new = []
        for i, x in enumerate(xcart):
            if i in imp_sub_positions: continue
            xcart_new.append(x)    



        if gbpos2:
            
            gbpos1 = gbpos2 - rprimd[0][0]/2.
            gbatoms = []
            
            for i, x in enumerate(xcart_new):
                # print i
                # if x[0] > gbpos1 - gbwidth/2. and x[0] < gbpos1 + gbwidth/2.:
                if abs(x[0] - gbpos1) < gbwidth/2.:
                    gbatoms.append(i)
                    # print i, x[0], abs(x[0] - gbpos1)
                if abs(x[0] - gbpos2) < gbwidth/2.:
                # if x[0] > gbpos2 - gbwidth/2. and x[0] < gbpos2 + gbwidth/2.:
                    # print i, x[0], abs(x[0] - gbpos2)
                    gbatoms.append(i)
            printlog( 'Atoms at GB:', gbatoms)
            atomselection = ''
            for i in gbatoms:
                atomselection +='Ti'+str(i+1+len(imp_positions))+','
            atomselection = atomselection[:-1]


        # elif withgb: # color half of cell
        # else: # color half of cell
        #     # pass
            # atomselection = 'atomno>'+str(0+len(imp_positions) )+' and atomno<'+str(( natom + len(imp_positions)  )/2-1)





        xyzfile = os.getcwd()+'/'+xyzfile
        scriptfile = basepath+name+".jmol"
        pngfile = os.getcwd()+'/'+basepath+name+".png"
        
        printlog( 'imp_positions = ',imp_positions)
        write_jmol(xyzfile, pngfile, scriptfile, atomselection, topview = topview, rprimd =rprimd, shift = shift, label = [(pos[3], pos[4]) for pos in imp_positions], 
            specialcommand = specialcommand, orientation = orientation, boundbox =boundbox, rotate = rotate)


    return



def write_lammps(st, filename = '', charges = None):
    """Writes structure in lammps format 

    charges (list of float) - list of charges for each atom type
    """

    rprimd = st.rprimd
    xcart  = st.xcart
    xred   = st.xred
    typat = st.typat
    ntypat = st.ntypat
    znucl = st.znucl
    name = st.name
    natom = st.natom

    if natom != len(xred) != len(xcart) != len(typat) or len(znucl) != max(typat): 
        print_and_log( "Error! write_xyz: check your structure"    )
    
    if name == '': 
        name = 'noname'
    if xcart == [] or len(xcart) != len(xred):
        print_and_log( "Warining! write_xyz: len(xcart) != len(xred) making xcart from xred.\n", imp = 'y')
        xcart = xred2xcart(xred, rprimd)
        #print xcart[1]

    if not filename:
        filename = 'lammps/'+name

    filename+='.inp'

    makedir(filename)



    """Write lammps structure file;  """
    if 1:
        """ My version; valid only for octahedral cells""" 
        printlog( "Warining! write_lammps(): this func supports only orthogonal cells", imp = 'Y')

        with open(filename+'','w') as f:
            f.write("Lammps format "+name+'\n')
            f.write(str(natom)+" atoms\n")
            f.write(str(ntypat)+" atom types\n")
            f.write("{:10.8f}  {:10.8f}  xlo xhi\n".format(0, rprimd[0][0]))
            f.write("{:10.8f}  {:10.8f}  ylo yhi\n".format(0, rprimd[1][1]))
            f.write("{:10.8f}  {:10.8f}  zlo zhi\n".format(0, rprimd[2][2]))
            f.write("0.00000000  0.00000000  0.00000000  xy xz yz\n")
            f.write("\nAtoms\n\n")

            for i, x in enumerate(xcart):
                f.write("{0:8d} {1:2d}".format(i+1, typat[i]))
                if charges:
                    f.write(" {:6.3f}".format(charges[typat[i]-1] ) )
                f.write(" {:12.6f}  {:12.6f}  {:12.6f}\n".format(x[0], x[1], x[2] ))
            


            f.write("\n")
    
        printlog('File', filename, 'was written', imp = 'y')



    else:
        """Write poscar and convert from poscar to lammps using external script; Valid for arbitary cells"""
        cl.write_structure('POSCAR', 'dir', path = 'voronoi_analysis/', state = state)
        runBash("voronoi_analysis/VASP-poscar2lammps.awk voronoi_analysis/POSCAR > "+filepath)
    


    if 0:
        """Write lammps.in file """
        with open('voronoi_analysis/voronoi.in','w') as f:
            f.write("""units           metal
                    atom_style atomic
                    boundary        p p p\n""")
            f.write("read_data   /home/aksenov/programs/Simulation_wrapper/siman1/voronoi_analysis/structure.lammps\n")
    #         f.write('lattice   custom 1 ')
    #         for i, a in enumerate(rprimd):
    #             f.write(' a'+str(i+1))
    #             for x in a:
    #                 f.write(' '+str(x))
            
    #         f.write(' &\n')
    #         for x in xred:
    #             f.write(' basis {0:f} {1:f} {2:f}&\n '.format(x[0], x[1], x[2]) )
    #         f.write("""\n
    # region 1 prism 0 1 0 1 0 1  1 0 0
    # create_box 1 prism
    # create_atoms 1 prism""")

            for i in range(ntypat):
                f.write('\nmass '+str(i+1)+' '+str(int(znucl[i]))+'\n')
            
            f.write('pair_style      lj/cut 2.0\n')
            for i in range(ntypat):
                for j in range(i, ntypat):
                    f.write('pair_coeff      '+str(i+1)+' '+str(j+1)+' 0.0 1.0\n')


            f.write("""compute v1 all voronoi/atom
                    dump    d1 all custom 1 /home/aksenov/programs/Simulation_wrapper/siman1/voronoi_analysis/dump.voro id type x y z c_v1[1] c_v1[2]
                    run 0
                    uncompute v1\n""")

    return
