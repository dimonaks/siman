#Copyright Aksyonov D.A
from __future__ import division, unicode_literals, absolute_import 
import os
import numpy  as np

from siman import header
from siman.header import printlog, runBash
from siman.functions import element_name_inv, unique_elements
from siman.small_functions import makedir, is_list_like, list2string
from siman.geo import local_surrounding, replic


def read_xyz(st, filename, rprimd = None):
    """
    Read xyz file into st

    rprimd (list of lists) - if None or [None, ] then Tv are read; if Tv does not exist then create automatically 

    """
    with open(filename,'r') as f:
        nlines = int(f.readline())
        st.name = f.readline().strip()
        
        # try:
        if 'SG' in st.name:
            printlog('Error! Space group record detected in xyz, please finish code', imp = 'Y')
            # st.name.split('SG')


        elements = []
        st.xcart = []
        st.rprimd = []
        for i in range(nlines):
            xc = f.readline().split()
            if len(xc) == 0:
                printlog('Warning! xyz file is broken, not enough lines')
                break

            if 'Tv' in xc[0]: 
                st.rprimd.append(np.asarray(xc[1:], dtype = float) )

            else:
                elements.append(xc[0])
                st.xcart.append(np.asarray(xc[1:], dtype = float) )

        st.natom = len(st.xcart)
        


        

    st.znucl = [element_name_inv(el) for el in unique_elements(elements)]
    

    elements_z = [element_name_inv(el) for el in elements]
    st.typat = []
    for z in elements_z:
        st.typat.append( st.znucl.index(z)+1 )

    st.ntypat = len(st.znucl)


    # print(st.rprimd)
    if rprimd == None or None in rprimd or 0 in rprimd or len(rprimd) != 3:
        printlog('None detected in *rprimd*, I use vectors from xyz file')
        if len(st.rprimd) != 3:
            printlog('Only these primitive vectors were found in xyz :\n', np.round(st.rprimd, 3), '\nI take rest from *rprimd*', imp ='y')
            if rprimd:
                for r in rprimd:
                    if is_list_like(r):
                        st.rprimd.append(r)
            else:
                printlog('Error! Please provide vector in *rprimd*')





    else:
        
        printlog('I use vectors from *rprimd*')
        st.rprimd = rprimd

    if len(st.rprimd) != 3:
        printlog('Error! Check *rprimd* or Tv in xyz')


    if st.get_volume() < 0:
        printlog('Warning! rprimd gives negative volume, I exchange vectors 2 and 3', imp = 'y')
        t = st.rprimd[1]
        st.rprimd[1] = st.rprimd[2]
        st.rprimd[2] = t
        if st.get_volume() < 0:
            printlog('Error! still negative volume, check your primitive vectors', imp = 'y')
        st.tmap = [1,3]
    else:
        st.tmap = [1,2]

    printlog('Final rprimd = \n', np.round(st.rprimd, 3), imp = 'y')



    st.nznucl = st.get_nznucl()

    st.recip = st.get_recip()

    st.update_xred()

    st.reorder_for_vasp(inplace = True)

    # print(st.perm)
    return st


def read_poscar(st, filename, new = True):
    # from classes import Structure
    #read poscar
    selective_dynamics = None

    elements_list = []

    # st = Structure()
    
    if new:
        st.name = os.path.basename(filename).replace('POSCAR', '').replace('CONTCAR', '')
    
    try:
        if '.' in st.name[-1]:
            st.name = st.name[0:-1]
    except:
        pass


    with open(filename,'r') as f:
        name = f.readline().strip()
        
        if new:
            st.des = name
        # print self.name, "self.name"


        # st.name = self.name
        # print(f.readline())
        mul = float( f.readline().split('!')[0] )
        # print 'mul', mul


        st.rprimd = []
        for i in 0, 1, 2:
            vec = f.readline().split()
            st.rprimd.append( np.asarray([float(vec[0])*mul, float(vec[1])*mul, float(vec[2])*mul]) )

        st.nznucl = []

        ilist = f.readline().split() #nznucl of elements?
        
        try:
            int(ilist[0])
            vasp5 = False
        except:
            vasp5 = True


        if vasp5:
            printlog('Vasp5 detected')
            for el in ilist:
                elements_list.append(el)
            printlog('elements_list:', elements_list)

            ilist = f.readline().split()
        else:
            printlog('Vasp4 detected')


        
        for z in ilist:
            st.nznucl.append( int(z)  )


        temp_line = f.readline()

        if temp_line[0] in ['s', 'S']:
            printlog('selective dynamics detected') 
            selective_dynamics = True
            temp_line = f.readline()

        type_of_coordinates = temp_line


        st.xred = []

        coordinates = []
        select = []


        if len(elements_list) > 0:
            read_elements = 0
        else:
            read_elements = 1

        for nz in st.nznucl:

            for i in range(nz):
                vec = f.readline().split()
                coordinates.append( np.asarray([float(vec[0]), float(vec[1]), float(vec[2])]) )

                if read_elements and len(vec) == 4: # elements may be added by pymatgen
                    # printlog("Probably elements names are added at the end of coordinates, trying to read them")
                    if vec[3] not in elements_list:
                        elements_list.append(vec[3])
                
                if selective_dynamics:
                    # convert 'T'/'F' to True/False
                    flagset = [True, True, True]
                    for fi, flag in enumerate(vec[3:6]):
                        if flag == 'F':
                            flagset[fi] = False
                    # print(flagset)
                    select.append(flagset)

        st.select = select

        if "Car" in type_of_coordinates or 'car' in type_of_coordinates:
            st.xcart  = coordinates
            st.update_xred()

            
        elif "dir" in type_of_coordinates or 'Dir' in type_of_coordinates:
            st.xred  = coordinates
            st.update_xcart()

        elif 'None' in type_of_coordinates:
            pass

        else:
            printlog("Error! The type of coordinates should be 'car' or 'dir' ")
            raise NameError



        if 'Species order:' in name:
            printlog('I detect that input file was generated by cif2cell\n')
            name = name.split(':')[-1]


        if not elements_list:
            elements_list = name.split('!')[0].strip().split()
            printlog('I take elements from the first line, The line is '+str(name.split('!'))+' you could use ! to add comment after name+\n')
            # print(elements_list)
            if 'i2a' in elements_list[0]:
                printlog('i2a list detected')
                el = elements_list[0].split('[')[-1].replace(']','')
                elements_list = el.split(',')
        else:
            printlog("Elements names have been taken from the end of coordinates, pymatgen file?\n")



        st.znucl = []
        for elname in elements_list:
            st.znucl.append( element_name_inv(elname) )
        # printlog('znucl is ')



        st.natom = len(st.xred)

        st.ntypat = len(st.znucl)

        st.typat = []
        for i, nz in enumerate(st.nznucl):
            for j in range(nz):
                st.typat.append(i+1)

        #Determine reciprocal vectors
        st.recip = st.get_recip()


            # if hasattr(self.init, 'vel'):
            #     print "I write to POSCAR velocity as well"
            #     f.write("Cartesian\n")
            #     for v in self.init.vel:
            #         f.write( '%.12f %.12f %.12f\n'%(v[0]*to_ang, v[1]*to_ang, v[2]*to_ang) )

    printlog('The following Z were read = '+ str(st.znucl)+'\n')


    printlog('VASP POSCAR format', filename, " was read\n")

    return st













def write_jmol(xyzfile, pngfile, scriptfile = None, atomselection = None, topview = 0, orientation = None,
    axis = False, bonds = True, rprimd = None, shift = None, rotate = None,
    label = None, high_contrast = None, specialcommand = None,
    boundbox = 2, atom_labels = None):
    """
    atomselection - string in gmol format with number of atoms to be nrotateSelected
    topview - additional top view, requires two models in xyz
    orientation - additional rotation
    axis - add axes
    rotate - rotation of all atoms around view axis in degrees
    label (tuple ()) - used for impurities, please decribe
    atom_labels (bool) - turn on atom labels
    """
    if not scriptfile:
        scriptfile = os.getcwd()+'/'+'temporary_jmol_script'
    with open(scriptfile,'w') as f:
        f.write('set frank off\n') #no jmol label
        if bonds:
            f.write('set autobond on\n')

        else:
            f.write('set autobond off\n set bonds off\n')

        f.write('load "'+xyzfile+'"\n')

        f.write('select all \n') #250
        if 0:
           f.write('cpk 250 \nwireframe 0.3\n') 

        f.write('background white \n')
        # f.write('select Ti* \ncolor [20,120,250] \nselect C* \ncolor [80,80,80]\n cpk 100\n')
        f.write('set perspectivedepth off\n')
        



        if boundbox:
            f.write('set boundbox ' +str(boundbox)+ ' \n')



        # f.write('set specular 85\n set specpower 85\n set diffuse 85\n')
        if high_contrast: #allows to make better view for black and white printing 
            f.write('set ambient 10 \nset specular 95\n set specpower 95\n set diffuse 95\n')



        
        if axis:
            f.write('set axes 10 \naxes scale 2.5 \n')
            f.write('axes labels "X" "Y" "Z" "" \n')
            f.write('color  axes  red \n')
            f.write('font axes 26 \n')


        if orientation:
            f.write(orientation+'\n')


        if atomselection:
            f.write('select '+atomselection+'\n')
            f.write('color purple    \n')

        if topview:
            f.write('select * /2  \ntranslateSelected 0 '+str(-rprimd[1][1]*shift)+' 0\nrotateSelected X 90\n')
        
            f.write('wireframe 0.1\ncpk 150\nmodel 0\n#zoom 60\n')



        if label:
            j = 1
            name_old = ''
            for i, el in enumerate(label):
                name  = el[0]+el[1]
                if name != name_old: j = 1
                label = str(j)+el[1]
                # print "label",label
                f.write('select '+el[0]+str(i+1)+'\ncpk 200\nset labeloffset 0 0\nset labelfront\ncolor label black\nlabel '+label+'\n font label 24 bold \n')
                j+=1
                name_old = name



        if atom_labels:
            f.write('select all\nset label "%e"\nset labeloffset 0 0\nset labelfront off\ncolor label black\nfont label 18 bold \n')




        if rotate:
            f.write('rotate z '+str(rotate)+'\n')

        if specialcommand:
            f.write(specialcommand+'\n')

        if 1:
            f.write('set displayCellParameters False ;\n')



        
        # f.write('write image 2800 2800 png "'+pngfile+'"')
        f.write('write image 1800 1800 png "'+pngfile+'"')
    
    # print(header.PATH2JMOL)
    # sys.exit()
    printlog( runBash(header.PATH2JMOL+' -ions '+scriptfile) )
    # print runBash('convert '+pngfile+' -shave 0x5% -trim '+pngfile) #cut by 5% from up and down (shave) and that trim left background
    printlog( pngfile )
    printlog( runBash('convert '+pngfile+' -trim '+pngfile)  ) # trim background
    printlog('png file by Jmol',pngfile, 'was written', imp = 'y' )
    # print(header.PATH2JMOL)
    return


def write_xyz(st = None, path = None, filename = None, file_name = None,
    include_vectors = True, repeat = 1, shift_2view = 1.0, replications = None, full_cell = False, 
    analysis = None, show_around = None, show_around_x = None,  nnumber = 6, only_elements = None,
    gbpos2 = None, gbwidth = 1, withgb = False, include_boundary = 2,
    imp_positions = [], imp_sub_positions = None,
    jmol = None, specialcommand = None, jmol_args = None, sts = None, mcif = 0, suf = ''
    ):
    """Writes st structure in xyz format in the folder xyz/path
    #void are visualized with Pu
    if repeat == 2: produces jmol script
    shift_2view - in rprimd[1][1] - shift of the second view
    gbpos2 - position of grain boundary in A
    gbwidth - atoms aroung gbpos2 will be colored differently

    imp_positions - (x1,x2,x3, element, label)- xcart and element name coordinates additionally to be added to structure; to visulaze all impurity positions: for jmol, additional key 's', 'i' can be added after element
    imp_sub_positions - list of atom numbers; the typat of these atoms is changed: not used now


    analysis - additional processing, allows to show only specifice atoms, 
        'imp_surrounding' - shows Ti atoms only around impurity
        nnumber - number of neighbours to show
        show_around - choose atom number around which to show, from 1
        show_around_x - show atoms around point, has higher priority
        only_elements - see local_surrounding

    replications - list of replications, (2,2,2) 

    full_cell - returns atoms to cell and replicate boundary atoms

    include_vectors (bool) - write primitive vectors to xyz

    jmol - 1,0 -  use jmol to produce png picture
    jmol_args - see write_jmol()
    mcif - write magnetic cif for jmol


    specialcommand - any command at the end of jmol script
    suf - additional suffix for name

    sts - list of Structure - write several structures to xyz file - other options are not working in this regime
    """

    if jmol_args == None:
        jmol_args = {}

    if st == None:
        st = sts[0]


    if replications:
        st = replic(st, mul = replications, inv = 1 )
  
    def update_var(st):
        if st.natom != len(st.xred) != len(st.xcart) != len(st.typat) or len(st.znucl) != max(st.typat): 
            printlog( "Error! write_xyz: check your arrays.\n\n"    )

        if st.xcart == [] or len(st.xcart) != len(st.xred):
            printlog( "Warining! write_xyz: len(xcart) != len(xred) making xcart from xred.\n")
            st.xcart = xred2xcart(st.xred, st.rprimd)
            #print xcart[1]

        return st.rprimd, st.xcart, st.xred, st.typat, st.znucl, len(st.xred)

    
    st = st.copy()
    rprimd, xcart, xred, typat, znucl, natom = update_var(st)




    if file_name:
        name = file_name
    elif filename:
        name = filename
    else:
        name = st.name+suf


    if sts:
        name+='_traj'

    
    printlog("write_xyz(): Name is", name, important = 'n')
    
    if name == '': 
        name = 'noname'



    
    if path:
        basepath = path
    else:
        basepath = 'xyz/'



    suf = ''





    """Processing section"""


    if analysis == 'imp_surrounding':
        printlog('analysis = imp_surrounding', imp = 'y')

        if show_around == 0:
            printlog('Error! number of atom *show_around* should start from 1')

        suf = '_loc'+str(show_around)
        lxcart = []
        ltypat = []
        i=0


        if is_list_like(show_around_x):
            x = show_around_x
            x_t = local_surrounding(x, st, nnumber, control = 'atoms', periodic = True, only_elements = only_elements)
            # print('write_xyz: local_surround:', x_t)
            lxcart+=x_t[0]
            ltypat+=x_t[1]            
        else:

            for t, x in zip(typat, xcart):
                
                condition = False
                # print show_around, 'show'
                if show_around:
                    # print i, condition
                    condition = (i + 1 == show_around)
                    # print i, condition

                else:
                    condition = (t > 1) # compat with prev behav, to show around any impurities (all atoms with typat more than one)
                
                # print 'se', condition

                if condition: 
                    # print('Atom at', x, 'used as central')
                    # lxcart.append(x)
                    # ltypat.append(t)
                    # print x, ' x'
                    x_t = local_surrounding(x, st, nnumber, control = 'atoms', periodic = True, only_elements = only_elements)
                    print (x_t)
                    lxcart+=x_t[0]
                    ltypat+=x_t[1]
                i+=1
        

        xcart = lxcart
        typat = ltypat
        natom = len(typat)
        # print natom, 'nat'
        # print('Number of neighbours', natom  )
        st.xcart = xcart
        st.typat = typat
        st.natom = natom
        st.update_xred()





    """Include atoms on the edge of cell"""
    if full_cell:
        # print xred
        # print natom
        # st = return_atoms_to_cell(st)
        # print xred
        st = replic(st, mul = (1,1,2), inv = 0, cut_one_cell = 1, include_boundary = include_boundary)
        # print natom, st.natom

        # print st.xred

        rprimd, xcart, xred, typat, znucl, natom = update_var(st)
        
    # asdegf

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


    if include_vectors:
        nvect = 3
    else:
        nvect = 0

    """Writing section"""   
    name+=suf
    xyzfile = os.path.join(basepath, name+".xyz")
    makedir(xyzfile)

    def write(st):
        rprimd, xcart, xred, typat, znucl, natom = update_var(st)

        f.write(str(natom + len(imp_positions)-nsub + nvect)+"\n") #+3 vectors
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
                if el == 'void':
                    el = 'Pu'

                f.write( el+" " )

            f.write( "%.5f %.5f %.5f \n"%( xcart[i][0], xcart[i][1], xcart[i][2] ) )

        if include_vectors:
            for r in st.rprimd:
                f.write('Tv {:.10f} {:.10f} {:.10f}\n'.format(*r)  )

    with open(xyzfile,'w') as f:
        if sts:
            for st in sts:
                write(st)
        else:
            for i in range(repeat):
                write(st)



    # os._exit(1)
    printlog('File', xyzfile, 'was written', imp = 'y')

    pngfile = None
    

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





        # xyzfile = os.getcwd()+'/'+xyzfile
        if mcif:
            xyzfile = st.write_cif(mcif = 1)
        else:
            xyzfile = st.write_poscar()
        
        scriptfile = basepath+name+".jmol"
        bn = (basepath+name).replace('.', '_')
        pngfile = os.getcwd()+'/'+bn+".png"
        
        printlog( 'imp_positions = ',imp_positions)
        write_jmol(xyzfile, pngfile, scriptfile, atomselection, rprimd =rprimd, shift = shift_2view, label = [(pos[3], pos[4]) for pos in imp_positions], 
            specialcommand = specialcommand, **jmol_args)


    return xyzfile, pngfile



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
        printlog( "Error! write_xyz: check your structure"    )
    
    if name == '': 
        name = 'noname'
    if xcart == [] or len(xcart) != len(xred):
        printlog( "Warining! write_xyz: len(xcart) != len(xred) making xcart from xred.\n", imp = 'y')
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


def write_occmatrix(occs, folder):
    #create OCCMATRIX 
    
    makedir(folder)
    printlog('I create OCCMATRIX in ', folder, imp = 'y')
    filename = folder+'/OCCMATRIX'
    with open(filename, 'w', newline = '') as f:
        numat = len(occs)
        f.write(str(numat)+'  #num of atoms to be specified\n')
        
        at_nums = occs.keys()
        at_spin = [] # # 2 or 1
        at_ltyp = [] # l - orbital type, 1 - s, 2 - d, 3 - f
        for key in occs: 
            occ = occs[key]
            if len(occ) == 10: # spin polarized, d orbital
                at_spin.append(2)
                at_ltyp.append(2)
            else:
                raise RuntimeError # please write by yourself for other cases


        for i, l, s in zip(at_nums, at_spin, at_ltyp):

            f.write(list2string([i+1, l, s])+'    #i, l, s\n')
            # for sp in range(s):
            f.write('spin 1\n')
            for row in occs[i][ 0:len(occs[i])//s ]:
                f.write(list2string(row)+'\n')
            if s == 2:
                f.write('spin 2\n')
                for row in occs[i][ len(occs[i])//s: ]:
                    f.write(list2string(row)+'\n')
            f.write('\n')
    return filename