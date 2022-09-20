#Copyright Aksyonov D.A
from __future__ import division, unicode_literals, absolute_import 
import os, io, re, math, sys
from csv import reader


import numpy  as np
try:
    import pandas as pd 
except:
    print('inout.py: Cant import pandas')

try:
    from tabulate import tabulate 
except:
    print('inout.py: Cant import tabulate')



from siman import header
from siman.header import printlog, runBash, plt
from siman.functions import element_name_inv, unique_elements, smoother
from siman.small_functions import makedir, is_list_like, list2string, red_prec
from siman.small_classes import empty_struct
from siman.geo import local_surrounding, replic



def determine_file_format(input_geo_file):

    supported_file_formats = {'abinit':'.geo',   'vasp':'POSCAR',   'cif':'.cif', 'xyz':'.xyz'} #format name:format specifier

    if ('POSCAR' in input_geo_file or 'CONTCAR' in input_geo_file) and not '.geo' in input_geo_file:
        input_geo_format = 'vasp'
    elif '.vasp' in input_geo_file:
        input_geo_format = 'vasp'
    elif '.cif' in input_geo_file:
        input_geo_format = 'cif'
    elif '.geo' in input_geo_file:
        input_geo_format = 'abinit'
    elif '.xyz' in input_geo_file:
        input_geo_format = 'xyz'
    elif '.gau' in input_geo_file:
        input_geo_format = 'gaussian'


    else:
        printlog("Attention! File format of",input_geo_file ,"is unknown, should be from ", supported_file_formats, 'skipping')
        input_geo_format = None
    return input_geo_format




def read_csv(filename, delimiter =','):
    with open(filename, 'r') as read_obj:
        # pass the file object to reader() to get the reader object
        csv_reader = reader(read_obj, delimiter = delimiter)
        # Pass reader object to list() to get a list of lists
        list_of_rows = list(csv_reader)
        # print(list_of_rows)
    return list_of_rows


def read_xyz(st, filename, rprimd = None):
    """
    Read xyz file into st

    rprimd (list of lists) - if None or [None, ] then Tv are read; if Tv does not exist then create automatically 

    """
    with open(filename,'r') as f:
        nlines = int(f.readline())
        st.name = f.readline().strip().split()[0]
        
        # try:
        if 'SG' in st.name:
            printlog('Warning! Space group record detected in xyz, please finish code', imp = 'Y')
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
    if rprimd is None or any(x is None for x in rprimd) or len(rprimd) != 3:
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
        mul_str = f.readline()
        # print(filename)
        # print(name)
        # print(mul_str)
        # sys.exit()
        mul = float( mul_str.split('!')[0] )
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


        velocities = []
        newline = f.readline() # if velocities are available - they are separated by one new line
        vel1 = f.readline()
        # print(vel1)
        if vel1:
            #usually vasp return zero velocities, which are not needed, skip them
            vec = vel1.split()
            vel_vec = np.asarray([float(vec[0]), float(vec[1]), float(vec[2])])
            vellen = abs(np.linalg.norm(vel_vec))
            # print(vellen)
        if vel1 and vellen > 1e-12:
            vec = vel1.split()
            velocities.append( np.asarray([float(vec[0]), float(vec[1]), float(vec[2])]) )
            printlog('Reading velocities from POS/CONT-CAR', imp = 'y')
            for i in range(len(coordinates)-1): # vel for first atom is already read
                vec = f.readline().split()
                velocities.append( np.asarray([float(vec[0]), float(vec[1]), float(vec[2])]) )
            st.vel = velocities

            newline = f.readline() # new empty line between blocks
            start = f.readline() # something appears after velocity
            if start:
                printlog('Reading predictor-corrector from POS/CONT-CAR', imp = 'y')
                predictor = start+f.read()
                # print(predictor)
                # sys.exit()
                st.predictor = predictor





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

    
    help:
    frame all - turn on all frames


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

        # f.write('frame all\n') #no jmol label
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
            # f.write('set axes 10 \naxes scale 2.5 \n')
            f.write('set axes 10 \naxes scale 2.0 \n')
            f.write('axes labels "X" "Y" "Z" "" \n')
            # f.write('color  axes  red \n')
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
    return pngfile


def write_xyz(st = None, path = None, filename = None, file_name = None,
    include_vectors = True, repeat = 1, shift_2view = 1.0, replications = None, full_cell = False, 
    analysis = None, show_around = None, show_around_x = None,  nnumber = 6, only_elements = None,
    gbpos2 = None, gbwidth = 1, withgb = False, include_boundary = 2,
    imp_positions = [], imp_sub_positions = None,
    jmol = None, 
    # specialcommand = None, # shoud be in jmol_args
    jmol_args = None, sts = None, mcif = 0, suf = ''
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
    jmol_args (dict) - arguments to write_jmol see write_jmol()
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
        if st.natom != len(st.xcart) != len(st.typat) or len(st.znucl) != max(st.typat): 
            printlog( "Error! write_xyz: check your arrays.\n\n"    )

        if st.xcart is [None] :
            printlog( "Warining! write_xyz: len(xcart) != len(xred) making xcart from xred.\n")
            st.xcart = xred2xcart(st.xred, st.rprimd)
            #print xcart[1]

        return st.rprimd, st.xcart, st.xred, st.typat, st.znucl, len(st.xcart)

    
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
        printlog('analysis == imp_surrounding', imp = '')

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
                    # print (x_t)
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

        # print(range(natom))

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
        elif sts:
            ''
        else:
            xyzfile = st.write_poscar()
        
        scriptfile = basepath+name+".jmol"
        bn = (basepath+name).replace('.', '_')
        pngfile = os.getcwd()+'/'+bn+".png"
        
        printlog( 'imp_positions = ',imp_positions)
        write_jmol(xyzfile, pngfile, scriptfile, atomselection, rprimd =rprimd, shift = shift_2view, label = [(pos[3], pos[4]) for pos in imp_positions], 
             **jmol_args)


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
            elif len(occ) == 6: # spin polarized, p orbital
                at_spin.append(2)
                at_ltyp.append(1)
            else:
                raise RuntimeError # please write by yourself for other cases


        for i, l, s in zip(at_nums, at_ltyp, at_spin):

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


def write_geometry_aims(st, filename, coord_type = 'cart', periodic = True):

    rprimd = st.rprimd
    with io.open(filename,'w', newline = '') as f:
        if periodic:
            for i in 0, 1, 2:
                f.write('lattice_vector   {:10.6f} {:10.6f} {:10.6f}\n'.format(rprimd[i][0],rprimd[i][1],rprimd[i][2]) )
            f.write("\n")

        if None in st.magmom:
            magmom = [None]*st.natom
        else:
            magmom = st.magmom

        for x, el, mag in zip(st.xcart, st.get_elements(), magmom ):
            f.write("atom  {:12.10f}  {:12.10f}  {:12.10f}  {:2s} \n".format(x[0], x[1], x[2], el) )
            if mag is not None:
                f.write('initial_moment '+str(mag)+'\n') 



def read_vasp_out(cl, load = '', out_type = '', show = '', voronoi = '', path_to_outcar = '', path_to_contcar = '', ):
    """Try to read xred from CONCAR and calculate xcart"""

    self = cl

    printlog('Path to CONTCAR', path_to_contcar)
    if os.path.exists(path_to_contcar):
        contcar_exist   = True
    else:
        contcar_exist   = False


    printlog('The status of CONTCAR file is', contcar_exist)
    # self.end.update_xred()

    if contcar_exist:
        try:
            self.end = read_poscar(self.end, path_to_contcar, new = False) # read from CONTCAR
            contcar_read = True
        except:
            contcar_read = False
        if not contcar_read:
            printlog('Attention!, error in CONTCAR:', path_to_contcar, '. I use data from outcar')

    else:
        printlog('Attention!, No CONTCAR:', path_to_contcar, '. I use data from outcar')
        contcar_read = False






    read = 1
    if read:
        if 0: #please use this only for linux or create cross-platform way
            nw = runBash('sed -n "/NPAR = approx SQRT( number of cores)/=" '+path_to_outcar) #remove warinig
            tmp = path_to_outcar+".tmp"
            if nw:
                nw = int(nw)
                runBash("sed '"+str(nw-11)+","+str(nw+8)+"d' "+path_to_outcar+">"+tmp+";mv "+tmp+" "+path_to_outcar)


        with open(path_to_outcar, 'rb') as outcar:
            
            printlog("Start reading from "+ path_to_outcar, imp = 'n')
            # outcarlines = outcar.readlines()
            text = outcar.read().decode(errors='replace')
            outcarlines = str(text).split('\n')




        re_lengths = re.compile("length of vectors")
        re_eltime = re.compile("Elapsed time")
        re_nkpts = re.compile("NKPTS")
        iterat = 0
        niter = 1
        i_line = 0
        mdstep_prev = 0
        dipol = None
        ediff = 0
        self.mdstep = 1
        warnings = 0#""
        self.time = 0
        self.memory = 0 # total memory per job
        nscflist = []; mdstep_old = 1; niter_old = 0
        maxforce = []; average = [];  gstress =[]
        # mforce = []
        self.list_e_sigma0 = []
        self.list_e_without_entr = []
        self.list_e_conv = [] # convergence of energy - all steps
        # try:
        #     self.end = copy.deepcopy(self.init) # below needed end values will be updated
        # except:
        if not contcar_read:
            self.end = self.init.new()

        # if not hasattr(self.end, "natom"): 
        #     self.end.natom = self.natom
        #Structure() #create structure object with end values after calculation
        #self.end.typat = self.typat
        #self.end.znucl = self.znucl
        self.end.name = self.name+'.end'
        self.end.list_xcart = []
        self.energy = empty_struct()
        self.end.zval = []
        de_each_md = 0 # to control convergence each md step
        de_each_md_list = []



        nsgroup = None
        magnitudes = []
        self.mag_sum = [] #toatal mag summed by atoms, +augmentation

        tot_mag_by_atoms = [] #magnetic moments by atoms on each step
        tot_chg_by_atoms = []
        tot_mag_by_mag_atoms = []


        ldauu = None
        e_sig0 = 0 #energy sigma 0 every scf iteration
        occ_matrices = {} # the number of atom is the key

        #dipole corrections
        self.dipole_min_pos = []  
        self.e_dipol_quadrupol_cor = []
        self.e_added_field_ion = []


        #detect neb, improve this 
        if hasattr(self.set, 'vasp_params'):
            images = self.set.vasp_params.get('IMAGES') or 1
        else:
            images = 1
        #which kind of forces to use
         # CHAIN + TOTAL  (eV/Angst)
        neb_flag = 0
        l = 0
        for line in outcarlines:
            # l+=1
            if 'energy of chain is' in line:
                neb_flag = 1
            if 'LOOP+' in line:
                break
        # print(l)

        if neb_flag:

            force_keyword = 'CHAIN + TOTAL  (eV/Angst)'
            ff  = (0, 1, 2)
            force_prefix = ' chain+tot '

        else:
            force_keyword = 'TOTAL-FORCE'
            ff  = (3, 4, 5)
            force_prefix = ' tot '
        # print(force_keyword)





        # try:
        #     spin_polarized = self.set.spin_polarized # again it will be better to determine this from outcar 
        # except:
        #     spin_polarized = None



        self.potcar_lines = []
        self.stress = None
        self.extpress = 0

        self.intstress = None
        spin_polarized = None
        ifmaglist = None
        for line in outcarlines:

            #Check bands

            # if 'band No.' in line:
            #     kpoint = float(outcarlines[i_line-1].split()[1])
            #     lastocc = float(outcarlines[i_line+self.nbands].split()[2])
            #     lastbandno = outcarlines[i_line+self.nbands].split()[0]
            #     if lastocc > 0:
            #         print "Warning!!! at kpoint ", kpoint, " last band No. ",lastbandno, " is not empty ", lastocc

            if 'TITEL' in line:
                self.potcar_lines.append( line.split()[2:] )

            if 'NGXF=' in line:
                # print(line)
                l = re.sub("[^\d\.\ ]", "", line) #remove everything except digits, dot, and space
                # sys.exit()
                self.ngxf = [int(n) for n in l.split()] # ngxf, ngyf, ngzf
                # print(self.ngxf ) 



            if 'LEXCH  =' in line:
                # print(line)
                self.xc_pot = line.split()[2].strip() #xc from potential 

            if 'GGA     =' in line:
                # print(line)
                self.xc_inc = line.split()[2].strip() #xc from incar

            if 'NELECT' in line:
                self.nelect = int(float(line.split()[2]))

            if '   ZVAL' in line:
                # zval  = line.split('')
                # print(line)
                self.end.zval = [int(float(n)) for n in line.split()[2:]]
                # print(self.end.zval)


            if 'ions per type =' in line:
                if not contcar_read:
                    self.end.nznucl = [int(n) for n in line.split()[4:]]
                    self.end.ntypat = len(self.end.nznucl)

                    self.end.natom  = sum(self.end.nznucl)

                    #correction of bug; Take into account that VASP changes typat by sorting impurities of the same type.
                    self.end.typat = []
                    for i, nz in enumerate(self.end.nznucl):
                        for j in range(nz):
                            self.end.typat.append(i+1)
                    #correction of bug



                    # print(self.potcar_lines)
                    elements = [t[1].split('_')[0] for t in self.potcar_lines]
                    # printlog('I read ',elements, 'from outcar')
                    self.end.znucl = [element_name_inv(el) for el in elements]
                    # print (self.end.znucl)

                ifmaglist, _ = self.end.get_maglist()


            if 'ISPIN' in line:
                if line.split()[2] == '2':
                    spin_polarized = True
                    self.spin_polarized = spin_polarized
                else:
                    spin_polarized = False
                    self.spin_polarized = False


            if "TOO FEW BANDS" in line:
                printlog("Warning! TOO FEW BANDS!!!\n\n\nWarning! TOO FEW BANDS!!!\n")



            #Check W(q)
            if 'operators is LMAX' in line:
                lmax = int(line.split()[7])
                # print 'lmax', lmax
            if "W(low)/X(q)" in line:
                kk = 1; 
                low = []; 
                high = [];
                
                while kk < 100:
                    if 'Optimization' in outcarlines[i_line + kk] or len(outcarlines[i_line + kk].split() ) != 7: 
                        break
                    if 'PSMAXN' in outcarlines[i_line + kk]:
                        # print(line)
                        printlog('Warning! PSMAXN for non-local potential too small')
                        break
                    # print( 'line', outcarlines[i_line + kk])

                    low.append(  float(outcarlines[i_line + kk].split()[4]) )
                    high.append( float(outcarlines[i_line + kk].split()[5]) )
                    kk+=1


                if any(v > 1e-3 for v in low+high):
                    printlog("W(q)/X(q) are too high, check output!\n", 'Y')
                    printlog('Low + high = ', low+high, imp = 'Y' )
                    printlog([v > 1e-3 for v in low+high], imp = 'Y' )
            
            if "direct lattice vectors" in line:
                if not contcar_read:
                    for v in 0,1,2:
                        line = outcarlines[i_line+1+v]
                        line = line.replace('-', ' -')
                        # print(line)
                        self.end.rprimd[v] = np.asarray( [float(ri) for ri in line.split()[0:3]   ] )


                #print self.end.rprimd
                #print self.rprimd
            if "POSITION" in line:
                # if not contcar_exist or out_type == 'xcarts':
                if not contcar_read or out_type == 'xcarts':
                    local_xcart = []
                    for i in range(self.end.natom):
                        #print outcarlines[i_line+1+i].split()[0:3] 
                        xcart = np.asarray ( 
                                    [   float(x) for x in outcarlines[i_line+2+i].split()[0:3]   ] 
                                )
                        
                        local_xcart.append( xcart )

                    self.end.xcart = local_xcart

            
                    if out_type == 'xcarts':
                        self.end.list_xcart.append(local_xcart) #xcart at each step only for dimer

                        #the change of typat is accounted below

            if "number of electron " in line:
                # print line
                # print line.split()[-1]
                try:
                    self.magn1 = float(line.split()[-1])
                except:
                    self.magn1 = 0

            if "augmentation part " in line:
                try:
                    self.magn2 = float(line.split()[-1])
                except:
                    self.magn2 = 0


            if force_keyword in line:
                # Calculate forces here...
                forces = []
                magnitudes = []

                # print(len(self.end.select), self.end.natom)
                try:
                    for j in range(self.end.natom):
                        parts = outcarlines[i_line+j+2].split()
                        # print ("parts", parts, j)
                        # sys.exit()
                        if hasattr(self.end, 'select') and self.end.select:
                            # print(float(parts[ff[0]]), self.end.select[j][0])
                            b = []
                            # print (self.end.select)
                            for kkk in 0,1,2:
                                cur = self.end.select[j][kkk]
                                # print(cur)
                                
                                if cur == False:# or 'F' in cur:
                                    b.append(0)
                                elif cur == True:# or 'T' in cur:
                                    b.append(1)
                                else:
                                    b.append(cur)
                            # if parts[ff[0]].isdigit():
                            x = float(parts[ff[0]]) * b[0]
                            y = float(parts[ff[1]]) * b[1]
                            z = float(parts[ff[2]]) * b[2]
                        else:
                            x = float(parts[ff[0]])
                            y = float(parts[ff[1]])
                            z = float(parts[ff[2]])
                        
                        
                        forces.append([x,y,z])
                        magnitudes.append(math.sqrt(x*x + y*y + z*z))
                except:
                    printlog('Warning! Problem with forces on ionic step', iterat)
                # print('new step:')
                # for f, s in zip(forces, self.end.select):
                #     print('{:5.2f} {:5.2f} {:5.2f} {}'.format(*f, s))
                # sys.exit()
                average.append( red_prec( sum(magnitudes)/self.end.natom * 1000 ) )
                imax = np.asarray(magnitudes).argmax()
                maxforce.append( [imax, round(magnitudes[imax] * 1000)]  )
                # mforce.append( round(magnitudes[imax] * 1000))
                
           

            #Check total drift
            if "total drift:" in line:
                #print line
                tdrift = [float(d) for d in line.split()[2:5]]
                #if any(d > 0.001 and d > max(magnitudes) for d in tdrift):
                    #printlog("Total drift is too high = "+str(tdrift)+", check output!\n")
                    #pass


            if "g(Stress)" in line:
                #print line
                gstress.append( round( float(line.split()[4])*1000 *100, 3 )  )
            #if "Total" in line:
                #gstress.append( red_prec(float(line.split()[4])*1000 *100, 1000 )  )
            if "volume of cell" in line:
                try:                     
                    self.end.vol = float(line.split()[4])
                except ValueError: 
                    printlog("Warning! Cant read volume in calc "+self.name+"\n")
                #print self.vol      

            if "generate k-points for:" in line: 
                self.ngkpt = tuple(  [int(n) for n in line.split()[3:]]  )
                #print self.set.ngkpt

              # Kohn-Sham hamiltonian: http://en.wikipedia.org/wiki/Kohn%E2%80%93Sham_equations
              #kinetic energy
              #+ the external potential + the exchange-correlation energy +
              #+ Hartree (or Coulomb) energy
            # print line
            
            if  "alpha Z        PSCENC" in line:
                # print line
                self.energy.alpha = float(line.split()[-1]) # the electrostatic interaction of the ions in a compensating electron gas.

            if  "Ewald energy   TEWEN" in line:
                self.energy.ewald = float(line.split()[-1]) # the electrostatic interaction of the ions in a compensating electron gas.
                # print self.energy.ewald
            if  "-1/2 Hartree   DENC" in line or "-Hartree energ DENC" in line:
                self.energy.hartree = float(line.split()[-1]) #Coulomb electron-electron energy
                # print self.energy.hartree
            if  "-V(xc)+E(xc)   XCENC" in line:
                self.energy.xc = float(line.split()[-1]) # Kohn-Sham exchange-correlation energy
            if  "PAW double counting" in line:
                self.energy.pawdc1 = float(line.split()[-2]) #
                self.energy.pawdc2 = float(line.split()[-1]) #
            if  "eigenvalues    EBANDS" in line:
                try:
                    self.energy.bands = float(line.split()[-1]) # - Kohn Sham eigenvalues - include kinetic energy , but not exactly
                except:
                    self.energy.bands = 0
            if  "atomic energy  EATOM" in line:
                self.energy.atomic = float(line.split()[-1]) #energy of atoms in the box



            if "energy  without entropy=" in line:
                #self.energy = float(line.split()[4])
                self.e_without_entr = float(line.split()[3]) #
                self.energy_sigma0 = float(line.split()[6]) #energy(sigma->0)
                self.e0 = self.energy_sigma0
                self.list_e_sigma0.append(  self.energy_sigma0  )
                self.list_e_without_entr.append(  self.e_without_entr  )

                de_each_md_list.append(de_each_md)


            if "energy without entropy =" in line:
                e_sig0_prev = e_sig0
                try:
                    e_sig0 = float(line.split()[7])
                except:
                    e_sig0 = 0
                de_each_md = e_sig0_prev - e_sig0
                self.list_e_conv.append(e_sig0)

            if "free  energy   TOTEN  =" in line:
                #self.energy = float(line.split()[4])
                self.energy_free = float(line.split()[4]) #F free energy
            




            if re_lengths.search(line):
                self.vlength = [red_prec( float(l),1000 ) for l in outcarlines[i_line + 1].split()[0:3]]
                #print self.vlength
            if "in kB" in line:
                # print(line)
                # try:
                line = line.replace('-', ' -')
                # print(line)
                lines_str = line.split()[2:]
                try:
                    self.stress = [float(i)*100 for i in lines_str]  # in MPa 
                except:
                    printlog('Warning! Some problem with *in kB* line of OUTCAR', imp = 'y')
                    printlog(line, imp = 'y')
                    self.stress = [0,0,0]
                # if '*' in line:
                    # self.stress = [0,0,0] # problem with stresses
                # else:
                # except:
            if "Total  " in line:
                # print(line)
                line = line.replace('-', ' -')
                try:
                    self.intstress = [int(float(i)*1000) for i in line.split()[1:]] #stress in internal units; can be regarded as forces
                except:
                    self.intstress = [0,0,0]
                    printlog('Warning! Some problem with *Total * line of OUTCAR')

            if "external pressure =" in line: 
                #print iterat
                self.extpress = float(line.split()[3]) * 100 # in MPa 
                if self.mdstep == 1 : self.extpress_init = self.extpress

            if "E-fermi :" in line: 
                # print line
                self.efermi = float(line.split()[2]) # in eV


            if "Elapsed time" in line:
                self.time = float(line.split()[3])
            
            if "Maximum memory used (kb):" in line:
                ''
                if hasattr(self, 'corenum'):
                    self.memory_max = float(line.split()[-1]) * self.corenum / 1024 / 1024 
                else:
                    self.memory_max = 0
            
            if "total amount of memory" in line:
                ''
                # self.memory = float(line.split()[-2])   * self.corenum / 1024 / 1024          



            if re_nkpts.search(line):
                self.NKPTS = int(line.split()[3])
            else:
                self.NKPTS = 0
            if "WARNING" in line:
                warnings += 1#line


            if "Subroutine DYNSYM returns" in line and not nsgroup:
                nsgroup = line.split()[4]#number of space group operations
            # if nsgroup == None:
            if "Subroutine GETGRP returns:" in line and not nsgroup:
                nsgroup = line.split()[4]    


            if "Iteration" in line:
                self.mdstep = int(line.split('(')[0].split()[2].strip())
                iterat +=1
                # print self.mdstep
                # print line
                if mdstep_old != self.mdstep:
                    nscflist.append( niter ) # add to list number of scf iterations during mdstep_old
                niter = int(line.split(')')[0].split('(')[-1].strip()) #number of scf iterations
                mdstep_old = self.mdstep


            if 'number of electron ' in line:
                # print (line)
                try:
                    self.mag_sum.append( [float(line.split()[5]), 0]) # magnetization at each relaxation step
                except:
                    pass

            if 'augmentation part' in line:
                # print (line)
                try:
                    self.mag_sum[-1][1]= float(line.split()[-1])
                except:
                    pass

            if 'total charge ' in line:
                chg = []
                for j in range(self.end.natom):
                    chg.append( float(outcarlines[i_line+j+4].split()[-1]) )
                
                tot_chg_by_atoms.append(np.array(chg))#[ifmaglist])                    


            if 'magnetization (x)' in line:
                # print(line)
                if ifmaglist is not None:
                    mags = []
                    for j in range(self.end.natom):
                        mags.append( float(outcarlines[i_line+j+4].split()[-1]) )
                    
                    tot_mag_by_atoms.append(np.array(mags))#[ifmaglist])
                    # print(ifmaglist)
                    tot_mag_by_mag_atoms.append(np.array(mags)[ifmaglist])
                # print tot_mag_by_atoms
                # magnetic_elements
                # ifmaglist
                # self.tot_mag_by_atoms = tot_mag_by_atoms



            if 'LDAUU' in line:
                ldauu = line


            if 'onsite density matrix' in line:
                i_at = int( outcarlines[i_line-2].split()[2]  ) #starting from one
                l_at = int( outcarlines[i_line-2].split()[8]  )
                # print (spin_polarized)
                spin1 = []
                spin2 = []
                nm = 2*l_at+1
                for i in range(nm):
                    try:
                        line = outcarlines[i_line+4+i]
                        spin1.append( np.array(line.split()).astype(float) )
                    except:
                        printlog('Warning! Somthing wrong with occ matrix:', line)
                if spin_polarized:
                    for i in range(nm):
                        try:
                            line = outcarlines[i_line+7+nm+i]
                            # print(line)
                            line = line.replace('-', ' -')
                            spin2.append( np.array(line.split()).astype(float) )
                        except:
                            printlog('Attention! Could not read spin2, probably not finishied')
                            spin2.append(0)        

                occ_matrices[i_at-1] = spin1+spin2
                # print (np.array(spin1) )


            if 'freq' in show:

                if 'Eigenvectors and eigenvalues of the dynamical matrix' in line:
                    freq = []

                    i = 0
                    while 'ELASTIC MODULI CONTR FROM IONIC RELAXATION' not in line:
                        i+=1
                        line = outcarlines[i_line+i]
                        if 'f  =' in line:
                            freq.append(float(line.split()[3]) ) #THz
                            # print(line)


            if 'TOTAL ELASTIC MODULI' in line:
                eltensor = []
                for i in range(9):
                    line = outcarlines[i_line+i]
                    # print(line.strip())
                    if i > 2:
                        eltensor.append([float(c)/10 for c in line.split()[1:]]) #GPa

                eltensor = np.asarray(eltensor)

                printlog('Elastic tensor, GPa:', imp = 'y')

                np.set_printoptions(formatter={'float_kind':"{:6.1f}".format})
                print(eltensor)
                w, v = np.linalg.eig(eltensor)
                printlog('Eigenvalues are:', w, imp = 'y')
                        # eltensor

            if 'average eigenvalue GAMMA=' in line:
                # print(line)
                gamma = float(line.split()[-1])
                if gamma > 1 and 'conv' in show:
                    printlog('average eigenvalue GAMMA >1', gamma, imp = 'y')
                # sys.exit()
            if 'EDIFF  =' in line:
                ediff = float(line.split()[2])


            # if 'DIPCOR: dipole corrections for dipol' in line:
            if self.mdstep > mdstep_prev:
                # print(self.mdstep, dipol)
                mdstep_prev = self.mdstep

            if 'dipolmoment' in line:
                dipol = line.split()[1:4]
                self.dipol = [float(d) for d in dipol]
                # print(line)


            if 'min pos' in line:
                dipole_min_pos = line.split()[-1][:-1]
                self.dipole_min_pos.append(int(dipole_min_pos))

            if 'dipol+quadrupol energy correction' in line:
                e_dq_cor = line.split()[3]
                # print(line)
                self.e_dipol_quadrupol_cor.append(float(e_dq_cor))
            
            if 'added-field ion interaction' in line:
                e_af_ii = line.split()[3]
                self.e_added_field_ion.append(float(e_af_ii))

                # print(line)


                # for i in range(1,4):
                #     line = outcarlines[i_line+i]
                #     print(line)



            # if 'irreducible k-points:': in line:
            #     self.nkpt = int(line.split()[1])




            i_line += 1
        # sys.exit()
        #Check total drift
        






    try:
        toldfe = self.set.toldfe  # eV
    except:
        if ediff:
            toldfe = ediff # 
        else:
            toldfe = 0.000000000001



    if len(magnitudes) > 0:
        max_magnitude = max(magnitudes)
    else:
        max_magnitude = 0
    
    try:
        max_tdrift    = max(tdrift)
    except:
        max_tdrift = 0

    self.force_prefix = force_prefix
    self.maxforce_list = maxforce
    self.average_list = average

    try:
        self.maxforce = maxforce[-1][1]
    except:
        self.maxforce = 0
    # if max_magnitude < self.set.toldff/10: max_magnitude = self.set.toldff
    # print 'magn', magnitudes
    # print 'totdr', tdrift
    # print 'max_magnitude', max_magnitude
    try: 
        
        if max_magnitude < self.set.tolmxf: 
            max_magnitude = self.set.tolmxf
    except:
        ''

    #if any(d > 0.001 and d > max_magnitude for d in tdrift):
    if max_tdrift > 0.001 and max_tdrift > max_magnitude:
        
        printlog( "Total drift is too high! At the end one component is {:2.1f} of the maximum force, check output!\n".format(max_tdrift)  )
        pass
    #else: maxdrift = 
    # print magn
    if tot_mag_by_atoms:
        self.end.magmom = tot_mag_by_atoms[-1].tolist()

    """update xred"""
    self.end.update_xred()

    # print(self.init.zval)
    # self.end.zval = self.init.zval








    #print "init pressure = ",self.extpress_init,"; final pressure =",self.extpress
    #print self.end.xred
    #self.vol = np.dot( self.rprimd[0], np.cross(self.rprimd[1], self.rprimd[2])  ); #volume
    nscflist.append( niter ) # add to list number of scf iterations during mdstep_old
    #print "Stress:", self.stress
    v = self.vlength
    self.end.vlength = self.vlength

    s = self.stress
    yznormal = np.cross(self.init.rprimd[1], self.init.rprimd[2])
    #print yznormal
    #print np.cross( yznormal, np.array([1,0,0]) )
    if not hasattr(self.init, 'gbpos'):
        self.init.gbpos = None#for compatability

    self.gbpos = self.init.gbpos #for compatability
    if self.gbpos:
        if any( np.cross( yznormal, np.array([1,0,0]) ) ) != 0: 
            printlog("Warning! The normal to yz is not parallel to x. Take care of gb area\n")
    self.end.yzarea = np.linalg.norm( yznormal )  #It is assumed, that boundary is perpendicular to x





    """Calculate voronoi volume"""
    # print hasattr(self, 'vorovol')
    voro = ''
    if voronoi:# and not hasattr(self, 'vorovol'):#out_type == 'e_seg':
        voro = calculate_voronoi(self)
        calculate_voronoi(self, state = 'init')


    #deal with ldauu
    u_hubbard = 0
    if ldauu: 
        ldauu = list(np.array(ldauu.split()[7:]).astype(float))
        # print (ldauu)
        #find first non-zero
        self.ldauu = ldauu
        u_hubbard = ( next((u for u in ldauu if u), 0) )
        # print ( np.unique(ldauu)  )
    else:
        self.ldauu = [0]

    #Check if energy is converged relative to relaxation
    e_diff_md = self.energy_sigma0
    if len(self.list_e_sigma0) > 2:
        e_diff_md = (self.list_e_sigma0[-1] - self.list_e_sigma0[-2])*1000 #meV

    e_diff = (e_sig0_prev - e_sig0)*1000 #meV
    # print(e_diff)
    self.e_diff = e_diff #convergence
    if abs(e_diff) > float(toldfe)*1000:
        toldfe_warning = '!'
        printlog("Attention!, SCF was not converged to desirable prec", 
            round(e_diff,3), '>', toldfe*1000, 'meV', imp = 'y')
    else:
        toldfe_warning = ''

    if 'conv' in show:
        for i, de in enumerate(de_each_md_list ):
            if de/toldfe > 1.01:
                printlog('Attention! bad SCF convergence {:6.1g} eV for MD step {:}; toldfe = {:6.0g} eV'.format(de, i+1, toldfe))
        self.plot_energy_conv()

    #  Construct beatifull table
    #self.a1 = float(v[0])/2 ; self.a2 = float(v[1])/2/math.sqrt(0.75); self.c = float(v[2])  # o1b
    
    try:
        self.a = self.hex_a ; self.c = self.hex_c  # c1b
    except AttributeError:
        self.a  = 0; self.c = 0 #calculations with full relaxation
    if self.a == None or self.a == [None]:
        self.a  = 0; self.c = 0

    j = (35,12,7,7,8,9,14,5,5,20,5,20,8,12,20,8,5,8,8,25,8,4,3)

    d = "|"
    name = ("%s.%s.%s" % (self.id[0],self.id[1], str(self.id[2]) )).ljust(j[0])
    etot = ("%.4f" % ( self.energy_sigma0 )).center(j[1])
    etot1 = ("%.4f" % ( self.energy_sigma0/self.end.natom )).center(j[1])
    # print self.a
    a = ("%.4f" %      ( self.a )      ).center(j[2])
    c = ("%.4f" %      ( self.c )      ).center(j[3])
    time = ("%.3f" % (self.time/3600.)    ).center(j[4])
    itertm = ("%.1f" % (self.time/1./iterat)    ).center(j[5])
    Nmd = ("%1i,%2i,%3i" % (self.mdstep, iterat/self.mdstep, iterat)    ).center(j[6])
    self.iterat = iterat
    War = ("%i" % (warnings)    ).center(j[7])
    #nbands = ("%i" % (self.set.vasp_params["NBANDS"])    ).center(j[8])
    #added = ("%.0f" % ( (self.set.add_nbands - 1) * 100 )    ).center(j[15])
    try:
        kmesh = ("%s" % (str(self.ngkpt) )    ).center(j[8])
        ks = self.calc_kspacings()
        kspacing = ("[%.2f,%.2f,%.2f]" % ( ks[0], ks[1], ks[2] )    ).center(j[9])
        ks1 = ("[%.2f]" % ( ks[0] )    ).center(j[9])
    except:
        kmesh = ''
        ks    = ''
        kspacing = ''
        ks1     = ''

    nkpt = ("%i" % ( self.NKPTS)     ).center(j[10])
    if self.stress:
        istrs = ("[%5i,%5i,%5i] " % ( self.intstress[0],self.intstress[1],self.intstress[2]  )    ).center(j[11])
        strs = ("%.0f,%.0f,%.0f " % ( self.stress[0],self.stress[1],self.stress[2]  )    ).center(j[11])   
        eprs = ("%.0f" % (self.extpress)).center(j[12])
    
    else:
        istrs = ''
        strs  = ''
        eprs =  ''
    try:
        tsm = ("%.0f" % (self.set.tsmear*1000)).center(j[13])
    except:
        tsm = ''

    entrr = ("%.3f" % (   (self.energy_free - self.energy_sigma0)/self.end.natom * 1000    )   ).center(j[14]) #entropy due to the use of smearing

    try:
        npar = ("%i" % (self.set.vasp_params["NPAR"])).center(j[16])
        lpl = ("%s" % (self.set.vasp_params["LPLANE"])).center(j[17])
        ecut = ("%s" % (self.set.vasp_params["ENCUT"]) ).center(j[18]) 
    except:
        npar = ''
        lpl  = ''
        ecut = ''

    # lens = ("%.2f;%.2f;%.2f" % (v[0],v[1],v[2] ) ).center(j[19])
    lens = "{:4.2f}, {:4.2f}, {:4.2f}".format(v[0],v[1],v[2] ) 
    r1 = ("%.2f" % ( v[0] ) ).center(j[19])            
    vol = ("%.1f" % ( self.end.vol ) ).center(j[20])
    nat = ("%i" % ( self.end.natom ) ).center(j[21])
    try:
        totd = ("%.0f" % (   max_tdrift/max_magnitude * 100      ) ).center(j[22])
    except:
        totd = ''

    nsg = ("%s" % (     nsgroup     ) ).center(j[22])
    Uhu   = " {:3.1f} ".format(u_hubbard)
    ed    = ' {:3.0f}'.format( e_diff)
    edg   = ' {:3.1f} '.format( e_diff_md)
    spg   = ' {:4s} '.format( self.end.sg(silent = 1)[0])
    """Warning! concentrations are calculated correctly only for cells with one impurity atom"""
    #gbcon = ("%.3f" % (     1./self.end.yzarea      ) ).center(j[23]) # surface concentation at GB A-2
    #bcon = ("%.1f" % (     1./self.natom * 100      ) ).center(j[24]) # volume atomic concentration, %
    #outstring_nbands = name+d+Etot+d+a+d+c+d+time+d+itertm+d+Nmd+d+War+d+nbands+d+added+"\\\\"
    #outstring_npar = name+d+Etot+d+a+d+c+d+time+d+itertm+d+Nmd+d+War+d+npar+d+lpl+"\\\\"

    #outstring_stress = name+d+Etot+d+a+d+c+d+time+d+itertm+d+Nmd+d+War+d+istrs+d+eprs

    #outstring_kp_ec = name+d+Etot+d+a+d+c+d+time+d+itertm+d+Nmd+d+War+d+strs+d+eprs+d+kmesh+d+ecut+"\\\\"
    outst_ecut= etot+d+a+d+c+                                            d+time+d+itertm+d+Nmd+d+War+d+ecut+"\\\\"
    outst_kp  = etot+d+a+d+c+                                            d+time+d+itertm+d+Nmd+d+War+d+kmesh+d+kspacing+d+nkpt+"\\\\"            
    outst_ts  = etot+d+a+d+c+                                            d+time+d+itertm+d+Nmd+d+War+d+kmesh+d+tsm+d+entrr+"\\\\"

    outst_all = voro+etot+d+a+d+c+d+lens+d+vol+d+kspacing+d+strs+d+eprs+d+nat+d+time+d+Nmd+d+War+d+totd+d+nsg+"\\\\"
    outst_seg = voro+etot+d+        lens+d+vol+d+ks1     +d+strs+d+eprs+d+nat+d+time+d+Nmd+d+War+d+totd+d+nsg+"\\\\" #for segregation
    outst_coseg=voro+etot+d+                                strs+d+eprs+d+nat+d+time+d+Nmd+d+War+d+totd+d+nsg+"\\\\" #for co-segregation; 
    outst_gbe = voro+etot+               d+vol+d+kspacing+d+strs+d+eprs+d+nat+d+time+d+Nmd+d+War+d+nsg+"\\\\" # For comparing gb energies and volume
    outst_imp = voro+etot+d+a+d+c+d+lens+d+vol+d+kspacing+d+       eprs+d+nat+d+time+d+Nmd+d+War+d+totd+d+nsg+"\\\\" # For comparing impurity energies
    
    outst_cathode = d.join([spg,etot, etot1, lens, vol,nkpt, strs, nat, time, Nmd, War, nsg, Uhu, ed, edg ])
    if 'help' in show:
        print('Name', 'Space group', 'Tot. energy, eV', 'Energy, eV/at', 'Vector lengths, A', 
        'Volume, A^3', 'number of k-points', 'Stresses, MPa', 'Number of atoms', 'time, h', 'Number of ionic steps; Number of electronic states per MD; Total number of electronic steps',
        'Number of warnings', 'Number of symmetry operations', 'U value', 'Electronic convergence, meV', 'Ionic steps convergence, meV' )
    # print self.end.xred[-1]
    #print outstring_kp_ec
    # print show
    # print 'force' in show


    if 'conv' in show:
        # print('asdf', de_each_md_list)
        # show achived convergence every step with respect to toldfe, should be less than 1
        # np.set_printoptions(linewidth=150, formatter={'float':lambda x: "%3.0f->" % x}) #precision=1,
        np.set_printoptions(precision=0, linewidth=150, )
        printlog('Conv each step, de/toldfe (toldfe = {:.0g} eV) =  \n{:};'.format(toldfe, np.array([de/toldfe for de in de_each_md_list ])), imp = 'Y')
    
    if 'time' in show:
        print('Time is {:.1f} s'.format(self.time))
        print('Time is {:.1f} s/it'.format(self.time/1./iterat))




    if 'sur' in show:
        self.sumAO = {}
        self.devAO = {}
        for el in 'Li', 'Na', 'Fe', 'O':
            if el in self.end.get_elements():
                pos  = determine_symmetry_positions(self.end, el)
                # print(pos)
                # sys.exit()
                # xc = self.end.xcart[pos[0][0]]
                for ps in pos:
                    print('position', ps[0])
                    xc = self.end.xcart[ps[0]]

                    if el == 'O':
                        neib = 6
                    else:
                        neib = 6
                    sumAO = local_surrounding(xc, self.end, neib, periodic = True, only_elements = [8, 9], control = 'av')#[0]
                    self.devAO[el+'-O'] = local_surrounding(xc, self.end, neib, periodic = True, only_elements = [8, 9], control = 'av_dev')[0]
                    print('d_av '+el+'-O:',sumAO )
                    print('dev_av '+el+'-O:',self.devAO[el+'-O'] )
                    AO = local_surrounding(xc, self.end, neib, periodic = True, only_elements = [8, 9], control = 'mavm')
                    print('d_min, d_avex, d_max: {:4.2f}, {:4.2f}, {:4.2f}'.format(*AO))




                    self.sumAO[el+'-O'] = sumAO
                    if self.id[2] in [1,12]:
                        self.end.write_xyz(show_around_x = xc, nnumber = neib, filename = self.end.name+'_'+el+'-OF'+str(neib), analysis = 'imp_surrounding', only_elements = [8, 9])

                    if el in ['Li', 'Na']:
                        neib = 2
                        sumAO = local_surrounding(xc, self.end, neib, periodic = True, only_elements = [26,], control = 'av')#[0]
                        if self.id[2] in [1,12]:
                            self.end.write_xyz(show_around_x = xc, nnumber = neib, filename = self.end.name+'_'+el+'-Fe'+str(neib), analysis = 'imp_surrounding', only_elements = [26])
                        
                        print(el+'-Fe',sumAO )
                        self.sumAO[el+'-Fe'] = sumAO






    if 'en' in show:
        # energy - max force

        self.plot_energy_force()

    if 'efav' in show:
        # energy - average force
        self.plot_energy_force(force_type = 'av')

    if 'est' in show: # e step
        self.plot_energy_step ()

    if 'smag' in show:
        # printlog('{:s}'.format([round(m) for m in self.mag_sum]), imp = 'Y' )
        printlog(np.array(self.mag_sum).round(2), imp = 'Y' )

    chosen_ion = None
    if 'maga' in show or 'occ' in show:
        from siman.analysis import around_alkali
        i_alk = self.end.get_specific_elements([3,11,19])
        numb, dist, chosen_ion = around_alkali(self.end, 4, i_alk[0])
        
        #probably not used anymore
        # dist_dic = {}
        # self.dist_numb = zip(dist, numb)
        # for d, n in self.dist_numb:
        #     dist_dic[n] = d 
        #probably not used anymore


    if 'mag' in show and tot_mag_by_atoms:
        # print ('\n\n\n')
        # printlog
        # print 'Final mag moments for atoms:'
        # print np.arange(self.end.natom)[ifmaglist]+1
        # print np.array(tot_mag_by_atoms)

        # print (tot_mag_by_atoms)
        # if tot_mag_by_atoms:
        # print ('first step ', tot_mag_by_atoms[0][numb].round(3) )
        # print ('first step all ', tot_mag_by_atoms[0][ifmaglist].round(3) )
        # for mag in tot_mag_by_atoms:
        #     print ('  -', mag[numb].round(3) )



        # print ('last  step ', tot_mag_by_atoms[-1][numb].round(3), tot_chg_by_atoms[-1][numb].round(3) )
        # mmm = tot_mag_by_atoms[-1][numb].round(3)

        # print ('atom:mag  = ', ', '.join('{}:{:4.2f}'.format(iat, m) for iat, m  in zip(  numb+1, mmm   )) )
        self.gmt()
        print ('\n')


        if 'a' in show:
            ''
            # print ('last  step all', tot_mag_by_atoms[-1][ifmaglist].round(3) )

            # sys.exit()
        if chosen_ion:
            printlog ('Dist from 1st found alkali ion ',element_name_inv( chosen_ion[1]),
                ' to sur. transition met atoms: (Use *alkali_ion_number* to choose ion manually)')
            print ('atom:dist = ', 
            ', '.join('{}:{:.2f}'.format(iat, d) for iat, d  in zip(  numb+1, dist   )  ) )

        # plt.plot(np.array(sur[3]).round(2), tot_mag_by_atoms[-1][numb]) mag vs dist for last step
        
        # print ('Moments on all mag atoms:\n', tot_mag_by_atoms[-1][ifmaglist].round(3))
        if 'magp' in show:
            plt.plot(np.array(tot_mag_by_mag_atoms)) # magnetization vs md step
            plt.show()
            plt.clf()

    if 'chg' in show:
        self.tot_chg_by_atoms = tot_chg_by_atoms[-1] #save for last step only
        # print(list(zip(self.end.get_elements(), self.tot_chg_by_atoms)))
        els  = self.end.get_elements()
        try:
            only_el = show.split('.')[1:]
        except:
            only_el = None
        print('\nMulliken charges are:')
        for el, ch in zip(els, self.tot_chg_by_atoms):
            if only_el == None or (only_el and el in only_el):
                print('{:s} {:4.2f};'.format(el, ch), end = ' ')
        print()
    if 'occ' in show:
        ''
        # print (matrices)
        # print (df)
        if chosen_ion:
            printlog('Distances (A) from alkali ion #',chosen_ion[0]+1,' to transition atoms:', 
                ',  '.join([ '({:}<->{:}): {:.2f}'.format(chosen_ion[0]+1, iat, d) for d, iat in zip(  dist, numb+1  )  ]), imp = 'Y'  )
        
        show_occ_for_atoms = [int(n) for n in re.findall(r'\d+', show)]
        # print (show_occ_for_atom)
        # sys.exit()
        if show_occ_for_atoms:
            iat = show_occ_for_atoms[0]
            # dist_toi = dist_dic[iat]
            i_mag_at = iat
        else:
            i = 0
            # dist_toi = dist[i]
            i_mag_at = numb[1:][i]
        # print (st.znucl[st.typat[i_mag_at]-1] )
        # print(numb, i_mag_at)
        l05 = len(occ_matrices[i_mag_at])//2

        df = pd.DataFrame(occ_matrices[i_mag_at]).round(5)
        els  = self.end.get_elements()
        printlog( 'Occ. matrix for atom ', i_mag_at, els[i_mag_at], end = '\n', imp = 'Y'  )
            # ':  ; dist to alk ion is ',  dist_toi, 'A', end = '\n' )
        printlog('Spin 1:',end = '\n', imp = 'Y'  )
        printlog(tabulate(df[0:l05], headers = ['dxy', 'dyz', 'dz2', 'dxz', 'dx2-y2'], floatfmt=".1f", tablefmt='psql'),end = '\n', imp = 'Y'  )
        # print(' & '.join(['d_{xy}', 'd_{yz}', 'd_{z^2}', 'd_{xz}', 'd_{x^2-y^2}']))
        # printlog(tabulate(occ_matrices[i_mag_at][l05:], headers = ['d_{xy}', 'd_{yz}', 'd_{z^2}', 'd_{xz}', 'd_{x^2-y^2}'], floatfmt=".1f", tablefmt='latex'),end = '\n' )
        # print(tabulate(a, tablefmt="latex", floatfmt=".2f"))
        printlog('Spin 2:',end = '\n', imp = 'Y'  )
        printlog(tabulate(df[l05:], floatfmt=".1f", tablefmt='psql'), imp = 'Y'  )
    self.occ_matrices = occ_matrices
    


    if 'freq' in show:
        dos = [1]*len(freq)
        # from scipy.ndimage.filters import gaussian_filter
        from scipy.signal import butter, lfilter, freqz
        # blurred = gaussian_filter(freq, sigma=7)
        fmin = min(freq)
        fmax = max(freq)
        fw   = fmax-fmin

        finefreq = np.linspace(fmin, fmax, 1000)
        dos = [0]*1000

        # for i in range(1000):
        # print(fw)
        for f in freq:
            # print(f)
            i = int( np.round( (f-fmin)/ fw * 999 ,0) )
            dos[i] += 1
            # print(i, finefreq[i], f)
        

        def butter_lowpass(cutoff, fs, order=5):
            nyq = 0.5 * fs
            normal_cutoff = cutoff / nyq
            b, a = butter(order, normal_cutoff, btype='low', analog=False)
            return b, a

        def butter_lowpass_filter(data, cutoff, fs, order=5):
            b, a = butter_lowpass(cutoff, fs, order=order)
            y = lfilter(b, a, data)
            return y

        order = 1
        fs = 500.0       # sample rate, Hz
        cutoff = 5  # desired cutoff frequency of the filter, Hz

        y = butter_lowpass_filter(dos, cutoff, fs, order)

        # plt.plot(finefreq, smoother(smoother(dos,10), 10), '-') 
        plt.plot(finefreq, y, '-') 
        # plt.plot(finefreq, dos, '-')
        plt.xlabel('Frequency, THz')
        plt.ylabel('DOS' )
        # plt.plot(finefreq, y, '-') 
        filename = 'figs/'+str(self.id)+'.pdf'
        plt.savefig(filename)
        printlog('Freq file saved to ', filename, imp = 'y')
        plt.show()
        plt.clf()

    # sys.exit()


    printlog("Reading of results completed\n\n", imp = 'n')
    self.end.outfile = path_to_outcar
    

    if header.pymatgen_flag:
        ''
        # self.end.write_cif(os.path.join(self.dir,self.name))
    

    # print(out_type)
    # sys.exit()
    if   out_type == 'gbe'  : outst = outst_gbe
    elif out_type == 'e_imp': outst = outst_imp
    elif out_type == 'e_seg': outst = outst_seg            
    elif out_type == 'coseg': outst = outst_coseg            
    elif 'ecut' in out_type : outst = outst_ecut
    elif 'kp' in out_type   : outst = outst_kp
    elif 'ts' in out_type   : outst = outst_ts
    
    elif not header.siman_run:
        # if not hasattr(cl, 'name'):
            # cl.name = 'noname'

        outst_simple = '|'.join([cl.name, etot, lens, strs, Nmd])
        # print("Bi2Se3.static.1               |  -20.1543  |    10.27;10.27;10.27    | -680,-680,-657 |   1,13, 13   |    ")
        if header.show_head :
            printlog("                          |  energy(eV)|    Vector lenghts (A)   | Stresses (MPa)     | N MD, N SCF   ", end = '\n', imp = 'Y')
            header.show_head = False
        
        outst = outst_simple
    else: 
        printlog('Output type: outst_cathode')
        outst = outst_cathode
    #else: printlog("Uknown type of outstring\n")


    #save cif file
    return outst




def read_aims_out(cl, load = '', out_type = '', show = ''):
    ''
    with open(cl.path['output'], 'r') as f:
        for line in f:
            if 'Total energy corrected' in line:
                cl.e0 = float(line.split()[5])
                cl.energy_sigma0 = cl.e0
                # print(cl.e0)
            if 'Number of self-consistency cycles' in line: 
                cl.iterat = float(line.split()[6])

            if '| Total time   ' in line:
                # print(line)
                cl.time = float(line.split()[4])


    cl.end = cl.init

    etot = ' {:5.4f} eV'.format(cl.e0)
    time = ' {:5.1f} h'.format(cl.time/3600)
    itrt = ' {:n},{:n},{:n} '.format(1, cl.iterat/1,cl.iterat)

    outstr = '|'.join([etot, time, itrt])

    return outstr


def read_atat_fit_out(filename, filter_names = None, i_energy = 1):
    """
    read fit.out of atat
    filter_names - do not read name numbers greater than filter_name
    i_energy - 1 for fit.out, 2 for predstr.out
    return - concentration list, energy list
    """
    X, E, nam = [] ,[], []
    count = 0
    with open(filename, 'r') as f:
        for line in f:


                # continue

            vals = line.split()
            
                # name = int(vals[-1])

            try:
                name = int(vals[-1])
            except:
                name = 0
            if filter_names and name > filter_names:
                continue

            e = round(float(vals[i_energy]),4)
            x = float(vals[0])
            # print(e)
            # if e in E:
            #     if X[E.index(e)] == x:
            #         continue
            
            # if E and abs(min(np.array(E)-e)) < 0.01:
            #     print('skipping point, too often')
            #     continue
            # if count > 10:
            #     count = 0
            X.append(x) # concentrations
            E.append(e) # formation energies
            nam.append(name) # names of structures
            count+=1
        # print(a)
    return X, E






def read_structure(filename = None, format = None, object_type = None, silent = 0):
    """
    Smart wrapper for reading files with atomic structure. New version of smart_structure_read 


    INPUT:
        - filename (str)
        - format (str) - file format. Normally the format is determined automatically from the name or extension
            - 'abinit'
            - 'vasp'
            - 'cif'
            - 'gaussian'
        - object_type (str)
            - 'molecule' - based on pymatgen
            - 'structure' - siman internal Structure() class


    returns Structure() if file is found and read successfully otherwise None
    """

    search_templates =       {'abinit':'*.geo*', 'vasp':'*POSCAR*', 'vasp-phonopy': 'POSCAR*', 'cif':'*.cif'}

    input_geo_format = determine_file_format(filename)
    printlog(input_geo_format,' format is detected')



    file_exists = os.path.exists(filename)

    if not file_exists:
        if not silent:
            printlog('Warning! file', filename, 'does not exist, return None')
        return None

    if input_geo_format   == 'abinit':
        # cl.read_geometry(input_geo_file)
        ''
    
    elif input_geo_format == 'vasp':
        # cl.read_poscar(input_geo_file, version = curver)
        ''

    elif input_geo_format == 'cif':
        ''
        # cif2poscar(input_geo_file, input_geo_file.replace('.cif', '.POSCAR'))
        # input_geo_file = input_geo_file.replace('.cif', '.POSCAR')
        # cl.read_poscar(input_geo_file)

    elif input_geo_format == 'xyz':
        from siman.core.structure import Structure
        from siman.core.molecule import Molecule

        if object_type == 'structure' or object_type is None:
            st = Structure()
            st = read_xyz(st, filename)
        elif object_type == 'molecule':
            # st = Molecule()
            st = Molecule.from_file(filename)
            st.filename = filename
            st.name = os.path.basename(filename).split('.')[0]



    elif input_geo_format == 'gaussian':
        ''
        from pymatgen.io.gaussian import GaussianInput
        from siman.core.molecule import Molecule
        gaus = GaussianInput(mol = None)
        st = gaus.from_file(filename).molecule

        # st = Molecule(pymatgen_molecule_object = st)
        st = Molecule.cast(st, filename)
        # print(st)
        st.name = os.path.basename(filename).split('.')[0]
    else:
        printlog("Error! smart_structure_read(): File format", input_geo_format, "is unknown")


        # printlog("Error! Input file was not properly read for some reason")
    

    return st