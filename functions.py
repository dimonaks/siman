# from classes import res_loop
import header
from header import *

def replic(structure, mul = (1,1,1), inv = 1, only_atoms = None, cut_one_cell = None, include_boundary = (1,1) ):
    """
    Replicate structure() according to: mul[i]*rprimd[i]
    
    Input:
    structure - structure() type
    mul[] - is tuple of three integer numbers
    Use from structure:
    xcart, typat, rprimd, natom, xred
    inv - 1 or -1 allows to replicate in different directions
    
    inv = 0 - cell is replicated in both directions by mul[i];  2 still gives -1 0 1 but 3 gives -2 -1 0 1 2; for 'only_matrix' may work not correctly


    only_atoms - allows to replicate only specific atoms; now 'only_matrix'

    cut_one_cell - allows to cut only one cell with replicated edge atoms
    include_boundary (A) - the width of region to include additional edge atoms (bottom, up)


    Return:
    replicated structure
    """
    st = copy.deepcopy(structure)
    # print 'Structure has before replication', st.natom,'atoms' 

    if hasattr(st, 'init_numbers'):
        numbers = st.init_numbers
    else:
        numbers = range(st.natom)


    #determine maximum and minimum values before replication
    xmax = [-1000000]*3
    xmin = [+1000000]*3
    for x in st.xcart:
        for j in 0,1,2:
            if xmax[j] < x[j]: xmax[j] = x[j] 
            if xmin[j] > x[j]: xmin[j] = x[j] 
    # print 'xmax, xmin', xmax, xmin



    inv_loc = inv
    for i in 0, 1, 2:
        
        axis_mul = range(mul[i])
        
        if inv == 0: # does not work propely; mul below is also should be accounted
            axis_mul = range(-mul[i]+1, mul[i])
            print 'axis_mul', axis_mul
            inv_loc = 1
        
        if only_atoms == 'only_matrix':
            st.xcart += [ x + inv_loc*k*st.rprimd[i] for x, t in zip(st.xcart[:], st.typat[:]) for k in axis_mul[1:] if t == 1] # fill by axis i by blocks
            
            st.typat += [ t                  for t in st.typat[:] for k in axis_mul[1:] if t == 1]
        else:
            st.xcart = [ x + inv_loc*k*st.rprimd[i] for x in st.xcart[:] for k in axis_mul ] # fill by axis i by blocks
            st.typat = [ t                  for t in st.typat[:] for k in axis_mul ]
            numbers = [n for n in numbers[:] for k in axis_mul]
            # print numbers
            # assert len(st.xcart) == abs(st.natom * reduce(lambda x, y: x*y, mul) )


        # print 'before', st.rprimd[i]
        
    if cut_one_cell:
        pass        
    else:

        for i in 0, 1, 2:
            if inv == 0:
                k = 2 * mul[i] - 1
            else:
                k = mul[i]
            st.rprimd[i] = st.rprimd[i] * k
        

    st.init_numbers = numbers

        # print st.init_numbers 
        # print 'after', st.rprimd[i]
    # print len(st.xcart)
    # print st.natom * reduce(lambda x, y: x*y, mul)
    #print st.xcart, 

    st.xred = xcart2xred(st.xcart, st.rprimd)


    if cut_one_cell:
        new_xred = []
        new_xcart = []
        new_typat = []

        precb = include_boundary[0]#/max(st.rprimd[2])
        precu = include_boundary[1]#/max(st.rprimd[2])
        


        bob = 0 - precb; upb = 1. + precu;
        print bob, upb
        

        n = 0 
        # print st.xred
        for t, xr in  zip(st.typat, st.xcart):
            for j in 0,1,2:
                if (xr[j]  < xmin[j] - precb) or (xr[j]  > xmax[j] + precu): break  
            else:
                new_xcart.append(xr)
                new_typat.append(t)

        st.typat = new_typat
        st.xcart = new_xcart
        print 'After removing, cell has ', len(st.xred)
        print st.xred
        # st.xcart = xred2xcart(st.xred, st.rprimd)
        st.xred = xcart2xred(st.xcart, st.rprimd)


    st.natom = len(st.xcart)
    # print 'Structure is replicated; now', st.natom,'atoms' 
    return st



def local_surrounding(x_central, st, n_neighbours, control = 'sum', periodic = False):
    """
    Return list of distances to n closest atoms around central atom. (By defauld sum of distances)
    
    Input:
    x_central - cartesian coordinates of central atom; vector
    st - structure with xcart list of coordinates of all atoms in system
    n_neighbours - number of needed closest neighbours

    control - type of output; sum - sum of distances, av - av distance, list - list of distances; 
              av_dev - average deviation, maximum deviation from average distance in mA.
              atoms  - coordinates of neighbours

    periodic - if True, then cell is additionaly replicated; needed for small cells

    """
    if periodic:
        st = replic(st, mul = (2,2,2), inv = 1 ) # to be sure that impurity is surrounded by atoms
        st = replic(st, mul = (2,2,2), inv = -1 )

    xcart = st.xcart
    typat = st.typat
    natom = st.natom
    # print x_central

    #print len(xcart)

    dlist = [np.linalg.norm(x_central - x) for x in xcart]# if all (x != x_central)] # list of all distances
    dlist_unsort = copy.deepcopy(dlist)
    dlist.sort()
    # write_xyz(st  )
    # print dlist
    dlistnn = dlist[1:n_neighbours+1] #without first impurity which is x_central
    # print dlistnn
    # os._exit(1)

    if control == 'list':
        output = dlistnn

    elif control == 'sum':
        output = round(sum(dlistnn), 2)
    
    elif control == 'av':
        n_neighbours = float(n_neighbours)
        dav = sum(dlistnn)/n_neighbours
        output = round(dav, 2)
       
    elif control == 'av_dev':
        n_neighbours = float(n_neighbours)
        dav = sum(dlistnn)/n_neighbours
        av_dev = sum( [abs(d-dav) for d in dlistnn] ) / n_neighbours
        max_dev = max([abs(d-dav) for d in dlistnn])
        output = (av_dev*1000, max_dev*1000)

    elif control == 'atoms':
        # print dlist_unsort
        if hasattr(st, 'init_numbers'):
            numbers = st.init_numbers
        else:
            numbers = range(natom)
        temp = zip(dlist_unsort, xcart, typat, numbers)
        # print temp
        temp.sort(key = itemgetter(0))
        temp2 = zip(*temp)
        dlist       = temp2[0][:n_neighbours+1]
        xcart_local = temp2[1][:n_neighbours+1]
        typat_local = temp2[2][:n_neighbours+1]
        numbers     = temp2[3][:n_neighbours+1]
        # print temp2[0][:n_neighbours]
        # print xcart_local[:n_neighbours]
        output =  (xcart_local, typat_local, numbers, dlist )
    return output

def salary_inflation():
    """Calculate salary growth in Russia taking into account inflation"""
    inflation2000_2014 = [
     5.34,
     6.45,
     6.58, 
     6.10,
     8.78, 
     8.80,
     13.28,
     11.87,
     9.00 ,
     10.91,
     11.74,
     11.99,
     15.06,
     18.8,
     20.1]
    init_salary = 1500 # in jan 2000; other sources 2000 - very important
    for i,  l in enumerate( reversed(inflation2000_2014)  ):
        init_salary = (1+l/100)*init_salary
        print init_salary, i+2000

    salary2014 = 30000
    increase = salary2014/init_salary
    print increase

# salary_inflation()

def element_name_inv(el):
    if type(el) == str:
        if el == "C": elinv = 6 
        elif el == "O": elinv = 8 
        elif el == "Ti": elinv = 22 
        elif el == "B": elinv = 5 
        elif el == "H": elinv = 1 
        elif el == "octa": elinv = 200 
        else:
            print_and_log("Unknown element\n")
            raise RuntimeError
    else:
        el = int(el)
        if el == 6: elinv = 'C'
        elif el == 5: elinv = 'B'
        elif el == 8: elinv = 'O'
        elif el == 22: elinv = 'Ti'
        elif el == 1: elinv = 'H'
        else:
            print_and_log("Unknown element\n")
            raise RuntimeError
    return elinv # inversed notion of element


def xcart2xred(xcart, rprimd):
    """Convert from cartesian coordinates xcart to
        dimensionless reduced coordinates 
        Input: xcart - list of numpy arrays, rprimd - list of numpy arrays
        Output: xred - list of numpy arrays"""
    xred = []
    gprimd = np.asarray( np.matrix(rprimd).I.T ) #Transpose of the inverse matrix of rprimd
    #print gprimd
    for xc in xcart:
        xred.append(  np.dot( gprimd , xc)  ) #dot product
    #print xred
    return xred

def xred2xcart(xred, rprimd):
    """Convert from dimensionless reduced coordinates to
    cartesian coordinates xcart;
        Input: xred - list of numpy arrays, rprimd - list of numpy arrays
        Output: xcart - list of numpy arrays"""
    xcart = []
    #print "rprimd ", rprimd
    for xr in xred:
        #for j in 0,1,2:
        #    print xr[0] * rprimd[0][j] + xr[1] * rprimd[1][j] + xr[2] * rprimd[2][j],
        #print ""
        #print np.dot( xr, rprimd)
        xcart.append(  np.dot( xr, rprimd)  ) #dot product

    #print xred
    return xcart

def return_atoms_to_cell(st):

    bob = 0; upb = 1;
    n = 0 
    # print st.xred
    for xr in st.xred:
        for j in 0,1,2:
            if xr[j]  < bob:  xr[j] = xr[j] - int(xr[j]) + 1 #allows to account that xr can be more than 2
            if xr[j]  > upb:  xr[j] = xr[j] - int(xr[j])
    n+=1
    # zmin = 100
    # for xr in st.xred:
    #     if xr[2]<zmin: zmin = xr[2]
    # if zmin < 0:
    #     for xr in st.xred:
    #         xr[2] = xr[2]-zmin



    st.xcart = xred2xcart(st.xred, st.rprimd)

    print_and_log(str(n)+" atoms were returned to cell.\n")
    #print st.xred
    return st



def calc_ac(a1, c1, a2, c2, a_b = 0.1, c_b = 0.1, type = "two_atoms"):
    """
    Calculate values of hexagonal lattice parameters for cell with two different atoms.
    The used assumption is:
    1. Provided lattice constants are for large enougth cells, in which excess volume (dV) of impurity does not depend on the size of cell.
    2. Two atoms do not interact with each other, which allows to use dV(CO) = dV(C) + dV(O)
    
    Two regimes:
    two_atoms - calculate cell sizes if additional atom was added
    double_cell - if cell was doubled; only first cell and second_cell are needed


    Input:
    a1, c1 - lattice constants of cell with first impurity atom (first cell)
    a2, c2 - lattice constants of cell with second impurity atom (second cell)
    a_b, c_b - lattice constants of cell with pure hexagonal metall
    
    Output:
    a, c - lattice constants of cell with two atoms
    """
    hstring = ("%s    #on %s"% (traceback.extract_stack(None, 2)[0][3],   datetime.date.today() ) )
    if hstring != header.history[-1]: header.history.append( hstring  )
    
    A = (a1**2 * c1) + (a2**2 * c2) - (a_b**2 * c_b)
    B = 0.5 * (c1/a1 + c2/a2)
    C = ( (a1**2 * c1) + (a2**2 * c2) ) * 0.5 #sum of cell volumes divided by 2 since during the construction of new cell we will use multiplication by 2
    print "A,B=",A,B
    a = (A/B)**(1./3)
    c = a * B
    a = round(a,5)
    c = round(c,5)
    print "a, c, c/a for cell with pure    hcp ", a_b, c_b, round(c_b/a_b,4)
    print "a, c, c/a for cell with first  atom ", a1, c1, round(c1/a1,4)
    print "a, c, c/a for cell with second atom ", a2, c2, round(c2/a2,4)

    #for double cell
    a3 = (C/B)**(1./3)
    c3 = a3 * B
    a3 = round(a3,5)
    c3 = round(c3,5)    

    if type == "two_atoms":
        print "a, c, c/a for cell with two   atoms ", a,  c, round(c/a,4), "# the same cell but with two atoms\n"
    elif type == "double_cell":
        print "a, c, c/a for new cell              ", a3,  c3, round(c3/a3,4), "# for cell with V = V(first_cell) + V(second cell), but only for the case if V(second cell) == V(first_cell)"

    return a, c



def write_jmol(xyzfile, pngfile, scriptfile = None, atomselection = None, topview = False, orientation = None,
    axis = False, bonds = True, rprimd = None, shift = None, rotate = None,
    label = None, high_contrast = None, specialcommand = None,
    boundbox = 2):
    """
    atomselection - string in gmol format with number of atoms to be nrotateSelected
    topview - additional top view, requires two models in xyz
    orientation - additional rotation
    axis - add axes
    rotate - rotation of all atoms around view axis in degrees
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

        f.write('select all \ncpk 250 \nwireframe 0.3\n') #250
        f.write('background white \nselect Ti* \ncolor [20,120,250] \nselect C* \ncolor [80,80,80]\n cpk 100\n')
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
            f.write('select Be*\ncpk 200\nset labeloffset 0 0\nset labelfront\ncolor label black\nlabel %i\n font label 24 bold \n')


        if rotate:
            f.write('rotate z '+str(rotate)+'\n')

        if specialcommand:
            f.write(specialcommand)

        
        f.write('write image 2800 2800 png "'+pngfile+'"')
    print runBash(header.path_to_jmol+' -ions '+scriptfile)
    # print runBash('convert '+pngfile+' -shave 0x5% -trim '+pngfile) #cut by 5% from up and down (shave) and that trim left background
    print pngfile
    print runBash('convert '+pngfile+' -trim '+pngfile) # trim background
    return




def write_xyz(st, path = '', repeat = 1, shift = 1.0,  gbpos2 = None, gbwidth = 1 , 
    imp_positions = [], specialcommand = None, analysis = None, replications = None, nnumber = 6, topview = True,
    file_name = None, full_cell = False, orientation = None, boundbox = 2, withgb = False,
    include_boundary = 2, rotate = None,
    ):
    """Writes st structure in xyz format in the folder xyz/pat

    if repeat == 2: produces jmol script
    shift - in rprimd[1][1] - shift of the second view
    gbpos2 - position of grain boundary in A
    gbwidth - atoms aroung gbpos2 will be colored differently

    imp_positions - type and xcart coordinates additionally to be added to structure; to visulaze all impurity positions 
    specialcommand - any command at the end of script

    analysis - additional processing, allows to show only specifice atoms, 'imp_surrounding' - shows Ti atoms only around impurity
    replications - list of replications, (2,2,2) 
    nnumber - number of neighbours to show

    full_cell - returns atoms to cell and replicate boundary atoms
    """
    if replications:
        st = replic(st, mul = replications, inv = 1 )
  
    def update_var():
        return st.rprimd, st.xcart, st.xred, st.typat, st.znucl, len(st.xred)

    rprimd, xcart, xred, typat, znucl, natom = update_var()


    if file_name:
        name = file_name
    else:
        name = st.name



    if natom != len(xred) != len(xcart) != len(typat) or len(znucl) != max(typat): 
        print "Error! write_xyz: check your arrays.\n\n"    
    # print st.natom, len(st.xred), len(st.xcart), len(st.typat), len(st.znucl), max(st.typat)
    
    print "Name is", name
    if name == '': name = 'noname'


    if xcart == [] or len(xcart) != len(xred):
        print "Warining! write_xyz: len(xcart) != len(xred) making xcart from xred.\n"
        xcart = xred2xcart(xred, rprimd)
        #print xcart[1]
    if path == None:
        basepath = './'
    else:
        basepath = 'xyz/'+path+'/'



    xyzfile = basepath+name+".xyz"
    if not os.path.exists(os.path.dirname(xyzfile)):
        os.makedirs(os.path.dirname(xyzfile))


    """Processing section"""


    if analysis == 'imp_surrounding':
        lxcart = []
        ltypat = []
        for t, x in zip(typat, xcart):
            if t > 1: 
                # lxcart.append(x)
                # ltypat.append(t)
                # print x
                x_t = local_surrounding(x, st, nnumber, control = 'atoms', periodic = True)
                # print x_t[1]
                lxcart+=x_t[0]
                ltypat+=x_t[1]
        xcart = lxcart
        typat = ltypat
        natom = len(typat)



    """Include atoms on the edge of cell"""
    if full_cell:
        # print xred
        print natom
        # st = return_atoms_to_cell(st)
        # print xred
        st = replic(st, mul = (1,1,2), inv = 0, cut_one_cell = 1, include_boundary = include_boundary)
        print natom, st.natom

        # print st.xred

        rprimd, xcart, xred, typat, znucl, natom = update_var()
        

    """Writing section"""

    print_and_log("Writing xyz "+xyzfile+" \n")

    with open(xyzfile,'w') as f:
        for i in range(repeat):
            f.write(str(natom + len(imp_positions ))+"\n")
            f.write(name+"\n")

            if imp_positions: 
                for i, el in enumerate(imp_positions):
                    f.write( "%s %.5f %.5f %.5f \n"%( el[3], el[0], el[1], el[2] ) )
                    # print 'composite -pointsize 60 label:{0:d} -geometry +{1:d}+{2:d} 1.png 2.png'.format(i, el[0], el[1])


            for i in range(natom):
                typ = typat[i] - 1
                
                z = int ( znucl[ typ ] )
                #print "typ", znucl
                #print "z", z
                if   z == 22:  f.write( "Ti " )
                elif z ==  6:  f.write( "C " )
                elif z ==  8:  f.write( "O " )
                elif z ==  1:  f.write( "H " )
                else:          f.write( "Pu " )
                f.write( "%.5f %.5f %.5f \n"%( xcart[i][0], xcart[i][1], xcart[i][2] ) )




    # os._exit(1)



    if repeat == 2:
        """
        script mode for jmol. Create script file as well for elobarate visualization
        """
        



        """Choose gb atoms to change their color"""
        print 'position of boundary 2', gbpos2
        atomselection = ''
        if gbpos2:
            
            gbpos1 = gbpos2 - rprimd[0][0]/2.
            gbatoms = []
            
            for i, x in enumerate(xcart):
                # print i
                if x[0] > gbpos1 - gbwidth/2. and x[0] < gbpos1 + gbwidth/2.:
                    gbatoms.append(i)
                if x[0] > gbpos2 - gbwidth/2. and x[0] < gbpos2 + gbwidth/2.:
                    gbatoms.append(i)
            
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
        
        write_jmol(xyzfile, pngfile, scriptfile, atomselection, topview = topview, rprimd =rprimd, shift = shift, label = True, 
            specialcommand = specialcommand, orientation = orientation, boundbox =boundbox, rotate = rotate)


    return


def write_lammps(cl, state, filepath = ''):
    """Writes structure in lammps format """
    if state == 'end':
        st = cl.end
    elif state == 'init':
        st = cl.init
    
    rprimd = st.rprimd
    xcart  = st.xcart
    xred   = st.xred
    typat = st.typat
    ntypat = st.ntypat
    znucl = st.znucl
    name = st.name
    natom = st.natom

    if natom != len(xred) != len(xcart) != len(typat) or len(znucl) != max(typat): 
        print "Error! write_xyz: check your arrays.\n\n"    
    
    if name == '': name = 'noname'
    if xcart == [] or len(xcart) != len(xred):
        print "Warining! write_xyz: len(xcart) != len(xred) making xcart from xred.\n"
        xcart = xred2xcart(xred, rprimd)
        #print xcart[1]

    if not filepath:
        filepath = 'lammps/'+name+".lammps"
    if not os.path.exists(os.path.dirname(filepath)):
        os.makedirs(os.path.dirname(filepath))

    """Write lammps structure file;  """
    if 0:
        """ My version; valid only for octahedral cells""" 
        with open(filepath+'','w') as f:
            f.write("Lammps format "+name+'\n')
            f.write(str(natom)+" atoms\n")
            f.write(str(ntypat)+" atom types\n")
            f.write("0 "+str(rprimd[0][0])+" xlo xhi\n")
            f.write("0 "+str(rprimd[1][1])+" ylo yhi\n")
            f.write("0 "+str(rprimd[2][2])+" zlo zhi\n")
            f.write("0.000000    0.000000    0.000000   xy xz yz\n")
            f.write("\nAtoms\n\n")

            for i, x in enumerate(xcart):
                f.write("{0:d} {1:d} {2:f} {3:f} {4:f}\n".format(i+1, typat[i], x[0], x[1], x[2] ))
            f.write("\n")
    else:
        """Write poscar and convert from poscar to lammps using external script; Valid for arbitary cells"""
        cl.write_structure('POSCAR', 'dir', path = 'voronoi_analysis/', state = state)
        runBash("voronoi_analysis/VASP-poscar2lammps.awk voronoi_analysis/POSCAR > "+filepath)
    



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



def image_distance(x1, x2, r, order = 1):
    """
    Calculate smallest distance and the next smallest distance between two atoms 
    correctly treating periodic boundary conditions and oblique cells.
    x1, x2 - vector[3] coordinates of two atoms
    r - rprimd of cell
    order - the order of periodic images which are accounted in the calcualtion of distances between atoms.
    for cubic cells, order = 1 always provide correct result.
    For highly oblique cell you should test and find the needed value of 'order' after which results are the same.

    return d1, d2 - the smallest and next smallest distances between atoms

    """
    d = [] # list of distances between 1st atom and images of 2nd atom
    for i in range(-order, order+1):
        for j in range(-order, order+1):
            for k in range(-order, order+1):
                x2i = x2 + (r[0] * i + r[1] * j + r[2] * k) #determine coordinates of image of atom 2 in corresponding image cells
                d.append(   np.linalg.norm(x1 - x2i)   )
    d.sort()
    #print d
    assert d[0] == min(d)
    return d[0], d[1] #, math.sqrt(dxl[0]**2 + dxl[1]**2 + dxl[2]**2)





def read_charge_den_vasp():
    """
    Read CHG vasp file and return ChargeDen object 
    """
    class ChargeDen():
        """docstring for ChargeDen"""
        def __init__(self, ):
            # self.arg = arg
            
            pass




def rotation_matrix(axis,theta):
    axis = axis/math.sqrt(np.dot(axis,axis))
    a = math.cos(theta/2)
    b,c,d = -axis*math.sin(theta/2)
    return np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                     [2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
                     [2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-b*b-c*c]])

def rotate():
    v = np.array([3,5,0])
    axis = np.array([4,4,1])
    theta = 1.2 

    print(np.dot(rotation_matrix(axis,theta),v))             
    # [ 2.74911638  4.77180932  1.91629719]

def plot_charge_den():
    """Test function; Was not used"""
    from mpl_toolkits.mplot3d import axes3d
    import matplotlib.pyplot as plt
    from matplotlib import cm

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    X, Y, Z = axes3d.get_test_data(0.05)
    print X
    print Y
    print Z

    ax.plot_surface(X, Y, Z, rstride=8, cstride=8, alpha=0.3)
    # cset = ax.contourf(X, Y, Z, zdir='z', offset=-100, cmap=cm.coolwarm)
    # cset = ax.contourf(X, Y, Z, zdir='x', offset=-40, cmap=cm.coolwarm)
    # cset = ax.contourf(X, Y, Z, zdir='y', offset=40, cmap=cm.coolwarm)

    ax.set_xlabel('X')
    ax.set_xlim(-40, 40)
    ax.set_ylabel('Y')
    ax.set_ylim(-40, 40)
    ax.set_zlabel('Z')
    ax.set_zlim(-100, 100)

    plt.show()

    return



def plot_interaction(calclist, calc):
    """
    For calculation of interaction parameter alpha;
    Take in mind that this parameter is obtained under aproximation of redular solution
    """
    e_seg = []
    dX = []
    for id in calclist:
        Xgb = calc[id].Xgb
        X = calc[id].X
        dX.append(Xgb/1 - X)
        e_seg.append(calc[id].e_seg)
        # print calc[id].e_seg
        # print calc[id].X
        #print dX
    coeffs1 = np.polyfit(dX, e_seg, 1)        
    
    fit_func1 = np.poly1d(coeffs1)
    print "list of seg energies: ", e_seg
    print "list of dX          : ", dX

    print "Fitting using linear function:"
    print fit_func1
    print "E_seg0 = {0:0.0f} meV, standart enthalpy of segregation".format(fit_func1[0])
    print "alpha  = {0:0.0f} meV, interaction coefficient".format(-fit_func1[1]/2)

    return

def calculate_voronoi(self, state = 'end'):
    # By default two quantities per atom are calculated by this compute. 
    # The first is the volume of the Voronoi cell around each atom. 
    # Any point in an atom's Voronoi cell is closer to that atom than any other. 
    # The second is the number of faces of the Voronoi cell, which 
    # is also the number of nearest neighbors of the atom in the middle of the cell. 
    # state - init or end; if init then saved in self.init.vorovol; if end than saved in self.vorovol

    write_lammps(self, state, filepath = 'voronoi_analysis/structure.lammps') #write structure for lammps
    runBash("rm voronoi_analysis/dump.voro; /home/aksenov/installed/lammps-1Feb14/src/lmp_serial < voronoi_analysis/voronoi.in > voronoi_analysis/log")

    if state == 'end':
        self.vorovol = []
        self.vorofaces = []
        vorovol = self.vorovol
        vorofaces = self.vorofaces
    elif state == 'init':
        self.init.vorovol = []
        self.init.vorofaces = []
        vorovol = self.init.vorovol
        vorofaces = self.init.vorofaces        

    vsum=0
    wlist = []    
    with open('voronoi_analysis/dump.voro','r') as volfile:  #analyze dump.voro
        for line in volfile:
            if 'ITEM: ATOMS ' in line:
                break
        for line in volfile:
            ll = line.split()
            if int(ll[1]) > 1:
                wlist.append( [ll[0], ll[5], ll[6], ll[2]] )
            # print 'Volume of atom ',ll[0],'is', ll[5]
            vsum= vsum+float(ll[5])
        print 'Check total volume ', vsum, self.end.vol

        wlist.sort(key = itemgetter(0)) #sort according to the position of atoms
        print "atom #, voronoi vol, voronoi faces, x coordinate: ", 
        print wlist
        for w in wlist:
            vorovol.append(float(w[1]))
            vorofaces.append(int(w[2]))
        # print 'Voro vol  ',self.end.vorovol
        # print 'Voro faces',self.end.vorofaces
        # print len(wlist)
    if hasattr(self, 'vorovol'): 
        voro = ''
        if len(vorovol) == 2: #C and O
            voro = " {0:5.2f} & {1:2d} & {2:5.2f} & {3:2d} ".format(vorovol[0], vorofaces[0], vorovol[1], vorofaces[1]  ).center(25)
        else: 
            voro = " {0:5.2f} & {1:2d} ".format(vorovol[0], vorofaces[0] ).center(25)
        voro+='&'
    else:
        voro = ""
    print "Voronoi volume", voro
    return voro

def log_history(hstring):
    try:
        if hstring != header.history[-1]: header.history.append( hstring  )
    except:
        header.history.append( hstring  )    
    return

