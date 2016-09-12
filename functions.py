# from classes import res_loop
import header
from header import *


def list2string(ilist):
    #' '.join(['{:}']*len(lis)).format(*lis)
    return ' '.join(np.array(ilist).astype(str))


def push_to_server(files = None, to = None,  addr = None):
    
    if not hasattr(files, '__iter__'):
        files = [files]    
    
    files_str = ' '.join(np.array(files ))
    

    return runBash('rsync -uaz  '+files_str+ ' '+addr+':'+to)

def get_from_server(files = None, to = None,  addr = None):
    """
    The zip file is checked only for the first file in list *files*
    """
    if not hasattr(files, '__iter__'):
        files = [files]
    
    files_str = ' :'.join(np.array(files ))
    print_and_log('Trying to download', files_str, 'from server')

    out = runBash('rsync -uaz  '+addr+':'+files_str+ ' '+to)
    # print 'out === ',out
    to_new = to+'/'+os.path.basename(files[0])
    if out:# and not os.path.exists(to_new):
        print_and_log('File', files[0], 'does not exist, trying gz', imp = 'n')
        files[0]+='.gz'
        # print files[0]
        out = runBash('rsync -uaz  '+addr+':'+files[0]+ ' '+to+';gunzip '+to_new)



    return out



def makedir(path):
    """
    Make directory at path if it not exist
    """
    if not os.path.exists(os.path.dirname(path)):
        os.makedirs(os.path.dirname(path))
        print_and_log("Directory "+os.path.dirname(path)+" was created\n")
    return



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

    if hasattr(st, 'init_numbers') and st.init_numbers:
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
        # print st.xred
        # st.xcart = xred2xcart(st.xred, st.rprimd)
        st.xred = xcart2xred(st.xcart, st.rprimd)


    st.natom = len(st.xcart)
    # print 'Structure is replicated; now', st.natom,'atoms' 
    return st



def local_surrounding(x_central, st, n_neighbours, control = 'sum', periodic = False, only_elements = None):
    """
    Return list of distances to n closest atoms around central atom. (By defauld sum of distances)
    
    Input:
    - x_central - cartesian coordinates of central atom; vector
    - st - structure with xcart list of coordinates of all atoms in system
    - n_neighbours - number of needed closest neighbours

    - control - type of output; sum - sum of distances, av - av distance, list - list of distances; 
              av_dev - average deviation, maximum deviation from average distance in mA.
              atoms  - coordinates of neighbours

    - periodic - if True, then cell is additionaly replicated; needed for small cells
    - *only_elements* - list of z of elements to which only the distances are needed
    """
    st_original = copy.deepcopy(st)
    st.init_numbers = None
    if periodic:
        st = replic(st, mul = (2,2,2), inv = 1 ) # to be sure that impurity is surrounded by atoms
        st = replic(st, mul = (2,2,2), inv = -1 )

    xcart = st.xcart
    typat = st.typat
    natom = st.natom
    # print x_central

    #print len(xcart)
    zlist = [int(st.znucl[t-1]) for t in st.typat]
    
    # if only_elements:
    #     dlist = [np.linalg.norm(x_central - x) for x, z in zip(xcart, zlist) if z in only_elements]
    # else:
    dlist = [np.linalg.norm(x_central - x) for x in xcart ]# if all (x != x_central)] # list of all distances

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
        if hasattr(st, 'init_numbers') and st.init_numbers:
            numbers = st.init_numbers
        else:
            numbers = range(natom)
        temp = zip(dlist_unsort, xcart, typat, numbers, zlist)
        if only_elements:
            # print temp[4]
            temp = [t for t in temp if t[4] in only_elements]
        # print temp
        temp.sort(key = itemgetter(0))
        temp2 = zip(*temp)
        dlist       = temp2[0][:n_neighbours+1]
        xcart_local = temp2[1][:n_neighbours+1]
        typat_local = temp2[2][:n_neighbours+1]
        numbers     = temp2[3][:n_neighbours+1]
        # print temp2[0][:n_neighbours]
        # print xcart_local[:n_neighbours]
        






        #check if atoms in output are from neighboring cells
        if 0:
            xred_local = xcart2xred(xcart_local, st_original.rprimd)
            # print 'xred_local', xred_local
            for x_l in xred_local:
                for i, x in enumerate(x_l):
                    if x > 1: 
                        x_l[i]-=1
                        # print 'returning to prim cell', x,x_l[i]
                    if x < 0: 
                        x_l[i]+=1
                        # print 'returning to prim cell', x,x_l[i]
            xcart_local = xred2xcart(xred_local, st_original.rprimd)

        # print 'Warning! local_surrounding() can return several atoms in one position due to incomplete PBC implementation; Improve please\n'

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
    el_dict = {'octa':200, 'n':0, 'H':1, 'He':2, 'Li':3, 'Be':4, 'B':5, 'C':6, 'N':7, 'O':8, 'F':9, 'Ne':10, 'Na':11, 'Mg':12, 'Al':13, 'Si':14, 'P':15, 'S':16, 'Cl':17, 'Ar':18, 'K':19, 'Ca':20, 'Sc':21, 'Ti':22, 'V':23, 'Cr':24, 'Mn':25, 'Fe':26, 'Co':27, 'Ni':28, 'Cu':29, 'Zn':30, 'Ga':31, 'Ge':32, 'As':33, 'Se':34, 'Br':35, 'Kr':36, 'Rb':37, 'Sr':38, 'Y':39, 'Zr':40, 'Nb':41, 'Mo':42, 'Tc':43, 'Ru':44, 'Rh':45, 'Pd':46, 'Ag':47, 'Cd':48, 'In':49, 'Sn':50, 'Sb':51, 'Te':52, 'I':53, 'Xe':54, 'Cs':55, 'Ba':56, 'La':57, 'Ce':58, 'Pr':59, 'Nd':60, 'Pm':61, 'Sm':62, 'Eu':63, 'Gd':64, 'Tb':65, 'Dy':66, 'Ho':67, 'Er':68, 'Tm':69, 'Yb':70, 'Lu':71, 'Hf':72, 'Ta':73, 'W':74, 'Re':75, 'Os':76, 'Ir':77, 'Pt':78, 'Au':79, 'Hg':80, 'Tl':81, 'Pb':82, 'Bi':83, 'Po':84, 'At':85, 'Rn':86, 'Fr':87, 'Ra':88, 'Ac':89, 'Th':90, 'Pa':91, 'U':92, 'Np':93, 'Pu':94, 'Am':95, 'Cm':96, 'Bk':97, 'Cf':98, 'Es':99, 'Fm':100, 'Md':101, 'No':102, 'Lr':103, 'Rf':104, 'Db':105, 'Sg':106, 'Bh':107, 'Hs':108, 'Mt':109, 'Ds':110, 'Rg':111, 'Cn':112, 'Uuq':114, 'Uuh':116, }
    nu_dict = { 200:'octa', 0:'n', 1:'H', 2:'He', 3:'Li', 4:'Be', 5:'B', 6:'C', 7:'N', 8:'O', 9:'F', 10:'Ne', 11:'Na', 12:'Mg', 13:'Al', 14:'Si', 15:'P', 16:'S', 17:'Cl', 18:'Ar', 19:'K', 20:'Ca', 21:'Sc', 22:'Ti', 23:'V', 24:'Cr', 25:'Mn', 26:'Fe', 27:'Co', 28:'Ni', 29:'Cu', 30:'Zn', 31:'Ga', 32:'Ge', 33:'As', 34:'Se', 35:'Br', 36:'Kr', 37:'Rb', 38:'Sr', 39:'Y', 40:'Zr', 41:'Nb', 42:'Mo', 43:'Tc', 44:'Ru', 45:'Rh', 46:'Pd', 47:'Ag', 48:'Cd', 49:'In', 50:'Sn', 51:'Sb', 52:'Te', 53:'I', 54:'Xe', 55:'Cs', 56:'Ba', 57:'La', 58:'Ce', 59:'Pr', 60:'Nd', 61:'Pm', 62:'Sm', 63:'Eu', 64:'Gd', 65:'Tb', 66:'Dy', 67:'Ho', 68:'Er', 69:'Tm', 70:'Yb', 71:'Lu', 72:'Hf', 73:'Ta', 74:'W', 75:'Re', 76:'Os', 77:'Ir', 78:'Pt', 79:'Au', 80:'Hg', 81:'Tl', 82:'Pb', 83:'Bi', 84:'Po', 85:'At', 86:'Rn', 87:'Fr', 88:'Ra', 89:'Ac', 90:'Th', 91:'Pa', 92:'U', 93:'Np', 94:'Pu', 95:'Am', 96:'Cm', 97:'Bk', 98:'Cf', 99:'Es', 100:'Fm', 101:'Md', 102:'No', 103:'Lr', 104:'Rf', 105:'Db', 106:'Sg', 107:'Bh', 108:'Hs', 109:'Mt', 110:'Ds', 111:'Rg', 112:'Cn', 114:'Uuq', 116:'Uuh', }
  
    if type(el) == str:
        try: 
            elinv = el_dict[el]
        except:
            print_and_log("Error! Unknown element: " +str(el))
            raise RuntimeError
    else:
        el = int(el)
        try:
            elinv = nu_dict[el]
        except:
            print_and_log("Error! Unknown element: "+str(el))
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
            j = 1
            name_old = ''
            for i, el in enumerate(label):
                name  = el[0]+el[1]
                if name != name_old: j = 1
                label = str(j)+el[1]
                print "label",label
                f.write('select '+el[0]+str(i+1)+'\ncpk 200\nset labeloffset 0 0\nset labelfront\ncolor label black\nlabel '+label+'\n font label 24 bold \n')
                j+=1
                name_old = name


        if rotate:
            f.write('rotate z '+str(rotate)+'\n')

        if specialcommand:
            f.write(specialcommand)

        
        # f.write('write image 2800 2800 png "'+pngfile+'"')
        f.write('write image 1800 1800 png "'+pngfile+'"')
    print runBash(header.path_to_jmol+' -ions '+scriptfile)
    # print runBash('convert '+pngfile+' -shave 0x5% -trim '+pngfile) #cut by 5% from up and down (shave) and that trim left background
    print pngfile
    print runBash('convert '+pngfile+' -trim '+pngfile) # trim background
    return




def write_xyz(st, path = '', repeat = 1, shift = 1.0,  gbpos2 = None, gbwidth = 1 , 
    imp_positions = [], specialcommand = None, analysis = None, show_around = None, replications = None, nnumber = 6, topview = True,
    file_name = None, full_cell = False, orientation = None, boundbox = 2, withgb = False,
    include_boundary = 2, rotate = None, imp_sub_positions = None, jmol = None
    ):
    """Writes st structure in xyz format in the folder xyz/path

    if repeat == 2: produces jmol script
    shift - in rprimd[1][1] - shift of the second view
    gbpos2 - position of grain boundary in A
    gbwidth - atoms aroung gbpos2 will be colored differently

    imp_positions - type and xcart coordinates additionally to be added to structure; to visulaze all impurity positions: for jmol
    imp_sub_positions - list of atom numbers; the typat of these atoms is changed: not used now

    specialcommand - any command at the end of script

    analysis - additional processing, allows to show only specifice atoms, 'imp_surrounding' - shows Ti atoms only around impurity
    replications - list of replications, (2,2,2) 
    nnumber - number of neighbours to show
    show_around - choose atom number around which to show

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
    else:
        name = st.name



    if natom != len(xred) != len(xcart) != len(typat) or len(znucl) != max(typat): 
        print "Error! write_xyz: check your arrays.\n\n"    
    # print st.natom, len(st.xred), len(st.xcart), len(st.typat), len(st.znucl), max(st.typat)
    
    print_and_log("write_xyz(): Name is", name, important = 'n')
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
    print_and_log("Writing xyz "+xyzfile+" \n")

    #analyze imp_positions
    if imp_sub_positions == None:
        imp_sub_positions = []
    nsub = 0
    for pos in imp_positions:
        if 's' not in pos[4]: continue # skip interstitial positions
        xs = np.asarray([pos[0],pos[1],pos[2]])
        nsub+=1
        # print xs
        for i, x in enumerate(xcart):
            # print np.linalg.norm( x-xs)
            if np.linalg.norm( x-xs) < 1:
                imp_sub_positions.append(i)

    if imp_sub_positions : 
        print imp_sub_positions, ': numbers of found atoms to be changed '


    # for i in sorted(indices, reverse=True):
    #     del somelist[i]




    with open(xyzfile,'w') as f:
        for i in range(repeat):
            f.write(str(natom + len(imp_positions)-nsub)+"\n")
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




    # os._exit(1)



    if jmol:
        """
        script mode for jmol. Create script file as well for elobarate visualization
        """
        
        """Choose gb atoms to change their color"""
        print 'position of boundary 2', gbpos2
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
            print 'Atoms at GB:', gbatoms
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
        
        print 'imp_positions = ',imp_positions
        write_jmol(xyzfile, pngfile, scriptfile, atomselection, topview = topview, rprimd =rprimd, shift = shift, label = [(pos[3], pos[4]) for pos in imp_positions], 
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





def latex_table(table, caption, label, header = None, fullpage = '', filename = None, writetype = 'w', header0 = None, size = None ):
    """
    If header is not provided, table[0] is used as a header

    header0 - additional header0 befor main header for complex tables
    
    path_to_paper should be provided

    """
    def myprint(string):
        if filename:
            f.write(string+"\n")
            print string
        else:
            print string


    if filename:
        path = path_to_paper+'/tab/'
        f = open(path+filename, writetype)
        print_and_log("Saving table to "+path+filename+'\n')


    n = len(table[0].split('&'))-2
    print 'Number of columns = ', n + 2
    
    myprint('\\begin{table'+fullpage+'}')
    myprint('\\center')
    if size: myprint('\\'+size)
        


    myprint('\\caption{'+caption+'}')
    myprint('\\label{'+label+'}')

    myprint('\\begin{tabular}{l'+ n*'c'+'r}')
    myprint('\\hline')

    if header0:
        myprint(header0+'\\\\')
        myprint('\\hline')

    if header:
        myprint(header+'\\\\')
        tabbeg = 0
    else:
        myprint(table[0]+'\\\\')
        tabbeg = 1



    myprint('\\hline')
    for r in table[tabbeg:] :
        if '&-' in r:
            r = r.replace('-','--')
        else:
            r = r.replace(' -','--') #to save beautiful columns 

        if 'hline' in r: 
            myprint(r)
        else:
            myprint(r + '\\\\')

    myprint('\\hline')
    myprint('\\end{tabular}')
    myprint('\\end{table'+fullpage+'}')

    if filename:
        f.close()
    return








def gb_energy_volume(gb,bulk):
    if (gb.end.rprimd[1] != bulk.end.rprimd[1]).any() or (gb.end.rprimd[2] != bulk.end.rprimd[2]).any():
        print_and_log("Warning! You are trying to calculate gb_energy from cells with different lateral sizes:"+str(gb.end.rprimd)+" "+str(bulk.end.rprimd)+"\n")
    #print bulk.vol
    V_1at = bulk.vol / bulk.natom #* to_ang**3

    E_1at = bulk.energy_sigma0 / bulk.natom 
    A = np.linalg.norm( np.cross(gb.end.rprimd[1], gb.end.rprimd[2])  ) #surface area of gb
    #print A
    gb.v_gb =      ( gb.vol              - V_1at * gb.natom) / A / 2. * 1000
    gb.e_gb =      ( gb.energy_sigma0    - E_1at * gb.natom) / A / 2. * eV_A_to_J_m * 1000
    gb.e_gb_init = ( gb.list_e_sigma0[0] - E_1at * gb.natom) / A / 2. * eV_A_to_J_m * 1000
    gb.bulk_extpress = bulk.extpress     
    #print "Calc %s; e_gb_init = %.3f J/m^2; e_gb = %.3f J/m; v_gb = %.3f angstrom "%(gb.name, gb.e_gb_init, gb.e_gb, gb.v_gb )
    outst = "%15s&%7.0f&%7.0f"%(gb.name, gb.e_gb, gb.v_gb)
    return outst




def headers():
    j = (7,12,14,7,8,9,9,5,5,20,5,20,8,12,20,8,5,8,8)
    d="&"
    header_for_bands= "Set".ljust(j[0])+d+"Etot".center(j[1])+d+"a1,a2".center(j[2])+d+"c".center(j[3])\
                +d+"time, m".center(j[4])+d+"ittime, s".center(j[5])+d+"Nmd,Avr.".rjust(j[6])+d\
                +"Warn!"+d+"nband"+d+"Added, \%"+"\\\\"

    header_for_ecut= "Set".ljust(j[0])+d+"Etot".center(j[1])+d+"a1,a2".center(j[2])+d+"c".center(j[3])\
                +d+"time, m".center(j[4])+d+"ittime, s".center(j[5])+d+"Nmd,Avr.".rjust(j[6])+d\
                +"Warn!"+d+"Ecut,eV"+"\\\\"

    header_for_npar= "Set".ljust(j[0])+d+"Etot".center(j[1])+d+"a1,a2".center(j[2])+d+"c".center(j[3])\
                +d+"time, m".center(j[4])+d+"ittime, s".center(j[5])+d+"Nmd,Avr.".rjust(j[6])+d\
                +"Warn!"+d+"NPAR".center(j[16])+d+"LPLANE".center(j[17])+"\\\\"

    header_for_kpoints= "Set".ljust(j[0])+d+"Etot".center(j[1])+d+"a1,a2".center(j[2])+d+"c".center(j[3])\
                +d+"time, m".center(j[4])+d+"ittime, s".center(j[5])+d+"Nmd,Avr.".rjust(j[6])+d\
                +"Warn!"+d+"k-mesh".center(j[8])+d+"k-spacings".center(j[9])+d+"nkpt".center(j[10])+"\\\\"
    header_for_tsmear= "Set".ljust(j[0])+d+"Etot".center(j[1])+d+"a1,a2".center(j[2])+d+"c".center(j[3])\
                +d+"time, m".center(j[4])+d+"ittime, s".center(j[5])+d+"Nmd,Avr.".rjust(j[6])+d\
                +"Warn!"+d+"k-mesh".center(j[8])+d+"tsmear, meV".center(j[13])+d+"Smearing error, meV/atom".center(j[14])+"\\\\"

    header_for_stress= "Set".ljust(j[0])+d+"Etot".center(j[1])+d+"a1,a2".center(j[2])+d+"c".center(j[3])\
                +d+"time, m".center(j[4])+d+"ittime, s".center(j[5])+d+"Nmd,Avr.".rjust(j[6])+d\
                +"Warn!"+d+"Stress, intr u.*1000".center(j[11])+d+"Pressure, MPa".center(j[12])
    #print "\\hline"
    return header_for_kpoints





def read_vectors(token, number_of_vectors, list_of_words):
    """Returns the list of numpy vectors for the last match"""

    number_of_matches = list_of_words.count( token )
    if number_of_matches == 0: 
        #print_and_log("Warning token '"+token+"' was not found! return empty\n")
        return [None]

    if number_of_matches > 1:
        print_and_log("Warning token '"+token+"' was found more than one times\n")
        raise RuntimeError


    index = list_of_words.index(token, number_of_matches - 1 )     #Return the index of the last match
    #print list_of_words[index]
    list_of_vectors = []
    vector = np.zeros((3))
    for i in range(number_of_vectors):
        vector[0] = float(list_of_words[index + 1])
        vector[1] = float(list_of_words[index + 2])
        vector[2] = float(list_of_words[index + 3])
        index+=3
        list_of_vectors.append(vector.copy())
    return list_of_vectors


def read_list(token, number_of_elements, ttype, list_of_words):
    """Input is token to find, number of elements to read, type of elements and list of words, 
    where to search
    Returns the list of elements for the last match"""
    

    number_of_matches = list_of_words.count( token )


    
    #if number_of_elements == 0:        raise RuntimeError
    
    if number_of_matches > 1:
        print_and_log("Warning token '"+token+"' was found more than one times\n")
        raise RuntimeError

    if number_of_matches == 0 or number_of_elements == 0: 
        #print_and_log("Warning token '"+token+"' was not found or asked number of elements is zero! set to [None]\n")
        #if ttype == str:
        #    return ['']*number_of_elements
        #else:
        #    return [0]*number_of_elements
        return [None]

    try:
        index = list_of_words.index(token, number_of_matches - 1 )     #Return the index of the last match

    except ValueError: 
        print_and_log("Warning!, token "+token+" was not found. I return [None]!\n")
        return [None]
    
    index+=1 #the position of token value
    list_of_elements = []
    
    #define function dependig on type:

    if   ttype == int  : 
        def convert(a): return int(a)
    elif ttype == float: 
        def convert(a): return float(a)
    elif ttype == str  : 
        def convert(a): return str(a)
    
    #print list_of_words[index], type(list_of_words[index])
    if list_of_words[index] == "None"  : 
        def convert(a): return [None]
    
    #Make convertion
    for i in range(number_of_elements):
        
        list_of_elements.append(    convert(  list_of_words[index]   )     )
        index+=1


    return list_of_elements


def words(fileobj):
    """Generator of words. However does not allow to use methods of list for returned"""
    for line in fileobj:
        for word in line.split():
            yield word

