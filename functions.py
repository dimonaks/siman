
from __future__ import division, unicode_literals, absolute_import 
import os, tempfile, copy, math, itertools, sys
import numpy as np
from operator import itemgetter


import header
from header import print_and_log, printlog, runBash
from small_functions import is_list_like, is_string_like, gunzip_file, makedir




def smoother(x, n, mul = 1):
    #mul - additionally multiplies values
     x_smooth = []
     L = len(x)
     store = np.zeros((n,1),float)
     for u in range(L-n):
          for v in range(n):
               store[v] = x[u+v]
          av = float(sum(store)) / n
          x_smooth.append(av*mul)
     
     for u in range(L-n,L):
          for v in range(L-u-1):
               store[v] = x[u+v]
          av = float(sum(store)) / n
          x_smooth.append(av*mul)
     return x_smooth




def list2string(ilist):
    #' '.join(['{:}']*len(lis)).format(*lis)
    return ' '.join(np.array(ilist).astype(str))


def run_on_server(command, addr):
    command = command.replace('\\', '/') # make sure is POSIX
    if header.ssh_object:
        out = header.ssh_object.run(command)
    else:
        out = runBash('ssh '+addr+' "'+command+'"')    
    return out



def push_to_server(files = None, to = None,  addr = None):
    """
    if header.ssh_object then use paramiko
    to (str)     - path to remote folder ! 
    """
    if not is_list_like(files):
        files = [files]    
    
    to = to.replace('\\', '/') # make sure is POSIX

    files_str = ' '.join(np.array(files ))

    printlog('push_to_server(): uploading files ', files, 'to', addr, to)
    
    command = ' mkdir -p {:}'.format( to )

    run_on_server(command, addr)


    if header.ssh_object:
        for file in files:
            # print(file, to)
            header.ssh_object.put(file,  to+'/'+os.path.basename(file) )
        out = ''
    else:
        out = runBash('rsync -uaz  '+files_str+ ' '+addr+':'+to)
    
    printlog(out)



    return out 


def file_exists_on_server(file, addr):

    file = file.replace('\\', '/') # make sure is POSIX

    printlog('Checking existence of file', file, 'on server', addr )
    if header.ssh_object:
        exist = header.ssh_object.fexists(file)
    else:
        exist = runBash('ssh '+addr+' ls '+file)
    if exist:
        res = True
    else:
        res = False

    printlog('File exist? ', res)

    return res



def get_from_server(files = None, to = None, to_file = None,  addr = None, trygz = True):
    """
    Download files using either  paramiko (higher priority) or rcync; 
    For paramiko header.ssh_object should be defined

    files (list of str)  - files on cluster to download 
    to (str)      - path to local folder ! 
    to_file (str) - path to local file (if name should be changed); in this case len(files) should be 1 

    The gz file is also checked

    RETURN
        result of download

    TODO:
    now for each file new connection is opened, 
    copy them in one connection
 


    """
    # print(addr)


    def download(file, to_file):

        if header.ssh_object:

            exist = file_exists_on_server(file, addr)
            # try:
            if exist:
                printlog('Using paramiko: ssh_object.get(): from  to ', file, to_file)
                header.ssh_object.get(file,  to_file  )
                out = ''
            # except FileNotFoundError:
            else:
                out = 'file not found'

        else:
            # print(addr,file,to_file)
            out = runBash('rsync -uaz  '+addr+':'+file+ ' '+to_file)

        if out:
            res = out
        else:
            res = 'OK'

        printlog('Download result is ', res)

        return out



    if not is_list_like(files):
        files = [files]
    
    files = [file.replace('\\', '/') for file in files] #make sure the path is POSIX




    files_str = ', '.join(np.array(files ))
    printlog('Trying to download', files_str, 'from server', imp = 'n')
    


    for file in files:

        if not to and not to_file: #use temporary file
            with tempfile.NamedTemporaryFile() as f:
                to_file_l = f.name #system independent filename

        elif not to_file: #obtain filename
            to_file_l = os.path.join(to, os.path.basename(file) )
        
        else:
            to_file_l = to_file

        makedir(to_file_l)

        out = download(file, to_file_l)

        if out and trygz:

            printlog('File', file, 'does not exist, trying gz', imp = 'n')
            file+='.gz'
            to_file_l+='.gz'
            out = download(file, to_file_l)

            if out:
                printlog('    No gz either!', imp = 'n')
            else:
                gunzip_file(to_file_l)


    return out








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
        print( init_salary, i+2000)

    salary2014 = 30000
    increase = salary2014/init_salary
    print( increase)

# salary_inflation()

def element_name_inv(el):
    el_dict = {'octa':200, 'n':0, 'H':1, 'He':2, 'Li':3, 'Be':4, 'B':5, 'C':6, 'N':7, 'O':8, 'F':9, 'Ne':10, 'Na':11, 'Mg':12, 'Al':13, 'Si':14, 'P':15, 'S':16, 'Cl':17, 'Ar':18, 'K':19, 'Ca':20, 'Sc':21, 'Ti':22, 'V':23, 'Cr':24, 'Mn':25, 'Fe':26, 'Co':27, 'Ni':28, 'Cu':29, 'Zn':30, 'Ga':31, 'Ge':32, 'As':33, 'Se':34, 'Br':35, 'Kr':36, 'Rb':37, 'Sr':38, 'Y':39, 'Zr':40, 'Nb':41, 'Mo':42, 'Tc':43, 'Ru':44, 'Rh':45, 'Pd':46, 'Ag':47, 'Cd':48, 'In':49, 'Sn':50, 'Sb':51, 'Te':52, 'I':53, 'Xe':54, 'Cs':55, 'Ba':56, 'La':57, 'Ce':58, 'Pr':59, 'Nd':60, 'Pm':61, 'Sm':62, 'Eu':63, 'Gd':64, 'Tb':65, 'Dy':66, 'Ho':67, 'Er':68, 'Tm':69, 'Yb':70, 'Lu':71, 'Hf':72, 'Ta':73, 'W':74, 'Re':75, 'Os':76, 'Ir':77, 'Pt':78, 'Au':79, 'Hg':80, 'Tl':81, 'Pb':82, 'Bi':83, 'Po':84, 'At':85, 'Rn':86, 'Fr':87, 'Ra':88, 'Ac':89, 'Th':90, 'Pa':91, 'U':92, 'Np':93, 'Pu':94, 'Am':95, 'Cm':96, 'Bk':97, 'Cf':98, 'Es':99, 'Fm':100, 'Md':101, 'No':102, 'Lr':103, 'Rf':104, 'Db':105, 'Sg':106, 'Bh':107, 'Hs':108, 'Mt':109, 'Ds':110, 'Rg':111, 'Cn':112, 'Uuq':114, 'Uuh':116, }
    nu_dict = { 200:'octa', 0:'n', 1:'H', 2:'He', 3:'Li', 4:'Be', 5:'B', 6:'C', 7:'N', 8:'O', 9:'F', 10:'Ne', 11:'Na', 12:'Mg', 13:'Al', 14:'Si', 15:'P', 16:'S', 17:'Cl', 18:'Ar', 19:'K', 20:'Ca', 21:'Sc', 22:'Ti', 23:'V', 24:'Cr', 25:'Mn', 26:'Fe', 27:'Co', 28:'Ni', 29:'Cu', 30:'Zn', 31:'Ga', 32:'Ge', 33:'As', 34:'Se', 35:'Br', 36:'Kr', 37:'Rb', 38:'Sr', 39:'Y', 40:'Zr', 41:'Nb', 42:'Mo', 43:'Tc', 44:'Ru', 45:'Rh', 46:'Pd', 47:'Ag', 48:'Cd', 49:'In', 50:'Sn', 51:'Sb', 52:'Te', 53:'I', 54:'Xe', 55:'Cs', 56:'Ba', 57:'La', 58:'Ce', 59:'Pr', 60:'Nd', 61:'Pm', 62:'Sm', 63:'Eu', 64:'Gd', 65:'Tb', 66:'Dy', 67:'Ho', 68:'Er', 69:'Tm', 70:'Yb', 71:'Lu', 72:'Hf', 73:'Ta', 74:'W', 75:'Re', 76:'Os', 77:'Ir', 78:'Pt', 79:'Au', 80:'Hg', 81:'Tl', 82:'Pb', 83:'Bi', 84:'Po', 85:'At', 86:'Rn', 87:'Fr', 88:'Ra', 89:'Ac', 90:'Th', 91:'Pa', 92:'U', 93:'Np', 94:'Pu', 95:'Am', 96:'Cm', 97:'Bk', 98:'Cf', 99:'Es', 100:'Fm', 101:'Md', 102:'No', 103:'Lr', 104:'Rf', 105:'Db', 106:'Sg', 107:'Bh', 108:'Hs', 109:'Mt', 110:'Ds', 111:'Rg', 112:'Cn', 114:'Uuq', 116:'Uuh', }
  
    # print type(el), el, type(str('sdf') )
    if is_string_like(el):
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




def return_atoms_to_cell(st):

    st = st.return_atoms_to_cell()
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
    # print "A,B=",A,B
    a = (A/B)**(1./3)
    c = a * B
    a = round(a,5)
    c = round(c,5)
    print_and_log( "a, c, c/a for cell with pure    hcp ", a_b, c_b, round(c_b/a_b,4), imp ='y' )
    print_and_log( "a, c, c/a for cell with first  atom ", a1, c1, round(c1/a1,4), imp ='y' )
    print_and_log( "a, c, c/a for cell with second atom ", a2, c2, round(c2/a2,4), imp ='y' )

    #for double cell
    a3 = (C/B)**(1./3)
    c3 = a3 * B
    a3 = round(a3,5)
    c3 = round(c3,5)    

    if type == "two_atoms":
        print_and_log( "a, c, c/a for cell with two   atoms ", a,  c, round(c/a,4), "# the same cell but with two atoms\n", imp ='y')
    elif type == "double_cell":
        print_and_log( "a, c, c/a for new cell              ", a3,  c3, round(c3/a3,4), "# for cell with V = V(first_cell) + V(second cell), but only for the case if V(second cell) == V(first_cell)", imp ='y')

    return a, c















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
    # print X
    # print Y
    # print Z

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
    print( "list of seg energies: ", e_seg  )
    print( "list of dX          : ", dX  )

    print( "Fitting using linear function:"  )
    print( fit_func1  )
    print( "E_seg0 = {0:0.0f} meV, standart enthalpy of segregation".format(fit_func1[0])  )
    print( "alpha  = {0:0.0f} meV, interaction coefficient".format(-fit_func1[1]/2)  )

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
        print_and_log( 'Check total volume ', vsum, self.end.vol)

        wlist.sort(key = itemgetter(0)) #sort according to the position of atoms
        print_and_log( "atom #, voronoi vol, voronoi faces, x coordinate: ", )
        print_and_log( wlist)
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
    print_and_log( "Voronoi volume = ", voro, imp = 'y')
    return voro

def log_history(hstring):
    try:
        if hstring != header.history[-1]: header.history.append( hstring  )
    except:
        header.history.append( hstring  )    
    return





def latex_table(table, caption, label, header = None, fullpage = '', filename = None, writetype = 'w', header0 = None, size = None,
    replace = None ):
    """
    If header is not provided, table[0] is used as a header

    header0 - additional header0 befor main header for complex tables
    
    path_to_paper should be provided

    replace - list of tuples for replacements

    """
    def myprint(string):
        if filename:
            f.write(string+"\n")
            print( string)
        else:
            print( string)


    if filename:
        # path = path_to_paper+'/tab/'
        path = ''
        f = open(path+filename, writetype)
        print_and_log("Saving table to "+path+filename+'\n')

    for i in range(len(table)):
        if is_list_like(table[i]):
            tab = ' & '.join([str(l) for l in table[i]])
            table[i] = tab


    n = len(table[0].split('&'))-2
    print( 'Number of columns = ', n + 2)
    
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
        myprint(table[0]+' \\\\')
        tabbeg = 1



    myprint('\\hline')
    for r in table[tabbeg:] :
        if '&-' in r:
            r = r.replace('-','--')
        else:
            r = r.replace(' -','--') #to save beautiful columns 
        r+=' '
        if '-- ' in r:
            r = r.replace('-- ',' - ')
        
        for rep in replace:
            # if rep[0] in r:

            r = r.replace(*rep)



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
        def convert(a): 
            return int(a)
    
    elif ttype == float: 
        def convert(a): 
            # print a
            return float(a)
    
    elif ttype == str  : 
        def convert(a): 
            return str(a)
    
    #print list_of_words[index], type(list_of_words[index])
    if list_of_words[index] == "None"  : 
        def convert(a): 
            return [None]
    
    #Make convertion
    for i in range(number_of_elements):
        
        if 'None' in list_of_words[index]:
            list_of_elements.append(None)
        else:
            list_of_elements.append(    convert(  list_of_words[index]   )     )
        

        index+=1


    return list_of_elements


def words(fileobj):
    """Generator of words. However does not allow to use methods of list for returned"""
    for line in fileobj:
        for word in line.split():
            yield word

