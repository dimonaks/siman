
from __future__ import division, unicode_literals, absolute_import 
import os, tempfile, copy, math, itertools, sys
import numpy as np
from operator import itemgetter
from itertools import product

try:
    import scipy

except:
    print('functions.py: no scipy, smoother() will not work()')

from siman import header
from siman.header import print_and_log, printlog, runBash, eV_A_to_J_m
from siman.small_functions import is_list_like, is_string_like, gunzip_file, makedir, grep_file, setting_sshpass


def unique_elements(seq, idfun=None): 
   # return only unique_elements order preserving
   if idfun is None:
       def idfun(x): return x
   seen = {}
   result = []
   for item in seq:
       marker = idfun(item)
       # in old Python versions:
       # if seen.has_key(marker)
       # but in new ones:
       if marker in seen: continue
       seen[marker] = 1
       result.append(item)
   return result



def smoother(x, n, mul = 1, align = 1):
    """
    mul - additionally multiplies values
    #align - find first non-zero point and return it to zero
    #n - smooth value, 
        if algo = 'gaus' than it is sigma
        use something like 0.8 
        if algo = 'my'
            n of 10-15 is good
    """
    algo = 'gaus'
    # algo = 'my'

    if algo == 'my':
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
    

    elif algo == 'gaus':
        x_smooth =x
        # x_smooth = scipy.ndimage.filters.median_filter(x,size =4)
        # print('sigma is ', n)
        x_smooth = scipy.ndimage.filters.gaussian_filter1d(x_smooth, n, order =0)
        # x_smooth = scipy.ndimage.interpolation.spline_filter1d(x, 4)

    else:
        x_smooth = x

    if align:
        # print(x_smooth[0])
        x_smooth[0] = 0
        # sys.exit()



    return np.asarray(x_smooth)







def run_on_server(command, addr = None):
    printlog('Running', command, 'on server ...')
    command = command.replace('\\', '/') # make sure is POSIX
    # sys.exit()
    
    # print(header.sshpass)
    # sys.exit()

    if addr is None:
        addr = header.cluster_address

    if header.ssh_object:
        # printlog('Using paramiko ...', imp = 'y')
        # if 'ne' in header.warnings:
        # sys.exit()

        out = header.ssh_object.run(command, noerror = True, printout = 'ne' in header.warnings)

    elif header.sshpass and header.sshpass == 'proxy':
        com = 'ssh -tt sdv sshpass -f '+ header.path2pass +' ssh '+addr+' "'+command+'"'
        # print(com)
        # sys.exit()

        out = runBash(com) 
        # print(out)
        out = out.split('Connection to')[0] # remove last message Connection to ipaddress closed
        # sys.exit()

    elif header.sshpass:
        com = 'sshpass -f '+header.path2pass+' ssh '+addr+' "'+command+'"'
        # print(com)
        # sys.exit()
        
        out = runBash(com)    

        # sys.exit()

    else:
        bash_comm = 'ssh '+addr+' "'+command+'"'
        # print(bash_comm)
        # sys.exit()
        out = runBash(bash_comm)    
    
    out = out.split('#')[-1].strip()

    printlog(out)
    # print(out)
    # sys.exit()
    


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

    
    command = ' mkdir -p {:}'.format( to )
    # print('asfsadfdsf', to)
    printlog('push_to_server():', command, run_on_server(command, addr))
    # sys.exit()

    printlog('push_to_server(): uploading files ', files, 'to', addr, to)

    if header.ssh_object:
        for file in files:
            #print(file, to+'/'+os.path.basename(file))
            header.ssh_object.put(file,  to+'/'+os.path.basename(file) )
        out = ''

    elif header.sshpass and header.sshpass == 'proxy':
        com = 'tar cf - '+ files_str + ' | ssh sdv "sshpass -f ~/.ssh/p ssh '+addr+' \\"cd '+header.cluster_home+' && tar xvf -\\"" '
        # print(com)
        # sys.exit()
        out = runBash(com)
    
        # print(out)
        # sys.exit()
    elif header.sshpass:
        # if '@' not in addr:
        #     printlog('Error! Please provide address in the form user@address')
        # l = addr.split('@')
        # print(l)
        # user = l[0]
        # ad   = l[1]
        # com = 'rsync --rsh='+"'sshpass -f /home/aksenov/.ssh/p ssh' "  +' -uaz  '+files_str+ ' '+addr+':'+to
        com = 'rsync --rsh='+"'sshpass -f "+header.path2pass+" ssh' "  +' -uaz  '+files_str+ ' '+addr+':'+to

        # print(com)
        # sys.exit()
        out = runBash(com)



    else:
        out = runBash('rsync -uaz  '+files_str+ ' '+addr+':'+to)
    
    printlog(out)



    return out 


def file_exists_on_server(file, addr):

    file = file.replace('\\', '/') # make sure is POSIX

    printlog('Checking existence of file', file, 'on server', addr )
    
    exist = run_on_server(' ls '+file, addr)

    # if header.ssh_object:
    #     exist = header.ssh_object.fexists(file)
    # else:
    #     exist = runBash('ssh '+addr+' ls '+file)
     

    if 'No such file' in exist:
        exist = ''
    else:
        exist = 'file exists'


    if exist:
        res = True
    else:
        res = False

    printlog('File exist? ', res)

    return res



def get_from_server(files = None, to = None, to_file = None,  addr = None, trygz = True):
    """
    Download files using either  paramiko (higher priority) or rsync; 
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
    # sys.exit()


    def download(file, to_file):

        # print(header.sshpass)
        if header.ssh_object:

            exist = file_exists_on_server(file, addr)
            # try:
            if exist:
                printlog('Using paramiko: ssh_object.get(): from  to ', file, to_file)
                header.ssh_object.get(file,  to_file  )
                out = ''
            # except FileNotFoundError:
            else:
                out = 'error, file not found'

        elif header.sshpass and header.sshpass == 'proxy':
            # com = 'ssh sdv "sshpass -f ~/.ssh/p ssh ' + addr + ' \\"tar zcf - '+ file +'\\"" | tar zxf - '+to_file # does not work?
            com = 'ssh sdv "sshpass -f ~/.ssh/p ssh ' + addr + ' \\"tar cf - '+ file +'\\"" > '+to_file
            # print('sshpass',com)
            # sys.exit()
            out = runBash(com)

        elif header.sshpass:
            #com = 'rsync --rsh='+"'sshpass -f /home/aksenov/.ssh/p ssh' "  +' -uaz  '+addr+':'+file+ ' '+to_file
            com = 'rsync --rsh='+"'sshpass -f "+header.path2pass+" ssh' "  +' -uaz  '+addr+':'+file+ ' '+to_file

            out = runBash(com)
            # print(addr)
            # sys.exit()

        else:
            # print(addr,file,to_file)
            out = runBash('rsync -uaz  '+addr+':'+file+ ' '+to_file)





        if 'error' in out:
            res = out
        else:
            res = 'OK'
            out = ''

        printlog('Download result is ', res)

        return out





    if '*' in files:
        printlog('get_from_server(): get by template')
        files = run_on_server('ls '+files, addr).splitlines()
        # print(files)
        # sys.exit()
        printlog('get_from_server(): I download', files)



    elif not is_list_like(files):
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
            # run_on_server
            files = run_on_server(' ls '+file+'*', addr)
            file = files.split()[-1]
            # print(file)
            nz = file.count('gz')
            ext = '.gz'*nz

            # file+='.gz'
            to_file_l+=ext

            if file:
                out = download(file, to_file_l)
                printlog('    gz found with multiplicity', ext, imp = 'n')

                for i in range(nz):
                    printlog('unzipping', to_file_l)
                    gunzip_file(to_file_l)
                    to_file_l = to_file_l[:-3]
            else:
                printlog('    No gz either!', imp = 'n')

            # if '5247' in file:
            #     sys.exit()



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
    return

def element_name_inv(el):
    el_dict = header.el_dict
    nu_dict = header.nu_dict
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

invert = element_name_inv


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




def rotation_matrix_from_vectors(vec1, vec2):
    """ Find the rotation matrix that aligns vec1 to vec2
    :param vec1: A 3d "source" vector
    :param vec2: A 3d "destination" vector
    :return mat: A transform matrix (3x3) which when applied to vec1, aligns it with vec2.
    """
    a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)
    kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
    return rotation_matrix






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





def read_vectors(token, number_of_vectors, list_of_words, type_func = None, lists = False):
    """Returns the list of numpy vectors for the last match"""
    # lists - return list of lists instead list of vectors

    if type_func is None:
        type_func = lambda a : float(a)

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
    list_of_lists = []
    vector = np.zeros((3))
    for i in range(number_of_vectors):
        vector[0] = type_func(list_of_words[index + 1])
        vector[1] = type_func(list_of_words[index + 2])
        vector[2] = type_func(list_of_words[index + 3])
        list3 = []
        for j in 1,2,3:
            list3.append(type_func(list_of_words[index + j]) )

        index+=3
        list_of_vectors.append(vector.copy())
        list_of_lists.append(list3)
    
    if lists:
        out = list_of_lists
    else:
        out = list_of_vectors

    return out


def read_string(token, length, string):
    sh = len(token)+1
    i = string.find(token)+sh
    # print('length', i, i+length)
    # sys.exit()
    if i == -1:
        return ''
    else:
        return string[i:i+length]



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
    # print('string 839 functions.py Blind Guardian! token', token)
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




def server_cp(copy_file, to, gz = True, scratch = False, new_filename = None):
    
    if scratch:
        if not header.PATH2ARCHIVE:
            printlog('Warning! PATH2ARCHIVE is empty! Please put path archive in ~/simanrc.py or ./project_conf.py ')

        copy_file = header.PATH2ARCHIVE + '/' + copy_file
    else:
        copy_file = header.project_path_cluster + '/' + copy_file

    filename = os.path.basename(copy_file)

    if new_filename is None:
        new_filename = filename


    if gz:
        command = 'cp '+copy_file + ' ' + to +'/'+new_filename + '.gz ; gunzip -f '+ to+ '/'+new_filename+'.gz'
    else:
        command = 'cp '+copy_file + ' ' + to +'/'+new_filename 



    printlog('Running on server', command, imp = '')
    if file_exists_on_server(copy_file, header.cluster_address):
        out = run_on_server(command, addr = header.cluster_address)
        printlog('Output of run_on_server', out, imp = '')
    else:
        out = 'error, file does not exist on server: '+copy_file                
    return out



def wrapper_cp_on_server(file, to, new_filename = None):
    """
    tries iterativly scratch and gz
    """
    copy_to   = to

    copy_file = file

    filename = os.path.basename(file)
    if new_filename:
        app = 'with new name '+new_filename
    else:
        app = ''

    for s, gz in product([0,1], ['', '.gz']):

        printlog('scratch, gz:', s, gz)

        out = server_cp(copy_file+gz, to = to, gz = gz, scratch = s, new_filename = new_filename)

        if out == '':
            printlog('File', filename, 'was succesfully copied to',to, app, imp = 'y')
            break
        # else:
    else:
        printlog('Warning! File was not copied, probably it does not exist. Try using header.warnings = "neyY" for more details', imp = 'y')


    return






def check_output(filename, check_string, load):
    """
    Check if file exist and it is finished by search for check_string
    """

    if filename and os.path.exists(filename):

        out = grep_file(check_string, filename, reverse = True)

        printlog('The grep result of',filename, 'is:', out)
        # sys.exit()
        if check_string in out or 'un' in load:
            state = '4. Finished'
        else:
            state = '5. Broken output file'

    else:
        state = '5. no output file'

    return state
