# Copyright (c) Siman Development Team.
# Distributed under the terms of the GNU License.

import os, subprocess, sys, shelve
try:
    import matplotlib as mpl
    """Global matplotlib control"""
    # size = 22 #for one coloumn figures
    size = 26 #for DOS
    # size = 16 #for two coloumn figures
    mpl.rc('font',family='Serif')
    # mpl.rc('xtick', labelsize= size) 
    # mpl.rc('ytick', labelsize= size) 
    # mpl.rc('axes', labelsize = size) 
    # mpl.rc('legend', fontsize= size) 
    # mpl.rc('Axes.annotate', fontsize= size) #does not work
    mpl.rcParams.update({'font.size': size})
    mpl.rcParams.update({'mathtext.fontset': "stix"})
    # plt.rcParams['mathtext.fontset'] = "stix"

    #paths to libraries needed by siman
    # sys.path.append('/home/dim/Simulation_wrapper/ase') #
    # plt = None
except:
    mpl = None
    plt = None
    print('Warning! matplotlib is not installed! Some functions will not work!')

history = []

#Defaults to some project_conf values
PBS_PROCS = False # if true than #PBS -l procs="+str(number_cores) is used
WALLTIME_LIMIT = False # now only for PBS if True 72 hours limit is used
CIF2CELL = False


# warnings = 'neyY'
warnings = 'yY' # level of warnings to show: n - all, e - normal, y - important, Y - very important

#Global constants, can be overriden in project_conf.py and simanrc.py
FROM_ONE = None # atom numbering convention; allows to override the default behaviour of functions where atomic numbers are provided as arguments
PATH2ARCHIVE_LOCAL  = None
SAVE_CONTCAR = 1 # adds 's' to savefile argument alowing to move CONTCAR to 1.CONTCAR

#Global variables
final_vasp_clean     = True 
clean_vasp_files     = []
clean_vasp_files_ignore = []
default_savefile = 'oc'
copy_to_cluster_flag = True
close_run = False # alows to control close run file automatically after each add_loop
first_run = True  # needed to write header of run script
ssh_object = None # paramiko ssh_object
sshpass = None # using sshpass wrapper for rsync; see functions.py/push_to_server()
show = None
corenum = 1
check_job = 1 # check job by additional ssh requests
reorganize = 0 # use this from time to time to optimize database file size
verbose_log = 0 # in addition to normal log write verbose log in any case by openning the log_verbose each time
cluster_address = ''
override_cluster_address = 0 # 1 or 0, override read calculations to header.CLUSTERS[cluster]['address'], usefull when switching between proxy and back of the same cluster
pymatgen_flag = None









#1. Read default global settings for siman package
from siman.default_project_conf import *

#2. Read user-related settings for siman
simanrc = os.path.expanduser("~/simanrc.py")
if os.path.exists(simanrc):
    # if 'n' in war
    print(simanrc, 'was read')
    sys.path.insert(0, os.path.dirname(simanrc))
    from simanrc import *

#3. Read project specific
if os.path.exists('./project_conf.py'):
    print('Reading ./project_conf.py')

    from project_conf import *
    siman_run = True
    log = open('log','a')
else:
    # print('Some module is used separately; default_project_conf.py is used')
    if mpl and not os.path.exists(simanrc):
        mpl.use('agg') #switch matplotlib on or off; for running script using ssh
    siman_run = False
    history.append('separate run')













try:
    import matplotlib.pyplot as plt
except:
    plt = None

calc_database = 'only_calc.gdbm3'

class CalcDict(dict):
    def __getitem__(self, key):
        # print(self)
        if type(key) == str:
            # print('String key detected', key)
            key_str = key
            l = key.split('.')
            if len(l) > 2 and l[-1].isdecimal():
                key = ('.'.join(l[0:-2]), l[-2], int(l[-1]))            
        # else:

        if dict.__contains__(self, key):
            # print('key', key, 'is  in self')
            val = dict.__getitem__(self, key)
        
        else:
            with shelve.open(calc_database, protocol = 3) as d:
                try:
                    val = d[str(key)]
                    dict.__setitem__(self, key, val)
                    # print('reading ',str(key), 'from db')
                except:
                    print('Problem with key',key,', return None')
                    val = None
                    val = d[str(key)]

                    # for k in d:
                    #     # print(k[0])
                    #     if str(key) in k:
                    #         print(k)
        

        return val

    def __contains__(self, key):

        if dict.__contains__(self, key):
            return True
        
        else:
            with shelve.open(calc_database, protocol = 3) as d:
                # print('checking if key',key, 'is in db:', str(key) in d)
                # print(self)
                return str(key) in d

    def items(self):

        keys = []

        with shelve.open(calc_database, protocol = 3) as d:
            for key in d:
                key = key.replace("'", "").replace('(', '').replace(')', '') #remove unnessary symbols from string
                l = key.split(',') # go back to tuple

                if len(l) == 3 and l[2].strip().isdigit():
                    keyt = (l[0].strip(), l[1].strip(), int(l[2]))
                    keys.append(keyt)

            return keys



db = CalcDict()
# global db
calc = db
conv = {};
varset = {};
struct_des = {};

sets = varset





#Constants
to_ang = 0.52917721092
to_eV = 27.21138386
Ha_Bohr_to_eV_A = 51.4220641868956
kB_to_GPa = 0.1
eV_A_to_J_m = 16.021765
eV_nm_to_J_m = 1.6021765
THz2eV = 0.00413566553853599
kJ_mol2eV = 1.0364e-2
J_mol_T2eV_T = 1.0364e-5
J2eV = 6.242e+18
cm_inv2eV = 1.23981e-4 
h = 4.135668e-15 # eV/Hz
F = 96485.3329 # sA/mol

kB = 8.617e-5 # eV/K
kB_SI = 1.380649e-23 # J/K
R = 8.3145 # J/mol/K
N_A = 6.02214e23
e_SI = 1.602176634e-19 #C
me_SI = 9.1093837015e-31 #kg

ANION_ELEMENTS = [6,7,8,9,16,17,34,35,53]
TRANSITION_ELEMENTS = [22, 23, 24, 25, 26, 27, 28, 29, 30, 40, 41,42,43,44, 74]
ALKALI_ION_ELEMENTS = [3, 11, 19, 37]
MAGNETIC_ELEMENTS = [26, 27, 28]
TM_MAG = {'NM':0.6,
          'Co2+hs':3, 'Co2+ls':1, 'Co3+hs':4, 'Co3+ls':0, 'Co4+hs':5, 'Co4+ls':1,
          'Ni2+':2, 'Ni3+hs':3, 'Ni3+ls':1, 'Ni4+hs':4, 'Ni4+ls':0,
          'Mn3+ls':0, 'Mn3+hs':4, 'Mn2+hs':5,'Mn2+ls':1, 'Mn4+':4,
          'Fe3+hs': 5, 'Fe3+ls': 1, 'Fe2+ls': 0,'Fe2+hs': 4, 'Fe4+': 4,
          'Cr3+':3, 'Cr4+':2, 'Cr2+hs':4, 'Cr2+ls':2,
          'V4+': 1, 'V2+': 3, 'V3+': 2,
          'Ti3+':1, 'Ti4+':0, 'Ti2+':2 ,
          'Al3+':0, 
          'W4+':0, 
          'Mg2+':0, 
          }
           # canonical magnetic moments of TM depending on oxidation state and spin state
dic_sol = {'zero': 0.0, 'elec':0.2, 'hole':-0.2, 'Co_Li':-0.1, 'Ni_Li': -0.05, 'Li_Co':0.15, 'Ni_Na': -0.18, 'Na_Ni': 0.1,
            'W_Co':0.1, 'Fels_Co':0.05, 'Fehs_Co':0.1, 'Ti_Co':0.15, 'Ti4+_Co':0.06, 'Cr3+_Co':0.7, 'Cr4+_Co':0.05,
            'Cr2+hs_Co':0.26, 'Cr2+ls_Co':0.19, 'Ni3+hs_Co':0.05,'Ni3+ls_Co':0, 'Ni4+hs_Co':-0.07, 'Ni2+_Co':0.14,
            'Mn3+hs_Co':0.1,'Mn3+ls_Co':0.03, 'Mn2+hs_Co':0.15,'Mg2+_Co':0.18, 'Mn4+_Co':-0.15, 'V4+_Co': 0.03,
            'Co3+hs_Ni': 0.04, 'Mg2+_Ni':0.15, 'V4+_Ni':0, 'Ti4+_Ni': 0.04, 'Mg2+_Ni4':-0.2,
            'Al_Co':0,}

ION_RADII_coord_6 = { 
                    'Co2+hs':0.745, 'Co2+ls':0.65, 'Co3+hs':0.61, 'Co3+ls':0.545, 'Co4+hs':0.53, 'Co4+ls': 0.53,
                    'Ni2+':0.69, 'Ni3+hs':0.6, 'Ni3+ls':0.56, 'Ni4+ls':0.48, 'Ni4+hs':0.48,
                    'Mn2+ls':0.67, 'Mn2+hs':0.83, 'Mn3+ls':0.58, 'Mn3+hs':0.645, 'Mn4+':0.53,
                    'Fe2+hs':0.78, 'Fe2+ls':0.61, 'Fe3+ls':0.55, 'Fe3+hs':0.645, 'Fe4+':0.585,
                    'Cr2+ls':0.73, 'Cr2+hs':0.8,  'Cr3+':0.615,  'Cr4+':0.55,
                    'V2+':0.79, 'V3+':0.64, 'V4+':0.58,
                    'Ti2+':0.86, 'Ti3+':0.67, 'Ti4+':0.605, 
                    'Mg2+':0.72,
                    'Al3+':0.535,
                    'Cu3+':0.54, 'Cu2+':0.73, 'Cu1+':0.77,
                    }

el_dict = {'void':300, 'octa':200, 'n':0, 'H':1, 'He':2, 'Li':3, 'Be':4, 'B':5, 'C':6, 'N':7, 'O':8, 'F':9, 'Ne':10, 'Na':11, 'Mg':12, 'Al':13, 'Si':14, 'P':15, 'S':16, 'Cl':17, 'Ar':18, 'K':19, 'Ca':20, 'Sc':21, 'Ti':22, 'V':23, 'Cr':24, 'Mn':25, 'Fe':26, 'Co':27, 'Ni':28, 'Cu':29, 'Zn':30, 'Ga':31, 'Ge':32, 'As':33, 'Se':34, 'Br':35, 'Kr':36, 'Rb':37, 'Sr':38, 'Y':39, 'Zr':40, 'Nb':41, 'Mo':42, 'Tc':43, 'Ru':44, 'Rh':45, 'Pd':46, 'Ag':47, 'Cd':48, 'In':49, 'Sn':50, 'Sb':51, 'Te':52, 'I':53, 'Xe':54, 'Cs':55, 'Ba':56, 'La':57, 'Ce':58, 'Pr':59, 'Nd':60, 'Pm':61, 'Sm':62, 'Eu':63, 'Gd':64, 'Tb':65, 'Dy':66, 'Ho':67, 'Er':68, 'Tm':69, 'Yb':70, 'Lu':71, 'Hf':72, 'Ta':73, 'W':74, 'Re':75, 'Os':76, 'Ir':77, 'Pt':78, 'Au':79, 'Hg':80, 'Tl':81, 'Pb':82, 'Bi':83, 'Po':84, 'At':85, 'Rn':86, 'Fr':87, 'Ra':88, 'Ac':89, 'Th':90, 'Pa':91, 'U':92, 'Np':93, 'Pu':94, 'Am':95, 'Cm':96, 'Bk':97, 'Cf':98, 'Es':99, 'Fm':100, 'Md':101, 'No':102, 'Lr':103, 'Rf':104, 'Db':105, 'Sg':106, 'Bh':107, 'Hs':108, 'Mt':109, 'Ds':110, 'Rg':111, 'Cn':112, 'Uuq':114, 'Uuh':116, }
nu_dict = {300:'void', 200:'octa', 0:'n', 1:'H', 2:'He', 3:'Li', 4:'Be', 5:'B', 6:'C', 7:'N', 8:'O', 9:'F', 10:'Ne', 11:'Na', 12:'Mg', 13:'Al', 14:'Si', 15:'P', 16:'S', 17:'Cl', 18:'Ar', 19:'K', 20:'Ca', 21:'Sc', 22:'Ti', 23:'V', 24:'Cr', 25:'Mn', 26:'Fe', 27:'Co', 28:'Ni', 29:'Cu', 30:'Zn', 31:'Ga', 32:'Ge', 33:'As', 34:'Se', 35:'Br', 36:'Kr', 37:'Rb', 38:'Sr', 39:'Y', 40:'Zr', 41:'Nb', 42:'Mo', 43:'Tc', 44:'Ru', 45:'Rh', 46:'Pd', 47:'Ag', 48:'Cd', 49:'In', 50:'Sn', 51:'Sb', 52:'Te', 53:'I', 54:'Xe', 55:'Cs', 56:'Ba', 57:'La', 58:'Ce', 59:'Pr', 60:'Nd', 61:'Pm', 62:'Sm', 63:'Eu', 64:'Gd', 65:'Tb', 66:'Dy', 67:'Ho', 68:'Er', 69:'Tm', 70:'Yb', 71:'Lu', 72:'Hf', 73:'Ta', 74:'W', 75:'Re', 76:'Os', 77:'Ir', 78:'Pt', 79:'Au', 80:'Hg', 81:'Tl', 82:'Pb', 83:'Bi', 84:'Po', 85:'At', 86:'Rn', 87:'Fr', 88:'Ra', 89:'Ac', 90:'Th', 91:'Pa', 92:'U', 93:'Np', 94:'Pu', 95:'Am', 96:'Cm', 97:'Bk', 98:'Cf', 99:'Es', 100:'Fm', 101:'Md', 102:'No', 103:'Lr', 104:'Rf', 105:'Db', 106:'Sg', 107:'Bh', 108:'Hs', 109:'Mt', 110:'Ds', 111:'Rg', 112:'Cn', 114:'Uuq', 116:'Uuh', }
#void 300 is not written to poscar
nme_list = ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Si', 'P', 'S', 'Cl', 'Ar', 'As', 'Se', 'Br', 'Kr', 'Te', 'I', 'Xe', 'At', 'Rn', 'Ts', 'Og']

EXCLUDE_NODES = False
MEM_CPU = None


def printlog(*logstrings, **argdic):
    """
    '' - silent
    e - errors and warnings
    a - attentions
    m - minimalistic output of scientific procedures - only obligatory mess are shown
    M - maximalistic output of scientific procedures  
    debug_level importance:
        'n' - not important at all - for debugging
        ''  - almost all actions, no flag is needed
        'y' - important - major actions
        'Y' - super important, or output asked by user
    """
    end = '\n\n'# no argument for end, make one separate line
    
    debug_level  = 'e' #empty
    for key in argdic:
        if 'imp' in key:
            debug_level = argdic[key]
        
        if 'end' in key:
            end = argdic[key]


    mystring = ''
    for m in logstrings:
        mystring+=str(m)+' '


    if len(mystring.splitlines()) == 1:
        mystring = '-- '+mystring
    else:
        ''
        # mystring = '    '+mystring.replace('\n', '\n    ') 
    
    mystring+=end


    if 'Error' in mystring:# or 'Warning' in mystring:
        mystring+='\n\n\n'
    
    if 'Warning' in mystring or 'Attention' in mystring:
        debug_level = 'Y'

    for level in 'neyY':
        # print(level, debug_level)
        if (level in warnings and level in debug_level):
            print (mystring,  end = "")

    # if warnings:
    #     ''
    #     # print(debug_level)
    #     if 'n' in debug_level and 'n' not in warnings:
    #         pass
    #     else:
    #         print (mystring,  end = "")

    if siman_run:
        log.write(mystring)
    
    if verbose_log:
        with open('verbose_log','a') as f:
            f.write(mystring)


    if 'Error!' in mystring:
        print (mystring)
        print ('Error! keyword was detected in message; invoking RuntimeError ')
        sys.exit()
        # raise RuntimeError

    return

print_and_log = printlog
open_terminal = False


def runBash(cmd, env = None, detached = False, cwd = None):
    """Input - string; Executes Bash commands and returns stdout
Need: import subprocess
    """
    global open_terminal

    if detached:
        stdout = None
        stderr = None
    else:
        stdout = subprocess.PIPE
        stderr = subprocess.STDOUT

    printlog('running in BASH:', cmd, '\n')
    # if 'sshpass' in cmd:
    #     sys.exit()

    my_env = os.environ.copy()
    # my_env["PATH"] = "/opt/local/bin:/opt/local/sbin:" + my_env["PATH"]
    
    if open_terminal:
        cmd = cmd[:-1]+'; exec bash;'+cmd[-1:]
        # cmd = cmd.replace('"', '\\"')
        cmd = cmd.replace('ssh ','ssh -t ' )
        cmd2 = "gnome-terminal -- bash -c ' "+cmd+" ; exec bash' "
        # print(cmd2)
        os.system(cmd2)

        open_terminal = False
    else:

        p = subprocess.Popen(cmd, 
            # executable='/bin/bash', 
            shell=True, stdout=stdout, stderr = stderr, stdin = None, env = my_env, cwd = cwd)
        

    # print (cmd)
    # print 'Bash output is\n'+out
    # print ( str(out, 'utf-8') ) 
    out = ''
    try:
        out = p.stdout.read().strip()

        out = str(out, 'utf-8')
    except:
        pass

    return out  #This is the stdout from the shell command




def run_win_cmd(cmd):
    #example for windows
    result = []
    process = subprocess.Popen(cmd,
                               shell=True,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    for line in process.stdout:
        result.append(line)
    errcode = process.returncode
    for line in result:
        print(line)
    if errcode is not None:
        raise Exception('cmd %s failed, see above for details', cmd)



def pickle_module_migration_script():
    """
    This script allows to update modules in database after moving them to siman package
    """
    import sys
    import siman
    from siman import header
    # sys.modules['set_functions'] = siman.set_functions
    # sys.modules['classes'] = siman.classes
    # sys.modules['siman.core.structure'] = siman.core.molecule
    for key in header.db.items():
        print(header.db[key])
