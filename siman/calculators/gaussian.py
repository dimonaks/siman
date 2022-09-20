# Copyright (c) Siman Development Team.
# Distributed under the terms of the GNU License.
import os, sys, glob
import numpy as np
from siman import header
from siman.header import printlog
from siman.small_functions import makedir
from siman.core.calculation import Calculation
from siman.core.molecule import Molecule

from siman.set_functions import InputSet, qe_keys
from siman.functions import push_to_server, check_output


from pymatgen.io.gaussian import GaussianOutput
from pymatgen.electronic_structure.core import Spin




class CalculationGaussian(Calculation):
    """object for Gaussian"""
    def __init__(self, inset = None, iid = None, output = None):
        super(CalculationGaussian, self).__init__(inset, iid, output)
        self.len_units = 'Angstrom'
        self.calculator = 'gaussian'

    def set_output_filenames(self, out_name, version):
        cl = self
        cl.path["output"] = cl.dir+cl.name+'.out'



    def actualize_set(self, curset = None, params = None):
        """
        """
        return
    def check_kpoints(self):

        return


    def write_structure(self, name_of_output_file, type_of_coordinates = 'dir', option = None, prevcalcver = None, path = None, state = 'init'):

        if path == None: 
            path = self.dir
        
        if state == 'init':
            st  = self.init
        elif state == 'end':
            st  = self.end
        else: 
            raise RuntimeError 
        
        version_str = name_of_output_file.split('.')[0]

        filename = os.path.join(path, version_str+'.xyz')

        makedir(filename)
        # st.write_xyz(filename = filename)
        st.to(filename=filename, fmt = 'xyz')

        # write_geometry_qe(st, filename, coord_type = type_of_coordinates, periodic = self.set.periodic)

        # sys.exit()
        return

    def add_potcar(self):
        'nothing is needed'

        self.path['potcar'] = None


        return

    def calculate_nbands(self, curset, path_to_potcar = None, params = None):
        ''
        return



    def make_incar(self):
        """
        Create input file for gaussian
        """
        from pymatgen.io.gaussian import GaussianInput
        d = self.dir
        
        incar = d+'input.gau' # 
        # with open(incar, 'r') as f:
        sp = self.set.params.copy()

        basis_set  = sp.get('basis_set')
        functional = sp.get('functional')
        job_type = sp.get('job_type')
        spin_multiplicity = sp.get('multiplicity')
        charge = sp.get('charge')
        SCRF = sp.get('SCRF')
        # route_parameters = None#{"SCF":"Tight"}
        route_parameters = {job_type:''}
        if SCRF:
            route_parameters[SCRF] = ''
        mem = self.cluster.get('memory')
        if not mem:
            printlog('Warning! no memory limit in cluster description, I set default of 24GB')
            mem = 24

        link0_parameters = {'%NProcShared':self.cluster['corenum'], '%mem':str(mem)+'GB'}

        inp = GaussianInput(self.init, basis_set = basis_set, charge = charge, spin_multiplicity = spin_multiplicity,
            functional = functional, route_parameters = route_parameters, 
            input_parameters = sp.get('optional'), link0_parameters = link0_parameters )

        inp.write_file(incar, cart_coords=True)
        # sys.exit()
        
        return [incar]

    def make_kpoints_file(self):
        'to be acomplished'
        return ['']

    def copy_to_cluster(self, list_to_copy, update):
        d = self.dir
        # list_to_copy.extend( glob.glob(   os.path.join(d, '*geometry*')  ) ) # 
        
        if "up" in update: #Copy to server
            printlog('Files to copy:', list_to_copy)

            push_to_server(list_to_copy,  self.project_path_cluster +'/'+ self.dir, self.cluster_address)

    def download(self, load):

        path_to_outcar  = self.path["output"]

        # print(path_to_outcar)
        # sys.exit()

        self.get_file(os.path.basename(path_to_outcar), up = load)

        return path_to_outcar

    def read_output_file(self):
        filename = self.path["output"]
        with open(filename, 'r') as f:
            outputlines = f.readlines()

        self.force_prefix = ''
        self.maxforce_list = []
        self.rmsforce_list = []
        self.low_frequencies = [] # first nine freqs
        f = float
        for i, line in enumerate(outputlines):

            if 'Dipole moment (field-independent basis, Debye)' in line:
                # print(outputlines[i+1])
                v = outputlines[i+1].split() # values
                # print(v)
                self.dipole = f(v[7]) # magnitude
                self.dipole_xyz = (f(v[1]),f(v[3]),f(v[5]))
                # print(self.dipole)

            if 'Cartesian Forces:  Max' in line: #Forces (Hartrees/Bohr)
                # print(line)
                maxforce = float(line.split()[3])*header.Ha_Bohr_to_eV_A*1000
                rmsforce = float(line.split()[5])*header.Ha_Bohr_to_eV_A*1000
                self.maxforce_list.append(maxforce) # meV/A
                self.rmsforce_list.append(rmsforce) # meV/A
            if 'Maximum Force' in line: # which one to use?
                ''
                # print(line)

            if 'Low frequencies' in line:
                # print(line)
                freqs = [f(a) for a in line.split()[3:]]
                self.low_frequencies.extend(freqs)


        return

    def read_results(self, load = '', out_type = '', voronoi = None, show = '', choose_outcar = None, alkali_ion_number = None):

        """
        Gaussian

        choose_outcar - for now is dummy
        alkali_ion_number - for now is dummy
        voronoi - dummy
        """
        to_eV = header.to_eV
        cl = self
        filename = cl.download(load) # wrapper for downloading output files


        # cl.state = check_output(filename, 'Normal termination of Gaussian', load)
        cl.state = check_output(filename, 'Elapsed time', load)
        cl.homo = 0
        cl.lumo = 0
        cl.dipole = 0

        if "4" in cl.state:

            cl.read_output_file()
            go = GaussianOutput(filename)
            
            # print(go.energies)
            # cl.e0 = go.energies[-1]*to_eV #the same as final
            cl.e0 = go.final_energy*to_eV

            nelec = go.electrons[0]

            # print(go.eigenvalues)
            if go.eigenvalues:
                eigen = go.eigenvalues[Spin.up]

            # print(eigen)

                cl.homo = eigen[nelec-1]*to_eV
                cl.lumo = eigen[nelec]*to_eV
                # print(homo, lumo, homo  - lumo)
            cl.gap = cl.lumo - cl.homo 
                # print(go.MO_coefficients)
                # printlog(outstr)

            head   = 'Energy | Dipole, D | HOMO, eV | LUMO, eV | Gap, eV | basis'
            basis = cl.set.params['basis_set']
            outstr = f'{cl.e0:7.5f} eV | {cl.dipole:5.2f} D | {cl.homo:5.2f} eV | {cl.lumo:5.2f} eV | {cl.gap:5.2f} eV | {basis:s}'
            
            cl.go = go
            cl.end =  Molecule.cast(go.final_structure)
            if hasattr(cl.init, 'name'):
                cl.end.name = cl.init.name + '_end'
            # print(go.cart_forces)
            # print(cl.end.name)

            # if hasattr()
            'Show frequencies'
            np.set_printoptions(precision=1, suppress=True)

            if go.frequencies:
                f = go.frequencies[-1] # can be more than one freq calculation
                # print( len(f), type(f) )
                # for a in f[0:1]:
                    # print(a)
                freq_list = np.asarray([a['frequency']*header.cm_inv2eV*1000 for a in f]) # meV
                low_freq_list = np.asarray([a*header.cm_inv2eV*1000  for a in cl.low_frequencies])

                print('Frequencies for rotation and translation (meV):',low_freq_list[0:6] )
                print('Frequencies for normal modes (meV):',freq_list )
                if sum(a < 0 for a in freq_list):
                    printlog('Warning! Imaginiry frequencies are detected')

                dif = low_freq_list[6:9]-freq_list[0:3]
                rms = np.sqrt(np.mean(np.square(dif)))*header.cm_inv2eV*1000 #meV
                # print( (rms) )
                if rms > 1e-3: # 1 microeV
                    printlog('Warning! The normal modes are contaminated by the rotational and translational modes, rms=', rms, 'meV; see https://gaussian.com/vib/')



        else:
            
            printlog('Status of calculation is', cl.state, 'continiue', imp = 'y')
            outstr = cl.state
        
        if 'help' in show:
            print(head)



        return outstr
