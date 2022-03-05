# Copyright (c) Siman Development Team.
# Distributed under the terms of the GNU License.
import os, sys, glob
from siman.header import printlog
from siman.small_functions import makedir
from siman.core.calculation import Calculation

from siman.set_functions import InputSet, qe_keys
from siman.functions import push_to_server

def write_geometry_qe(st, filename, coord_type, periodic):
    """

    """

    



def read_qe_out(cl, out_type, show):
    """

    """



class CalculationGaussian(Calculation):
    """object for Gaussian"""
    def __init__(self, inset = None, iid = None, output = None):
        super(CalculationGaussian, self).__init__(inset, iid, output)
        self.len_units = 'Angstrom'
        self.calculator = 'gaussian'

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
        # route_parameters = None#{"SCF":"Tight"}
        route_parameters = {job_type:''}

        inp = GaussianInput(self.init, basis_set = basis_set, 
            functional = functional, route_parameters = route_parameters, 
            input_parameters = sp.get('optional') )

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

    def read_results(self, load = '', out_type = '', voronoi = None, show = '', choose_outcar = None, alkali_ion_number = None):

        """
        Gaussian

        choose_outcar - for now is dummy
        alkali_ion_number - for now is dummy
        voronoi - dummy
        """
        cl = self
        filename = cl.download(load) # wrapper for downloading output files


        cl.state = check_output(filename, 'string with succesfull ', load)
        
        if "4" in cl.state:

            outstr = read_qe_out(cl, out_type, show)
            
            printlog(outstr)

        else:
            
            printlog('Status of calculation is', cl.state, 'continiue', imp = 'y')
            outstr = cl.state
        

        return outstr
