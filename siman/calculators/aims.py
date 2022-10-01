
# Copyright (c) Siman Development Team.
# Distributed under the terms of the GNU License.
from siman.core.calculation import Calculation

class CalculationAims(Calculation):
    """object for Aims code """
    def __init__(self, inset = None, iid = None, output = None):
        super(CalculationAims, self).__init__(inset, iid, output)
        self.len_units = 'Angstrom'
        self.calculator = 'aims'
        self.init = Structure()
        self.end = Structure()

    def set_output_filenames(self, out_name, version):
        cl = self

        cl.path["output"] = cl.dir+cl.name+'.log'



    def write_structure(self, name_of_output_file, type_of_coordinates = 'dir', option = None, prevcalcver = None, path = None, state = 'init'):

        if path == None: 
            path = self.dir
        
        if state == 'init':
            st  = self.init
        elif state == 'end':
            st  = self.end
        else: 
            raise RuntimeError 
        
        filename = os.path.join(path, 'geometry.in')

        makedir(filename)

        write_geometry_aims(st, filename, coord_type = type_of_coordinates, periodic = self.set.periodic)


    def add_potcar(self):

        d = self.dir

        incar = d+'control.in'

        with open(self.set.path_to_potcar, 'r') as f:
            fil = f.read()

        with open(incar, 'w') as f:
            f.write(fil)

        self.path['potcar'] = self.set.path_to_potcar

    def make_incar(self):
        d = self.dir
        
        incar = d+'control.in'
        with open(incar, 'r') as f:
            fil = f.read()
        vp = self.set.params
        
        N = self.check_kpoints()
        # print(N)
        # self.exit()
        if N:
            vp['k_grid'] = list2string(N)

        with open(incar, 'w') as f:
            f.write(vp['universal'])
            f.write('\n')
            for key in vp:
                if key in aims_keys:
                    # print(key, self.set.params[key])
                    if vp[key] is not None:
                        f.write(key+' '+str(vp[key])+'\n')
            f.write(fil)
        
        return [incar]

    def make_kpoints_file(self):
        printlog( "Attention! ngkpt for kpoints file are created from kspacing\n")
        N = self.check_kpoints()
        self.set.ngkpt = N
        return ['']


    def copy_to_cluster(self, list_to_copy, update):
        d = self.dir
        list_to_copy.extend( glob.glob(   os.path.join(d, '*geometry*')  ) )
        
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
        Aims

        choose_outcar - for now is dummy
        alkali_ion_number - for now is dummy
        voronoi - dummy
        """
        cl = self
        filename = cl.download(load) # wrapper for downloading output files


        cl.state = check_output(filename, 'Have a nice day', load)
        
        if "4" in cl.state:

            outstr = read_aims_out(cl, out_type, show)
            
            printlog(outstr)

        else:
            
            printlog('Status of calculation is', cl.state, 'continiue', imp = 'y')
            outstr = cl.state
        

        return outstr

