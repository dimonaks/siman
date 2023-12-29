# Copyright (c) Siman Development Team.
# Distributed under the terms of the GNU License.
from siman.core.calculation import Calculation
from siman.set_functions import InputSet, qe_keys
from siman.small_functions import makedir, list2string
from siman.core.structure import Structure
from siman.header import printlog
from siman import header
from siman.functions import (read_vectors,
                             read_list,
                             words,
                             read_string,
                             element_name_inv,
                             invert,
                             calculate_voronoi,
                             get_from_server,
                             push_to_server,
                             run_on_server,
                             smoother,
                             file_exists_on_server,
                             check_output,
                             )
import os
import copy
import glob
import shutil
import sys
import gzip
# Some functions are left as placeholders and will need to be implemented as needed.


class CalculationQE(Calculation):
    """CalculationQE backend class to run Quantum Espresso DFT code within Siman 
            prototype  - CalculationVasp
    """

    def __init__(self, inset=None, iid=None, output=None):
        """
        Args:
            inset (_type_, optional): input set. Defaults to None.
            iid (_type_, optional): index name. Defaults to None.
            output (_type_, optional): . Defaults to None.
        """
        super(CalculationQE, self).__init__(inset, iid, output)

        self.len_units = "Angstrom"
        self.calculator = "qe"
        # tmp solution
        self.input_instance = {}
        self.input_params = {}
        self.list_pot = []
        self.list_tmp = []

        printlog("\n Attention! This calculator is in a test mode!\n")

        self.list_e_sigma0 = []
        self.maxforce_list = []
        self.ry2ev = 13.6057039763
        self.init = Structure()
        self.end = Structure()

    def write_structure(self, name_of_output_file, type_of_coordinates="dir", path=None, option=None,
                        prevcalcver=None, state="init",):
        """Saves relaxed structure

        Args:
            name_of_output_file (string): output file
            type_of_coordinates (str, optional): Dircet or else. Defaults to "dir".
            state (str, optional): State of the structure relaxation. Defaults to "init".
            path (str, optional): path to the output file. Defaults to None.

            Legacy from the same method in CalculationVasp:
            option (_type_, optional): _description_. Defaults to None.
            prevcalcver (_type_, optional): _description_. Defaults to None.


        Raises:
            RuntimeError:

        Returns:
            _type_: _description_
        """
        if path is None:
            path = self.dir
        if state == "init":
            self.st = self.init
        elif state == "end":
            self.st = self.end
        else:
            raise RuntimeError

        path2poscar = os.path.join(path, name_of_output_file)
        self.st.write_espresso(filename=path2poscar)

        return path2poscar

    def set_output_filenames(self, out_name, version):
        """ Defines the name of the ouput file

        Args:
            out_name (string): custom name of the output file
            version(string): legacy
        """
        if version:
            self.path["output"] = self.dir+out_name
        else:
            self.path["output"] = self.dir + self.name + ".out"  # qe
        # tmp solution
        # self.path["output"] = self.dir + self.name + ".out"
        #
        return self.path["output"]

    def calculate_nbands(self, curset, path_to_potcar=None, params=None):
        """Placeholder
        """
        pass

    def actualize_set(self, curset=None, params=None):
        """Updates current input set / placeholder

        Args:
            curset (dictionary, optional): Current input set. Defaults to None.
            params (dictionary, optional): Dictionary with parameters to update. Defaults to None.
        """
        for section in self.input_params:
            for param in self.input_params[section]:
                if param.lower() in [item for item in self.set.params]:
                    self.input_instance.input_params[section][param] = self.set.params[
                        param
                    ]
                    print(
                        f"Updated: {section} {param} => {self.set.params[param]}")
        pass  # Placeholder

    def add_potcar(self):
        """Generate POTCAR and creates list of psp to copy"""

        path2potcar = f"{self.dir}/qe_input.potcar.in"
        path2pot = header.PATH2POTENTIALS

        with open(path2potcar, "w") as f:
            f.write("ATOMIC_SPECIES \n")
            for z, el in zip(self.init.znucl, self.init.typat):
                f.write(
                    f"{element_name_inv(z)} {z} {element_name_inv(z)}.upf \n")
                self.list_pot.append(f"{path2pot}/{element_name_inv(z)}.upf")
        self.list_tmp.append(path2potcar)

        # legacy
        self.path["potcar"] = header.PATH2POTENTIALS

    def make_kpoints_file(self):
        """Generate file with kpoints"""

        path2kcar = f"{self.dir}/qe_input.kcar.in"
        # self.init.ngkpt = [3,3,3]
        with open(path2kcar, "w") as f:
            f.writelines("K_POINTS automatic \n")
            f.writelines(f" ".join(str(nk) for nk in self.init.ngkpt))
            # f.writelines(f' '.join(str(nk) for nk in self.init.ngkpt))
            f.writelines(" 0 0 0 \n")
            f.writelines("\n")

        self.list_tmp.append(path2kcar)
        # legacy
        return []

    def copy_to_cluster(self, list_to_copy, update):
        """Method to upload files to the server

        Args:
            list_to_copy (list): list of files to copy
            update (str): mode
        """
        list_to_copy += self.list_pot
        list_to_copy.extend(glob.glob(os.path.join(self.dir, "*POSCAR*")))
        if "up" in update:
            push_to_server(
                list_to_copy,
                self.project_path_cluster + "/" + self.dir,
                self.cluster_address,
            )

    def download(self, load):
        """Method to download output files

        Args:
            load (str): mode

        Returns:
            path_to_outcar: path to output card 
        """
        self.get_file(os.path.basename(self.path["output"]), up=load)

        return path_to_outcar

    def make_incar(self):
        """Method to generate part of the QE input

        Returns:
                incar_list : list "incar" files
        """
        incar_list = []
        setseq = [self.set]
        if hasattr(self.set, 'set_sequence') and self.set.set_sequence:
            for s in self.set.set_sequence:
                setseq.append(s)
        nsets = len(setseq)
        for i, curset in enumerate(setseq):
            if nsets == 1:
                name_mod = ''
            else:
                name_mod = curset.ise+'.'
            path2input = f"{self.dir}/{name_mod}INCAR"
            path2incar = f"{self.dir}/{name_mod}qe_input.incar.in"
            self.input_params = curset.params
            self.input_params["system"]["nat"] = self.st.natom
            self.input_params["system"]["ntyp"] = len(
                [z for z in self.st.znucl])
            self.write_input_file(output_filename=path2incar)
            self.list_tmp.append(path2incar)

            # Generates INPUT
            with open(path2input, "w") as outfile:
                for fname in self.list_tmp[::-1]:
                    with open(fname) as infile:
                        outfile.write(infile.read())
                    os.remove(fname)
            incar_list.append(path2input)
        return incar_list

    def write_input_file(self, output_filename="./scf.in"):
        """Writes down the the input file in qe format

        Args:
            output_filename (str, optional): name of the ouput file. Defaults to "./scf.in".
        """
        # tmp solution
        if "KPOINTS" in self.input_params.keys():
            del self.input_params["KPOINTS"]
        elif "KSPACING" in self.input_params.keys():
            del self.input_params["KSPACING"]

        with open(output_filename, "w") as f:
            for section, params in self.input_params.items():
                f.write(f"&{section}\n")
                for key, value in params.items():
                    f.write(f"    {key} = {value}\n")
                f.write("/\n")

    def read_input_file(self, input_filename="./scf.in"):
        """Reads QE input file

        Args:
            input_filename (str, optional): Name of the input file. Defaults to "./scf.in".
        """
        with open(input_filename, "r") as f:
            lines = f.readlines()
            section = None  # To keep track of the current section
            for line in lines:
                line = line.strip()
                # Check for the beginning of a new section
                if line.startswith("&"):
                    section = line.strip("&").strip()
                    self.input_params[section] = {}
                elif line.startswith("/"):
                    section = None
                elif section:
                    key, value = map(str.strip, line.split("="))
                    self.input_params[section][key] = value

    def read_results(self, load="", out_type="", voronoi=False, show="", choose_outcar=None,
                     alkali_ion_number=None, only_load=False):
        """Method to pull and read QE inputfiles
                loosely based on same method from CalculationsVasp
                requires fix of associated_outcars

        Args:
         ---- require fix

        Returns:
            state:  str
        """
        join = os.path.join
        dirname = os.path.dirname

        if header.show:
            show += header.show

        if not hasattr(self, "dir"):
            self.dir = os.path.dirname(self.path["output"])

        if (choose_outcar and hasattr(self, "associated_outcars") and self.associated_outcars
                and len(self.associated_outcars) >= choose_outcar and len(self.associated_outcars) > 1):
            # print('associated outcars = ', self.associated_outcars)
            printlog("read_results(): choose_outcar", choose_outcar)
            path_to_outcar = join(
                dirname(self.path["output"]), self.associated_outcars[choose_outcar - 1])
            printlog(self.associated_outcars)
        else:
            path_to_outcar = self.path["output"]

        printlog("read_results() path to outcar",
                 self.project_path_cluster + path_to_outcar)
        # sys.exit()

        if not os.path.exists(path_to_outcar):
            load = load + "o"

        """Copy from server """

        printlog("The load flag is ", load)
        # print('debug',load, self.cluster_address)

        if "o" in load and hasattr(self, "cluster_address"):
            # print('debug')
            if "un2" in load:
                out_name = os.path.basename(path_to_outcar)
                path_to_outcar = path_to_outcar.replace(out_name, "scf.out")
            files = [self.project_path_cluster + "/" + path_to_outcar]
            # print(load)
            # print(files)
            for file in files:
                debug = self.get_file(os.path.basename(file), up=load)
                # print(debug)
        if os.path.exists(path_to_outcar):
            outcar_exist = True
        else:
            outcar_exist = False
            path_to_zip = path_to_outcar + ".gz"
            if os.path.exists(path_to_zip):
                with gzip.open(path_to_zip, "rb") as f_in:  # unzip OUTCAR
                    with open(path_to_outcar, "wb") as f_out:
                        shutil.copyfileobj(f_in, f_out)
                if os.path.exists(path_to_outcar):
                    outcar_exist = True

        """Start reading """

        self.state = self.check_state_qe(path2output=path_to_outcar)

        # print(path_to_outcar, self.state)
        if "4" in self.state:
            outst = self.output_paraser(load=load, out_type=out_type,
                                        show=show, voronoi=voronoi, path_to_outcar=path_to_outcar,
                                        )
        else:
            try:
                os.rename(self.path["output"],
                          self.path["output"] + "_unfinished")
                printlog(
                    "read_results():", cl.id, "is unfinished, continue:", self.dir, self.cluster["address"], imp="y",)
                self.state = "4. Unfinished!"
                outst = self.state
            except:
                printlog("read_results():", self.id,
                         "probably was not submitted:", self.dir, imp="y",)
                outst = "6. not submitted"

        return outst

    def check_state_qe(self, path2output=""):
        """Checks the state of quantum espresso calculations

        Args:
            path2output (str, optional): path to the output file. Defaults to "".

        Returns:
            state:  str
        """
        try:
            with open(path2output) as f:
                for lines in f.readlines()[::-1]:
                    if "JOB DONE." in lines:
                        state = "4. Finished"
                        break
                    else:
                        state = "5. Unfinished!"
        except:
            printlog("No output file found")
            state = "5. no outputfile"
        return state

    def output_paraser(self, load="", out_type="", show="", voronoi="", path_to_outcar=""):
        """QE output specific parser. For now based on ase.io.read() method.
                Magentic moments are not working.
                Have to be rewritten latter.

        Args:
            load (str, optional): mode. Defaults to "".
            out_type (str, optional): legacy. Defaults to "".
            show (str, optional): mode show output. Defaults to "".
            voronoi (str, optional): legacy. Defaults to "".
            path_to_outcar (str, optional): path to the ouput file. Defaults to "".


        Returns:
            outst: str (require fix)
        """

        self.maxforce_list = []
        self.force_prefix = ''
        self.ry2ev = 13.6057039763
        self.au2angst = 0.529177249
        self.ry_per_au2ev_per_angst = header.Ha_Bohr_to_eV_A/2
        self.end = Structure()
        try:
            from ase.io import read as aseread
        except:
            printlog(
                'No module ase is found, switching to parser with limited functionallity')
            outst = self._output_paraser(load=load, out_type=out_type,
                                         show=show, voronoi=voronoi,
                                         path_to_outcar=path_to_outcar,
                                         path_to_contcar=path_to_contcar,)
            return outst

        output = aseread(path_to_outcar, format='espresso-out')

        self.force = output.get_forces()
        self.maxforce_list = [max(abs(item)) for item in self.force]
        self.end = self.qe_structure(self.end, output)
        self.e0 = output.get_total_energy()
        self.energy_sigma0 = self.e0
        self.list_e_sigma0.append(self.e0)

        # this is a placeholder (require fix)
        try:
            self.magn1 = output.get_magnetic_moments()
        except:
            self.magn1 = 0
        try:
            self.magn2 = output.get_magnetic_moments()
        except:
            self.magn2 = 0

        outst = f"{self.e0} eV | {self.end.vol} A^3"

        return outst

    def qe_structure(self, st, qe):
        """Reads structure from QE output

        Args:
            st (SimanObject): Structure object.
            qe (aseObject): readout of QE output file.

        Returns:
            st: SimanObject
        """
        st.name = qe.symbols
        st.rprimd = qe.cell
        st.natom = len(qe.arrays['numbers'])
        st.znucl = set(qe.arrays['numbers'])
        st.xcart = qe.arrays['positions']
        st.ntypat = len(st.znucl)
        st.vol = qe.get_volume()

        return st

    def _output_paraser(self, load="", out_type="", show="", voronoi="",
                        path_to_outcar="", path_to_contcar="",):
        """ In case ase module is not found code reffers to this method to read some parameters from
                QE output.

        Args:
            load (str, optional): mode. Defaults to "".
            out_type (str, optional): legacy. Defaults to "".
            show (str, optional): mode show output. Defaults to "".
            voronoi (str, optional): legacy. Defaults to "".
            path_to_outcar (str, optional): path to the ouput file. Defaults to "".

        Returns:
            outstate:  str (require fix of the format)
        """

        with open(path_to_outcar, "rb") as outcar:
            printlog("Start reading from " + path_to_outcar, imp="n")
            # outcarlines = outcar.readlines()
            text = outcar.read().decode(errors="replace")
            outcarlines = str(text).split("\n")
        self.maxforce_list = []
        self.force_prefix = ''
        self.ry2ev = 13.6057039763
        self.au2angst = 0.529177249
        self.ry_per_au2ev_per_angst = header.Ha_Bohr_to_eV_A/2

        for il, line in enumerate(outcarlines):
            if 'Total force' in line:
                self.force = float(line.split()[3])
            if 'Forces acting on atoms' in line:
                for tmpline in range(2, self.st.natom+2):
                    max_force = max(
                        [float(outcarlines[il+tmpline].split()[j])*self.ry_per_au2ev_per_angst*1000 for j in [6, 7, 8]])
                    self.maxforce_list.append(max_force)
            if "unit-cell volume" in line:
                self.st.vol = float(line.split()[3])

            if "!" in line:
                self.energy_sigma0 = float(line.split()[4])*self.ry2ev
                self.e0 = self.energy_sigma0
                self.list_e_sigma0.append(self.energy_sigma0)
                outst = str(float(line.split()[4])*self.ry2ev)
        outst += " eV "
        return outst
