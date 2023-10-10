# Copyright (c) Siman Development Team.
# Distributed under the terms of the GNU License.
from siman.core.calculation import Calculation
from siman.set_functions import InputSet, qe_keys
from siman.small_functions import makedir, list2string
from siman.header import printlog
from siman import header
from siman.functions import (
    read_vectors,
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
import os, copy, glob, shutil, sys

# Some functions are left as placeholders and will need to be implemented as needed.


class CalculationQE(Calculation):
    def __init__(self, inset=None, iid=None, output=None):
        super(CalculationQE, self).__init__(inset, iid, output)
        self.len_units = "Angstrom"
        self.calculator = "qe"
        self.input_instance = {}
        self.input_params = {}
        self.list_pot = []
        self.list_tmp = []
        printlog("Attention! This calculator is in a test mode\n")

        self.list_e_sigma0 = []

    def write_structure(
        self,
        name_of_output_file,
        type_of_coordinates="dir",
        option=None,
        prevcalcver=None,
        path=None,
        state="init",
    ):
        if path is None:
            path = self.dir
        if state == "init":
            self.st = self.init
        elif state == "end":
            self.st = self.end
        else:
            raise RuntimeError

        path2poscar = os.path.join(path, name_of_output_file)

        # makedir(filename)
        self.st.write_espresso(filename=path2poscar)
        # self.list_pos.append(path2poscar)
        return path2poscar

    def set_output_filenames(self, out_name, version):
        self.path["output"] = self.dir + self.name + ".log"

    def calculate_nbands(self, curset, path_to_potcar=None, params=None):
        pass  # Placeholder

    def actualize_set(self, curset=None, params=None):
        for section in self.input_params:
            for param in self.input_params[section]:
                if param.lower() in [item for item in self.set.params]:
                    self.input_instance.input_params[section][param] = self.set.params[
                        param
                    ]
                    print(f"Updated: {section} {param} => {self.set.params[param]}")
        pass  # Placeholder

    def add_potcar(self):
        """Generate POTCAR and creates list of psp to copy"""
        path2potcar = f"{self.dir}/qe_input.potcar.in"
        path2pot = header.PATH2POTENTIALS

        with open(path2potcar, "w") as f:
            f.write("ATOMIC_SPECIES \n")
            for z, el in zip(self.init.znucl, self.init.typat):
                f.write(f"{element_name_inv(z)} {z} {element_name_inv(z)}.upf \n")
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
        list_to_copy += self.list_pot
        list_to_copy.extend(glob.glob(os.path.join(self.dir, "*POSCAR*")))
        # list_to_copy+=self.list_pos
        # print(self.list_pos)
        if "up" in update:
            # printlog('Files to copy:', list_to_copy)
            push_to_server(
                list_to_copy,
                self.project_path_cluster + "/" + self.dir,
                self.cluster_address,
            )

    def download(self, load):
        path_to_outcar = self.path["output"]
        self.get_file(os.path.basename(path_to_outcar), up=load)
        return path_to_outcar

    def make_incar(self):
        path2input = f"{self.dir}/INCAR"
        path2incar = f"{self.dir}/qe_input.incar.in"
        # Generate Head (system/electron/ion)
        try:
            self.read_input_file(input_filename="./INCAR")
            printlog("QE input  found")
        except:
            printlog("QE input not found; attempting to create one from scratch")
            self._init_defualt()
        self.input_params["system"]["nat"] = self.st.natom
        self.input_params["system"]["ntyp"] = len([z for z in self.st.znucl])
        self.write_input_file(output_filename=path2incar)
        self.list_tmp.append(path2incar)

        # Generate INPUT
        print(self.list_tmp[::-1])
        with open(path2input, "w") as outfile:
            for fname in self.list_tmp[::-1]:
                with open(fname) as infile:
                    outfile.write(infile.read())
                os.remove(fname)
        return [path2input]

    def _init_defualt(self):
        default_input_params = {
            "control": {
                "calculation": '"scf"',
                "restart_mode": '"from_scratch"',
                "prefix": '"lno"',
                "pseudo_dir": '"./"',
                "outdir": '"./outdir"',
                "tprnfor": ".TRUE.",
            },
            "system": {
                "ibrav": "0",
                "nat": "1",
                "ntyp": "1",
                "ecutwfc": "30",
                "occupations": "'smearing'",
                "smearing": "'marzari-vanderbilt'",
                "degauss": "0.01",
                # 'nspin':'2',
                # 'starting_magnetization(2)' : '0.5',
            },
            "electrons": {},
            "ions": {},
            "cell": {"cell_dofree": "'ibrav'"},
        }
        self.input_params = default_input_params.copy()

    def update_params(self, section, key, value):
        if section not in self.input_params:
            raise ValueError(
                f"Section '{section}' does not exist in the input parameters."
            )
        self.input_params[section][key] = value

    def write_input_file(self, output_filename="./scf.in"):
        with open(output_filename, "w") as f:
            for section, params in self.input_params.items():
                f.write(f"&{section}\n")
                for key, value in params.items():
                    f.write(f"    {key} = {value}\n")
                f.write("/\n")

    def read_input_file(self, input_filename="./scf.in"):
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
        print(self.input_params)

    def read_results(
        self,
        load="",
        out_type="",
        voronoi=False,
        show="",
        choose_outcar=None,
        alkali_ion_number=None,
        only_load=False,
    ):
        """
        Download and read QE output file

        ###INPUT:
            - load (str) - 'x' - download xml, o - download outcar and contcar, un - read unfinished
            - show (str) - print additional information
                alkali_ion_number - show mag around this ion
            - choose_outcar - see description in res_loop(), from 1
            # - out_type - controls the return string
            #     see in code, add here
            #     also controls reading of OUTCAR
                'xcarts' read xcart every relaxation step and write into self.end.list_xcart

            - only_load (bool) - if true - only load the files (used for database)

        ###RETURN:


        ###DEPENDS:
        TODO:
        please split into outcar parser, downloader, and checker-formatter

        """

        # print (choose_outcar, hasattr(self, 'associated_outcars'), self.associated_outcars)
        join = os.path.join
        dirname = os.path.dirname

        if header.show:
            show += header.show
        # print('debug', self.path["output"])
        if not hasattr(self, "dir"):
            self.dir = os.path.dirname(self.path["output"])

        if (
            choose_outcar
            and hasattr(self, "associated_outcars")
            and self.associated_outcars
            and len(self.associated_outcars) >= choose_outcar
            and len(self.associated_outcars) > 1
        ):
            # print ('associated outcars = ',self.associated_outcars)
            printlog("read_results(): choose_outcar", choose_outcar)
            path_to_outcar = join(
                dirname(self.path["output"]), self.associated_outcars[choose_outcar - 1]
            )
            printlog(self.associated_outcars)
        else:
            path_to_outcar = self.path["output"]

        # print('debug', self.project_path_cluster+path_to_outcar)
        printlog(
            "read_results() path to outcar", self.project_path_cluster + path_to_outcar
        )
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
        # print('      printlog')
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

        self.state = self.check_state_qe(path_to_outcar)

        print(self.state)
        if "f" in self.state:
            outst = self.output_paraser(
                load=load,
                out_type=out_type,
                show=show,
                voronoi=voronoi,
                path_to_outcar=path_to_outcar,
            )
        else:
            cl = self
            try:
                os.rename(cl.path["output"], cl.path["output"] + "_unfinished")
                printlog(
                    "read_results():",
                    cl.id,
                    "is unfinished, continue:",
                    cl.dir,
                    cl.cluster["address"],
                    imp="y",
                )
                cl.state = "5. Unfinished"
            except:
                printlog(
                    "read_results():",
                    cl.id,
                    "probably was not submitted:",
                    cl.dir,
                    imp="y",
                )
            outst = "here"
        return outst

    def output_paraser(
        self,
        load="",
        out_type="",
        show="",
        voronoi="",
        path_to_outcar="",
        path_to_contcar="",
    ):
        # self - CalculationQE class
        with open(path_to_outcar, "rb") as outcar:
            printlog("Start reading from " + path_to_outcar, imp="n")
            # outcarlines = outcar.readlines()
            text = outcar.read().decode(errors="replace")
            outcarlines = str(text).split("\n")
        for line in outcarlines:
            if "!" in line:
                self.energy_sigma0 = float(line.split()[4])
                self.e0 = self.energy_sigma0
                self.list_e_sigma0.append(self.energy_sigma0)
                outst = line.split()[4]
        outst += " Ry "
        return outst

    def check_state_qe(self, path2output=""):
        try:
            with open(path2output) as f:
                for lines in f.readlines()[::-1]:
                    if "JOB DONE." in lines:
                        state = "finished"
                        break
                    else:
                        state = "running"
        except:
            printlog("No output file found")
            state = "no ouputfile"
        return state
