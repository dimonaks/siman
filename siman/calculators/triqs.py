
from typing import Any
from pymatgen.core import Element


from typing import Any

from pymatgen.core import Element

def dc_value(n, U, J):
    """
    Fully Localized Limit (FLL) double-counting potential.

    Parameters
    ----------
    n : float
        Total d-shell occupation.
    U : float
        Hubbard U (eV).
    J : float
        Hund's coupling (eV).

    Returns
    -------
    float
        Double-counting potential V_DC (eV).
    """
    return U * (n - 0.5) - J * (n / 2 - 0.5)


def write_parameters_for_vasp2triqs(cl, params=None):
    """
    Write a plo.cfg file for the VASP-to-TRIQS converter.

    Parameters
    ----------
    cl
        SIMAN calculation object. The atomic structure is taken from
        ``cl.init``.

    params : dict, optional
        Parameters overriding the default settings.

        The projected elements can be specified as a list:

            params = {
                "el": ["Ni", "O"],
            }

        In this case, the default settings are:

        - d-block elements: LSHELL = 2, CORR = True
        - p-block elements: LSHELL = 1, CORR = False

        The CORR value can be specified explicitly using a dictionary:

            params = {
                "el": {
                    "Ni": True,
                    "O": False,
                },
            }

        Shell-specific parameters can also be supplied as dictionaries:

            params = {
                "el": {"Ni": True, "O": False},
                "LSHELL": {"Ni": 2, "O": 1},
                "IONS": {"Ni": [1], "O": [2, 3]},
                "TRANSFORM": {
                    "Ni": [...],
                    "O": [...],
                },
            }

    Returns
    -------
    str
        Path to the generated plo.cfg file.
    """
    params = {} if params is None else params.copy()

    convert_config_filename = cl.dir + "/plo.cfg"
    python_script_filename = cl.dir + "/convert.py"
    shell_script_filename = cl.dir + "/convert.sh"
    plot_script_filename = cl.dir + "/plot.py"

    cfg = {
        "DOSMESH": "-10 10 3001",
        "EWINDOW": "-10 10",
        "BANDS": "5 14",
    }

    # Parameters that may differ between individual shells.
    elements_param = params.pop("el", ["Ni"])
    lshell_param = params.pop("LSHELL", {})
    ions_param = params.pop("IONS", {})
    corr_param = params.pop("CORR", {})
    transform_param = params.pop("TRANSFORM", {})

    # All remaining parameters are treated as general/group parameters.
    cfg.update(params)

    # A dictionary passed through "el" defines both the elements and CORR.
    if isinstance(elements_param, dict):
        elements = list(elements_param)
        corr_from_el = elements_param
    else:
        elements = list(elements_param)
        corr_from_el = {}

    def get_element_value(value, symbol, default=None):
        """Return an element-specific value or a common value."""
        if isinstance(value, dict):
            return value.get(symbol, default)

        if value is None or value == {}:
            return default

        return value

    def get_default_lshell(symbol):
        """Select the default shell from the element block."""
        element = Element(symbol)

        if element.block == "d":
            return 2

        if element.block == "p":
            return 1

        raise ValueError(
            f"Cannot automatically select LSHELL for {symbol}. "
            f"Specify it explicitly, for example "
            f"'LSHELL': {{'{symbol}': 0}}."
        )

    def get_default_corr(lshell):
        """Enable correlations for d shells and disable them otherwise."""
        return lshell == 2

    def make_identity_transform(lshell):
        """Create an identity transformation matrix for a given shell."""
        size = 2 * lshell + 1

        return [
            "  ".join(
                "1.0" if row == column else "0.0"
                for column in range(size)
            )
            for row in range(size)
        ]

    def get_ion_numbers(symbol):
        """Return one-based atom indices for the selected element."""
        atomic_number = Element(symbol).Z

        return [
            index + 1
            for index in cl.init.get_specific_elements([atomic_number])
        ]

    shells = []

    for symbol in elements:
        lshell = get_element_value(
            lshell_param,
            symbol,
            get_default_lshell(symbol),
        )

        corr = corr_from_el.get(
            symbol,
            get_element_value(
                corr_param,
                symbol,
                get_default_corr(lshell),
            ),
        )

        ions = get_element_value(ions_param, symbol)

        if ions is None:
            ions = get_ion_numbers(symbol)

        if isinstance(ions, (list, tuple, set)):
            ions = " ".join(map(str, ions))
        else:
            ions = str(ions)

        transform = get_element_value(transform_param, symbol)

        if transform is None:
            transform = make_identity_transform(lshell)

        shells.append(
            {
                "element": symbol,
                "lshell": lshell,
                "ions": ions,
                "corr": corr,
                "transform": transform,
            }
        )

    # Group SHELLS must contain shell numbers, e.g. "1 2", not their count.
    shell_numbers = " ".join(
        map(str, range(1, len(shells) + 1))
    )

    with open(convert_config_filename, "w") as file:
        file.write(
            "[General]\n"
            f"DOSMESH = {cfg['DOSMESH']}\n\n"
            "[Group 1]\n"
            f"SHELLS = {shell_numbers}\n"
            f"EWINDOW = {cfg['EWINDOW']}\n"
            f"BANDS = {cfg['BANDS']}\n"
        )

        for shell_number, shell in enumerate(shells, start=1):
            orbital_name = {
                0: "s",
                1: "p",
                2: "d",
                3: "f",
            }.get(
                shell["lshell"],
                f"l={shell['lshell']}",
            )

            file.write(
                f"\n[Shell {shell_number}] "
                f"# {shell['element']} {orbital_name} shell\n"
                f"LSHELL = {shell['lshell']}\n"
                f"IONS = {shell['ions']}\n"
                f"CORR = {shell['corr']}\n"
            )

            transform = shell["transform"]

            if transform:
                file.write(f"TRANSFORM = {transform[0]}\n")

                indent = " " * len("TRANSFORM = ")

                for row in transform[1:]:
                    file.write(f"{indent}{row}\n")




    with open(python_script_filename, "w") as f:
        f.write(f"""from triqs_dft_tools.converters import VaspConverter
import triqs_dftkit.vasp.plovasp.converter as plo_converter

plo_converter.generate_and_output_as_text('plo.cfg', vasp_dir='./')

Converter = VaspConverter(filename='vasp')
Converter.convert_dft_input()
""")

    with open(plot_script_filename, "w") as f:
        f.write(f"""import py4vasp
import matplotlib.pyplot as plt
import numpy as np

# get VASP DOS
calc = py4vasp.Calculation.from_file('vaspout.h5')
dft_dos = calc.dos.to_dict()

# load projector (PLO) DOS from TRIQS
plo_dos_Ni = np.loadtxt('pdos_0_0.dat')

fig = plt.figure(figsize=(7,4), dpi=150)

plt.plot(dft_dos['energies'], dft_dos['total'], label=r'VASP total DOS')

# sum PLO DOS for all orbitals
plt.plot(plo_dos_Ni[:,0],np.sum(plo_dos_Ni[:,1:],axis=1), label=r'PLO Ni_d$')

plt.xlim(-10,4)
plt.ylim(0,12)
plt.xlabel(r'epsilon - E_F$ (eV)')
plt.ylabel(r'states/eV')
plt.legend()
#plt.show()
plt.savefig('dos.png')
""")


    with open(shell_script_filename, "w") as f:
        f.write(f"""module load q-ch/triqs/4.0.0
~/miniconda3/envs/triqs_4_build/bin/python convert.py
~/miniconda3/envs/triqs_4_build/bin/python plot.py
""")



    return convert_config_filename, python_script_filename, plot_script_filename, shell_script_filename



from pathlib import Path
from typing import Any


def write_dmft_config_file(filename, params = None, skip_dft = 1):

    config = {
        "general": {
            "seedname": "vasp",
            "jobname": "new_job",
            "prec_mu": 0.001,
            "n_iw": 501,
            "n_tau": 10001,
            "enforce_off_diag": False,
            "block_threshold": 0.001,
            "h_int_type": "density_density",
            "h_int_basis": "vasp",
            "U": 8.0,
            "J": 1.0,
            "beta": 20,
            "n_iter_dmft": 12,
            "dc_type": 0,
            "dc": True,
            "dc_dmft": False,
            "h5_save_freq": 1,
            "set_rot": None,
            "calc_energies": True,
            "magnetic": None,
            "magmom": None,
            "csc": None,
            "n_iter_dmft_first": None,
            "n_iter_dmft_per": None,
            "load_sigma": None,
            "path_to_sigma": None,
            "sigma_mix":None,

        },
        "solver": {
            "type": "ctseg",
            "length_cycle": 500,
            "n_warmup_cycles": 1e4,
            "n_cycles_tot": 1e7,
            "perform_tail_fit": True,
            "fit_max_moment": 4,
            "fit_min_w": 20,
            "fit_max_w": 30,
        },

        "dft": {
            "plo_cfg": "plo.cfg",
            "dft_code": "vasp",
            "projector_type": "plo",
            "n_iter": 4,
            "n_cores": 16,
            "mpi_env": "default",
            "mpi_exe": "/opt/intel/oneapi/mpi/2021.12/bin/mpirun",
            "dft_exec": "/opt/cms/q-ch/vasp/6.6.0_w90/bin/vasp_std",
        }, 
        "advanced": {
            "dc_fixed_value": None,
        },
    }

    for key, value in (params or {}).items():
        # if value is None:
            # continue
        # print(key, value)
        if "." in key:
            section, parameter = key.split(".", 1)
        else:
            matches = [
                section for section, params2 in config.items()
                if key in params2
            ]
            if len(matches) != 1:
                raise KeyError(
                    f"Parameter {key!r} is unknown or ambiguous; "
                    "use 'section.parameter'"
                )
            section, parameter = matches[0], key

        if section not in config or parameter not in config[section]:
            raise KeyError(f"Unknown parameter: {section}.{parameter}")

        config[section][parameter] = value

    if "jobname" not in (params or {}) and "general.jobname" not in (params or {}):
        g = config["general"]
        s = config["solver"]
        g["jobname"] = (
            f'U{g["U"]:.1f}-'
            f'J{g["J"]:.1f}-'
            f'beta{g["beta"]}-'
            f'qmc{s["n_cycles_tot"]:.0e}'
        )

    def toml_value(value: Any) -> str:
        if isinstance(value, bool):
            return str(value).lower()
        if isinstance(value, str):
            return f'"{value}"'
        if isinstance(value, (list, tuple)):
            # return ", ".join(map(str, value))+", "
            return str(value)
        return str(value)

    lines = []
    for section, parameters in config.items():
        if section == 'dft' and skip_dft:
            continue
        lines.append(f"[{section}]")
        lines.extend(
            f"{key} = {toml_value(value)}"
            for key, value in parameters.items()
            if value is not None
        )
        # print(lines)
        lines.append("")

    Path(filename).write_text("\n".join(lines), encoding="utf-8")

    return filename, g["jobname"]


def write_dmft_config(cl,    params = None, up = 'up2'):
    """Write dmft_config.toml, overriding only supplied parameters."""

    os_config_filename = cl.dir+"/dmft_config.toml"
    csc_config_filename = cl.dir+"/dmft_config_csc.toml"
    slurm_script_filename = cl.dir+"/triqs.run"
    slurm_csc_script_filename = cl.dir+"/csc_triqs.run"
    plot_self_script_filename = cl.dir+"/self_plot.py"
    maxent_script_filename = cl.dir+"/gimp_maxent.py"
    spectral_plot_script_filename = cl.dir+"/spec_plot.py"

    # print(params)
    cl.path['triqs_out'] = cl.dir+"/out"
    cl.path['csc_triqs_out'] = cl.dir+"/csc/out"

    _, jobname = write_dmft_config_file(os_config_filename, params = params)

    params.update(
        csc=True,
        dc_dmft=True, #added
        n_iter_dmft_first=2,
        n_iter_dmft_per=2,
        n_iter_dmft=20,
        load_sigma=True,
        sigma_mix = 0.5, #added
        path_to_sigma=jobname+"/vasp.h5"
    )

    write_dmft_config_file(csc_config_filename, params = params, skip_dft = 0)

    if 'csc' in up:
        vasph5 = 'vasp.h5'

    else:
        vasph5 = jobname+'/vasp.h5'

    cl.jobname = jobname
    with open(maxent_script_filename, "w") as f:
        f.write(f"""import solid_dmft.postprocessing.maxent_gf_imp as gimp_maxent

Aimp_csc_w = gimp_maxent.main(external_path='{vasph5}',
                              omega_min=-30, omega_max=15,
                              maxent_error=0.03,
                              sum_spins=True,
                              n_points_maxent=400,
                              n_points_alpha=20)
""")



    with open(slurm_script_filename, "w") as f:
        f.write(f"""#!/bin/bash   
#SBATCH -J triqs_{cl.id[0]} 
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -o {cl.project_path_cluster+'/'+cl.dir}/sbatch.out
#SBATCH -e {cl.project_path_cluster+'/'+cl.dir}/sbatch.err

cd {cl.project_path_cluster+'/'+cl.dir}
source /etc/profile.d/modules.sh; module load q-ch/vasp/6.6.0_w90 
module load q-ch/triqs/4.0.0 
export PATH=$PATH:/home/aksenov//tools
touch RUNNING
#Basic run:
./convert.sh
mpirun solid_dmft 1>out 2>err
mpirun ~/miniconda3/envs/triqs_4_build/bin/python gimp_maxent.py 1>max_out 2>err_max
~/miniconda3/envs/triqs_4_build/bin/python self_plot.py
~/miniconda3/envs/triqs_4_build/bin/python spec_plot.py
rm RUNNING
""")


    with open(slurm_csc_script_filename, "w") as f:
        f.write(f"""#!/bin/bash   
#SBATCH -J csc_triqs_{cl.id[0]} 
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -o {cl.project_path_cluster+'/'+cl.dir}/csc/sbatch.out
#SBATCH -e {cl.project_path_cluster+'/'+cl.dir}/csc/sbatch.err

export I_MPI_FABRICS=shm
cd {cl.project_path_cluster+'/'+cl.dir + '/csc'}
source /etc/profile.d/modules.sh; module load q-ch/vasp/6.6.0_w90 
module load q-ch/triqs/4.0.0-intel
export PATH=$PATH:/home/aksenov//tools
export HDF5_USE_FILE_LOCKING=FALSE
touch RUNNING

#cp ../1.POSCAR POSCAR
#cp ../1.CHGCAR CHGCAR
#cp ../1.WAVECAR WAVECAR
#cp -r ../{jobname} .
cp ../dmft_config_csc.toml dmft_config.toml
cp ../INCAR ../KPOINTS ../POTCAR ../plo.cfg .
cp ../self_plot.py ../spec_plot.py ../gimp_maxent.py .
sed -i \
-e '/^[[:space:]]*ICHARG[[:space:]]*=/Id' \
-e '$aICHARG = 5' \
-e '/^[[:space:]]*AMIX[[:space:]]*=/Id' \
-e '$aAMIX = 0.08' \
-e '/^[[:space:]]*LSYNCH5[[:space:]]*=/Id' \
-e '$aLSYNCH5 = .TRUE.' \
INCAR

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

#mpirun -np 32 -genv I_MPI_PIN 0  solid_dmft 1>out 2>err
mpirun python gimp_maxent.py >max_out 2>err_max
python self_plot.py
python spec_plot.py
rm RUNNING
""")




    with open(spectral_plot_script_filename, "w") as f:
        f.write(rf"""# load the results
import numpy as np
import matplotlib.pyplot as plt

from h5 import HDFArchive


h5_file = '{vasph5}'

with HDFArchive(h5_file, 'r') as ar:
    A_imp_csc_w = ar['DMFT_results/last_iter/Aimp_maxent_0']


omega = np.asarray(A_imp_csc_w['mesh']).real
spectra = A_imp_csc_w['Aimp_w_line_fit']
block_names = list(spectra.keys())

print('Available MaxEnt blocks:', block_names)


# NiO-like case:
# t2g and eg are stored as separate blocks total_0 and total_2
if 'total_0' in block_names and 'total_2' in block_names:

    A_t2g = np.asarray(spectra['total_0'])[0, 0, :].real
    A_eg = np.asarray(spectra['total_2'])[0, 0, :].real

    print('Using separate t2g and eg blocks')


# LiNiO2-like case:
# all five d orbitals are stored in one 5x5 block
elif 'total_0' in block_names:

    A_full = np.asarray(spectra['total_0']).real

    if A_full.ndim != 3:
        raise ValueError(
            f'Expected a three-dimensional spectral array, '
            f'but got shape {{A_full.shape}}'
        )

    n_orb = A_full.shape[0]

    if A_full.shape[1] != n_orb:
        raise ValueError(
            f'Expected a square orbital matrix, but got shape {{A_full.shape}}'
        )

    if n_orb != 5:
        raise ValueError(
            f'Only one MaxEnt block was found, but it contains {{n_orb}} orbitals. '
            f'Automatic t2g/eg splitting is currently defined only for a '
            f'five-orbital d shell.'
        )

    # Orbital order for the present LiNiO2 projector basis
    t2g_indices = [0, 1, 3]
    eg_indices = [2, 4]

    # Average spectral function per orbital
    A_t2g = np.mean(
        [A_full[i, i, :] for i in t2g_indices],
        axis=0,
    )

    A_eg = np.mean(
        [A_full[i, i, :] for i in eg_indices],
        axis=0,
    )

    print('Using one 5x5 block')
    print('t2g orbitals:', t2g_indices)
    print('eg orbitals:', eg_indices)


else:
    raise KeyError(
        f'Could not identify spectral blocks. '
        f'Available blocks: {{block_names}}'
    )


fig, ax = plt.subplots(1, 1, figsize=(10, 5), dpi=100)

ax.plot(
    omega,
    A_t2g,
    label='Ni $t_{{2g}}$',
)

ax.plot(
    omega,
    A_eg,
    label='Ni $e_{{g}}$',
)

ax.set_xlim(-10, 10)
ax.set_ylim(bottom=0)

ax.set_ylabel(
    r'A$_{{\mathrm{{proj}}}}$ states/eV'
)
ax.set_xlabel(
    r'$\omega$ (eV)'
)

ax.legend()

fig.tight_layout()
fig.savefig(
    'spec.png',
    dpi=200,
)

            """)




    with open(plot_self_script_filename, "w") as f:
        f.write(rf"""import numpy as np

# plotting
import matplotlib.pyplot as plt

# TRIQS modules
from h5 import HDFArchive
from triqs.gf import *
from triqs.plot.mpl_interface import plt,oplot
  
def plot_sigma_iw(S_iw, ax, color, label='', marker='-o', subtract=True):
    # convert the mesh to a numpy array
    
    # Plot Matsubara self-energy Σ(iωn), characterizing correlation effects
    # and quasiparticle renormalization.
    
    mesh = mesh_to_np_arr(S_iw.mesh)
    mid = len(mesh)//2
    if subtract:
        ax[0].plot(mesh[mid:], S_iw.data[mid:,0,0].real-S_iw.data[-1,0,0].real, marker, color = color, label=label)
    else:
        ax[0].plot(mesh[mid:], S_iw.data[mid:,0,0].real, marker, color = color, label=label)
    ax[1].plot(mesh[mid:], -1*S_iw.data[mid:,0,0].imag, marker, color = color, label=label)

    ax[0].set_xlabel(r"$i \omega_n$ (eV)")
    ax[1].set_xlabel(r"$i \omega_n$ (eV)")
    ax[0].set_ylabel(r"$Re \Sigma (i \omega_n)$  (eV)")
    ax[1].set_ylabel(r"$- Im \Sigma (i \omega_n)$  (eV)")
  
    return

def plot_sigma_w(S_w, ax, color, label='', marker='-', subtract=True):

    #Plot real-frequency self-energy Σ(ω), describing energy renormalization
    #and finite quasiparticle lifetime.

    mesh = mesh_to_np_arr(S_w.mesh)
    mid = len(mesh)//2
    if subtract:
        ax[0].plot(mesh, S_w.data[:,0,0].real-S_w.data[mid,0,0].real, marker, color = color, label=label)
    else:
        ax[0].plot(mesh, S_w.data[:,0,0].real, marker, color = color, label=label)
    ax[1].plot(mesh, -1*S_w.data[:,0,0].imag, marker, color = color, label=label)
 
    ax[0].set_xlabel(r"$\omega$ (eV)")
    ax[1].set_xlabel(r"$\omega$ (eV)")
    ax[0].set_ylabel(r"$Re \Sigma (\omega)$  (eV)")
    ax[1].set_ylabel(r"$- Im \Sigma (\omega)$  (eV)")
    return



def plot_sigma_iw_orb(
    S_iw,
    ax,
    color,
    label='',
    marker='-',
    subtract=True,
    orbital=0,
):
    mesh = mesh_to_np_arr(S_iw.mesh)
    mid = len(mesh) // 2

    sigma = S_iw.data[:, orbital, orbital]

    if subtract:
        ax[0].plot(
            mesh[mid:],
            sigma[mid:].real - sigma[-1].real,
            marker,
            color=color,
            label=label,
        )
    else:
        ax[0].plot(
            mesh[mid:],
            sigma[mid:].real,
            marker,
            color=color,
            label=label,
        )

    ax[1].plot(
        mesh[mid:],
        -sigma[mid:].imag,
        marker,
        color=color,
        label=label,
    )


def plot_sigma_iw_orbitals(
    S_iw,
    orbitals,
    ax,
    color,
    label='',
    marker='-',
    subtract=True,
):
    mesh = mesh_to_np_arr(S_iw.mesh)
    mid = len(mesh) // 2

    sigma = sum(
        S_iw.data[:, i, i]
        for i in orbitals
    ) / len(orbitals)

    if subtract:
        real_part = sigma.real - sigma[-1].real
    else:
        real_part = sigma.real

    ax[0].plot(
        mesh[mid:],
        real_part[mid:],
        marker,
        color=color,
        label=label,
    )

    ax[1].plot(
        mesh[mid:],
        -sigma[mid:].imag,
        marker,
        color=color,
        label=label,
    )

    ax[0].set_xlabel(r'$i\omega_n$ (eV)')
    ax[1].set_xlabel(r'$i\omega_n$ (eV)')
    ax[0].set_ylabel(r'$\mathrm{{Re}}\,\Sigma(i\omega_n)$ (eV)')
    ax[1].set_ylabel(r'$-\mathrm{{Im}}\,\Sigma(i\omega_n)$ (eV)')





def mesh_to_np_arr(mesh):
    from triqs.gf import MeshImTime, MeshReFreq, MeshImFreq

    if isinstance(mesh, MeshReFreq):
        mesh_arr = np.linspace(mesh.w_min, mesh.w_max, len(mesh))
    elif isinstance(mesh, MeshImFreq):
        mesh_arr = np.linspace(mesh(mesh.first_index).imag, mesh(mesh.last_index).imag, len(mesh))
    elif isinstance(mesh, MeshImTime):
        mesh_arr = np.linspace(0, mesh.beta, len(mesh))
    else:
        raise AttributeError('input mesh must be either MeshReFreq, MeshImFreq, or MeshImTime')
 
    return mesh_arr



with HDFArchive('{vasph5}','r') as ar:
    S_os_iw = ar['DMFT_results/last_iter/Sigma_freq_0']

fig, ax = plt.subplots(1,2, figsize=(14,5), dpi=100, sharex=True)
ax = ax.reshape(-1)

#plot_sigma_iw(S_os_iw['up_0'],ax,color='C0',label='$t_{{2g}}$',marker='-')
#plot_sigma_iw(S_os_iw['up_2'],ax,color='C1',label='$e_{{g}}$',marker='-')



#plot_sigma_iw(S_os_iw['up_0'], ax, color='C0', label=r'one_orbital', marker='-',    orbital=0,)

sigma = S_os_iw['up_0']

plot_sigma_iw_orbitals(sigma, orbitals=[0, 1, 3], ax=ax, color='C0', label=r'$t_{{2g}}$', )

plot_sigma_iw_orbitals(sigma, orbitals=[2, 4], ax=ax, color='C1', label=r'$e_g$', )



ax[0].set_xlim(0,30)
ax[0].set_ylim(-1,2)
ax[1].set_ylim(0,5)

plt.legend()
#plt.show()
plt.savefig('self.png')
""")













    return os_config_filename, csc_config_filename, plot_self_script_filename, spectral_plot_script_filename, maxent_script_filename, slurm_script_filename, slurm_csc_script_filename



from siman.calculators.solid_dmft_out_parser import TriqsOutput



from math import sqrt


def relative_dmft_phase_energies(
    cl_r3m,
    cl_jt,
    cl_bd,
    last_n=10, out = ''
):
    """
    Calculate JT and BD energies relative to R-3m
    in eV per LiNiO2 formula unit.

    The mean energies and standard deviations are obtained using
    ``cl.dmft_conv.energy_statistics(last_n=last_n)``.

    Standard deviations are propagated assuming independent energy
    estimates.
    """

    atoms_per_fu = 4  # LiNiO2

    calculations = {
        "R-3m": cl_r3m,
        "JT": cl_jt,
        "BD": cl_bd,
    }

    energies = {}

    for phase, cl in calculations.items():
        if 'csc' in out:
            stats = cl.csc_dmft_conv.energy_statistics(last_n=last_n)
        else:
            stats = cl.dmft_conv.energy_statistics(last_n=last_n)

        total_energy = stats["Total energy"]

        n_fu = cl.end.natom / atoms_per_fu

        energies[phase] = {
            "mean": total_energy["mean"] / n_fu,
            "std": total_energy["std"] / n_fu,
            "n_fu": n_fu,
        }

    reference = energies["R-3m"]
    result = {}

    for phase in ("JT", "BD"):
        relative_mean = (
            energies[phase]["mean"]
            - reference["mean"]
        )

        relative_std = sqrt(
            energies[phase]["std"] ** 2
            + reference["std"] ** 2
        )

        result[phase] = {
            "mean": relative_mean,
            "std": relative_std,
        }
    fmt = ".2f"
    energies = result 
    print(
        out+" DMFT energies: "
        f"JT: {energies['JT']['mean']:{fmt}} ± "
        f"{energies['JT']['std']:{fmt}} eV/f.u.; "
        f"BD: {energies['BD']['mean']:{fmt}} ± "
        f"{energies['BD']['std']:{fmt}} eV/f.u."
    )



    return result