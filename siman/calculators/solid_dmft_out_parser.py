from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import numpy as np


Number = Union[int, float]


def _float(value: str) -> float:
    """Convert Fortran/Python scientific notation to float."""
    return float(value.replace("D", "E").replace("d", "e"))


def _numbers(text: str) -> List[float]:
    """Extract all floating-point numbers from a string."""
    pattern = r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?"
    return [_float(x) for x in re.findall(pattern, text)]


def _parse_matrix(text: str) -> np.ndarray:
    """
    Parse a matrix printed as:

    [[1.0 0.0]
     [0.0 1.0]]
    """
    rows = re.findall(r"\[([^\[\]]+)\]", text)

    matrix = []
    for row in rows:
        values = _numbers(row)
        if values:
            matrix.append(values)

    return np.array(matrix, dtype=float)


@dataclass
class ImpurityData:
    index: int

    dft_occupancy: Optional[float] = None
    occupancy: Optional[float] = None
    g_beta_half: Optional[float] = None

    dc_potential: Optional[float] = None
    dc_energy: Optional[float] = None
    dc_method: Optional[str] = None
    dc_u_mean: Optional[float] = None

    average_sign: Optional[float] = None
    perturbation_order: Optional[float] = None
    number_of_measures: Optional[int] = None
    cycles_per_second: Optional[float] = None
    solver_time: Optional[float] = None

    density_matrix: Dict[str, np.ndarray] = field(default_factory=dict)
    density_eigenvalues: Dict[str, np.ndarray] = field(default_factory=dict)

    Z: Dict[str, List[float]] = field(default_factory=dict)

    convergence: Dict[str, Any] = field(default_factory=dict)


@dataclass
class TriqsOutput:
    path: Path
    text: str

    version_solid_dmft: Optional[str] = None
    version_triqs: Optional[str] = None
    solver: Optional[str] = None

    U: Optional[float] = None
    J: Optional[float] = None
    beta: Optional[float] = None

    iteration: Optional[int] = None
    requested_iterations: Optional[int] = None

    chemical_potential: Optional[float] = None
    total_density: Optional[float] = None
    total_correlated_charge: Optional[float] = None

    energy: Dict[str, Optional[float]] = field(default_factory=dict)
    convergence: Dict[str, Optional[float]] = field(default_factory=dict)

    overall_elapsed_time: Optional[float] = None

    impurities: Dict[int, ImpurityData] = field(default_factory=dict)

    @classmethod
    def from_file(
        cls,
        filename: Union[str, Path],
    ) -> "TriqsOutput":
        path = Path(filename)
        text = path.read_text(encoding="utf-8", errors="replace")

        obj = cls(path=path, text=text)
        obj._parse()
        return obj

    def _get_impurity(self, index: int) -> ImpurityData:
        if index not in self.impurities:
            self.impurities[index] = ImpurityData(index=index)
        return self.impurities[index]

    def _last(
        self,
        pattern: str,
        group: int = 1,
        flags: int = 0,
    ) -> Optional[str]:
        matches = list(re.finditer(pattern, self.text, flags))
        if not matches:
            return None
        return matches[-1].group(group)

    def _last_float(
        self,
        pattern: str,
        group: int = 1,
        flags: int = 0,
    ) -> Optional[float]:
        value = self._last(pattern, group=group, flags=flags)
        return _float(value) if value is not None else None

    def _last_int(
        self,
        pattern: str,
        group: int = 1,
        flags: int = 0,
    ) -> Optional[int]:
        value = self._last(pattern, group=group, flags=flags)
        return int(value) if value is not None else None

    def _parse(self) -> None:
        self._parse_general()
        self._parse_dft_occupancies()
        self._parse_dc()
        self._parse_solver_blocks()
        self._parse_impurity_summary()
        self._parse_convergence()
        self._parse_energies()

    def _parse_general(self) -> None:
        self.version_solid_dmft = self._last(
            r"solid_dmft version:\s*([^\s]+)"
        )
        self.version_triqs = self._last(
            r"TRIQS version:\s*([^\s]+)"
        )

        self.solver = self._last(
            r"Using\s+([A-Za-z0-9_-]+)\s+solver"
        )

        self.U = self._last_float(
            r"Assuming U\s*=\s*([-+0-9.EeDd]+)"
        )
        self.J = self._last_float(
            r"Assuming J\s*=\s*([-+0-9.EeDd]+)"
        )
        self.beta = self._last_float(
            r"(?:general\.)?beta\s*[=:]\s*([-+0-9.EeDd]+)"
        )

        iteration_match = self._last(
            r"Running iteration:\s*(\d+)\s*/\s*(\d+)",
            group=0,
        )
        if iteration_match:
            match = re.search(
                r"Running iteration:\s*(\d+)\s*/\s*(\d+)",
                iteration_match,
            )
            if match:
                self.iteration = int(match.group(1))
                self.requested_iterations = int(match.group(2))

        self.chemical_potential = self._last_float(
            r"Total Density\s*=\s*[-+0-9.EeDd]+"
            r"\s*;\s*Chemical Potential\s*=\s*([-+0-9.EeDd]+)"
        )

        self.total_density = self._last_float(
            r"Total Density\s*=\s*([-+0-9.EeDd]+)"
            r"\s*;\s*Chemical Potential"
        )

        self.total_correlated_charge = self._last_float(
            r"Total charge of all correlated shells\s*:\s*"
            r"([-+0-9.EeDd]+)"
        )

        self.overall_elapsed_time = self._last_float(
            r"overall elapsed time:\s*([-+0-9.EeDd]+)\s*seconds"
        )

    def _parse_dft_occupancies(self) -> None:
        pattern = (
            r"Total density for imp\s+(\d+)\s+from DFT:\s*"
            r"([-+0-9.EeDd]+)"
        )

        for match in re.finditer(pattern, self.text):
            index = int(match.group(1))
            impurity = self._get_impurity(index)
            impurity.dft_occupancy = _float(match.group(2))

    def _parse_dc(self) -> None:
        """
        Parse blocks such as:

        Calculating standard DC for impurity 0 ...
        Detected use_dc_formula=1, changing to cHeld
        DC potential computed using the cHeld method, V_DC = ...
        E_DC using the cHeld method, E_DC = ...
        Held method for 5 orbitals, computed U_mean=...
        """
        starts = list(
            re.finditer(
                r"Calculating standard DC for impurity\s+(\d+)",
                self.text,
            )
        )

        for position, start in enumerate(starts):
            index = int(start.group(1))

            end = (
                starts[position + 1].start()
                if position + 1 < len(starts)
                else len(self.text)
            )
            block = self.text[start.start():end]
            impurity = self._get_impurity(index)

            method_match = re.search(
                r"DC potential computed using the\s+(\S+)\s+method",
                block,
            )
            if method_match:
                impurity.dc_method = method_match.group(1)

            potential_matches = re.findall(
                r"V_DC\s*=\s*([-+0-9.EeDd]+)",
                block,
            )
            if potential_matches:
                impurity.dc_potential = _float(potential_matches[-1])

            energy_matches = re.findall(
                r"E_DC\s*=\s*([-+0-9.EeDd]+)",
                block,
            )
            if energy_matches:
                impurity.dc_energy = _float(energy_matches[-1])

            u_mean_matches = re.findall(
                r"computed U_mean\s*=\s*([-+0-9.EeDd]+)",
                block,
            )
            if u_mean_matches:
                impurity.dc_u_mean = _float(u_mean_matches[-1])

    def _parse_solver_blocks(self) -> None:
        """
        Parse output produced for each successive impurity solver.

        The impurity number is inferred from:
            Solving the impurity problem for shell N ...

        If that heading is absent, blocks are assigned sequentially.
        """
        solver_starts = list(
            re.finditer(
                r"Solving the impurity problem for shell\s+(\d+)",
                self.text,
            )
        )

        if solver_starts:
            for position, start in enumerate(solver_starts):
                index = int(start.group(1))

                end = (
                    solver_starts[position + 1].start()
                    if position + 1 < len(solver_starts)
                    else len(self.text)
                )

                block = self.text[start.start():end]
                self._parse_single_solver_block(index, block)
        else:
            # Fallback for outputs without explicit shell headings.
            blocks = re.split(
                r"(?=Total number of measures:)",
                self.text,
            )

            sequential_index = 0
            for block in blocks:
                if "Total number of measures:" not in block:
                    continue

                self._parse_single_solver_block(
                    sequential_index,
                    block,
                )
                sequential_index += 1

    def _parse_single_solver_block(
        self,
        index: int,
        block: str,
    ) -> None:
        impurity = self._get_impurity(index)

        value = self._block_last_float(
            block,
            r"Total number of measures:\s*(\d+)",
        )
        if value is not None:
            impurity.number_of_measures = int(value)

        impurity.cycles_per_second = self._block_last_float(
            block,
            r"Total cycles \(measures\) / second:\s*"
            r"([-+0-9.EeDd]+)",
        )

        impurity.average_sign = self._block_last_float(
            block,
            r"Average sign:\s*([-+0-9.EeDd]+)",
        )

        impurity.perturbation_order = self._block_last_float(
            block,
            r"Average perturbation order in Delta:\s*"
            r"([-+0-9.EeDd]+)",
        )

        impurity.solver_time = self._block_last_float(
            block,
            r"Actual time for solver:\s*"
            r"([-+0-9.EeDd]+)\s*s",
        )

        qmc_occupancy = self._block_last_float(
            block,
            r"Total charge of impurity problem:\s*"
            r"([-+0-9.EeDd]+)",
        )
        if qmc_occupancy is not None:
            impurity.occupancy = qmc_occupancy

        convergence = self._block_last_float(
            block,
            r"Total charge convergency of impurity problem:\s*"
            r"([-+0-9.EeDd]+)",
        )
        if convergence is not None:
            impurity.convergence["solver_charge"] = convergence

        self._parse_density_matrix(block, impurity)

    @staticmethod
    def _block_last_float(
        block: str,
        pattern: str,
    ) -> Optional[float]:
        matches = re.findall(pattern, block, flags=re.MULTILINE)
        if not matches:
            return None
        return _float(matches[-1])

    def _parse_density_matrix(
        self,
        block: str,
        impurity: ImpurityData,
    ) -> None:
        """
        Parse the measured QMC density matrix.

        Only a section explicitly preceded by:
            Density matrix:
        is selected. This avoids confusing it with:
            Estimated density matrix:
            DFT density matrix
        """
        density_sections = list(
            re.finditer(
                r"(?<!Estimated )Density matrix:\s*",
                block,
            )
        )

        if not density_sections:
            return

        start = density_sections[-1].end()
        density_text = block[start:]

        stop_patterns = [
            r"\n\s*\*\*\*",
            r"\n\s*Post-processing",
            r"\n\s*Total charge of all correlated shells",
            r"\n\s*sumk calc_mu",
            r"\n\s*Solving the impurity problem",
        ]

        stop = len(density_text)
        for pattern in stop_patterns:
            match = re.search(pattern, density_text)
            if match:
                stop = min(stop, match.start())

        density_text = density_text[:stop]

        spin_pattern = re.compile(
            r"(up_0|down_0),\s*real part\s*"
            r"(\[\[.*?\]\])"
            r"\s*eigenvalues:\s*\[([^\]]+)\]",
            flags=re.DOTALL,
        )

        for match in spin_pattern.finditer(density_text):
            spin = match.group(1).replace("_0", "")
            matrix_text = match.group(2)
            eigenvalues_text = match.group(3)

            impurity.density_matrix[spin] = _parse_matrix(matrix_text)
            impurity.density_eigenvalues[spin] = np.array(
                _numbers(eigenvalues_text),
                dtype=float,
            )

    def _parse_impurity_summary(self) -> None:
        occupancy_pattern = (
            r"total occupan(?:c)?y of impurity\s+(\d+):\s*"
            r"([-+0-9.EeDd]+)"
        )
        for match in re.finditer(
            occupancy_pattern,
            self.text,
            flags=re.IGNORECASE,
        ):
            index = int(match.group(1))
            self._get_impurity(index).occupancy = _float(match.group(2))

        g_pattern = (
            r"G\(beta/2\) occ of impurity\s+(\d+):\s*"
            r"([-+0-9.EeDd]+)"
        )
        for match in re.finditer(g_pattern, self.text):
            index = int(match.group(1))
            self._get_impurity(index).g_beta_half = _float(match.group(2))

        z_pattern = re.compile(
            r"Z \(simple estimate\) of impurity\s+(\d+)\s+per orb:\s*"
            r"\n\s*up:\s*([^\n]+)"
            r"\n\s*down:\s*([^\n]+)",
            flags=re.MULTILINE,
        )

        for match in z_pattern.finditer(self.text):
            index = int(match.group(1))
            impurity = self._get_impurity(index)

            impurity.Z["up"] = _numbers(match.group(2))
            impurity.Z["down"] = _numbers(match.group(3))

    def _parse_convergence(self) -> None:
        self.convergence["mu"] = self._last_float(
            r"δμ:\s*([-+0-9.EeDd]+)"
        )
        self.convergence["total_energy"] = self._last_float(
            r"δE_tot:\s*([-+0-9.EeDd]+)"
        )

        pattern = re.compile(
            r"Impurity\s+(\d+):\s*"
            r"\nδn imp\s*:\s*([-+0-9.EeDd]+)"
            r"\s*\nδn orb\s*:\s*([^\n]+)"
            r"\s*\nδ Gimp\s*:\s*([-+0-9.EeDd]+)"
            r"\s*\nδ G0\s*:\s*([-+0-9.EeDd]+)"
            r"\s*\nδ Σ\s*:\s*([-+0-9.EeDd]+)",
            flags=re.MULTILINE,
        )

        for match in pattern.finditer(self.text):
            index = int(match.group(1))
            impurity = self._get_impurity(index)

            impurity.convergence.update(
                {
                    "n_imp": _float(match.group(2)),
                    "n_orb": _numbers(match.group(3)),
                    "Gimp": _float(match.group(4)),
                    "G0": _float(match.group(5)),
                    "Sigma": _float(match.group(6)),
                }
            )

    def _parse_energies(self) -> None:
        self.energy = {
            "total": self._last_float(
                r"total energy:\s*([-+0-9.EeDd]+)"
            ),
            "dft": self._last_float(
                r"DFT energy:\s*([-+0-9.EeDd]+)"
            ),
            "correlation": self._last_float(
                r"correlation energy:\s*([-+0-9.EeDd]+)"
            ),
            "dft_band_correction": self._last_float(
                r"DFT band correction:\s*([-+0-9.EeDd]+)"
            ),
            "kinetic_dft": self._last_float(
                r"Kinetic energy contribution dft part:\s*"
                r"([-+0-9.EeDd]+)"
            ),
            "kinetic_dmft": self._last_float(
                r"Kinetic energy contribution dmft part:\s*"
                r"([-+0-9.EeDd]+)"
            ),
        }

    @property
    def results(self) -> Dict[str, Any]:
        return {
            "version_solid_dmft": self.version_solid_dmft,
            "version_triqs": self.version_triqs,
            "solver": self.solver,
            "U": self.U,
            "J": self.J,
            "beta": self.beta,
            "iteration": self.iteration,
            "requested_iterations": self.requested_iterations,
            "chemical_potential": self.chemical_potential,
            "total_density": self.total_density,
            "total_correlated_charge": self.total_correlated_charge,
            "energy": self.energy,
            "convergence": self.convergence,
            "overall_elapsed_time": self.overall_elapsed_time,
            "impurities": {
                index: {
                    "dft_occupancy": impurity.dft_occupancy,
                    "occupancy": impurity.occupancy,
                    "g_beta_half": impurity.g_beta_half,
                    "dc_potential": impurity.dc_potential,
                    "dc_energy": impurity.dc_energy,
                    "dc_method": impurity.dc_method,
                    "dc_u_mean": impurity.dc_u_mean,
                    "average_sign": impurity.average_sign,
                    "perturbation_order": impurity.perturbation_order,
                    "number_of_measures": impurity.number_of_measures,
                    "cycles_per_second": impurity.cycles_per_second,
                    "solver_time": impurity.solver_time,
                    "density_matrix": impurity.density_matrix,
                    "density_eigenvalues": impurity.density_eigenvalues,
                    "Z": impurity.Z,
                    "convergence": impurity.convergence,
                }
                for index, impurity in sorted(self.impurities.items())
            },
        }

    def print_summary(self) -> None:
        fmt = ".2f"
        exp = ".1e"

        print("=" * 72)
        print("solid_dmft summary")
        print("=" * 72)

        if self.version_solid_dmft:
            print(f"solid_dmft version : {self.version_solid_dmft}")
        if self.version_triqs:
            print(f"TRIQS version      : {self.version_triqs}")
        if self.solver:
            print(f"Solver             : {self.solver}")

        if self.iteration is not None:
            if self.requested_iterations is not None:
                print(
                    f"Iteration          : "
                    f"{self.iteration} / {self.requested_iterations}"
                )
            else:
                print(f"Iteration          : {self.iteration}")

        if self.U is not None:
            print(f"U                  : {self.U:{fmt}} eV")
        if self.J is not None:
            print(f"J                  : {self.J:{fmt}} eV")
        if self.beta is not None:
            print(f"beta               : {self.beta:{fmt}} eV^-1")

        if self.chemical_potential is not None:
            print(
                f"Chemical potential : "
                f"{self.chemical_potential:{fmt}} eV"
            )

        if self.total_density is not None:
            print(f"Total density      : {self.total_density:{fmt}}")

        if self.total_correlated_charge is not None:
            print(
                f"Correlated charge  : "
                f"{self.total_correlated_charge:{fmt}}"
            )

        print()
        print("Impurity observables")
        print("-" * 72)

        for index, impurity in sorted(self.impurities.items()):
            print(f"Impurity {index}")

            if impurity.dft_occupancy is not None:
                print(
                    f"  DFT occupancy       : "
                    f"{impurity.dft_occupancy:{fmt}}"
                )

            if impurity.occupancy is not None:
                print(
                    f"  DMFT occupancy      : "
                    f"{impurity.occupancy:{fmt}}"
                )

            if impurity.g_beta_half is not None:
                print(
                    f"  G(beta/2)           : "
                    f"{impurity.g_beta_half:{fmt}}"
                )

            if impurity.dc_method is not None:
                print(f"  DC method           : {impurity.dc_method}")

            if impurity.dc_potential is not None:
                print(
                    f"  V_DC                : "
                    f"{impurity.dc_potential:{fmt}} eV"
                )

            if impurity.dc_energy is not None:
                print(
                    f"  E_DC                : "
                    f"{impurity.dc_energy:{fmt}} eV"
                )

            if impurity.dc_u_mean is not None:
                print(
                    f"  DC U_mean           : "
                    f"{impurity.dc_u_mean:{fmt}} eV"
                )

            if impurity.average_sign is not None:
                print(
                    f"  Average sign        : "
                    f"{impurity.average_sign:{fmt}}"
                )

            if impurity.perturbation_order is not None:
                print(
                    f"  Perturbation order  : "
                    f"{impurity.perturbation_order:{fmt}}"
                )

            if impurity.number_of_measures is not None:
                print(
                    f"  Measures            : "
                    f"{impurity.number_of_measures}"
                )

            if impurity.cycles_per_second is not None:
                print(
                    f"  Measures/s          : "
                    f"{impurity.cycles_per_second:{fmt}}"
                )

            if impurity.solver_time is not None:
                print(
                    f"  Solver time         : "
                    f"{impurity.solver_time:{fmt}} s"
                )

            for spin in ("up", "down"):
                if spin in impurity.Z:
                    values = " ".join(
                        f"{value:{fmt}}"
                        for value in impurity.Z[spin]
                    )
                    print(f"  Z {spin:4s}             : {values}")

            if impurity.convergence:
                convergence = impurity.convergence

                if "solver_charge" in convergence:
                    print(
                        f"  Solver charge delta : "
                        f"{convergence['solver_charge']:{exp}}"
                    )

                if "n_imp" in convergence:
                    print(
                        f"  delta n_imp         : "
                        f"{convergence['n_imp']:{exp}}"
                    )

                if "n_orb" in convergence:
                    values = " ".join(
                        f"{value:{exp}}"
                        for value in convergence["n_orb"]
                    )
                    print(f"  delta n_orb         : {values}")

                if "Gimp" in convergence:
                    print(
                        f"  delta Gimp          : "
                        f"{convergence['Gimp']:{exp}}"
                    )

                if "G0" in convergence:
                    print(
                        f"  delta G0            : "
                        f"{convergence['G0']:{exp}}"
                    )

                if "Sigma" in convergence:
                    print(
                        f"  delta Sigma         : "
                        f"{convergence['Sigma']:{exp}}"
                    )

            if impurity.density_matrix:
                print("  Density matrices:")

                for spin in ("up", "down"):
                    matrix = impurity.density_matrix.get(spin)

                    if matrix is not None:
                        print(f"    {spin}:")

                        matrix_text = np.array2string(
                            matrix,
                            formatter={
                                "float_kind": lambda value: f"{value:{fmt}}"
                            },
                            suppress_small=True,
                        )

                        for line in matrix_text.splitlines():
                            print("      " + line)

            if impurity.density_eigenvalues:
                print("  Density-matrix eigenvalues:")

                for spin in ("up", "down"):
                    eigenvalues = impurity.density_eigenvalues.get(spin)

                    if eigenvalues is not None:
                        values = " ".join(
                            f"{value:{fmt}}"
                            for value in eigenvalues
                        )
                        print(f"    {spin:4s}: {values}")

            print()

        print("Energetics")
        print("-" * 72)

        labels = {
            "total": "Total energy",
            "dft": "DFT energy",
            "correlation": "Correlation energy",
            "dft_band_correction": "DFT band correction",
            "kinetic_dft": "Kinetic energy, DFT part",
            "kinetic_dmft": "Kinetic energy, DMFT part",
        }

        for key, label in labels.items():
            value = self.energy.get(key)

            if value is not None:
                print(f"{label:<28s}: {value:{fmt}} eV")

        print()
        print("Global convergence")
        print("-" * 72)

        if self.convergence.get("mu") is not None:
            print(
                f"delta mu            : "
                f"{self.convergence['mu']:{exp}}"
            )

        if self.convergence.get("total_energy") is not None:
            print(
                f"delta E_tot         : "
                f"{self.convergence['total_energy']:{exp}}"
            )

        if self.overall_elapsed_time is not None:
            print()
            print(
                f"Overall elapsed time: "
                f"{self.overall_elapsed_time:{fmt}} s "
                f"({self.overall_elapsed_time / 3600:{fmt}} h)"
            )

        print("=" * 72)


        return

    def print_summary_comp(self) -> None:
        fmt = ".1f"
        exp = ".1e"

        def number(value, format_spec=fmt, missing="-"):
            """Format a scalar value."""
            if value is None:
                return missing

            try:
                return f"{float(value):{format_spec}}"
            except (TypeError, ValueError):
                return str(value)

        def vector(values, format_spec=fmt, missing="-"):
            """Format a vector as comma-separated values."""
            if values is None:
                return missing

            try:
                values = np.asarray(values).ravel()
            except (TypeError, ValueError):
                return missing

            if values.size == 0:
                return missing

            return ",".join(
                f"{float(np.real(value)):{format_spec}}"
                for value in values
            )

        def orbital_quantities(impurity):
            """
            Calculate orbital occupancies and spin polarization from the
            diagonal elements of the density matrices.

            Returns
            -------
            occ_up : ndarray or None
                Orbital-resolved spin-up occupancy.
            occ_down : ndarray or None
                Orbital-resolved spin-down occupancy.
            occ_total : ndarray or None
                n_up + n_down for each orbital.
            spin_orbital : ndarray or None
                n_up - n_down for each orbital.
            total_moment : float or None
                Sum_i(n_up_i - n_down_i), in mu_B.
            """
            matrix_up = impurity.density_matrix.get("up")
            matrix_down = impurity.density_matrix.get("down")

            occ_up = None
            occ_down = None
            occ_total = None
            spin_orbital = None
            total_moment = None

            if matrix_up is not None:
                matrix_up = np.asarray(matrix_up)
                occ_up = np.real(np.diag(matrix_up))

            if matrix_down is not None:
                matrix_down = np.asarray(matrix_down)
                occ_down = np.real(np.diag(matrix_down))

            if occ_up is not None and occ_down is not None:
                if len(occ_up) == len(occ_down):
                    occ_total = occ_up + occ_down
                    spin_orbital = occ_up - occ_down
                    total_moment = float(np.sum(spin_orbital))

            elif occ_up is not None:
                occ_total = occ_up.copy()

            elif occ_down is not None:
                occ_total = occ_down.copy()

            return (
                occ_up,
                occ_down,
                occ_total,
                spin_orbital,
                total_moment,
            )

        print("=" * 150)
        print("solid_dmft compact summary")
        print("=" * 150)

        # ------------------------------------------------------------
        # General calculation information
        # ------------------------------------------------------------

        general = []

        if self.iteration is not None:
            if self.requested_iterations is not None:
                general.append(
                    f"iteration: {self.iteration}/"
                    f"{self.requested_iterations}"
                )
            else:
                general.append(f"iteration: {self.iteration}")

        if self.U is not None:
            general.append(f"U: {self.U:{fmt}} eV")

        if self.J is not None:
            general.append(f"J: {self.J:{fmt}} eV")

        if self.beta is not None:
            general.append(f"beta: {self.beta:{fmt}} eV^-1")

        if self.chemical_potential is not None:
            general.append(
                f"mu: {self.chemical_potential:{fmt}} eV"
            )

        if self.total_density is not None:
            general.append(f"N: {self.total_density:{fmt}}")

        if general:
            print(" | ".join(general))

        print("-" * 150)

        # ------------------------------------------------------------
        # Main impurity table
        # ------------------------------------------------------------

        header = (
            f"{'Imp':>3s} "
            f"{'n_DFT':>7s} "
            f"{'n_DMFT':>7s} "
            f"{'G(beta/2)':>11s} "
            f"{'V_DC':>8s} "
            f"{'sign':>6s} "
            f"{'<k>':>7s} "
            f"{'moment':>8s} "
            f"{'Z up':>23s} "
            f"{'Z down':>23s} "
            f"{'dSigma':>10s}"
        )

        print(header)
        print("-" * len(header))

        orbital_data = {}

        for index, impurity in sorted(self.impurities.items()):
            (
                occ_up,
                occ_down,
                occ_total,
                spin_orbital,
                total_moment,
            ) = orbital_quantities(impurity)

            orbital_data[index] = (
                occ_up,
                occ_down,
                occ_total,
                spin_orbital,
                total_moment,
            )

            z_up = vector(impurity.Z.get("up"))
            z_down = vector(impurity.Z.get("down"))

            delta_sigma = impurity.convergence.get("Sigma")

            row = (
                f"{index:>3d} "
                f"{number(impurity.dft_occupancy):>7s} "
                f"{number(impurity.occupancy):>7s} "
                f"{number(impurity.g_beta_half):>11s} "
                f"{number(impurity.dc_potential):>8s} "
                f"{number(impurity.average_sign):>6s} "
                f"{number(impurity.perturbation_order):>7s} "
                f"{number(total_moment):>8s} "
                f"{z_up:>23s} "
                f"{z_down:>23s} "
                f"{number(delta_sigma, exp):>10s}"
            )

            print(row)

        print("-" * len(header))
        print("moment = sum_i(n_up,i - n_down,i), in mu_B")

        # ------------------------------------------------------------
        # Orbital-resolved occupancies and spin polarization
        # ------------------------------------------------------------

        print()
        print("Orbital-resolved occupancies and spin polarization")
        print("-" * 150)

        orbital_header = (
            f"{'Imp':>3s} "
            f"{'n_up':>31s} "
            f"{'n_down':>31s} "
            f"{'n_total':>31s} "
            f"{'n_up - n_down':>31s} "
            f"{'M, mu_B':>8s}"
        )

        print(orbital_header)
        print("-" * len(orbital_header))

        for index, impurity in sorted(self.impurities.items()):
            (
                occ_up,
                occ_down,
                occ_total,
                spin_orbital,
                total_moment,
            ) = orbital_data[index]

            print(
                f"{index:>3d} "
                f"{vector(occ_up):>31s} "
                f"{vector(occ_down):>31s} "
                f"{vector(occ_total):>31s} "
                f"{vector(spin_orbital):>31s} "
                f"{number(total_moment):>8s}"
            )

        print("-" * len(orbital_header))

        # ------------------------------------------------------------
        # Energies and global convergence
        # ------------------------------------------------------------

        footer = []

        if self.energy.get("total") is not None:
            footer.append(
                f"E_tot: {self.energy['total']:{fmt}} eV"
            )

        if self.energy.get("correlation") is not None:
            footer.append(
                f"E_corr: {self.energy['correlation']:{fmt}} eV"
            )

        if self.convergence.get("mu") is not None:
            footer.append(
                f"delta mu: {self.convergence['mu']:{exp}}"
            )

        if self.convergence.get("total_energy") is not None:
            footer.append(
                f"delta E_tot: "
                f"{self.convergence['total_energy']:{exp}}"
            )

        if self.overall_elapsed_time is not None:
            footer.append(
                f"time: {self.overall_elapsed_time / 3600:{fmt}} h"
            )

        if footer:
            print(" | ".join(footer))

        print("=" * 150)


from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np


@dataclass
class DmftConvergenceHistory:
    """History of observables from a solid_dmft convergence table."""

    iteration: np.ndarray
    chemical_potential: np.ndarray

    g_beta_half: np.ndarray
    orbital_occupancy: np.ndarray

    impurity_occupancy: np.ndarray

    total_energy: np.ndarray
    dft_energy: np.ndarray
    band_correction_energy: np.ndarray
    impurity_interaction_energy: np.ndarray
    double_counting_energy: np.ndarray

    # ------------------------------------------------------------
    # Reading
    # ------------------------------------------------------------

    @classmethod
    def from_file(
        cls,
        filename: Union[str, Path],
    ) -> "DmftConvergenceHistory":
        """
        Read a solid_dmft convergence table.

        Expected row format:

        it | mu | G(beta/2) per orbital | orbital occs up+down |
        impurity occ | E_tot |
        E_DFT E_bandcorr E_int_imp E_DC
        """
        filename = Path(filename)

        if not filename.is_file():
            raise FileNotFoundError(
                f"Convergence file not found: {filename}"
            )

        rows = []

        with filename.open("r", encoding="utf-8") as file:
            for line_number, line in enumerate(file, start=1):
                line = line.strip()

                if not line:
                    continue

                # Skip the header and separator lines.
                if line.lower().startswith("it"):
                    continue

                if set(line) <= {"-", "+", "|", " "}:
                    continue

                parts = [part.strip() for part in line.split("|")]

                if len(parts) != 7:
                    raise ValueError(
                        f"Invalid convergence-table row at line "
                        f"{line_number}: expected 7 pipe-separated "
                        f"columns, got {len(parts)}.\n"
                        f"Line: {line}"
                    )

                try:
                    iteration = int(parts[0])
                    chemical_potential = float(parts[1])

                    g_beta_half = [
                        float(value)
                        for value in parts[2].split()
                    ]

                    orbital_occupancy = [
                        float(value)
                        for value in parts[3].split()
                    ]

                    impurity_occupancy = float(parts[4])
                    total_energy = float(parts[5])

                    energy_terms = [
                        float(value)
                        for value in parts[6].split()
                    ]

                except ValueError as error:
                    raise ValueError(
                        f"Could not parse line {line_number}:\n"
                        f"{line}"
                    ) from error

                if len(energy_terms) != 4:
                    raise ValueError(
                        f"Expected four energy components at line "
                        f"{line_number}, got {len(energy_terms)}:\n"
                        f"{line}"
                    )

                if len(g_beta_half) != len(orbital_occupancy):
                    raise ValueError(
                        f"Different number of orbitals in G(beta/2) "
                        f"and occupancy at line {line_number}: "
                        f"{len(g_beta_half)} and "
                        f"{len(orbital_occupancy)}."
                    )

                rows.append(
                    {
                        "iteration": iteration,
                        "chemical_potential": chemical_potential,
                        "g_beta_half": g_beta_half,
                        "orbital_occupancy": orbital_occupancy,
                        "impurity_occupancy": impurity_occupancy,
                        "total_energy": total_energy,
                        "dft_energy": energy_terms[0],
                        "band_correction_energy": energy_terms[1],
                        "impurity_interaction_energy": energy_terms[2],
                        "double_counting_energy": energy_terms[3],
                    }
                )

        if not rows:
            raise ValueError(
                f"No convergence data found in {filename}"
            )

        orbital_counts = {
            len(row["g_beta_half"])
            for row in rows
        }

        if len(orbital_counts) != 1:
            raise ValueError(
                "The number of orbitals changes between iterations: "
                f"{sorted(orbital_counts)}"
            )

        return cls(
            iteration=np.asarray(
                [row["iteration"] for row in rows],
                dtype=int,
            ),
            chemical_potential=np.asarray(
                [row["chemical_potential"] for row in rows],
                dtype=float,
            ),
            g_beta_half=np.asarray(
                [row["g_beta_half"] for row in rows],
                dtype=float,
            ),
            orbital_occupancy=np.asarray(
                [row["orbital_occupancy"] for row in rows],
                dtype=float,
            ),
            impurity_occupancy=np.asarray(
                [row["impurity_occupancy"] for row in rows],
                dtype=float,
            ),
            total_energy=np.asarray(
                [row["total_energy"] for row in rows],
                dtype=float,
            ),
            dft_energy=np.asarray(
                [row["dft_energy"] for row in rows],
                dtype=float,
            ),
            band_correction_energy=np.asarray(
                [
                    row["band_correction_energy"]
                    for row in rows
                ],
                dtype=float,
            ),
            impurity_interaction_energy=np.asarray(
                [
                    row["impurity_interaction_energy"]
                    for row in rows
                ],
                dtype=float,
            ),
            double_counting_energy=np.asarray(
                [
                    row["double_counting_energy"]
                    for row in rows
                ],
                dtype=float,
            ),
        )

    # ------------------------------------------------------------
    # Basic properties
    # ------------------------------------------------------------

    @property
    def number_of_iterations(self) -> int:
        return len(self.iteration)

    @property
    def number_of_orbitals(self) -> int:
        return self.g_beta_half.shape[1]

    # ------------------------------------------------------------
    # Observable names
    # ------------------------------------------------------------

    @property
    def observable_names(self) -> Dict[str, str]:
        """
        Internal name -> intuitive display name.
        """
        return {
            "chemical_potential": "Chemical potential",
            "g_beta_half": "G(beta/2) per orbital",
            "orbital_occupancy": "Orbital occupancy",
            "impurity_occupancy": "Total impurity occupancy",
            "total_energy": "Total energy",
            "dft_energy": "DFT energy",
            "band_correction_energy": "Band correction energy",
            "impurity_interaction_energy":
                "Impurity interaction energy",
            "double_counting_energy":
                "Double-counting energy",
        }

    def print_observables(self) -> None:
        """Print names accepted by plot()."""
        print("Available observables:")

        for internal_name, display_name in self.observable_names.items():
            if internal_name in {
                "g_beta_half",
                "orbital_occupancy",
            }:
                print(
                    f"  {display_name:<32s} "
                    f"all orbitals or orbital=0..."
                    f"{self.number_of_orbitals - 1}"
                )
            else:
                print(f"  {display_name}")

    def _resolve_observable_name(self, name: str) -> str:
        """
        Accept either an internal name or an intuitive display name.
        """
        normalized = (
            name.strip()
            .lower()
            .replace("-", " ")
            .replace("_", " ")
        )

        aliases = {
            "mu": "chemical_potential",
            "chemical potential": "chemical_potential",

            "g": "g_beta_half",
            "g beta 2": "g_beta_half",
            "g(beta/2)": "g_beta_half",
            "g beta half": "g_beta_half",
            "g(beta/2) per orbital": "g_beta_half",

            "occupancy": "orbital_occupancy",
            "orbital occupancy": "orbital_occupancy",
            "orbital occupancies": "orbital_occupancy",

            "impurity occupancy": "impurity_occupancy",
            "total impurity occupancy": "impurity_occupancy",
            "n impurity": "impurity_occupancy",

            "total energy": "total_energy",
            "e tot": "total_energy",
            "etot": "total_energy",

            "dft energy": "dft_energy",
            "e dft": "dft_energy",

            "band correction": "band_correction_energy",
            "band correction energy": "band_correction_energy",
            "bandcorr": "band_correction_energy",

            "interaction energy": "impurity_interaction_energy",
            "impurity interaction energy":
                "impurity_interaction_energy",
            "e int imp": "impurity_interaction_energy",

            "double counting": "double_counting_energy",
            "double counting energy": "double_counting_energy",
            "dc energy": "double_counting_energy",
            "e dc": "double_counting_energy",
        }

        # Internal names are also accepted.
        internal_normalized = {
            key.replace("_", " "): key
            for key in self.observable_names
        }

        if normalized in internal_normalized:
            return internal_normalized[normalized]

        if normalized in aliases:
            return aliases[normalized]

        available = ", ".join(self.observable_names.values())

        raise ValueError(
            f"Unknown observable: {name!r}.\n"
            f"Available observables: {available}"
        )

    def get_observable(
        self,
        observable: str,
        orbital: Optional[int] = None,
    ) -> Tuple[np.ndarray, str]:
        """
        Return observable values and a display label.

        For orbital-resolved quantities:
          orbital=None  -> all orbitals;
          orbital=i     -> selected orbital.
        """
        internal_name = self._resolve_observable_name(observable)
        values = getattr(self, internal_name)
        label = self.observable_names[internal_name]

        if values.ndim == 2 and orbital is not None:
            if not 0 <= orbital < self.number_of_orbitals:
                raise IndexError(
                    f"Orbital index {orbital} is outside the valid "
                    f"range 0...{self.number_of_orbitals - 1}."
                )

            values = values[:, orbital]
            label = f"{label}, orbital {orbital}"

        elif values.ndim == 1 and orbital is not None:
            raise ValueError(
                f"{label} is not orbital resolved; "
                f"the orbital argument must be None."
            )

        return values, label

    # ------------------------------------------------------------
    # Plotting
    # ------------------------------------------------------------

    def plot(
        self,
        observable: str,
        orbital: Optional[int] = None,
        orbital_labels: Optional[Sequence[str]] = None,
        last_n: Optional[int] = None,
        show: bool = True,
        ax=None,
    ):
        """
        Plot an observable versus DMFT iteration.

        Parameters
        ----------
        observable
            Intuitive or internal observable name, for example:
            "Chemical potential", "Total energy",
            "G(beta/2)", "Orbital occupancy".
        orbital
            Orbital index for an orbital-resolved observable.
            If None, all orbitals are plotted.
        orbital_labels
            Optional labels such as:
            ["d_xy", "d_yz", "d_z2", "d_xz", "d_x2-y2"].
        last_n
            Plot only the last n iterations.
        show
            Call plt.show().
        ax
            Existing matplotlib axis. A new one is created by default.
        """
        values, label = self.get_observable(
            observable,
            orbital=orbital,
        )

        iterations = self.iteration

        if last_n is not None:
            if last_n <= 0:
                raise ValueError("last_n must be positive.")

            iterations = iterations[-last_n:]
            values = values[-last_n:]

        if ax is None:
            _, ax = plt.subplots()

        if values.ndim == 1:
            ax.plot(
                iterations,
                values,
                marker="o",
                label=label,
            )

        elif values.ndim == 2:
            if orbital_labels is not None:
                if len(orbital_labels) != self.number_of_orbitals:
                    raise ValueError(
                        "orbital_labels must contain exactly "
                        f"{self.number_of_orbitals} labels."
                    )

                labels = list(orbital_labels)
            else:
                labels = [
                    f"orbital {index}"
                    for index in range(self.number_of_orbitals)
                ]

            for orbital_index, orbital_label in enumerate(labels):
                ax.plot(
                    iterations,
                    values[:, orbital_index],
                    marker="o",
                    label=orbital_label,
                )

            ax.legend()

        ax.set_xlabel("DMFT iteration")
        ax.set_ylabel(label)
        ax.set_title(f"Convergence of {label}")
        ax.grid(True)

        if values.ndim == 1:
            ax.legend()

        if show:
            plt.show()

        return ax

    # ------------------------------------------------------------
    # Energy statistics
    # ------------------------------------------------------------

    @property
    def energy_observables(self) -> Dict[str, np.ndarray]:
        return {
            "Total energy": self.total_energy,
            "DFT energy": self.dft_energy,
            "Band correction energy":
                self.band_correction_energy,
            "Impurity interaction energy":
                self.impurity_interaction_energy,
            "Double-counting energy":
                self.double_counting_energy,
        }

    def energy_statistics(
        self,
        last_n: int = 5,
        ddof: int = 1,
    ) -> Dict[str, Dict[str, float]]:
        """
        Calculate energy mean and standard deviation over the last n
        iterations.

        Parameters
        ----------
        last_n
            Number of final iterations included.
        ddof
            Delta degrees of freedom for numpy.std.
            ddof=1 gives the sample standard deviation.
            ddof=0 gives the population standard deviation.
        """
        if last_n <= 0:
            raise ValueError("last_n must be positive.")

        used_n = min(last_n, self.number_of_iterations)

        if used_n <= ddof:
            effective_ddof = 0
        else:
            effective_ddof = ddof

        statistics = {}

        for name, values in self.energy_observables.items():
            selected = values[-used_n:]

            statistics[name] = {
                "mean": float(np.mean(selected)),
                "std": float(
                    np.std(selected, ddof=effective_ddof)
                ),
                "last": float(selected[-1]),
                "minimum": float(np.min(selected)),
                "maximum": float(np.max(selected)),
                "n": used_n,
            }

        return statistics

    def print_energy_statistics(
        self,
        last_n: int = 5,
        ddof: int = 1,
        fmt: str = ".6f",
        std_fmt: str = ".3e",
    ) -> None:
        """
        Print energy value, mean and standard deviation over the last n
        iterations.
        """
        statistics = self.energy_statistics(
            last_n=last_n,
            ddof=ddof,
        )

        used_n = next(iter(statistics.values()))["n"]

        print(
            f"Energy statistics over the last {used_n} "
            f"DMFT iterations"
        )

        header = (
            f"{'Energy':<30s}"
            f"{'Last, eV':>15s}"
            f"{'Mean, eV':>15s}"
            f"{'Std, eV':>15s}"
        )

        print(header)
        print("-" * len(header))

        for name, values in statistics.items():
            print(
                f"{name:<30s}"
                f"{values['last']:>{15}{fmt}}"
                f"{values['mean']:>{15}{fmt}}"
                f"{values['std']:>{15}{std_fmt}}"
            )

