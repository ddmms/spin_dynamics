import numpy as np
from ase.calculators.calculator import Calculator, all_changes

# hbar in eV * fs
HBAR = 0.6582119569


class SpinLatticeHeisenberg(Calculator):
    """Coupled spin-lattice Heisenberg calculator.

    Implements equations from Tranchida et al. (2018).
    """

    implemented_properties = ["energy", "forces", "magnetic_forces", "mag_energy"]
    default_parameters = {
        "r0": 2.5,
        "K": 10.0,  # eV/A^2
        "J0": 0.1,  # eV
        "alpha": 2.0,  # dimensionless
        "hbar": HBAR,  # eV * fs
        "H_ext": (0.0, 0.0, 0.0),  # Tesla
        "g_i": 2.0,  # Lande factor (scalar or array)
    }

    def __init__(self, **kwargs):
        Calculator.__init__(self, **kwargs)

    def calculate(self, atoms=None, properties=None, system_changes=all_changes):
        properties = properties or ["energy"]
        Calculator.calculate(self, atoms, properties, system_changes)

        natoms = len(self.atoms)
        gi = self._get_gi_array(natoms)
        mu_B = 5.7883818012e-5  # Bohr magneton in eV/T

        # Zeeman contribution
        zeeman_e, zeeman_w = self._calculate_zeeman(natoms, gi, mu_B)

        # Pairwise (Exchange + Lattice)
        pair_e, pair_me, pair_f, pair_w = self._calculate_pairwise(natoms)

        self.results["mag_energy"] = zeeman_e + pair_me
        self.results["energy"] = self.results["mag_energy"] + pair_e
        self.results["forces"] = pair_f
        self.results["magnetic_forces"] = zeeman_w + pair_w

    def _get_gi_array(self, natoms):
        gi_param = self.parameters.g_i
        if np.isscalar(gi_param):
            return np.full(natoms, gi_param)
        gi = np.array(gi_param)
        if len(gi) != natoms:
            raise ValueError(f"g_i array length {len(gi)} != natoms {natoms}")
        return gi

    def _calculate_zeeman(self, natoms, gi, mu_B):
        H_ext = np.array(self.parameters.H_ext)
        spins = self.atoms.get_array("spins")
        mag_energy = 0.0
        mag_forces = np.zeros((natoms, 3))
        for i in range(natoms):
            e_z = -mu_B * gi[i] * np.dot(spins[i], H_ext)
            mag_energy += e_z
            mag_forces[i] += (gi[i] * mu_B / self.parameters.hbar) * H_ext
        return mag_energy, mag_forces

    def _calculate_pairwise(self, natoms):
        energy, mag_energy = 0.0, 0.0
        forces = np.zeros((natoms, 3))
        mag_forces = np.zeros((natoms, 3))
        pos = self.atoms.get_positions()
        spins = self.atoms.get_array("spins")

        for i in range(natoms):
            for j in range(i + 1, natoms):
                e, me, f, w_i, w_j = self._compute_pair(i, j, pos, spins)
                energy += e
                mag_energy += me
                forces[i] += f
                forces[j] -= f
                mag_forces[i] += w_i
                mag_forces[j] += w_j
        return energy, mag_energy, forces, mag_forces

    def _compute_pair(self, i, j, pos, spins):
        r0, K, J0, alpha = (
            self.parameters.r0,
            self.parameters.K,
            self.parameters.J0,
            self.parameters.alpha,
        )
        diff = pos[j] - pos[i]
        r = np.linalg.norm(diff)
        if r < 1e-10:
            return 0.0, 0.0, np.zeros(3), np.zeros(3), np.zeros(3)

        unit_vec = diff / r
        # Mechanical
        v_ij = 0.5 * K * (r - r0) ** 2
        f_mech = -K * (r - r0) * (-unit_vec)
        # Exchange
        j_r = J0 * np.exp(-alpha * (r / r0 - 1.0))
        dot_s = np.dot(spins[i], spins[j])
        # Force & Field (Eq 7 & 8)
        hbar = self.parameters.hbar
        f_sl = (-(alpha / r0) * j_r * dot_s) * (-unit_vec)
        f_tot = f_mech + f_sl
        w_i = (j_r / hbar) * spins[j]
        w_j = (j_r / hbar) * spins[i]
        return v_ij, -j_r * dot_s, f_tot, w_i, w_j

    def get_magnetic_forces(self, atoms=None):
        return self.get_property("magnetic_forces", atoms)
