import numpy as np
from ase.calculators.calculator import Calculator, all_changes

# hbar in eV * fs
# hbar_eV_s = 6.582119569e-16
# hbar_eV_fs = 0.6582119569
HBAR = 0.6582119569


class SpinLatticeHeisenberg(Calculator):
    """Coupled spin-lattice Heisenberg calculator.

    Implements the total spin-lattice Hamiltonian (Equation 4):
    H_sl = H_mag + sum_i p_i^2 / (2m_i) + sum_{i<j} V(r_ij)

    Where the magnetic Hamiltonian H_mag is (Equation 3):
    H_mag = -mu_B * mu_0 * sum_i g_i s_i . H_ext - sum_{i<j} J(r_ij) s_i . s_j

    And the mechanical potential V(r) is a harmonic well:
    V(r) = 0.5 * K * (r - r0)^2

    The exchange J(r) follows an exponential decay:
    J(r) = J0 * exp(-alpha * (r/r0 - 1))
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

    def calculate(self, atoms=None, properties=["energy"], system_changes=all_changes):
        Calculator.calculate(self, atoms, properties, system_changes)

        r0 = self.parameters.r0
        K = self.parameters.K
        J0 = self.parameters.J0
        alpha = self.parameters.alpha
        hbar = self.parameters.hbar
        H_ext = np.array(self.parameters.H_ext)
        gi_param = self.parameters.g_i

        # Bohr magneton in eV/T (CODATA 2018)
        mu_B = 5.7883818012e-5

        positions = self.atoms.get_positions()
        spins = self.atoms.get_array("spins")
        natoms = len(self.atoms)

        # Handle per-atom Lande factors (g_i in Eq 3)
        if np.isscalar(gi_param):
            gi = np.full(natoms, gi_param)
        else:
            gi = np.array(gi_param)
            if len(gi) != natoms:
                raise ValueError(f"g_i array length {len(gi)} != natoms {natoms}")

        energy = 0.0  # Mechanical/Lattice potential energy (sum V)
        mag_energy = 0.0  # Magnetic Hamiltonian H_mag (Eq 3)
        forces = np.zeros((natoms, 3))
        mag_forces = np.zeros((natoms, 3))  # Effective fields omega_i (Eq 8)

        # Zeeman contribution (First term of Equation 3)
        # H_zeeman = -mu_B * mu_0 * sum_i gi * si . H_ext
        for i in range(natoms):
            e_zeeman = -mu_B * gi[i] * np.dot(spins[i], H_ext)
            mag_energy += e_zeeman
            # omega_i_zeeman = gi * mu_B / hbar * H_ext
            mag_forces[i] += (gi[i] * mu_B / hbar) * H_ext

        # Exchange and Lattice terms
        for i in range(natoms):
            for j in range(i + 1, natoms):
                diff = positions[j] - positions[i]
                r = np.linalg.norm(diff)
                if r < 1e-10:
                    continue
                unit_vec = diff / r

                # Lattice mechanical part (Last term of Equation 4)
                v_ij = 0.5 * K * (r - r0) ** 2
                f_mag = -K * (r - r0)  # -dV/dr
                forces[i] += f_mag * (-unit_vec)
                forces[j] -= f_mag * (-unit_vec)
                energy += v_ij

                # Exchange part (Second term of Equation 3)
                # J(r)
                j_r = J0 * np.exp(-alpha * (r / r0 - 1.0))
                dot_si_sj = np.dot(spins[i], spins[j])
                e_exch = -j_r * dot_si_sj
                mag_energy += e_exch

                # Magnetic effective fields (Equation 8)
                # omega_i = 1/hbar * sum_j J(r_ij) s_j
                mag_forces[i] += (j_r / hbar) * spins[j]
                mag_forces[j] += (j_r / hbar) * spins[i]

                # Spin-lattice force on atoms (Equation 7)
                # F_i = sum_j [ -dV/dr + dJ/dr (si.sj) ] e_ij
                dj_dr = -(alpha / r0) * j_r
                f_sl_mag = dj_dr * dot_si_sj
                forces[i] += f_sl_mag * (-unit_vec)
                forces[j] -= f_sl_mag * (-unit_vec)

        self.results["mag_energy"] = mag_energy
        self.results["energy"] = energy + mag_energy
        self.results["forces"] = forces
        self.results["magnetic_forces"] = mag_forces

    def get_magnetic_forces(self, atoms=None):
        return self.get_property("magnetic_forces", atoms)
