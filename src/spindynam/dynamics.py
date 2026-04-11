import numpy as np
from ase.md.md import MolecularDynamics
from ase.units import fs


class SpinLatticeDynamics(MolecularDynamics):
    """Coupled Spin-Lattice Dynamics integrator.

    Following Tranchida et al., J. Comput. Phys. 372 (2018) 406-425.
    """

    def __init__(
        self,
        atoms,
        timestep,
        trajectory=None,
        logfile=None,
        loginterval=1,
        append_trajectory=False,
        **kwargs,
    ):
        MolecularDynamics.__init__(
            self,
            atoms,
            timestep,
            trajectory,
            logfile,
            loginterval,
            append_trajectory=append_trajectory,
            **kwargs,
        )
        self._ensure_spins()

    def _ensure_spins(self):
        if "spins" in self.atoms.arrays:
            return
        magmoms = self.atoms.get_initial_magnetic_moments()
        spins = np.zeros((len(self.atoms), 3))
        if magmoms.ndim == 1:
            spins[:, 2] = magmoms
        else:
            spins = magmoms.copy()
        self.atoms.set_array("spins", spins)

    def step(self, forces=None):
        dt = self.dt
        # 1. Update momenta (t+dt/2)
        forces = forces if forces is not None else self.atoms.get_forces(md=True)
        p = self.atoms.get_momenta() + 0.5 * dt * forces
        self.atoms.set_momenta(p)

        # 2. Update spins (t+dt/2)
        self.update_spins(0.5 * dt)

        # 3. Update positions (t+dt)
        r = self.atoms.get_positions()
        m = self.atoms.get_masses()[:, np.newaxis]
        self.atoms.set_positions(r + dt * p / m)

        # 4. Update spins (t+dt)
        self.update_spins(0.5 * dt)

        # 5. Update momenta (t+dt)
        forces = self.atoms.get_forces(md=True)
        self.atoms.set_momenta(self.atoms.get_momenta() + 0.5 * dt * forces)
        return forces

    def update_spins(self, dt):
        """Update spins using a symmetric sequential update (Suzuki-Trotter)."""
        if not hasattr(self.atoms.calc, "get_magnetic_forces"):
            raise RuntimeError("Calculator must implement get_magnetic_forces")

        natoms = len(self.atoms)
        half_dt_fs = 0.5 * dt / fs

        # Forward pass (1 to N)
        for i in range(natoms):
            self._rotate_single_spin(i, half_dt_fs)

        # Backward pass (N to 1)
        for i in range(natoms - 1, -1, -1):
            self._rotate_single_spin(i, half_dt_fs)

        self.atoms.calc.results.clear()

    def _rotate_single_spin(self, idx, dt_fs):
        self.atoms.calc.results.clear()
        omega = self.atoms.calc.get_magnetic_forces(self.atoms)
        spins = self.atoms.get_array("spins")
        w, s = omega[idx], spins[idx]
        w_norm_sq = np.dot(w, w)

        if w_norm_sq > 1e-30:
            w_cross_s = np.cross(w, s)
            w_c_w_c_s = np.cross(w, w_cross_s)
            denom = 1.0 + 0.25 * dt_fs**2 * w_norm_sq
            delta_s = (dt_fs * w_cross_s + 0.5 * dt_fs**2 * w_c_w_c_s) / denom
            s_new = s + delta_s
            spins[idx] = s_new / np.linalg.norm(s_new)
            self.atoms.set_array("spins", spins)
