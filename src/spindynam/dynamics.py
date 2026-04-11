import numpy as np
from ase.md.md import MolecularDynamics
from ase.units import fs


class SpinLatticeDynamics(MolecularDynamics):
    """Coupled Spin-Lattice Dynamics integrator.

    Following Tranchida et al., J. Comput. Phys. 372 (2018) 406-425.

    The integrator solves the Equations of Motion (Equations 6, 7, and 8)
    derived from the spin-lattice Hamiltonian (Equation 4) using the
    specified Poisson bracket (Equation 5).

    The algorithm uses a second-order symmetric Suzuki-Trotter decomposition
    (Equation 13):
    exp(L_p*dt/2) exp(L_s*dt/2) exp(L_r*dt) exp(L_s*dt/2) exp(L_p*dt/2)
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

        # Check if spins are available
        if "spins" not in self.atoms.arrays:
            # Initialize with magmoms if available, otherwise zeros
            magmoms = self.atoms.get_initial_magnetic_moments()
            if magmoms.ndim == 1:
                # Assuming z-aligned spins if only scalars provided
                spins = np.zeros((len(self.atoms), 3))
                spins[:, 2] = magmoms
            else:
                spins = magmoms.copy()
            self.atoms.set_array("spins", spins)

    def step(self, forces=None):
        atoms = self.atoms
        dt = self.dt

        # 1. Apply exp(L_p * dt/2): update momenta
        if forces is None:
            forces = atoms.get_forces(md=True)
        p = atoms.get_momenta()
        p += 0.5 * dt * forces
        atoms.set_momenta(p)

        # 2. Apply exp(L_s * dt/2): update spins
        # We need effective fields omega_i
        # Calculator should provide 'magnetic_forces' or 'effective_fields'
        self.update_spins(0.5 * dt)

        # 3. Apply exp(L_r * dt): update positions
        r = atoms.get_positions()
        m = atoms.get_masses()[:, np.newaxis]
        # p is already updated to t + dt/2
        r_new = r + dt * p / m
        atoms.set_positions(r_new)

        # 4. Apply exp(L_s * dt/2): update spins (with new positions)
        self.update_spins(0.5 * dt)

        # 5. Apply exp(L_p * dt/2): update momenta
        # Need new forces at new positions and spins
        forces = atoms.get_forces(md=True)
        p = atoms.get_momenta()
        p += 0.5 * dt * forces
        atoms.set_momenta(p)

        return forces

    def update_spins(self, dt):
        """Update spins using a symmetric sequential update (Suzuki-Trotter)."""
        if not hasattr(self.atoms.calc, "get_magnetic_forces"):
            raise RuntimeError(
                "Calculator must implement get_magnetic_forces for spin dynamics"
            )

        natoms = len(self.atoms)
        dt_fs = dt / fs

        # We perform a symmetric sequence of updates to preserve
        # symplecticity/reversibility.
        # i.e., exp(L_s1*dt/2) ... exp(L_sN*dt) ... exp(L_s1*dt/2)
        # For simplicity, we'll do:
        # 1. Forward pass: rotate each spin by dt/2
        # 2. Backward pass: rotate each spin by dt/2
        # We must recalculate effective fields when a neighbor changes.

        def single_spin_update(idx, time_step):
            # Recalculate omega for this atom (neighbor dependency)
            # For this simple implementation, we recalculate all
            self.atoms.calc.results.clear()
            omega = self.atoms.calc.get_magnetic_forces(self.atoms)
            spins = self.atoms.get_array("spins")

            w = omega[idx]
            s = spins[idx]
            w_norm_sq = np.dot(w, w)
            if w_norm_sq > 1e-30:
                # Equation 15: Rational approximation for single-spin propagation
                # s(t+dt) = s(t) + [ dt*(w x s) + (dt**2/2)*w_x_w_x_s ] / (1 + ...)

                # Precompute terms
                w_cross_s = np.cross(w, s)
                w_cross_w_cross_s = np.cross(w, w_cross_s)

                denom = 1.0 + 0.25 * time_step**2 * w_norm_sq

                dt2_05 = 0.5 * time_step**2
                delta_s = (time_step * w_cross_s + dt2_05 * w_cross_w_cross_s) / denom

                s_new = s + delta_s

                # Numerical stability: although Eq 15 is norm-preserving,
                # we keep the spin normalized to machine precision
                # (The paper says it eliminates the need for rescaling,
                # but for long trajectories it's safer)
                spins[idx] = s_new / np.linalg.norm(s_new)
                self.atoms.set_array("spins", spins)

        # Forward pass (1 to N)
        for i in range(natoms):
            single_spin_update(i, 0.5 * dt_fs)

        # Backward pass (N to 1)
        for i in range(natoms - 1, -1, -1):
            single_spin_update(i, 0.5 * dt_fs)

        if self.atoms.calc is not None:
            self.atoms.calc.results.clear()
