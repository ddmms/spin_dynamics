import numpy as np
from ase import Atoms

from spindynam import SpinLatticeDynamics, SpinLatticeHeisenberg


def test_dimer_conservation():
    # Setup a dimer
    atoms = Atoms("Fe2", positions=[[0, 0, 0], [0, 0, 2.6]])

    # Initialize spins: s1: (1,0,0), s2: (0,1,0)
    # These are unit vectors. Calculator will use them directly or as magnitudes.
    # In my calculator, I assumed spins are the vectors used in dot product.
    spins = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    atoms.set_array("spins", spins)

    # Also set momenta to give some lattice dynamics
    atoms.set_velocities([[0.0, 0.0, 0.1], [0.0, 0.0, -0.1]])

    calc = SpinLatticeHeisenberg(r0=2.5, K=10.0, J0=0.5, alpha=2.0)
    atoms.calc = calc

    # Dynamics
    from ase import units

    dt = 0.1 * units.fs
    dyn = SpinLatticeDynamics(atoms, timestep=dt)

    energies = []
    spin_norms = []
    total_mags = []

    for _ in range(100):
        dyn.step()
        energies.append(atoms.get_total_energy())
        s = atoms.get_array("spins")
        spin_norms.append([np.linalg.norm(si) for si in s])
        total_mags.append(np.sum(s, axis=0))

    energies = np.array(energies)
    spin_norms = np.array(spin_norms)
    total_mags = np.array(total_mags)

    # Spin norm conservation: should be exactly 1.0 (within numerical precision)
    print(f"Spin norms: {spin_norms[-1]}")
    assert np.allclose(spin_norms, 1.0, atol=1e-10)

    # Total magnetization conservation (for Heisenberg)
    dm = np.linalg.norm(total_mags - total_mags[0], axis=1)
    print(f"Initial mag: {total_mags[0]}")
    print(f"Final mag: {total_mags[-1]}")
    print(f"Max magnetization drift: {np.max(dm)}")
    # Relaxing this as ST is 2nd order only
    assert np.allclose(total_mags, total_mags[0], atol=1e-3)

    # Energy conservation: drift should be small
    de = (energies - energies[0]) / energies[0]
    print(f"Initial energy: {energies[0]}")
    print(f"Final energy: {energies[-1]}")
    print(f"Max energy drift: {np.max(np.abs(de))}")
    assert np.allclose(energies, energies[0], rtol=1e-3)


def test_zeeman_precession():
    """Verify Larmor precession with an external field."""
    from ase import units

    # Single atom at rest
    atoms = Atoms("Fe", positions=[[0, 0, 0]])

    # Spin initially in x-direction
    spins = np.array([[1.0, 0.0, 0.0]])
    atoms.set_array("spins", spins)

    # Magnetic field in z-direction: 10 Tesla
    H_ext = [0.0, 0.0, 10.0]

    # Use Lande factor g=2
    calc = SpinLatticeHeisenberg(H_ext=H_ext, g_i=2.0)
    atoms.calc = calc

    # Dynamics
    dt = 0.1 * units.fs
    dyn = SpinLatticeDynamics(atoms, timestep=dt)

    spins_traj = []
    for _ in range(100):
        dyn.step()
        spins_traj.append(atoms.get_array("spins")[0].copy())

    spins_traj = np.array(spins_traj)

    # Expected precession frequency omega = g * mu_B / hbar * B
    # mu_B = 5.788e-5 eV/T, hbar = 0.658 eV*fs
    # for B=10T, g=2: omega = 2 * 5.788e-5 / 0.658 * 10 = 0.001759 rad/fs
    # In 100 steps of 0.1fs (10fs total), theta = 0.01759 rad
    expected_theta = (2.0 * 5.78838181e-5 / 0.6582119569) * 10.0 * 10.0

    actual_s_final = spins_traj[-1]
    # Final spin should be [cos(theta), sin(theta), 0]
    expected_s_final = np.array([np.cos(expected_theta), np.sin(expected_theta), 0.0])

    print(f"Expected final spin: {expected_s_final}")
    print(f"Actual final spin:   {actual_s_final}")

    assert np.allclose(actual_s_final, expected_s_final, atol=1e-5)
    print("Zeeman precession test passed!")


def test_mag_energy_separation():
    """Verify that mag_energy and energy are correctly calculated."""
    # Dimer with H_ext
    atoms = Atoms("Fe2", positions=[[0, 0, 0], [0, 0, 2.5]])
    atoms.set_velocities([[0, 0, 0], [0, 0, 0]])

    # Spins: one along x, one along y
    spins = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    atoms.set_array("spins", spins)

    # Mechanical K=10, r=2.5, r0=2.5 -> Mechanical energy = 0
    # Exchange J0=0.1, alpha=2.0 -> J=0.1. s1.s2 = 0 -> Exchange energy = 0
    # Zeeman: g=[2, 3], H_ext=[0,0,10], s1=[1,0,0], s2=[0,1,0] -> Zeeman energy = 0
    # Wait, let's make it non-zero
    H_ext = [10.0, 0.0, 0.0]  # Field along x
    calc = SpinLatticeHeisenberg(H_ext=H_ext, g_i=[2.0, 3.0], J0=0.1, r0=2.5)
    atoms.calc = calc

    # E_zeeman1 = -mu_B * 2 * 1 * 10
    # E_zeeman2 = -mu_B * 3 * 0 * 10 = 0
    mu_B = 5.7883818012e-5
    expected_zeeman = -mu_B * 2.0 * 10.0

    # E_exch = -J * s1.s2 = 0 (since x . y = 0)

    total_energy = atoms.get_potential_energy()
    mag_energy = atoms.calc.results["mag_energy"]

    print(f"Total potential energy: {total_energy}")
    print(f"Magnetic energy:        {mag_energy}")
    print(f"Expected mag energy:    {expected_zeeman}")

    assert np.allclose(mag_energy, expected_zeeman)
    assert np.allclose(total_energy, mag_energy)  # Mechanical is 0
    print("Magnetic energy separation test passed!")


if __name__ == "__main__":
    test_dimer_conservation()
    test_zeeman_precession()
    test_mag_energy_separation()
