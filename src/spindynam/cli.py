#!/usr/bin/env python3
import argparse

import numpy as np
from ase.io import read
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution

from .calculator import SpinLatticeHeisenberg
from .dynamics import SpinLatticeDynamics


def main():
    parser = argparse.ArgumentParser(
        description="Coupled Spin-Lattice Dynamics Simulation CLI"
    )

    # Input/Output
    parser.add_argument(
        "--input",
        "-i",
        type=str,
        required=True,
        help="Input structure file (ASE readable)",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        default="simulation.traj",
        help="Trajectory output file",
    )
    parser.add_argument(
        "--steps", "-s", type=int, default=1000, help="Number of MD steps"
    )
    parser.add_argument("--dt", type=float, default=0.1, help="Timestep in fs")

    # Hamiltonian Parameters
    parser.add_argument(
        "--j0", type=float, default=0.1, help="Exchange constant J0 (eV)"
    )
    parser.add_argument(
        "--alpha", type=float, default=2.0, help="Exchange decay constant alpha"
    )
    parser.add_argument(
        "--k", type=float, default=10.0, help="Mechanical spring constant K (eV/A^2)"
    )
    parser.add_argument(
        "--r0",
        type=float,
        help="Equilibrium distance (defaults to initial mean distance)",
    )
    parser.add_argument(
        "--hext",
        type=float,
        nargs=3,
        default=[0.0, 0.0, 0.0],
        help="External magnetic field in Tesla (Bx By Bz)",
    )
    parser.add_argument("--g", type=float, default=2.0, help="Lande g-factor")

    # Simulation Parameters
    parser.add_argument(
        "--temp",
        "-t",
        type=float,
        default=0.0,
        help="Initial temperature for velocities (K)",
    )
    parser.add_argument(
        "--initial-spins",
        choices=["random", "ferro", "file"],
        default="ferro",
        help="Initial spin configuration",
    )
    parser.add_argument(
        "--loginterval",
        type=int,
        default=100,
        help="Interval for logging and traj saving",
    )

    args = parser.parse_args()

    # 1. Load system
    print(f"Loading structure from {args.input}...")
    atoms = read(args.input)
    natoms = len(atoms)

    # 2. Initialize Velocities
    if args.temp > 0:
        print(f"Initializing velocities for T = {args.temp} K...")
        MaxwellBoltzmannDistribution(atoms, temperature_K=args.temp)

    # 3. Initialize Spins
    if args.initial_spins == "random":
        print("Initializing spins randomly...")
        spins = np.random.uniform(-1, 1, (natoms, 3))
        norms = np.linalg.norm(spins, axis=1)[:, np.newaxis]
        spins /= norms
        atoms.set_array("spins", spins)
    elif args.initial_spins == "file" and "spins" in atoms.arrays:
        print("Using spins from input file.")
    else:
        print("Initializing spins ferromagnetically (along z)...")
        spins = np.zeros((natoms, 3))
        spins[:, 2] = 1.0
        atoms.set_array("spins", spins)

    # 4. Setup Calculator
    r0 = args.r0
    if r0 is None:
        if natoms > 1:
            r0 = atoms.get_all_distances().mean()
        else:
            r0 = 2.5

    print(
        f"Setting up SpinLatticeHeisenberg calculator "
        f"(J0={args.j0}, alpha={args.alpha}, K={args.k}, r0={r0:.3f})..."
    )
    calc = SpinLatticeHeisenberg(
        J0=args.j0, alpha=args.alpha, K=args.k, r0=r0, H_ext=args.hext, g_i=args.g
    )
    atoms.calc = calc

    # 5. Setup Dynamics
    print(f"Initializing SpinLatticeDynamics (dt={args.dt} fs)...")
    dyn = SpinLatticeDynamics(
        atoms, timestep=args.dt, trajectory=args.output, loginterval=args.loginterval
    )

    # 6. Logging Function
    def print_status():
        epot = atoms.get_potential_energy()
        ekin = atoms.get_kinetic_energy()
        spins = atoms.get_array("spins")
        avg_mag = np.mean(spins, axis=0)
        mag_norm = np.linalg.norm(avg_mag)
        msg = (
            f"Step: {dyn.get_number_of_steps():>6d} | "
            f"E_tot: {epot + ekin:>10.4f} eV | "
            f"Mag: [{avg_mag[0]:>6.3f}, {avg_mag[1]:>6.3f}, {avg_mag[2]:>6.3f}] | "
            f"|M|: {mag_norm:>6.3f}"
        )
        print(msg)

    dyn.attach(print_status, interval=args.loginterval)

    # 7. Run
    print(f"Running simulation for {args.steps} steps...")
    print_status()  # Initial status
    dyn.run(args.steps)
    print("Simulation complete.")


if __name__ == "__main__":
    main()
