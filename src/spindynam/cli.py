#!/usr/bin/env python3
import argparse

import numpy as np
from ase.io import read
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution

from .calculator import SpinLatticeHeisenberg
from .dynamics import SpinLatticeDynamics


def get_parser():
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
    parser.add_argument("--j0", type=float, default=0.1, help="Exchange J0 (eV)")
    parser.add_argument("--alpha", type=float, default=2.0, help="Decay alpha")
    parser.add_argument("--k", type=float, default=10.0, help="Spring K (eV/A^2)")
    parser.add_argument("--r0", type=float, help="Equilibrium r0")
    parser.add_argument(
        "--hext",
        type=float,
        nargs=3,
        default=[0, 0, 0],
        help="External field B (Tesla)",
    )
    parser.add_argument("--g", type=float, default=2.0, help="Lande factor")

    # Simulation Parameters
    parser.add_argument("--temp", "-t", type=float, default=0.0, help="Temp (K)")
    parser.add_argument(
        "--initial-spins",
        choices=["random", "ferro", "file"],
        default="ferro",
        help="Spin config",
    )
    parser.add_argument("--loginterval", type=int, default=100, help="Log interval")
    return parser


def init_atoms(args):
    atoms = read(args.input)
    if args.temp > 0:
        MaxwellBoltzmannDistribution(atoms, temperature_K=args.temp)

    if args.initial_spins == "random":
        spins = np.random.uniform(-1, 1, (len(atoms), 3))
        atoms.set_array("spins", spins / np.linalg.norm(spins, axis=1)[:, np.newaxis])
    elif not (args.initial_spins == "file" and "spins" in atoms.arrays):
        spins = np.zeros((len(atoms), 3))
        spins[:, 2] = 1.0
        atoms.set_array("spins", spins)
    return atoms


def setup_calc(args, atoms):
    r0 = args.r0 or (atoms.get_all_distances().mean() if len(atoms) > 1 else 2.5)
    calc = SpinLatticeHeisenberg(
        J0=args.j0, alpha=args.alpha, K=args.k, r0=r0, H_ext=args.hext, g_i=args.g
    )
    atoms.calc = calc
    return atoms


def run_sim(args, atoms):
    dyn = SpinLatticeDynamics(
        atoms, timestep=args.dt, trajectory=args.output, loginterval=args.loginterval
    )

    def print_status():
        epot, ekin = atoms.get_potential_energy(), atoms.get_kinetic_energy()
        avg_mag = np.mean(atoms.get_array("spins"), axis=0)
        print(
            f"Step: {dyn.get_number_of_steps():>6d} | E: {epot + ekin:>10.4f} | "
            f"|M|: {np.linalg.norm(avg_mag):>6.3f}"
        )

    dyn.attach(print_status, interval=args.loginterval)
    print_status()
    dyn.run(args.steps)


def main():
    args = get_parser().parse_args()
    atoms = init_atoms(args)
    atoms = setup_calc(args, atoms)
    run_sim(args, atoms)


if __name__ == "__main__":
    main()
