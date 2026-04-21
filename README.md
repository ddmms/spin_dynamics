# spindynam

[![CI](https://github.com/ddmms/spin_dybamics/actions/workflows/ci.yml/badge.svg)](https://github.com/ddmms/spin_dynamics/actions/workflows/ci.yml)
[![Documentation](https://github.com/ddmms/spin_dybamics/actions/workflows/docs_deploy.yml/badge.svg)](https://ddmms.github.io/spin_dynamics/)
[![Coverage Status](https://coveralls.io/repos/github/ddmms/spin_dynamics/badge.svg?branch=main)](https://coveralls.io/github/ddmms/spin_dynamics?branch=main)
[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

Coupled Spin-Lattice Dynamics for ASE based on Tranchida et al. (2018).

**this is vibe coding experiment, do not use it for production**

[**Read the Documentation**](https://ddmms.github.io/spin_dynamics/)

## Installation

```bash
uv sync
```

## Usage

```bash
uv run spin-dynamics --input structure.xyz --steps 1000
```

