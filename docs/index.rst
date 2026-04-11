Welcome to spindynam's documentation!
=====================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   modules
   physics

Overview
--------

** this is vibe coding experiment, do not use it for production **

spindynam is a Python package for coupled spin-lattice dynamics simulations
using the Atomic Simulation Environment (ASE).

**Reference Paper**: :download:`Tranchida et al., J. Comput. Phys. 372 (2018) 406-425  <spin_dynamics.pdf>`

It implements the symplectic integration algorithm and equations of motion
described in **Tranchida et al., J. Comput. Phys. 372 (2018) 406-425**.

Key Features:
- Coupled SD-MD with symplectic ST decomposition.
- Magneto-elastic coupling via distance-dependent exchange interactions.
- Optimized geometric spin updates (Equation 15).
- Seamless integration with the ASE calculator interface.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
