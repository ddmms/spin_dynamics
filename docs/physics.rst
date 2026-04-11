Physics and Equations
=====================

The ``spindynam`` package implements the coupled spin-lattice dynamics formalism
described in **Tranchida et al., J. Comput. Phys. 372 (2018) 406-425**. 

Hamiltonian
-----------

The total Hamiltonian of the coupled spin-lattice system is given by Equation 4:

.. math::

   \mathcal{H}_{sl} = \mathcal{H}_{mag} + \sum_{i=1}^N \frac{|\mathbf{p}_i|^2}{2m_i} + \sum_{i<j} V(r_{ij})

Where:
- :math:`\mathbf{p}_i` is the momentum of atom :math:`i`.
- :math:`m_i` is the mass of atom :math:`i`.
- :math:`V(r_{ij})` is the interatomic potential (mechanical part).

The magnetic Hamiltonian :math:`\mathcal{H}_{mag}` includes Zeeman and Exchange interactions (Equation 3):

.. math::

   \mathcal{H}_{mag} = - \mu_B \sum_{i=1}^N g_i \mathbf{s}_i \cdot \mathbf{H}_{ext} - \sum_{i<j} J(r_{ij}) \mathbf{s}_i \cdot \mathbf{s}_j

Where:
- :math:`\mathbf{s}_i` is the atomic spin (unit vector).
- :math:`J(r_{ij})` is the distance-dependent exchange interaction.
- :math:`\mathbf{H}_{ext}` is the external magnetic field.
- :math:`\mu_B` is the Bohr magneton and :math:`g_i` is the Landé g-factor.

Functional Forms
----------------

The ``SpinLatticeHeisenberg`` calculator utilizes the following functional forms for the interactions:

Exchange Interaction :math:`J(r)`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The exchange coupling follows an exponential decay with distance:

.. math::

   J(r) = J_0 \exp\left[-\alpha \left(\frac{r}{r_0} - 1.0\right)\right]

Where:
- :math:`J_0` is the exchange constant at the equilibrium distance :math:`r_0`.
- :math:`\alpha` is the dimensionless decay parameter.

Atomic Potential :math:`V(r)`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The interatomic mechanical potential is modeled as a harmonic well:

.. math::

   V(r) = \frac{1}{2} K (r - r_0)^2

Where:
- :math:`K` is the spring constant.
- :math:`r_0` is the equilibrium distance.

Equations of Motion
-------------------

The time evolution of any observable :math:`A` in the coupled spin-lattice phase 
space is governed by the Poisson bracket (Equation 5):

.. math::

   \{F, G\} = \sum_{i=1}^N \left( \frac{\partial F}{\partial \mathbf{r}_i} \cdot \frac{\partial G}{\partial \mathbf{p}_i} - \frac{\partial F}{\partial \mathbf{p}_i} \cdot \frac{\partial G}{\partial \mathbf{r}_i} \right) + \sum_{i=1}^N \frac{\partial F}{\partial \mathbf{s}_i} \cdot \left( \frac{\mathbf{s}_i}{\hbar} \times \frac{\partial G}{\partial \mathbf{s}_i} \right)

The resulting equations of motion for the positions, momenta, and spins are 
given by (Equation 12):

.. math::

   \dot{\mathbf{r}}_i &= \{\mathbf{r}_i, \mathcal{H}_{sl}\} = \frac{\mathbf{p}_i}{m_i} \\
   \dot{\mathbf{p}}_i &= \{\mathbf{p}_i, \mathcal{H}_{sl}\} = \mathbf{F}_i \\
   \dot{\mathbf{s}}_i &= \{\mathbf{s}_i, \mathcal{H}_{sl}\} = \boldsymbol{\omega}_i \times \mathbf{s}_i

Forces and Effective Fields
---------------------------

The total force :math:`\mathbf{F}_i` acting on atom :math:`i` consists of mechanical 
and magnetic contributions (Equation 7):

.. math::

   \mathbf{F}_i = - \frac{\partial \mathcal{H}_{sl}}{\partial \mathbf{r}_i} = - \nabla_{\mathbf{r}_i} \sum_{j < k} V(r_{jk}) + \sum_{j \neq i} \frac{\partial J(r_{ij})}{\partial \mathbf{r}_i} (\mathbf{s}_i \cdot \mathbf{s}_j)

The magnetic effective field :math:`\boldsymbol{\omega}_i` (or torque) acting on 
spin :math:`i` is defined as (Equation 8):

.. math::

   \boldsymbol{\omega}_i = \frac{1}{\hbar} \left( g_i \mu_B \mathbf{H}_{ext} + \sum_{j \neq i} J(r_{ij}) \mathbf{s}_j \right)

Symplectic Integrator
---------------------

To preserve the geometric structure of the phase space (including the spin norm 
and total energy for NVE), we use a second-order symmetric Suzuki-Trotter 
decomposition (Equation 13):

.. math::

   e^{\mathcal{L} \Delta t} \approx e^{\mathcal{L}_p \frac{\Delta t}{2}} e^{\mathcal{L}_s \frac{\Delta t}{2}} e^{\mathcal{L}_r \Delta t} e^{\mathcal{L}_s \frac{\Delta t}{2}} e^{\mathcal{L}_p \frac{\Delta t}{2}}

Where:
- :math:`\mathcal{L}_p` updates momenta using forces.
- :math:`\mathcal{L}_r` updates positions using momenta.
- :math:`\mathcal{L}_s` is the magnetic Liouvillian, updated via a sequential rotation of each spin.

Spin Rotation (Equation 15)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Individual spin updates are performed using a norm-preserving rational propagator:

.. math::

   \mathbf{s}_i(t+\Delta t) = \mathbf{s}_i(t) + \frac{\Delta t (\boldsymbol{\omega}_i \times \mathbf{s}_i) + \frac{\Delta t^2}{2} \boldsymbol{\omega}_i \times (\boldsymbol{\omega}_i \times \mathbf{s}_i)}{1 + \frac{\Delta t^2}{4} |\boldsymbol{\omega}_i|^2}

This formulation avoids the need for ad-hoc normalization and strictly preserves 
the spin magnitude throughout the simulation.

References
----------

* **Tranchida et al.**, *"A symplectic algorithm for coupled atomistic spin-lattice dynamics"*, **Journal of Computational Physics**, 372 (2018) 406-425. `DOI: 10.1016/j.jcp.2018.06.042 <https://doi.org/10.1016/j.jcp.2018.06.042>`_
