Welcome to WaveBlocksND's documentation!
========================================

.. image:: _static/waveblocks.png

Reusable building blocks for simulations with semiclassical wavepackets for
solving the time-dependent Schrödinger equation.

Contents:

.. toctree::
   :maxdepth: 2

Source code documentation
=========================


WaveBlocks Classes
-------------------

Basic numerics
^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 1

   waveblocks_classes/ComplexMath

   waveblocks_classes/Grid
   waveblocks_classes/DenseGrid
   waveblocks_classes/TensorProductGrid
   waveblocks_classes/GridWrapper

   waveblocks_classes/QuadratureRule
   waveblocks_classes/GaussHermiteQR
   waveblocks_classes/TensorProductQR

   waveblocks_classes/MatrixExponential

   waveblocks_classes/Utils

Basic quantum mechanics
^^^^^^^^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 1

   waveblocks_classes/WaveFunction

   waveblocks_classes/MatrixPotential
   waveblocks_classes/MatrixPotential1S
   waveblocks_classes/MatrixPotential2S
   waveblocks_classes/MatrixPotentialMS

   waveblocks_classes/PotentialFactory

Wavepackets
^^^^^^^^^^^

.. toctree::
   :maxdepth: 1

   waveblocks_classes/BasisShape
   waveblocks_classes/HyperCubicShape
   waveblocks_classes/HyperbolicCutShape
   waveblocks_classes/LimitedHyperbolicCutShape

   waveblocks_classes/Wavepacket
   waveblocks_classes/HagedornWavepacketBase
   waveblocks_classes/HagedornWavepacket
   waveblocks_classes/HagedornWavepacketInhomogeneous

   waveblocks_classes/Quadrature
   waveblocks_classes/HomogeneousQuadrature
   waveblocks_classes/InhomogeneousQuadrature

   waveblocks_classes/BasisTransformation
   waveblocks_classes/BasisTransformationWF
   waveblocks_classes/BasisTransformationHAWP


Observables
^^^^^^^^^^^

.. toctree::
   :maxdepth: 1

   waveblocks_classes/Observables
   waveblocks_classes/ObservablesHAWP


Time propagation
^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 1

   waveblocks_classes/KineticOperator
   waveblocks_classes/Propagator
   waveblocks_classes/FourierPropagator
   waveblocks_classes/HagedornPropagator
   waveblocks_classes/HagedornPropagatorInhomogeneous
   waveblocks_classes/SimulationLoop
   waveblocks_classes/SimulationLoopFourier
   waveblocks_classes/SimulationLoopHagedorn
   waveblocks_classes/SimulationLoopHagedornInhomogeneous


Simulation result storage I/O
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. toctree::
   :maxdepth: 1

   waveblocks_classes/IOManager

   waveblocks_classes/IOM_plugin_parameters
   waveblocks_classes/IOM_plugin_grid
   waveblocks_classes/IOM_plugin_wavefunction
   waveblocks_classes/IOM_plugin_fourieroperators
   waveblocks_classes/IOM_plugin_wavepacket
   waveblocks_classes/IOM_plugin_inhomogwavepacket
   waveblocks_classes/IOM_plugin_norm
   waveblocks_classes/IOM_plugin_energy


Other classes
^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 1

   waveblocks_classes/GlobalDefaults
   waveblocks_classes/ParameterLoader
   waveblocks_classes/ParameterProvider
   waveblocks_classes/FileTools

Etc
===

.. toctree::
   :maxdepth: 2

   citation


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
