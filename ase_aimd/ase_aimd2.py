import os
import subprocess
from ase import Atoms
from ase.io import read, write, Trajectory
from ase.md import Langevin
from ase.md.nvtberendsen import NVTBerendsen
import ase.units as units
from ase.calculators.vasp import Vasp
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from bader import attach_charges

VASP_FLAGS = {
    "nsw": 0,
    "nelm": 300,
    "isif": 0,
    "isym": 0,
    "prec": "Normal",
    "lreal": "Auto",
    "ediff": 1e-5,
    "encut": 400.0,
    "ismear": 0,
    "sigma": 0.05,
    "algo": "Fast",
    "laechg": True,
    "lwave": True,
    "lcharg": True,
    "istart": 1,
    "icharg": 1,
    "ivdw": 11,
    "xc": "PBE",
    "kpts": (1, 1, 1),
    "gamma": True,
}

# Define the initial atomic structure
atoms = read('initial.vasp', format='vasp')

# Set up the VASP calculator
calc = Vasp(directory='vasp_sp', **VASP_FLAGS)
atoms.set_calculator(calc)

# Set the momenta corresponding to T=100K
MaxwellBoltzmannDistribution(atoms, temperature_K=500)

# Create Langevin object
aimd = NVTBerendsen(atoms, timestep=1 * units.fs,
                    temperature=500 * units.kB,
                    taut=0.5*600*units.fs,
                    logfile='aimd_with_bader.log')

# Create a trajectory file to save the simulation snapshots
traj = Trajectory('aimd_with_bader.traj', 'w', atoms)
# aimd.attach(traj.write, interval=1)

# Run the AIMD simulation
total_step = 5000
for step in range(total_step):
    aimd.run(1)  # Run 1 timesteps per iteration

    # Capture energy and force during the simulation
    energy = atoms.get_potential_energy()
    forces = atoms.get_forces()

    # Calculate Bader charges using subprocess
    os.chdir('vasp_sp')
    subprocess.run(['chgsum.pl', 'AECCAR0', 'AECCAR2'])
    subprocess.run(['bader', 'CHGCAR', '-ref', 'CHGCAR_sum'])
    attach_charges(atoms, 'ACF.dat')
    os.chdir("..")

    traj.write(atoms)

# Close.traj file
traj.close()  
print('AIMD simulation with Bader charge calculation completed.')
