import os
import subprocess
from ase import Atoms
from ase.io import read, write, Trajectory
from ase.calculators.vasp import Vasp
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.verlet import VelocityVerlet
import ase.units as units
from bader import attach_charges
from ase.md import Langevin

VASP_FLAGS = {
    "nsw": 0,
    "nelm": 300,
    "isif": 2,
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

# # Set the temperature of the system
# temperature_initial = 500  # the initial temperature
# temperature_final = 400  # the final temperature
# MaxwellBoltzmannDistribution(atoms, temperature_initial)

# Set up the VASP calculator
calc = Vasp(directory="vasp_sp", **VASP_FLAGS)
atoms.set_calculator(calc)

# Create .traj file to save trajectory
trajectory_file = 'aimd_with_bader.traj'
trajectory = Trajectory(trajectory_file, 'w')

# Define tloghe simulation parameters
total_time = 10  # the total simulation time500K-2ps
time_step = 1.0  # time_step（unit：fs）

# Define the simulation parameters of the thermostatic process
dyn = Langevin(atoms, timestep=1 * units.fs,
                temperature=500 * units.kB,
                friction=0.01,
                logfile='aimd_with_bader.log')
for step in range(int(total_time / time_step)):
    dyn.run(1)

    # Capture energy and force during the simulation
    energy = atoms.get_potential_energy()
    forces = atoms.get_forces()

    # Calculate Bader charges using subprocess
    os.chdir("vasp_sp")
    subprocess.run(['chgsum.pl', 'AECCAR0', 'AECCAR2'])
    subprocess.run(['bader', 'CHGCAR', '-ref', 'CHGCAR_sum'])
    # TODO: modify attach_charges func for vasp
    attach_charges(atoms, 'ACF.dat')
    os.chdir("..")

    # Save the trajectory to the .traj file
    trajectory.write(atoms,energy=energy, forces=forces)

# Define the simulation parameters in the process of changing temperature
total_steps_ramp = 5  # the total steps in the changing temperature process
temperature_initial = 500
temperature_final = 400

# Temperature_step
temperature_step = (temperature_final - temperature_initial) / total_steps_ramp

dyn_ramp = Langevin(atoms, 1.0, temperature_initial,logfile='aimd_with_bader.log2')

# Conduct the changing temperature process and save the trajectory
for step in range(total_steps_ramp):

    # Set the temperature of the current step
    current_temperature = temperature_initial + step * temperature_step
    dyn_ramp.set_temperature(current_temperature)

    # Run one step
    dyn_ramp.run(1)

    # Capture energy and force during the simulation
    energy = atoms.get_potential_energy()
    forces = atoms.get_forces()

    # Calculate Bader charges using subprocess
    os.chdir("vasp_sp")
    subprocess.run(['chgsum.pl', 'AECCAR0', 'AECCAR2'])
    subprocess.run(['bader', 'CHGCAR', '-ref', 'CHGCAR_sum'])
    # TODO: modify attach_charges func for vasp
    attach_charges(atoms, 'ACF.dat')
    os.chdir("..")

    # Save the trajectory to .traj file
    trajectory.write(atoms,energy=energy, forces=forces)
 
 # Close.traj file
trajectory.close()  

