# preparation stage
## initialize molecule system and minimize energy for MD simulation

# from openmm.app import * #read pdb files, define topology and force field
from openmm.app import PDBFile, ForceField, Simulation, NoCutoff, HBonds
from openmm import Platform, LangevinIntegrator #contain core classes
from openmm.unit import kelvin, picosecond, picoseconds #contain units
import sys
import os
import time

# if len(sys.argv) < 2:
#     print("❌ Error: Please provide a PDB filename as an argument.")
#     print("✅ Usage: python setup_system.py your_protein.pdb")
#     sys.exit(1)
#
# pdb_filename = sys.argv[1]
# pdb = PDBFile(pdb_filename)
# print(f"📥 Loading PDB file: {pdb_filename}")
# basename = os.path.splitext(os.path.basename(pdb_filename))[0]
# output_dir = input("Where do yo want to save the minimized pdb? Enter full path, or leave empty for current director: ").strip()
# if output_dir == '':
#     output_dir = '.'

input_pdb_path = "target_protein.pdb"  # <--- Change this to your actual PDB file
# Provide the output directory (use '.' for current directory)
output_dir = "./results/minimized"

if not os.path.isfile(input_pdb_path):
    raise FileNotFoundError(f"❌ Input file not found: {input_pdb_path}")

basename = os.path.splitext(os.path.basename(input_pdb_path))[0]
output_filename = os.path.join(output_dir, f"{basename}_minimized.pdb")

os.makedirs(output_dir, exist_ok=True)

forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
#load force field parameters,
# define the interaction principles between molecules
# amber14-all.xml - standard AMBER14 force field for proteins and nucleic acids
#amber14/tip3pfb.xml - water model parameters (TIP3P with improved bond lengths and angles)

system = forcefield.createSystem(pdb.topology, nonbondedMethod=NoCutoff, constraints=HBonds)
#create a system object, containing all forces and particles for simulation,
# contains all force information that need to be simulated (potential energy function)
#nonbondedMethod=NoCutoff means nocutoff van der waals forces (suitable for vacumm simulations) (no long-range electrostatics)
#constraints=HBonds constraints bonds involving hydrogen atoms, allowing larger timestep (2 fs)

integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
#set MD integrator, controlling how the system proceed as times goes
#use Langevin Dynamics, brownian motion at simulated temp, containing friction and noise terms
#300*kelvin: target temp, default set
#1/picosecond: friction coefficient (damping )
#0.002*picoseconds: timestep

platform = Platform.getPlatformByName('CPU')
# specify running platform

simulation = Simulation(pdb.topology, system, integrator, platform)
simulation.context.setPositions(pdb.positions)

print("🛠 Starting energy minimization...")
start_time = time.time()
simulation.minimizeEnergy()
end_time = time.time()
print(f"✔️energy minimization complete in {end_time - start_time:.2f} seconds.")


with open(output_filename, 'w') as f:
    PDBFile.writeFile(simulation.topology,
                      simulation.context.getState(getPositions=True).getPositions(),
                      f)
print(f"Minimized structure saved to: {output_filename}")

simulation.context.getState(getPositions=True, getVelocities=False)

## ref: https://openmm.org/documentation.html chap3&4
##