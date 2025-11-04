import sys
import argparse
import random
from pyrosetta import *
from rosetta.core.pack.task import operation
from rosetta.core.pack.task import TaskFactory
from rosetta.core.pack import pack_rotamers
from rosetta.core import kinematics, optimization
from pyrosetta.rosetta.protocols.moves import AddPyMOLObserver, MonteCarlo
'''
init(extra_options="-ignore_unrecognized_res")

# Initialize parser
parser = argparse.ArgumentParser()
parser.add_argument("-s","--structure", required=True, help="PDB file to load")
args = parser.parse_args()

# Load structure 1UBQ
mypose = pose_from_pdb(args.structure)

print(f"Loaded pose with {mypose.total_residue()} residues from: {args.structure}")

# Set up score function 
sfxn = rosetta.core.scoring.get_score_function()
scoring = sfxn(mypose)
print(f"The score of{mypose} is {scoring}")

# Set up task factory 
tf = TaskFactory() 
task = tf.create_task_and_apply_taskoperations(mypose)
task.restrict_to_repacking()

# Set up MoveMap 
movemap = kinematics.MoveMap()
movemap.set_bb(True)
movemap.set_chi(True)

# Set up Minimizer
min_opts = optimization.MinimizerOptions("lbfgs_armijo_atol", 0.01, True)
minimizer = optimization.AtomTreeMinimizer()

temp = 1.0
monte_carlo_initial = MonteCarlo(mypose,sfxn,temp)

# PyMol Observer 
the_observer = AddPyMOLObserver(mypose)
the_observer.pymol().apply(mypose)

# Monte Carlo Loop 
counter_True = 0
counter_False = 0
total =  mypose.total_residue()

for i in range(10):
    pack_rotamers(mypose,sfxn,task)
    minimizer.run(mypose,movemap,sfxn,min_opts)

    randres = random.randrange(1, total)

    phi_pert = random.random()
    psi_pert = random.random()

    orig_phi = mypose.phi(randres)
    orig_psi = mypose.psi(randres)

    mypose.set_phi(randres, orig_phi + phi_pert)
    mypose.set_psi(randres, orig_psi + psi_pert)
    
    accepted_mc = monte_carlo_initial.boltzmann(mypose)
    print(f'The Monte Carlo Acceptance is {accepted_mc}')


    if accepted_mc == True:
        counter_True += 1
    elif accepted_mc == False:
        counter_False += 1

    if (i+1) % 10 == 0:
        accep_rate = counter_True/10
        print(f'The acceptance rate is{accep_rate}')
        avg_energy = mypose.energies().total_energy()/10
        print(f'The average evergy is{avg_energy}')
    else:
        continue 

    print(f"The lowest score is",monte_carlo_initial.lowest_score())
    monte_carlo_initial.lowest_score_pose().dump_pdb('Monte_Carlo.pdb')


print(f'The amount of accepted Monte Carlo Structure is{counter_True}')
print(f'The amount of rejected Monte Carlo Structure is{counter_False}')
'''
def identify_secondary_structure_spans(ss):
    blocks = []
    current = None
    for i,j in enumerate(ss, start = 1):
        if j in ('E','H'):
            if current is None:
                current = j
                start = i
            elif j != current:
                blocks.append((start, i-1))
                current = j
                start = i
        else:
            if current is not None:
                blocks.append((start, i-1))
                current = None
                start = None
    if current is not None:
        blocks.append((start,len(ss)))
    return blocks
