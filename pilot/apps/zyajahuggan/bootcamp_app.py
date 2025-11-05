import sys
import argparse
import random
from pyrosetta import *
from rosetta.core.pack.task import operation
from rosetta.core.pack.task import TaskFactory
from rosetta.core.pack import pack_rotamers
from rosetta.core import kinematics, optimization
from pyrosetta.rosetta.protocols.moves import AddPyMOLObserver, MonteCarlo
from pyrosetta.rosetta.core.scoring.dssp import Dssp
from pyrosetta.rosetta.core.kinematics import FoldTree
from pyrosetta.rosetta.core.pose import correctly_add_cutpoint_variants


def identify_secondary_structure_spans(ss):
    blocks = []
    current = None
    start = None

    for i, j in enumerate(ss, start=1):
        # If we are entering a structured segment
        if j in ('E', 'H'):
            if current is None:
                current = j
                start = i
            elif j != current:
                blocks.append((start, i-1))
                current = j
                start = i
        else:  # j == 'L'
            if current is not None:
                blocks.append((start, i-1))
                current = None
                start = None
    # Close final block if necessary
    if current is not None:
        blocks.append((start, len(ss)))
    return blocks

def fold_tree_from_dssp_string(ss_string):
    N = len(ss_string)

    # --- Identify SSE blocks (runs of E/H) ---
    sse = []
    current = None
    start = None
    for i, ch in enumerate(ss_string, start=1):
        if ch in ('E','H'):
            if current is None:
                current = ch
                start = i
            elif ch != current:
                sse.append((start, i-1))
                current = ch
                start = i
        else:
            if current is not None:
                sse.append((start, i-1))
                current = None
                start = None
    if current is not None:
        sse.append((start, N))

    if not sse:
        raise ValueError("No helices or strands detected.")

    # --- Compute midpoints ---
    def mid(a,b): return (a+b)//2
    sse_mid = [mid(a,b) for (a,b) in sse]

    # --- Identify loops between SSEs ---
    loops = []
    loop_mid = []
    for (a1,b1),(a2,b2) in zip(sse, sse[1:]):
        if a2 - b1 > 1:
            loops.append((b1+1, a2-1))
            loop_mid.append(mid(b1+1, a2-1))

    root = sse_mid[0]  # anchor midpoint

    ft = FoldTree()
    jump_id = 1

    # --- First SSE: back to 1, forward to SSE1 end ---
    a0,b0 = sse[0]
    if root > 1:
        ft.add_edge(root, 1, -1)
    if root < b0:
        ft.add_edge(root, b0, -1)

    # --- Interleave each loop and the next SSE ---
    rows = max(len(loops), len(sse)-1)
    for k in range(rows):

        # Loop branch
        if k < len(loops):
            lm = loop_mid[k]; la,lb = loops[k]
            ft.add_edge(root, lm, jump_id); jump_id += 1
            if lm > la: ft.add_edge(lm, la, -1)
            if lm < lb: ft.add_edge(lm, lb, -1)

        # Next SSE branch
        if (k+1) < len(sse):
            sm = sse_mid[k+1]; sa,sb = sse[k+1]
            ft.add_edge(root, sm, jump_id); jump_id += 1
            forward_end = N if (k+1)==len(sse)-1 else sb
            if sm > sa: ft.add_edge(sm, sa, -1)
            if sm < forward_end: ft.add_edge(sm, forward_end, -1)

    try:
        ft.reorder(root)
    except Exception:
        pass

    return ft

def fold_tree_from_ss(pose):
    ss = Dssp(pose).get_dssp_secstruct()
    ss = ss[:pose.total_residue()].ljust(pose.total_residue(), 'L')
    return fold_tree_from_dssp_string(ss)


if __name__ == "__main__":

    import argparse 
    from pyrosetta import pose_from_pdb
    from pyrosetta.rosetta.core.scoring import ScoreType

    init(extra_options="-ignore_unrecognized_res")

    # Initialize parser
    parser = argparse.ArgumentParser()
    parser.add_argument("-s","--structure", required=False, help="PDB file to load")
    args = parser.parse_args()

    # Load structure 1UBQ
    mypose = pose_from_pdb(args.structure)
    print(f"Loaded pose with {mypose.total_residue()} residues from: {args.structure}")

    correctly_add_cutpoint_variants(mypose)
    # Set up score function 
    sfxn = rosetta.core.scoring.get_score_function()
    sfxn.set_weight(ScoreType.linear_chainbreak, 1.0)
    scoring = sfxn(mypose)
    print(f"The score of{mypose} is {scoring}")


    # Initialize parser
    parser = argparse.ArgumentParser()
    parser.add_argument("-s","--structure", required=False, help="PDB file to load")
    args = parser.parse_args()


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

    for i in range(100):
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

        if (i+1) % 100 == 0:
            accep_rate = counter_True/10
            print(f'The acceptance rate is{accep_rate}')
            avg_energy = mypose.energies().total_energy()/100
            print(f'The average evergy is{avg_energy}')
        else:
            continue 

        print(f"The lowest score is",monte_carlo_initial.lowest_score())
        monte_carlo_initial.lowest_score_pose().dump_pdb('Monte_Carlo.pdb')


    print(f'The amount of accepted Monte Carlo Structure is{counter_True}')
    print(f'The amount of rejected Monte Carlo Structure is{counter_False}')

    #ft = fold_tree_from_ss(mypose)
    #print(ft)


    #dssp = Dssp(mypose)
    #ss_raw = dssp.get_dssp_secstruct()
    #ss_clean = ''.join(ch if ch in ('H','E','L','C','-') else 'L' for ch in ss_raw)
    ##print(dssp)
    #N = mypose.total_residue()
    #ss = ss_clean[:N].ljust(N, 'L')
    #print(ss)

    #print( identify_secondary_structure_spans(ss) )
    #print(get_midpoints_with_loops(ss))
    #[(2, 7), (12, 16), (22, 22), (23, 34), (38, 40), (41, 45), (48, 49), (55, 55), (57, 59), (66, 71)]
    #[(2, 7), (12, 16), (22, 22), (23, 34), (38, 40), (41, 45), (48, 49), (55, 55), (57, 59), (66, 71)]
    #print("len(ss) =", len(ss), "pose length =", mypose.total_residue())
    #dssp = Dssp(mypose)
    #ss = dssp.get_dssp_secstruct()
    #print(fold_tree_from_dssp_string(ss))
    #ft = fold_tree_from_ss(mypose)
    #print(ft)

    # My acceptance rate was 42% which was worse than when I didnt have jumps. Energy is lower but the acceptance rate is lower