from pyrosetta import *
from pyrosetta.rosetta.core.scoring.dssp import Dssp
from pyrosetta.rosetta.core.kinematics import FoldTree
from dataclasses import dataclass
from typing import List
from pyrosetta.rosetta.protocols.loops import Loop
from bootcamp_app import fold_tree_from_ss, fold_tree_from_dssp_string

@dataclass
class _Built:
    ft: FoldTree
    loops: List[Loop]
    loop_for_residue: List[int]  # 1..N, values in 0..len(loops)

class FoldTreeFromSS:

    def __init__(self, pose: rosetta.core.pose.Pose, loop_left: int = 2, loop_right: int = 3):
        self.pose = pose
        self.loop_left = loop_left
        self.loop_right = loop_right 
   
    def fold_tree(self,ss_string) -> FoldTree:
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

    def loop(self, index: int) -> Loop:
        '''Find start, stop and cutpoints of loops that need to be closed'''
        ft = self.fold_tree_from_ss(self.pose) # Get connected fold tree 
        N = self.pose.total_residue() # Total length of pose 
        
        cutpoints = [] #store all cutpoints

        for i in range(1,N):
            if ft.is_cutpoint(i):
                cutpoints.append(i)

        if index < 1 or index > len(cutpoints):
            raise IndexError(f"Loop index {index} out of range (1..{len(cutpoints)})")
        
        cut = cutpoints[index -1]
        start = max(1, cut - self.loop_left)
        stop = min(N,cut + self.loop_right)

        return Loop(start,stop,cut)

    def loop_for_residue(self, seqpos: int) -> int:
        '''Finds if a particular residue is a part of a loop'''
        ft = fold_tree_from_ss(self.pose)
        N = self.pose.total_residue()

        cutpoints = []

        for i in range(1,N):
            if ft.is_cutpoint(i):
                cutpoints.append(i)

        for i, cut in enumerate(cutpoints, start = 1):
            start = max(1,cut - self.loop_left)
            stop = min(N, cut + self.loop_right)
            if start <= seqpos <= stop:
                return i
        return 0 


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

    def fold_tree_from_ss(self,pose):
        ss = Dssp(pose).get_dssp_secstruct()
        ss = ss[:pose.total_residue()].ljust(pose.total_residue(), 'L')
        return self.fold_tree(ss)
