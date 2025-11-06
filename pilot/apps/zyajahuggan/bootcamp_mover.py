from pyrosetta import *
from rosetta.core.pack.task import operation
from rosetta.core.pack.task import TaskFactory
from rosetta.core.pack import pack_rotamers
from rosetta.core import kinematics, optimization
from pyrosetta.rosetta.protocols.moves import AddPyMOLObserver, MonteCarlo
from pyrosetta.rosetta.core.scoring.dssp import Dssp
from pyrosetta.rosetta.core.kinematics import FoldTree
from pyrosetta.rosetta.core.pose import correctly_add_cutpoint_variants
from pyrosetta.rosetta.core.scoring import ScoreType, get_score_function


class BootCampMover:
    def __init__(self):
        super().__init__(self) # or rosetta.protocols.moves.Mover.__init__(self)

    def apply(self,pose):
        ss = Dssp(pose).get_dssp_secstruct()

        ft = fold_tree_from_dssp_string(ss)  
        pose.fold_tree(ft)

        correctly_add_cutpoint_variants(pose)

        sfxn = get_score_function()
        sfxn.set_weight(ScoreType.linear_chainbreak, 1.0)

        score = sfxn(pose)
        print("Score after fold tree + cutpoints:", score)
                
    def get_name(self):
        self.__class__.__name__