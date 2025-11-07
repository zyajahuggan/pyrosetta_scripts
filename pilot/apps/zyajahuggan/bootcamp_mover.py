import random
from pyrosetta import *
from pyrosetta.rosetta.core.pack.task import operation
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack import pack_rotamers
from pyrosetta.rosetta.core import kinematics, optimization
from pyrosetta.rosetta.protocols.moves import AddPyMOLObserver, MonteCarlo
from pyrosetta.rosetta.core.scoring.dssp import Dssp
from pyrosetta.rosetta.core.kinematics import FoldTree
from pyrosetta.rosetta.core.pose import correctly_add_cutpoint_variants
from pyrosetta.rosetta.core.scoring import ScoreType, get_score_function, parse_score_function, attributes_for_parse_score_function_w_description
from pyrosetta.rosetta.utility.tag import XMLSchemaAttribute, XMLSchemaComplexTypeGenerator, XMLSchemaDataType, XMLSchemaCommonType, XMLSchemaType
from pyrosetta.rosetta.protocols.moves import xsd_type_definition_w_attributes
from bootcamp_app import fold_tree_from_dssp_string
from pyrosetta.rosetta.protocols.moves import Mover
from FoldTreeFromSS import FoldTreeFromSS
from dataclasses import dataclass
from typing import List
from pyrosetta.rosetta.std import list_utility_tag_XMLSchemaAttribute_t
from pyrosetta.rosetta.core.scoring import parse_score_function


class BootCampMover(Mover):
    _clones = list()
    def __init__(self, sfxn  = None, num_iterations: int = 10):
        super().__init__() # or rosetta.protocols.moves.Mover.__init__(self)
        if sfxn is None:
            from pyrosetta.rosetta.core.scoring import get_score_function
            sfxn = get_score_function()
        self._sfxn = sfxn
        self._num_iterations = num_iterations 
    
    def get_score_function(self):
        return self._sfxn

    def set_score_function(self,new_sfxn):
        self._sfxn = new_sfxn

    def get_num_iterations(self):
        return self._num_iterations
    
    def set_num_iterations(self,new_value):
        self._num_iterations = new_value

    def apply(self,pose):
        ss = Dssp(pose).get_dssp_secstruct()

        ftfss = FoldTreeFromSS(pose)
        ft = ftfss.fold_tree_from_ss(pose)
        pose.fold_tree(ft)

        ft = fold_tree_from_dssp_string(ss)  
        pose.fold_tree(ft)

        correctly_add_cutpoint_variants(pose)

        sfxn = get_score_function()
        sfxn.set_weight(ScoreType.linear_chainbreak, 1.0)

        score = sfxn(pose)
        print("Score after fold tree + cutpoints:", score)

        # Set up task factory 
        tf = TaskFactory() 
        task = tf.create_task_and_apply_taskoperations(pose)
        task.restrict_to_repacking()

        # Set up MoveMap 
        movemap = kinematics.MoveMap()
        movemap.set_bb(True)
        movemap.set_chi(True)

        # Set up Minimizer
        min_opts = optimization.MinimizerOptions("lbfgs_armijo_atol", 0.01, True)
        minimizer = optimization.AtomTreeMinimizer()

        temp = 1.0
        monte_carlo_initial = MonteCarlo(pose,sfxn,temp)

        # PyMol Observer 
        the_observer = AddPyMOLObserver(pose)
        the_observer.pymol().apply(pose)

        # Monte Carlo Loop 
        counter_True = 0
        counter_False = 0
        total =  pose.total_residue()

        for i in range(5):
            pack_rotamers(pose,sfxn,task)
            minimizer.run(pose,movemap,sfxn,min_opts)

            randres = random.randrange(1, total)

            phi_pert = random.random()
            psi_pert = random.random()

            orig_phi = pose.phi(randres)
            orig_psi = pose.psi(randres)

            pose.set_phi(randres, orig_phi + phi_pert)
            pose.set_psi(randres, orig_psi + psi_pert)
            
            idx = ftfss.loop_for_residue(randres)
            if idx > 0:
                ranloop = ftfss.loop(idx)
                print(f"Closing loop: start={ranloop.start()} stop={ranloop.stop()} cut={ranloop.cut()}")
                ccd = pyrosetta.rosetta.protocols.loops.loop_closure.ccd.CCDLoopClosureMover(ranloop, movemap)
                ccd.apply(pose)
            accepted_mc = monte_carlo_initial.boltzmann(pose)
            print(f'The Monte Carlo Acceptance is {accepted_mc}')


            if accepted_mc == True:
                counter_True += 1
            elif accepted_mc == False:
                counter_False += 1

            if (i+1) % 5== 0:
                accep_rate = counter_True/5
                print(f'The acceptance rate is {accep_rate}')
                avg_energy = pose.energies().total_energy()/100
                print(f'The average evergy is{avg_energy}')
            else:
                continue 

            print(f"The lowest score is",monte_carlo_initial.lowest_score())
            monte_carlo_initial.lowest_score_pose().dump_pdb('Monte_Carlo.pdb')
                
    def get_name(self):
        self.__class__.__name__
        return "BootCampMover"

    @staticmethod
    def mover_name():
        return "BootCampMover"

    def parse_my_tag(self,tag,datamap):
        self.set_score_function(parse_score_function(tag, datamap))
        if tag.hasOption("num_iterations"):
            iters = tag.get_option_int("num_iterations", 1)    

    @staticmethod
    def provide_xml_schema(xsd):
        attrs = list_utility_tag_XMLSchemaAttribute_t()
        # ---- num_iterations attribute ----
        attrs.append(
            XMLSchemaAttribute.attribute_w_default(
                "num_iterations",
                XMLSchemaType(XMLSchemaCommonType.xsct_positive_integer),
                "Number of iteration",
                "10"
            )
        )

        # ---- ScoreFunction attribute helper ----
        attributes_for_parse_score_function_w_description(
            attrs,
            "ScoreFunction to use for sampling"
        )

        # Register this mover’s schema
        xsd_type_definition_w_attributes(
            xsd,
            BootCampMover.mover_name(),
            "BootCampMover constructs a FoldTree from DSSP and performs refinement using Monte Carlo.",
            attrs
        )

    def clone(self):
        copy = BootCampMover(self._sfxn, self._num_iterations)
        BootCampMover._clones.append(copy)
        return copy

    def fresh_instance(self):
        return BootCampMover()


