from pyrosetta import *
from pyrosetta.rosetta.protocols.moves import MoverRegistrator
from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects

from bootcamp_mover import BootCampMover 

EMBEDDED_XML = """
<ROSETTASCRIPTS>
  <SCOREFXNS>
    <ScoreFunction name="sfxn" weights="ref2015"/>
  </SCOREFXNS>
  <MOVERS>
    <BootCampMover name="bcm" scorefxn="sfxn" num_iterations="100"/>
  </MOVERS>
  <PROTOCOLS>
    <Add mover_name="bcm"/>
  </PROTOCOLS>
</ROSETTASCRIPTS>
"""

def main():
    MoverRegistrator.register_mover("BootCampMover", BootCampMover)
    init()

    pose = pose_from_pdb("input.pdb")

    xmlobj = XmlObjects.create_from_string(EMBEDDED_XML)
    protocol = xmlobj.get_mover("ParsedProtocol")

    protocol.apply(pose)


    pose.dump_pdb("refined_output.pdb")


if __name__ == "__main__":
    main()