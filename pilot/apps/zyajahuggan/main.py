from pyrosetta import *
from pyrosetta.rosetta.protocols.moves import MoverFactory
from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects

from bootcamp_mover import BootCampMover 
from register_mover import BootCampMoverCreator

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
    BootCampMoverCreator.register()
    init()
    pose = pose_from_pdb("1UBQ.pdb.gz")

    xmlobj = XmlObjects.create_from_string(EMBEDDED_XML)
    protocol = xmlobj.get_mover("ParsedProtocol")
    protocol.apply(pose)

    pose.dump_pdb("refined_output.pdb")
    print("✅ Finished! Output written to refined_output.pdb")



if __name__ == "__main__":
    main()