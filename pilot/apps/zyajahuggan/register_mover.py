from pyrosetta.rosetta import protocols 
from bootcamp_mover import BootCampMover

class BootCamMoverCreator(protocols.moves.MoverCreator):
    _instances = list()

    def __init__(self):
        protocols.moves.MoverCreator.__init__(self)

    def create_mover(self):
        mover = BootCampMover()
        self._instances.append(mover)
        return mover 

    def keyname(self):
        return BootCampMover.mover_name()

    def provide_xml_schema(self, xsd):
        print("creator provide_xml_schema is called")
        BootCampMover.provide_xml_schema(xsd)