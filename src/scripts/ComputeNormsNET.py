"""The WaveBlocks Project

Compute the norms of the different wavepackets or wavefunctions.
This script does not do an eigentransformation of the data.

@author: R. Bourquin
@copyright: Copyright (C) 2012 R. Bourquin
@license: Modified BSD License
"""

import sys

from WaveBlocksND import IOManager


if __name__ == "__main__":

    iom = IOManager()

    # Read file with simulation data
    try:
        iom.open_file(filename=sys.argv[1])
    except IndexError:
        iom.open_file()

    # Iterate over all blocks
    for blockid in iom.get_block_ids():
        print("Computing the norms in data block '"+str(blockid)+"'")

        if iom.has_norm(blockid=blockid):
            print("Datablock '"+str(blockid)+"' already contains norm data, silent skip.")
            continue

        # TODO: Add new algorithms here

        # We test for an inhomogeneous wavepacket next
        if iom.has_inhomogwavepacket(blockid=blockid):
            import NormWavepacket
            NormWavepacket.compute_norm_inhawp(iom, blockid=blockid, eigentrafo=False)
        # We test for a homogeneous wavepacket next
        elif iom.has_wavepacket(blockid=blockid):
            import NormWavepacket
            NormWavepacket.compute_norm_hawp(iom, blockid=blockid, eigentrafo=False)
        # We have no wavepacket, then we try for a wavefunction
        elif iom.has_wavefunction(blockid=blockid):
            import NormWavefunction
            NormWavefunction.compute_norm(iom, blockid=blockid, eigentrafo=False)
        # If there is also no wavefunction, then there is nothing to compute the norm
        else:
            print("Warning: Not computing any norm in block '"+str(blockid)+"'!")

    iom.finalize()
