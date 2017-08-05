# project a dynamic property on the structure using the B-factor (tempfactor) field

import numpy as np
import MDAnalysis
import MDAnalysis.analysis.align

from MDAnalysis.tests.datafiles import PSF, DCD


u = MDAnalysis.Universe(PSF, DCD)
ref = MDAnalysis.Universe(PSF, DCD)  # copy of u

CORE_selection = "resid 1-29 or resid 60-121 or resid 160-214"
pdbtrj = "adk_distance_bfac.pdb"


# workaround for Issue #1359 https://github.com/MDAnalysis/mdanalysis/issues/1359
# which will be fixed in 0.17.0
u.add_TopologyAttr(MDAnalysis.core.topologyattrs.Tempfactors(np.zeros(len(u.atoms))))


with MDAnalysis.Writer(pdbtrj, multiframe=True, bonds=None, n_atoms=u.atoms.n_atoms) as PDB:
    # reference coordinates: set to first frame
    ref.trajectory[0]
    # iterate through our trajectory
    for ts in u.trajectory:
        # superimpose on the reference CORE (at t=0)
        rmsd = MDAnalysis.analysis.align.alignto(u.atoms, ref.atoms, select=CORE_selection)
        distances = np.sqrt(np.sum((u.atoms.positions - ref.atoms.positions)**2, axis=1))
        # project displacement on structure via bfactor ("tempfactor") field
        u.atoms.tempfactors = distances
        PDB.write(u.atoms)
        print("Frame {0}: CORE RMSD before/after superposition: {1[0]:.1f} / {1[1]:.1f} A. "
              "min-max displacement: {2:.1f}...{3:.1f} A".format(ts.frame, rmsd, distances.min(), distances.max()))

print("Wrote PDB trajectory {0} with distances in bfactor field".format(pdbtrj))
