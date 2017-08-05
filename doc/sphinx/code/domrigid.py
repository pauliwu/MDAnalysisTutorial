import numpy as np
from MDAnalysis.analysis.rms import rmsd


if __name__ == "__main__":
    import MDAnalysis
    import matplotlib
    import matplotlib.pyplot as plt

    # load AdK DIMS trajectory
    from MDAnalysis.tests.datafiles import PSF, DCD
    u = MDAnalysis.Universe(PSF, DCD)

    # one AtomGroup per domain
    domains = {
        'CORE': u.select_atoms("(resid 1-29 or resid 60-121 or resid 160-214) and name CA"),
        'LID': u.select_atoms("resid 122-159 and name CA"),
        'NMP': u.select_atoms("resid 30-59 and name CA"),
        }
    colors = {'CORE': 'black', 'NMP': 'blue', 'LID': 'red'}

    u.trajectory[0]   # rewind trajectory
    xref0 = dict((name, g.positions - g.center_of_mass()) for name, g in domains.items())

    nframes = len(u.trajectory)
    results = dict((name, np.zeros((nframes, 2), dtype=np.float64)) for name in domains)

    for iframe, ts in enumerate(u.trajectory):
        for name, g in domains.items():
            results[name][iframe, :] = (u.trajectory.time,
                                        rmsd(g.positions, xref0[name],
                                             center=True, superposition=True))


    # plot
    fig = plt.figure(figsize=(5,5))
    ax = fig.add_subplot(111)
    for name in "CORE", "NMP", "LID":
        data = results[name]
        ax.plot(data[:,0], data[:,1], linestyle="-", color=colors[name], lw=2, label=name)
    ax.legend(loc="best")
    ax.set_xlabel(r"time  $t$ (ps)")
    ax.set_ylabel(r"C$_\alpha$ RMSD from $t=0$, $\rho_{\mathrm{C}_\alpha}$ ($\AA$)")

    for ext in ('svg', 'pdf', 'png'):
        fig.savefig("AdK_domain_rigidity.{0}".format(ext))
