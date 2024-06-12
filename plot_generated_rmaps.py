#!/usr/bin/python
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
def plot_rmap(fa_phi,fa_psi,pa_phi,pa_psi,phi,psi,pdb_id,chain_id,residue,resname):
    plt.figure(figsize=(7,7))
    if len(fa_phi)!=0:
        plt.scatter(fa_phi,fa_psi,s=5,alpha=1,c="mediumseagreen",color="mediumseagreen",zorder=1)
    if len(pa_phi)!=0:
        plt.scatter(pa_phi,pa_psi,s=5,alpha=1,c="antiquewhite",color="antiquewhite",zorder=1)
    plt.scatter(phi,psi,s=20,alpha=1,c='red',zorder=1)
    plt.hlines(0,-180,180)
    plt.vlines(0,-180,180)
    plt.xlabel('phi (in $\degree$)',fontsize=12)
    plt.ylabel('psi (in $\degree$)',fontsize=12)
    plt.xlim(-180, 180)
    plt.ylim(-180, 180)
    plt.savefig("results/"+pdb_id+"_"+chain_id+"_"+resname+str(residue)+"_rmap.png",format="png",dpi=200, bbox_inches="tight")
    plt.close()



