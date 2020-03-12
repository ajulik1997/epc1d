import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.animation as animation
from matplotlib import cm
from tqdm import tqdm

def ani_update(num, loaded, phase_plot, density_plot, vel_plot, progress_bar):
    phase_plot.set_offsets(np.c_[loaded["positions"][num], loaded["velocities"][num]])
    density_plot.set_data(loaded["cells"], loaded["densities"][num])
    vel_plot.set_data(loaded["vhist"][num], loaded["vbins"][num])
    progress_bar.update()
    
    return phase_plot, density_plot, vel_plot
    

with np.load("part=22500;steps=10k.npz", mmap_mode="r") as loaded:
    plt.figure()
    plt.plot(loaded["output_times"], loaded["first_harmonic"], "k-")
    plt.xlabel("Time [Normalised]")
    plt.ylabel("First harmonic amplitude [Normalised]")
    plt.yscale('log')
    plt.tight_layout()
    plt.show()
    
    #####
    
    colourmap = cm.get_cmap("plasma", len(loaded["positions"][0]))
    fig = plt.figure(figsize=[16,9])
    gs = gridspec.GridSpec(4, 4)
    ax1 = fig.add_subplot(gs[0:3,0:3])
    phase_plot = ax1.scatter(loaded["positions"][0], loaded["velocities"][0], s=32, c=colourmap.colors, marker=".", alpha=0.5)
    ax1.set_title("Phase space")
    ax1.set_xlim((np.min(loaded["positions"]), np.max(loaded["positions"])))
    ax1.set_ylim((np.min(loaded["velocities"]), np.max(loaded["velocities"])))

            
    ax2 = fig.add_subplot(gs[3,0:3])
    density_plot = ax2.plot([], [], 'k-')[0]
    ax2.set_xlim((np.min(loaded["cells"]), np.max(loaded["cells"])))
    ax2.set_ylim((np.min(loaded["densities"]), np.max(loaded["densities"])))
            
    ax3 = fig.add_subplot(gs[0:3,3])
    vel_plot = ax3.plot([], [], 'k-')[0]
    ax3.set_xlim((np.min(loaded["vhist"]), np.max(loaded["vhist"])))
    ax3.set_ylim((np.min(loaded["vbins"]), np.max(loaded["vbins"])))
    
    progress_bar = tqdm(total=len(loaded["output_times"]))
    ani = animation.FuncAnimation(fig, ani_update, len(loaded["output_times"]), fargs=(loaded, phase_plot, density_plot, vel_plot, progress_bar), interval=50, blit=True)
    
    plt.tight_layout()
    #plt.show()
    ani.save("part=22500;steps=10k.mp4", fps=30, dpi=96)
