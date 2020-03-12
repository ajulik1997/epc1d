import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.animation as animation
from matplotlib import cm
from tqdm import tqdm
from sys import argv, exit

def draw_all(num, data, phase_plot, density_plot, vel_plot, progress_bar):
    phase_plot.set_offsets(np.c_[data["positions"][num], data["velocities"][num]])
    density_plot.set_data(data["cells"], data["densities"][num])
    vel_plot.set_data(data["vhist"][num], data["vbins"][num])
    progress_bar.update()
    return phase_plot, density_plot, vel_plot

def draw_pv(num, data, phase_plot, progress_bar):
    phase_plot.set_offsets(np.c_[data["positions"][num], data["velocities"][num]])
    progress_bar.update()
    return phase_plot,

def harmonic_plot(data):
    plt.figure()
    plt.plot(data["output_times"], data["first_harmonic"], "k-")
    plt.xlabel("Time [Normalised]")
    plt.ylabel("First harmonic amplitude [Normalised]")
    plt.yscale('log')
    plt.tight_layout()
    plt.show()
    
def render_pv(data):
    colourmap = cm.get_cmap("plasma", len(data["positions"][0]))
    fig = plt.figure(figsize=[16,9])
    ax = fig.add_axes([0, 0, 1, 1]); ax.axis('off')
    phase_plot = ax.scatter(data["positions"][0], data["velocities"][0], s=32, c=colourmap.colors, marker=".", alpha=0.5)
    ax.set_xlim((np.min(data["positions"]), np.max(data["positions"])))
    ax.set_ylim((np.min(data["velocities"]), np.max(data["velocities"])))
    
    progress_bar = tqdm(total=len(data["output_times"])+1)
    ani = animation.FuncAnimation(fig, draw_pv, len(data["output_times"]), fargs=(data, phase_plot, progress_bar), interval=50, blit=True)
    ani.save("part={};steps={}.mp4".format(len(data["positions"][0]), len(data["output_times"])), fps=30, dpi=96, savefig_kwargs={'facecolor':'black'})
    
def render_all(data):
    colourmap = cm.get_cmap("plasma", len(data["positions"][0]))
    fig = plt.figure(figsize=[16,9])
    gs = gridspec.GridSpec(4, 4)
    ax1 = fig.add_subplot(gs[0:3,0:3])
    phase_plot = ax1.scatter(data["positions"][0], data["velocities"][0], s=32, c=colourmap.colors, marker=".", alpha=0.5)
    ax1.set_title("Phase space")
    ax1.set_xlim((np.min(data["positions"]), np.max(data["positions"])))
    ax1.set_ylim((np.min(data["velocities"]), np.max(data["velocities"])))

    ax2 = fig.add_subplot(gs[3,0:3])
    density_plot = ax2.plot([], [], 'k-')[0]
    ax2.set_xlim((np.min(data["cells"]), np.max(data["cells"])))
    ax2.set_ylim((np.min(data["densities"]), np.max(data["densities"])))
            
    ax3 = fig.add_subplot(gs[0:3,3])
    vel_plot = ax3.plot([], [], 'k-')[0]
    ax3.set_xlim((np.min(data["vhist"]), np.max(data["vhist"])))
    ax3.set_ylim((np.min(data["vbins"]), np.max(data["vbins"])))
    
    plt.tight_layout()
    
    progress_bar = tqdm(total=len(data["output_times"])+1)
    ani = animation.FuncAnimation(fig, draw_all, len(data["output_times"]), fargs=(data, phase_plot, density_plot, vel_plot, progress_bar), interval=50, blit=True)
    ani.save("part={};steps={}.mp4".format(len(data["positions"][0]), len(data["output_times"])), fps=30, dpi=96)
    
if __name__ == "__main__":
    if len(argv) != 2: print("1 argument expected"); exit(1)
    with np.load(argv[-1], mmap_mode="r") as data:
        if "first_harmonic" in data:
            harmonic_plot(data)
        if len(data) == 3:
            render_pv(data)
        else:
            render_all(data)
