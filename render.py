################################################################################
#																			   #
#				Code for rendering data created with epc1d					   #
#																			   #
#	Author:		Alexander Liptak											   #
#	Contact:	alexander.liptak.jr@googlemail.com							   #
#	Date:		March 2020													   #
#																			   #
################################################################################

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.animation as animation
from matplotlib import cm
from os import listdir, mkdir, path, rmdir
from PIL import Image
from shutil import rmtree
from subprocess import call, DEVNULL
from sys import argv, exit
from tqdm import tqdm

################################################################################

size = 4	 # size of each particle
alpha = 0.1	 # opacity of each particle
split = 2	 # split frame into chunks (higher is slower but uses less memory)

################################################################################

def interrupted(name):
	"""Attempts to check whether render was interrupted and can be continued
	
	Arguments:
	----------
	name	 : str
		Name of render job
		
	Returns:
	----------
	index	 : int
		Index at which to resume render"""
		
	if not path.exists(name): return 0	# no directory exists, start from 0
	try:  # find highest number in folder, start from one higher than that
		return max([int(f.split(".")[0]) for f in listdir(name)]) + 1
	except:				# the above will fail if folder is empty
		rmdir(name)		# remove empty folder
		return 0		# start from the beginning

def png_write(fig, fname):
	"""Custom png writer that doesn't call draw() like matplotlib
	
	Arguments:
	----------
	fig      : matplotlib.figure.Figure
		Figure to be saved
	fname    : str
		Filename of image to be saved, including file extenstion
	
	Returns:
	----------
	None"""
	
	# get pixel width and height of figure from size and DPI
	width, height = map(int, fig.get_size_inches() * fig.get_dpi())
	
	image = np.frombuffer(fig.canvas.tostring_argb(),  # ARGB str from canvas
	                      dtype='uint8')  # convert binary to uint8 array
	image = image.reshape(height, width, 4)  # shape to image with ARGB per px
	image = np.roll(image, -1, 2)         # convert ARGB to RGBA format
	Image.fromarray(image, 'RGBA').save(fname)  # create PIL image and save
	
def gen_pv_frames(data, colourmap, name, index, padding):
	"""Generates and saves postions-velocity frames to temporary folder
	
	Arguments:
	----------
	data      : dict
		Named arrays loaded from file
	colourmap : list
		Colours from colourmap with one colour for each particle
	name	  : str
		Name of render job
	index     : int
		Frame index at which to start drawing
	padding   : int
		Zero-padded width of each saved frame
		
	Returns:
	----------
	None"""
	
	# generate array of indices in each chunk to be drawn per frame
	print(f"Generating {split} chunks...")
	chunks = np.split(np.arange(data["positions"].shape[1], dtype=int), split)
	
	print("Calculating limits...")
	xlim = (np.min(data["positions"]), np.max(data["positions"]))
	ylim = (np.min(data["velocities"]), np.max(data["velocities"]))
	
	# render and save each frame, drawing each chunk per frame
	for frame in tqdm(range(index, data["output_times"].shape[0]),
	                  desc="Drawing frame:"):
					  
		# set up figure and axis
		fig = plt.figure(figsize=[16,9], facecolor="k")
		ax = fig.add_axes([0, 0, 1, 1]); ax.axis('off')
		ax.set_xlim(xlim); ax.set_ylim(ylim)
		
		# create "empty" arist for drawing on canvas
		phase_plot = ax.scatter([0], [0],     # dummy point
		                        s=size,       # set size from global variable
								c=["#00000000"],  # dummy transparent colour
								marker='.',   # draw spheres
								alpha=alpha)  # set particle transparency
		fig.canvas.draw()  # draw canvas before any real artist draw events
		
		# process particles in each chunk, setting their position and colour
		for chunk in tqdm(chunks, desc="Processing chunk:", leave=False):
			phase_plot.set_offsets(np.c_[data["positions"][frame][chunk],
			                             data["velocities"][frame][chunk]])
			phase_plot.set_color(colourmap[chunk])
			ax.draw_artist(phase_plot)	# draw computed particles to image
		
		# save generated image for stitching
		png_write(fig, f"{name}/{frame:0{padding}}.png")

def gen_all_frames(data, colourmap, name, index, padding):
	"""Generates and saves frames with all computed data to temporary folder
	(assumes that as there is more data to be drawn, there is not enough
	particles to warrant chunking; this makes things slightly faster)
	
	Arguments:
	----------
	data      : dict
		Named arrays loaded from file
	colourmap : list
		Colours from colourmap with one colour for each particle
	name	  : str
		Name of render job
	index     : int
		Frame index at which to start drawing
	padding   : int
		Zero-padded width of each saved frame
		
	Returns:
	----------
	None"""
	
	print("Calculating limits...")
	main_xlim = (np.min(data["positions"]), np.max(data["positions"]))
	main_ylim = (np.min(data["velocities"]), np.max(data["velocities"]))
	dens_xlim = (np.min(data["cells"]), np.max(data["cells"]))
	dens_ylim = (np.min(data["densities"]), np.max(data["densities"]))
	hist_xlim = (np.min(data["vhist"]), np.max(data["vhist"]))
	hist_ylim = (np.min(data["vbins"]), np.max(data["vbins"]))
	
	# render and save each frame, drawing each chunk per frame
	for frame in tqdm(range(index, data["output_times"].shape[0]),
	                  desc="Drawing frame:"):
					  
		# set up figure with gridspec
		fig = plt.figure(figsize=[16,9])
		gs = gridspec.GridSpec(4, 4)
		
		# set up axes
		ax1 = fig.add_subplot(gs[0:3,0:3]); ax1.set_title("Phase space")
		ax2 = fig.add_subplot(gs[3,0:3])
		ax3 = fig.add_subplot(gs[0:3,3])

		# set up axes limits
		ax1.set_xlim(main_xlim); ax1.set_ylim(main_ylim)
		ax2.set_xlim(dens_xlim); ax2.set_ylim(dens_xlim)
		ax3.set_xlim(hist_xlim); ax3.set_ylim(hist_xlim)

		# draw all artists
		ax1.scatter(data["positions"][frame], data["velocities"][frame],
		            s=size, c=colourmap, marker=".", alpha=alpha)
		ax2.plot(data["cells"], data["densities"][frame], 'k-')
		ax3.plot(data["vhist"][frame], data["vbins"][frame], 'k-')
		
		plt.tight_layout()  # apply tight layout formatting
		
		# save generated image for stitching
		png_write(fig, f"{name}/{frame:0{padding}}.png")

def harmonic_plot(data):
	"""Draws the preliminary first-harmonic plot to screen
	
	Arguments:
	----------
	data      : dict
		Named arrays loaded from file

	Returns:
	----------
	None"""
	
	plt.figure()
	plt.plot(data["output_times"], data["first_harmonic"], "k-")
	plt.xlabel("Time [Normalised]")
	plt.ylabel("First harmonic amplitude [Normalised]")
	plt.yscale('log')
	plt.tight_layout()
	plt.show()
	
if __name__ == "__main__":

	# check if argv has two items (one argument), quit otherwise
	if len(argv) != 2: print("1 argument expected"); exit(1)
	
	# check whether script has access to ffmpeg
	try:
		call("ffmpeg", stdout=DEVNULL, stderr=DEVNULL)
	except:
		print("ffmpeg not found"); exit(1)
			
	name = argv[-1].split(".")[0]	# extract filename from argument
	
	# open requested array and begin rendering
	with np.load(argv[-1], mmap_mode="r") as data:
	
		# data contains arrays for plotting first harmonic
		if "first_harmonic" in data:
			harmonic_plot(data)
			
		index = interrupted(name)  # get index at which to start drawing
		print(f"Starting from frame {index}")
		
		# calculate length required for consistent zero-padding
		padding = len(str(data["output_times"].shape[0]))
		
		if index == 0: mkdir(name)  # if starting from scratch, create wkdir
		
		# generate colourmap of correct length, extract colours and order them
		print("Generating colourmap...")
		colourmap = cm.get_cmap("plasma", data["velocities"].shape[1])
		colourmap = colourmap.colors[np.argsort(data["velocities"][0])]
		
		if len(data) == 3:	# only data for particle-velocity plot available
			gen_pv_frames(data, colourmap, name, index, padding)
		else:               # data for all plots is available
			gen_all_frames(data, colourmap, name, index, padding)
		
		call(["ffmpeg",		# call ffmpeg to stitch all images together
		      "-i", f"./{name}/%04d.png",	# input images
			  "-r", "30",                   # framerate
			  "-s", "1920x1080",            # resolution
			  "-vcodec", "libx264",         # video codec
			  "-crf", "15",                 # constant-rate factor
			  "-pix_fmt", "yuv420p",        # pixel format
			  f"{name}.mp4"])               # output video

		rmtree(name)		# clean up images after render
