#!/usr/bin/env cython
# -*- coding: utf-8 -*-
# cython: infer_types=True
# cython: boundscheck=False
# cython: wraparound=False
# cython: nonecheck=False

############################ TOP LEVEL ERROR HANDLER ###########################

cdef void error_handler(unsigned short errno,
                        str info):
    """Handles error by printing relevant message and exiting
    
    Parameters:
    ----------
    errno     : unsigned short
        Error number
    info      : str
        Extra information regarding error
    
    Returns:
    ----------
    void"""
    
    cdef dict errors = {1: "Failed to import module",
                        2: "Failed to parse arguments"}
    print("ERRNO {} ({}): {}".format(errno, errors[errno], info))
    from sys import exit; exit(errno)

try:
################################ PYTHON IMPORTS ################################
    from argparse import ArgumentParser, ArgumentTypeError, ArgumentError
    from tqdm import tqdm
    from numpy.fft import fft, ifft
    import numpy as np

################################ CYTHON IMPORTS ################################
    from libc.math cimport sqrt, M_PI
    from cpython cimport bool
    cimport numpy as cnp
    
except ImportError as e:
    error_handler(1, str(e))

########################## GLOBAL AUXILIARY FUNCTIONS ##########################

cdef void rk4step(float dt,
                  float length,
                  unsigned int cells,
                  unsigned int particles):
    """Takes a single step using the RK4 method
    
    Parameters:
    ----------
    dt          : float
        Single time step
    length      : float
        Length of system
    cells       : unsigned int
        Number of cells
    particles   : unsigned int
        Number of particles
        
    Returns:
    ----------
    void"""
    
    temp["rk4"][0] = pic(temp["concat"],                         length, cells, particles)
    temp["rk4"][1] = pic(temp["concat"] + 0.5*dt*temp["rk4"][0], length, cells, particles)
    temp["rk4"][2] = pic(temp["concat"] + 0.5*dt*temp["rk4"][1], length, cells, particles)
    temp["rk4"][3] = pic(temp["concat"] + dt*temp["rk4"][2],     length, cells, particles)
    temp["concat"] += (temp["rk4"][0] + 2*temp["rk4"][1] + 2*temp["rk4"][2] + temp["rk4"][3])*dt / 6
    
cdef void calc_density(cnp.ndarray[cnp.double_t, ndim=1] position,
                       float length,
                       unsigned int cells,
                       unsigned int particles,
                       float dx):
    """Calculates charge density given particle positions
    
    Parameters:
    ----------
    positions   : np.ndarray
        Array of positions
    length      : float
        Length of system
    cells       : unsigned int
        Number of cells
    particles   : unsigned int
        Number of particles
    dx          : float
        Cell spacing
        
    Returns:
    ----------
    void"""
    
    temp["frac_pos"] = position / dx                     # positions wrt cells
    temp["int_pos"] = temp["frac_pos"].astype(int)       # cells to left of particles
    temp["offset"] = temp["frac_pos"] - temp["int_pos"]  # offset from left cells
        
    temp["density"].fill(0)                              # reset array to all-zero  
    # add each particle to left and right of cell to density count
    np.add.at(temp["density"], temp["int_pos"], 1 - temp["offset"])
    np.add.at(temp["density"], (temp["int_pos"] + 1) % cells, temp["offset"])
    temp["density"] *= cells/particles                   # normalise to 1

cdef cnp.ndarray[cnp.double_t, ndim=1] periodic_interp(cnp.ndarray[cnp.double_t, ndim=1] interp, 
                                                       cnp.ndarray[cnp.double_t, ndim=1] idx,
                                                       unsigned int cells):
    """Carries out linear interpolation of a periodic array interp at index idx
    
    Parameters:
    ----------
    interp      : np.ndarray
        Periodic array to be interpolated
    idx         : np.ndarray
        Array of indexes at which ti interpolate
    cells       : unsigned int
        Number of cells
        
    Returns:
    ----------
    y[x]        : np.ndarray
        Interpolated y with non-integer x"""
        
    temp["interp_xl"] = np.floor(idx).astype(int)   # Left index
    temp["interp_dx"] = idx - temp["interp_xl"]     # Offset
    # ensure values of array lie between 0 and cells-1 inclusive
    temp["interp_xl"] = ((temp["interp_xl"] % cells) + cells) % cells
    
    return (interp[temp["interp_xl"]] * (1. - temp["interp_dx"]) 
            + interp[(temp["interp_xl"] + 1) % cells] * temp["interp_dx"])

cdef cnp.ndarray[cnp.double_t, ndim=1] fft_integrate(cnp.ndarray[cnp.double_t, ndim=1] y,
                                                     unsigned int n):
    """Integrates a period function using FFTs
    
    Parameters:
    ----------
    y           : np.ndarray
        Points of periodic function
    n           : unsigned int
        Number of points
        
    Returns:
    ----------
    arr         : np.ndarray
        FFT integrated periodic function"""

    temp["fft"][0] = fft(y) # Take FFT
    # Result is in standard layout with positive frequencies first then negative
    # n even: [ f(0), f(1), ... f(n/2), f(1-n/2) ... f(-1) ]
    # n odd:  [ f(0), f(1), ... f((n-1)/2), f(-(n-1)/2) ... f(-1) ]
    
    if n % 2 == 0:      # even number of points
        temp["fft"][1] = np.concatenate((np.arange(0, n/2+1), np.arange(1-n/2, 0)))
    else:               # odd number of points
        temp["fft"][1] = np.concatenate((np.arange(0, (n-1)/2+1), np.arange( -(n-1)/2, 0)))
    temp["fft"][1] = (2 * M_PI * temp["fft"][1]) / n
    
    # modify frequencies by dividing by ik
    temp["fft"][0][1:] /= (1j * temp["fft"][1][1:]) 
    temp["fft"][0][0] = 0 # set the arbitrary zero-frequency term to zero
    
    return ifft(temp["fft"][0]).real # return Reverse Fourier Transform

cdef cnp.ndarray[cnp.double_t, ndim=1] pic(cnp.ndarray[cnp.double_t, ndim=1] posvel,
                                           float length,
                                           unsigned int cells,
                                           unsigned int particles):
    """Particle-in-cell (PIC) function
    
    Parameters:
    ----------
    posvel   : np.ndarray
        Array of positions and velocities
    length      : float
        Length of system
    cells       : unsigned int
        Number of cells
    particles   : unsigned int
        Number of particles
        
    Returns:
    ----------
    arr         : np.ndarray
        Array of updated positions and velocities"""

    cdef float dx = length / cells      # cell spacing

    temp["pic"] = posvel[:particles]    # recover position of each particle
    temp["pic"] = ((temp["pic"] % length) + length) % length    # ensure pos is between 0 and length    
    calc_density(temp["pic"], length, cells, particles, dx)     # update number density
    
    return np.concatenate((posvel[particles:],
                           -periodic_interp(-fft_integrate(temp["density"] - 1,
                                                           cells) * dx,
                                            temp["pic"] / dx, cells)))

############################# DATA SAVING FUNCTION #############################

cdef void save_data(float length,
                    unsigned int cells,
                    unsigned int particles,
                    str name):
    """Saves data dictionary to disk
    
    Parameters:
    ----------
    length      : float
        Length of system
    cells       : unsigned int
        Number of cells
    particles   : unsigned int
        Number of particles
    name        : str
        Name of file to be saved
        
    Returns:
    ----------
    void"""
    
    cdef Py_ssize_t idx                                     # index of current output time step
    cdef cnp.ndarray[cnp.double_t, ndim=1] post_t, vel_t    # positions and velocities at time step
    cdef cnp.ndarray[cnp.double_t, ndim=1] bins             # edges of histogram bins
    cdef float dx = length / cells                          # length of a cell
    cdef unsigned int sqrt_part = int(sqrt(particles))    # floor of square root of particle number
    
    if len(data) != 3:
        for idx, (pos_t, vel_t) in tqdm(enumerate(zip(data["positions"], data["velocities"])),
                                        desc="Saving...",
                                        total=len(data["output_times"])):
            calc_density(pos_t, length, cells, particles, dx)
            data["densities"][idx] = temp["density"]
            data["first_harmonic"][idx] = 2 * np.abs(fft(data["densities"][idx])[1]) / cells
            data["vhist"][idx], bins = np.histogram(vel_t, sqrt_part)
            data["vbins"][idx] = 0.5 * (bins[1:] + bins[:-1])
    np.savez("{}.npz".format(name), **data)

############################### MAIN RUN FUNCTION ##############################

cdef void run(float length,
              unsigned int cells,
              unsigned int particles):
    """Main function that runs the simulation
    
    Parameters:
    ----------
    length      : float
        Length of system
    cells       : unsigned int
        Number of cells
    particles   : unsigned int
        Number of particles
        
    Returns:
    ----------
    void"""
        
    cdef float hdx = 0.5 * (length / cells)  # half of distance covered by each cell
    cdef float dt = 0                        # single time step
    cdef float time = 0                      # current time
    
    cdef Py_ssize_t idx                      # index of current output time step
    cdef double tnext                        # next output time from output_times
    cdef bool stepping                       # whether to keep advancing or save data
    
    # define starting state
    temp["concat"] = np.concatenate((data["positions"][0], data["velocities"][0]))

    # begin simulation, saving state at every defined output time
    for idx, tnext in tqdm(enumerate(data["output_times"]),
                           desc="Simulating...",
                           total=len(data["output_times"])):
        stepping = True                      # advancing to next output time
        while stepping:         # maximum distance a particle can move is one cell
            dt = hdx / np.max(np.abs(data["velocities"][idx if idx==0 else idx-1]))
            if time + dt >= tnext:           # next step will reach or exceed tnext
                stepping = False             # in that case, stop stepping
                dt = tnext - time            # update dt in reference to current time
            # calculate next position and velocity using rk4 function
            rk4step(dt, length, cells, particles)
            time += dt                       # increment time
            
        # extract position and velocities to be saved
        data["positions"][idx] = ((temp["concat"][0:particles] % length) + length) % length
        data["velocities"][idx] = temp["concat"][particles:]

###################### INITIAL CONDITION SETUP FUNCTIONS #######################

cdef void landau(unsigned int n_particles,
                 float length,
                 float alpha,
                 unsigned short newton):
    """Creates initial conditions for Landau damping
    
    Parameters:
    ----------
    n_particles : unsigned int
        Number of particles
    length      : float
        Length of system
    alpha       : float
        Perturbation amplitude
    newton      : unsigned short
        Number of Newton iterations
    
    Returns:
    ----------
    void
    """
    
    cdef float k = (2*M_PI) / length
    cdef Py_ssize_t _
    
    # start with uniform distribution of positions
    data["positions"][0] = np.random.uniform(0, length, n_particles)
    data["positions"][1] = data["positions"][0].copy()
    
    # adjust distribution using Newton iteration
    for _ in range(newton):
        data["positions"][0] -= ((data["positions"][0] - data["positions"][1] 
                                  + (alpha / k) * np.sin(k * data["positions"][0])) 
                                 / (1 + (alpha * np.cos(k * data["positions"][0]))))
    
    # normal velocity distribution
    data["velocities"][0] = np.random.normal(0, 1, n_particles)

cdef void twostream(unsigned int n_particles,
                    float length,
                    float v_beam):
    """Creates initial conditions for a two-stream instability
    
    Parameters:
    ----------
    n_particles : unsigned int
        Number of particles
    length      : float
        Length of system
    v_beam      : float
        Beam velocity
        
    Returns:
    ----------
    void"""
    
    cdef int half = int(n_particles / 2)
    
    # Start with a uniform distribution of positions and normal velocity distribution
    data["positions"][0] = np.random.uniform(0, length, n_particles)
    data["velocities"][0] = np.random.normal(0, 1, n_particles)
    
    data["velocities"][0][:half] += v_beam  # half of particles moving one way
    data["velocities"][0][half:] -= v_beam  # other half moving the opposite way

######################## MEMORY PRE-ALLOCATION ROUTINES ########################

cdef dict allocate_temp(args):
    """Allocates memory for all temporary arrays
    
    Parameters:
    ----------
    args      : argparse.Namespace
    
    Returns:
    ----------
    temp      : dict
        Dictionary of temporary named arrays for intermediate calculations"""

    cdef cnp.ndarray[cnp.double_t,   ndim=1] concat, density, frac_pos, offset, pic, interp_dx
    cdef cnp.ndarray[cnp.double_t,   ndim=2] rk4
    cdef cnp.ndarray[double complex, ndim=2] fft
    cdef cnp.ndarray[cnp.uint_t,     ndim=1] int_pos, interp_xl
        
    concat    = np.empty((args.particles*2),    dtype=np.double)   # concatenated arrays
    density   = np.empty((args.cells),          dtype=np.double)   # used in calc_density
    frac_pos  = np.empty((args.particles),      dtype=np.double)   # fractional positions
    offset    = np.empty((args.particles),      dtype=np.double)   # position offsets
    int_pos   = np.empty((args.particles),      dtype=np.uint)     # floor frac positions
    rk4       = np.empty((4, args.particles*2), dtype=np.double)   # used in RK4 routine
    pic       = np.empty((args.particles),      dtype=np.double)   # used in PIC routine
    fft       = np.empty((2, args.cells),       dtype=np.cdouble)  # storing FFT data
    interp_xl = np.empty((args.particles),      dtype=np.uint)     # used in periodic_interp
    interp_dx = np.empty((args.particles),      dtype=np.double)   # used in periodic_interp
    
    return {"concat"   : concat,
            "density"  : density,
            "frac_pos" : frac_pos,
            "offset"   : offset,
            "int_pos"  : int_pos,
            "rk4"      : rk4,
            "pic"      : pic,
            "fft"      : fft,
            "interp_xl": interp_xl,
            "interp_dx": interp_dx}

cdef dict allocate_data(args):
    """Allocates memory for all arrays to be saved
    
    Parameters:
    ----------
    args      : argparse.Namespace
    
    Returns:
    ----------
    data      : dict
        Dictionary of named arrays to be saved"""
    
    cdef int sqrt_npart = int(sqrt(args.particles))
    
    cdef cnp.ndarray[cnp.double_t, ndim=1] output_times, first_harmonic, cells
    cdef cnp.ndarray[cnp.double_t, ndim=2] positions, velocities, densities, vbins
    cdef cnp.ndarray[cnp.uint_t,   ndim=2] vhist
    
    output_times   = np.linspace(0, args.time, args.steps,    dtype=np.double)
    positions      = np.empty((args.steps, args.particles),   dtype=np.double)
    velocities     = np.empty((args.steps, args.particles),   dtype=np.double)
    
    if args.render:
        return {"output_times"   : output_times,
                "positions"      : positions,
                "velocities"     : velocities}
    
    else:
        first_harmonic = np.empty((args.steps),                   dtype=np.double)
        densities      = np.empty((args.steps, args.cells),       dtype=np.double)
        cells          = np.linspace(0, args.length, args.cells,  dtype=np.double)
        vhist          = np.empty((args.steps, sqrt_npart),       dtype=np.uint)
        vbins          = np.empty((args.steps, sqrt_npart),       dtype=np.double)
        
        return {"output_times"   : output_times,
                "first_harmonic" : first_harmonic,
                "positions"      : positions,
                "velocities"     : velocities,
                "densities"      : densities,
                "cells"          : cells,
                "vhist"          : vhist,
                "vbins"          : vbins}

############################### ARGUMENT PARSER ################################

cdef setup_argparse():
    """Sets up argument parser
    
    Parameters:
    ----------
    None
    
    Returns:
    ----------
    parser    : argparse.ArgumentParser
        Argument parser"""
        
    parser = ArgumentParser(description="epc1d: Electrostatic PIC code in a 1D cyclic domain")
    
    parser.add_argument('-l', '--landau',    action="store_true",       help="Generate initial conditions for Landau damping")
    parser.add_argument('-t', '--twostream', action="store_true",       help="Generate initial conditions for a two-stream instability")
    parser.add_argument('-r', '--render',    action="store_true",       help="If simulating for rendering, only positions and velocities are saved")
    
    parser.add_argument('-L', '--length',    type=float, default=0,     help="Length of the domain")
    parser.add_argument('-c', '--cells',     type=int,   default=0,     help="Number of cells")
    parser.add_argument('-p', '--particles', type=int,   default=0,     help="Number of particles")
    parser.add_argument('-a', '--alpha',     type=float, default=0.2,   help="Landau damping perturbation amplitude")
    parser.add_argument('-N', '--newton',    type=int,   default=10,    help="Number of Newton iterations with which to adjust position distribution in Landau damping")
    parser.add_argument('-v', '--vbeam',     type=float, default=2,     help="Two-stream instability beam velocity")
    
    parser.add_argument('-T', '--time',      type=float, default=20,    help="Simulation time")
    parser.add_argument('-s', '--steps',     type=int,   default=50,    help="Number of time steps at which to save data")
    parser.add_argument('-n', '--name',      type=str,   default="tmp", help="Name of saved array (will overwrite existing file)")
    
    return parser

cdef parse_args(parser):
    """Parses arguments and verifies them for validity
    
    Parameters:
    ----------
    parser    : argparse.ArgumentParser
         Argument parser
    
    Returns:
    ----------
    args      : argparse.Namespace
         Namespace object containing parsed arguments"""

    args = parser.parse_args()
        
    # check for initial conditions
    if not (args.landau or args.twostream): error_handler(2, "No initial conditions selected")
    if (args.landau and args.twostream): error_handler(2, "Two initial conditions selected")
    
    # set defaults for unset parameters in case of landau damping
    if args.landau:
        if args.length == 0: args.length = 4*M_PI
        if args.cells == 0: args.cells = 20
        if args.particles == 0: args.particles = 1000
    
    # set defaults for unset parameters in case of two-stream instability
    if args.twostream:
        if args.length == 0: args.length = 100
        if args.cells == 0: args.cells = 20
        if args.particles == 0: args.particles = 10000
    
    # finally, check all arguments are positive and nonzero
    for arg in ["length", "cells", "particles", "newton", "time", "steps"]:
        if getattr(args, arg) <= 0:
            error_handler(2, "Argument {} must be positive and nonzero".format(arg))
    
    return args

##################################### MAIN #####################################

if __name__ == '__main__':
    parser = setup_argparse()
    args = parse_args(parser)
    
    # set up global dictionary of temporary and data arrays
    global temp, data
    temp = allocate_temp(args)
    data = allocate_data(args)
    
    # set up initial conditions for landau damping
    if args.landau:
        landau(args.particles, args.length, args.alpha, args.newton)
    
    # set up initial conditions for two-stream instability
    if args.twostream:
        twostream(args.particles, args.length, args.vbeam)
        
    # run simulation and save data
    run(args.length, args.cells, args.particles)
    save_data(args.length, args.cells, args.particles, args.name)
