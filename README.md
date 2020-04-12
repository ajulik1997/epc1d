# epc1d - Electrostatic PIC code in a 1D cyclic domain

One-dimensional particle-in-cell code written in Cython for simulating Landau damping and the two-stream instability.

<img src="https://github.com/ajulik1997/epc1d/blob/master/example.gif" alt="Two-stream instability (full view)">

## Quick-start

Compile Cython code with: `./compile.sh`

Generate default two-stream instability data: `./epc1d --twostream`

Render to video: `python render.py tmp.npz` (ensure `ffmpeg` is available)

## Arguments

|    Argument   	| Abbr. 	|  Type 	|                  Default                  	|       Units       	|               Required?               	|
|:-------------:	|:-----:	|:-----:	|:-----------------------------------------:	|:-----------------:	|:-------------------------------------:	|
|   `--landau`  	|  `-l` 	|  bool 	|                   False                   	|                   	|       Yes, if not `--twostream`       	|
| `--twostream` 	|  `-t` 	|  bool 	|                   False                   	|                   	|         Yes, if not `--landau`        	|
|   `--render`  	|  `-r` 	|  bool 	|                   False                   	|                   	|                Optional               	|
|   `--length`  	|  `-L` 	| float 	|  4π if `--landau`<br>100 if `--twostream` 	|   Debye lengths   	|                Optional               	|
|   `--cells`   	|  `-c` 	|  int  	|                     20                    	|       cells       	|                Optional               	|
| `--particles` 	|  `-p` 	|  int  	| 1e3 if `--landau`<br>1e4 if `--twostream` 	|     particles     	|                Optional               	|
|   `--alpha`   	|  `-a` 	| float 	|                    0.2                    	|      unitless     	|   Optional,<br>used with `--landau`   	|
|   `--newton`  	|  `-N` 	|  int  	|                     10                    	| Newton iterations 	|   Optional,<br>used with `--landau`   	|
|   `--vbeam`   	|  `-v` 	| float 	|                     2                     	|   thermal speeds  	| Optional,<br> used with `--twostream` 	|
|    `--time`   	|  `-T` 	| float 	|                     20                    	|  normalised time  	|                Optional               	|
|   `--steps`   	|  `-s` 	|  int  	|                     50                    	|     time steps    	|                Optional               	|
|    `--name`   	|  `-n` 	|  str  	|                    tmp                    	|                   	|                Optional               	|

## Examples

- 22500 particles, 10000 steps, full view: [link](https://mega.nz/file/t1pD3JbQ#NHD7XyDFssapQbd9zTxML9my6iaDc2K5EctsXqsPpl4)
- 500000 particles, 1000 steps, full view: [link](https://mega.nz/file/5sx1Qb4A#B-bHV7i5HWAP1s-aSH9-cidjtT1f4snYspV2qPvX1Po)
- 1000000 particles, 1400 steps, part-view (shown below): [link](https://mega.nz/file/I0xD2b7C#E5DfHJJ58yqeyOLJSvPMc8shwWS30MUZsMDDsaa5O8M)

<img src="https://github.com/ajulik1997/epc1d/blob/master/example2.gif" alt="Two-stream instability (pv view)">
