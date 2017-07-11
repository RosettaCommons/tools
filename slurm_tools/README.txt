README: Submit Script Generator
This script was written in 2016 for LK peptide studies in Jeff Gray's lab, by 
Joseph Lubin. The structure is posted here because it might be useful to adapt 
for other projects on SLURM systems that require the creation of multiple 
different simulation conditions simultaneously.

Many of the functions are specific to my project, and you might strip them out
for your own purposes, and pare it down to options for writing flags files and 
sbatch scripts. For my own purposes, I varied the starting PDB file I used, 
the score functions (and here, I created talaris variants with modified 
weights), and specific atom properties. The script could distinguish between 
inputs of score terms with varied weight and atom properties automatically.


More widely useful, the function creates a set of subfolders, one for each 
simulation condition, and in each creates a flags file, an sbatch script, 
and folders for the resulting models and cluster outout. I included an option 
to adjust the number of decoys requested. Another option allows for the 
submission of all created sbatch files to the cluster for calculation. The 
calculations for wall time were specific to my simulations, but can easily be 
adapted to other simulations.