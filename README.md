# nprunstats

Measures depth and statistics from CP processed imaging, SExtractored source catalogs and PSFEx psf files.
So step zero is to create the source catalogs and psf files as you see fit. I recommend using the rapala
reduction pipeline https://github.com/legacysurvey/rapalaource.

These files should be saved in a common directory structure bass/CPYYYYMMDD, where YYYMMDD refers to the
observing night. Then set the environment variable CPCATS to point to bass/.

Set the environment variable BASSDATA to point to the CP processed NinetyPrime images. At NERSC these are
located in /project/projectdirs/cosmo/staging/bok/BOK_CP/

Set the environment variable PS1CHUNKS to point to PS1 catalogs. At NERSC these are located at
/global/project/projectdirs/cosmo/work/ps1/cats/chunks-qz-star-v2/

Customize the SLURM script runstats.sl to your specific needs.

Run ./runrunstats YYYYMMDD for each night you'd like process. Output, ouput reports and error reports will 
be written to the working directory. 
