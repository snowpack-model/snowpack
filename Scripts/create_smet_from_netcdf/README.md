# SNOWPACK_forcing
Digest climate model output into SNOWPACK atmospheric forcing.

## Current functionality
This workflow currently digests MERRA-2 atmospheric reanalysis, CESM or RACMO-2 into suitable atmospheric forcing for SNOWPACK, which are .smet files. The MERRRA-2 and CESM example demonstrate GRID_EXTRACT functionality, where the closest model grid point is extracted. The RACMO-2 example demonstrates GRID_SMART functionality, where data from the four model grid points cornering the location of interest is downscaled using interpolation techniques while taking into account the local topography defined by the digital elevation model.

## Usage
1. Define target sites: Update `station_list.lst` to include the latitude and longitude of the sites you would like SNOWPACK forcing files for. See `station_list.lst` for an example which would create .smet forcing files at 5 locations on the Antarctic ice sheet. 
2. Workflow setup: Run workflow setup to print bash commands, which will be executed at runtime, to `to_exec.lst`. Edit the settings in `setup.sh`, including atmospheric model, output temporal resolution, and starting and ending year. 
```
bash setup.sh
```
3. Runtime: Submit workflow as Slurm job by executing `job.sbatch`. Make sure to edit appropriate Slurm settings! We will soon add an example using the PBS job scheduler, but for now we only support Slurm.
```
sbatch job.sbatch
```
4. Postprocessing: Concatenate and zip individual yearly atmospheric forcing files into a continous time series located at `output/smet_forcing.zip`. `postprocess.sh` requires the atmospheric model as an argument, for example:
```
bash postprocess.sh MERRA-2
```

## To do list
1. Update workflow to provide CESM support. 
