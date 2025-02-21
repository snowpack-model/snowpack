# SNOWPACK_forcing
Digest climate model output in NetCDF files into SNOWPACK atmospheric forcing files \*.smet.

## Current functionality
This workflow currently provides templates for MERRA-2, CESM, ERA5, COSMO and RACMO-2 to create atmospheric forcing for SNOWPACK, which are \*.smet files. The MERRRA-2, ERA5, COSMO and CESM examples demonstrate GRID_EXTRACT functionality, where the closest model grid point is extracted. The RACMO-2 example demonstrates GRID_SMART functionality, where data from the four model grid points cornering the location of interest is downscaled using interpolation techniques while taking into account the local topography defined by the digital elevation model.

## Usage
1. Define target sites: Update `station_list.lst` to include the latitude and longitude of the sites you would like SNOWPACK forcing files for. See `station_list.lst` for an example which would create .smet forcing files at 5 locations on the Antarctic ice sheet.
2. Modify the `io_files/<model>.ini` to have DEMFILE and SOURCE_DEM point to the DEM, and METEOPATH and GRID2DPATH to the paths containing NetCDF files with the parameters
3. Make sure the `meteoio_timeseries` executable is available in `./create_smet_from_netcdf/` folder (either make a copy, or a soft-link). Test if the executable executes (`./meteoio_timeseries`). Check for example if modules need to be loaded (and make sure they are loaded in the sbatch script too!).
4. Workflow setup: Run workflow setup to print bash commands, which will be executed at runtime, to the file `to_exec.lst`. Edit the settings in `setup.sh`, including atmospheric model, output temporal resolution, and starting and ending year. 
```
bash setup.sh
```
5. Runtime: Now, there are three ways to generate the \*.smet files:
    - Submit workflow as a Slurm or PBS job by executing:

        ```sbatch job.sbatch```

        Make sure to edit appropriate Slurm/PBS settings! **Particularly, count the number of lines *n* in `to_exec.lst`**, for example using `wc -l to_exec.lst`, **and modify the job array specification: `#SBATCH --array=1-n`**
    - Use the tool `parallel`. For example, when having *n* cores available, execute:

        ```parallel -j n < to_exec.lst```

    - When only one CPU available, run sequential using:

        ```bash to_exec.lst```

6. Postprocessing: Concatenate and zip individual yearly atmospheric forcing files into a continous time series located at `output/smet_forcing.zip`. `postprocess.sh` requires the atmospheric model as an argument, for example:
```
bash postprocess.sh MERRA-2
```
