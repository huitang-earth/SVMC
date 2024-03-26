# SVMC
[![FMI](docs/img/FMI_logo.jpg "Finnish Meteorological Institute")](https://en.ilmatieteenlaitos.fi/)
[![LUKE](docs/img/Luke_logo.jpg "Natural Resources Institute Finland")](https://www.luke.fi/en)

## Compiling:

- Go to `src` folder. Run the following command:
```
# For Puhti:  

module load gcc/11.3.0 openmpi/4.1.4 netcdf-fortran/4.5.4
gmake -f Makefile_puhti

# For containers (docker or openshift)
gmake -f Makefile_container

```
- The program **SVMC** will appear in the `src` folder.
- Need to adapt the [library path](https://github.com/huitang-earth/SVMC/blob/f177abfb135c5111fdfd56436410c0f973d79f2b/src/Makefile_puhti#L16) and [include path](https://github.com/huitang-earth/SVMC/blob/f177abfb135c5111fdfd56436410c0f973d79f2b/src/Makefile_puhti#L8) of netcdf library when needed. 
- Only **gfortran** is tested at the moment!

## Running the test (Qvidja, 2021)

`./SVMC`

- To set the control parameters (e.g., input/output folders, running periods, time steps), please check `readctrl_mod.f90` at [here](https://github.com/huitang-earth/SVMC/blob/9738c7a5dc576ad472d9b5fc4664d79b576b7891/src/readctrl_mod.f90#L60-L72)
- To set the parameters of p-hydro, please check `readvegpara_mod.f90` at [here](https://github.com/huitang-earth/SVMC/blob/9738c7a5dc576ad472d9b5fc4664d79b576b7891/src/readvegpara_mod.f90#L147-L161)
- To set the parameters related spafhy, please check `readsoilpara_mod.f90` at [here](https://github.com/huitang-earth/SVMC/blob/9738c7a5dc576ad472d9b5fc4664d79b576b7891/src/readsoilpara_mod.f90#L138-L155).
- The input and output subroutines with netcdf are kept in `io_mod.f90`. If you want to add more output data into the output netcdf file, please define the variable [here](https://github.com/huitang-earth/SVMC/blob/9738c7a5dc576ad472d9b5fc4664d79b576b7891/src/io_mod.f90#L73-L198) in `io_mod.f90` and call the writing variable subroutine in `SVMC.f90`[here](https://github.com/huitang-earth/SVMC/blob/9738c7a5dc576ad472d9b5fc4664d79b576b7891/src/SVMC.f90#L316-L328). 
- Input and Output are both kept in the `data` folder [here](https://github.com/huitang-earth/SVMC/tree/main/data)

## Additional notes:
- Namelist files are implemented in the code, but have not be tested yet. So, **you need to change parameters in the code and re-compile the model** when doing testing. 
- Some parameters are still kept in p-hydro or spafhy code at the moment. **The plan is to move all the parameters into parameter module**, which will be easy for parameter calibration in the future.  
