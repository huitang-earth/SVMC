#!/bin/bash

#SBATCH --time=02:00:00
#SBATCH --partition=fmi
#SBATCH --job-name=svm
#SBATCH --account=project_2006422
#SBATCH --mem-per-cpu=4000
#SBATCH --ntasks=1


######### Set "run_***" to "T" to run the site simulation
run_case="F"            # T or F, swtich for running long experiments
post_proc="T"         # T or F, swtich for archiving experiments. 



#my_variable=$(python my_script.py)


params1=(0.5 5 10)
params2=(-1 -2 -3)

expids=(km0.5_p1_c5 km0.5_p2_c5 km0.5_p3_c5 km5_p1_c5 km5_p2_c5 km5_p3_c5 km10_p1_c5 km10_p2_c5 km10_p3_c5)

for i in {0..8}; do


if [ ${run_case} == "T" ]
then

module load gcc/11.3.0 openmpi/4.1.4 netcdf-fortran/4.5.4

cd /users/tanghui1/SVMC/src

cat << EOF > ctrl_namelist
&ctrl_namelist
start_date_day=20180510
start_date_hour=000000
time_step=1
end_date_day=20230101
end_date_hour=000000
num_sites=1
sites_name='Qvidja'
input_dir='../data/input/'  
output_dir='../data/output/'     
time_step_output=1.0
obs_lai=.true.
obs_soilmoist=.false.
obs_snowdepth=.false.
obs_manage=.true.
yasso_year=.false.
phydro_debug=.false.
yasso_debug=.true.
water_debug=.true.
log_level=1
experiment_id='0.6lr0.3tl_invert0_${expids[i]}'
/
EOF

cat << EOF > veg_namelist
&veg_namelist
num_pft=1
pft_type="grass"
conductivity=5e-17
psi50=-1
b=2
alpha=0.08
gamma=1
opt_hypothesis='PM'
rdark=0.015
/
EOF

cat << EOF > alloc_namelist
&alloc_namelist
 cratio_resp=5.0e-8
 cratio_leaf=0.6
 cratio_root=0.4
 cratio_biomass=0.42
 harvest_index=0.8
 turnover_cleaf=0.03 
 turnover_croot=0.004
 sla=10
 q10=2
 invert_option=0
/
EOF

cat << EOF > soilhydro_namelist
&soilhydro_namelist
 soil_depth=0.6
 max_poros=0.54
 fc=0.29
 wp=0.09
 ksat=2.0e-6
 org_depth=0.04 
 org_poros=0.9
 org_fc=0.3 
 maxpond=0.0
 org_sat=1.0   
 n_van=1.18     
 watres=0.0
 alpha_van=3.35
 watsat=0.54
 wmax=0.5
 wmaxsnow=4.5
 hc=0.6
 w_leaf=0.01
 rw=0.20
 rwmin=0.02
 gsoil=5e-3
 kmelt=2.8934e-05
 kfreeze=5.79e-6
 frac_snowliq=0.05
 zmeas=2.0
 zground=0.1
 zo_ground=0.01
/
EOF

cat << EOF > soilyasso_namelist
&soilyasso_namelist
 ! param_y20_map=0.51, 5.19, 0.13, 0.1, 0.5, 0., 1., 1., 0.99, 0., 0., 0., 0., 0., 0.163, 0., 0., 0., 0., 0., 0., 0.158, -0.002, 0.17, -0.005, 0.067, 0.0, -1.44, -2.0, -6.9, 0.0042, 0.0015, -2.55, 1.24, 0.20
 ! nc_mb=0.2
 ! cue_min=0.1
 ! nc_h_max=0.1
 awenh_fineroot=0.46, 0.32, 0.04, 0.18, 0.0
 awenh_leaf=0.46, 0.32, 0.04, 0.18, 0.0
 awenh_soluble=0.0, 1.0, 0.0, 0.0, 0.0
 awenh_compost=0.69, 0.09, 0.02, 0.20, 0.0
 tempr_c=5.4
 tempr_ampl=20
 precip_day=1.87
 totc=16.0
 !cn_input=50
 fract_root_input=0.6
 fract_legacy_soc=0.0
/ 
EOF

srun ./SVMC

fi

if [ ${post_proc} == "T" ]
then
 
 module load geoconda
 module load nco
 cd /users/tanghui1/SVMC/data/output

 for year in {2018..2022}; do
    ncks --quench --mk_rec_dmn time SVM_Qvidja.${year}.day_0.6lr0.3tl_invert0_${expids[i]}.nc test_day_${year}_${expids[i]}.nc
    ncks --quench --mk_rec_dmn time SVM_Qvidja.${year}.hr_0.6lr0.3tl_invert0_${expids[i]}.nc test_hr_${year}_${expids[i]}.nc
 done
 ncrcat test_day_*_${expids[i]}.nc  SVM_Qvidja_day_0.6lr0.3tl_invert0_${expids[i]}.nc
 ncrcat test_hr_*_${expids[i]}.nc  SVM_Qvidja_hr_0.6lr0.3tl_invert0_${expids[i]}.nc
 
 rm test_day_*_${expids[i]}.nc
 rm test_hr_*_${expids[i]}.nc

fi

done
