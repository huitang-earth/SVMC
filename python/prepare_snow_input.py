import cdsapi
import argparse

parser = argparse.ArgumentParser()
# set start and end time
parser.add_argument('-out', '--output', type=str, default='lai.nc', help="netcdf file of lai.")
parser.add_argument('-y', '--year', type=str, default='1980', help="year for downloading the data.")
parser.add_argument('-lat', '--latitude', type=float, default='60', help="center latitude of a site.")
parser.add_argument('-lon', '--longitude', type=float, default='20', help="ceter longitude of a site.")
parser.add_argument('-site', '--sitename', type=str, default='test', help="site name.")

args = parser.parse_args()

year=args.year
lat=args.latitude
lon=args.longitude


c = cdsapi.Client()

c.retrieve(
    'reanalysis-era5-single-levels',
    {
        'product_type': 'reanalysis',
        'format': 'netcdf',
        'variable': [
            '10m_u_component_of_wind', '10m_v_component_of_wind', '2m_dewpoint_temperature',
            '2m_temperature', 'surface_pressure', 'mean_surface_downward_long_wave_radiation_flux', 
            'mean_surface_downward_short_wave_radiation_flux', 'total_precipitation',
        ],
        'year': [year],
        'month': [
            '01', '02', '03',
            '04', '05', '06',
            '07', '08',
        ],
        'day': [
            '01', '02', '03',
            '04', '05', '06',
            '07', '08', '09',
            '10', '11', '12',
            '13', '14', '15',
            '16', '17', '18',
            '19', '20', '21',
            '22', '23', '24',
            '25', '26', '27',
            '28', '29', '30',
            '31',
        ],
        'time': [
            '00:00', '01:00', '02:00',
            '03:00', '04:00', '05:00',
            '06:00', '07:00', '08:00',
            '09:00', '10:00', '11:00',
            '12:00', '13:00', '14:00',
            '15:00', '16:00', '17:00',
            '18:00', '19:00', '20:00',
            '21:00', '22:00', '23:00',
        ],
        'area': [
            ulat, llon, llat, rlon,
        ],
    },
    'era5_ctsm_'+year+'.nc')