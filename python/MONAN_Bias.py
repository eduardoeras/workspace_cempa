"""
MONAN — Bias and MAE Analysis Against GPM, GSMAP, and MSWEP

Purpose
-------
This script evaluates the performance of MONAN 24-hour accumulated precipitation
forecasts by computing Bias and Mean Absolute Error (MAE) relative to multiple
observational precipitation products (GPM IMERG, GSMAP, and MSWEP).

The workflow is adapted from an existing operational Python script developed by
Dr. André Lyra.

Scope
-----
- Read MONAN 24-hour accumulated precipitation NetCDF files
- Read observational 24-hour accumulated precipitation datasets (GPM, GSMAP, MSWEP)
- Remap observational datasets to the MONAN grid using CDO
- Compute Bias and Mean Absolute Error (MAE) for each forecast cycle
- Compute global and regional statistics (Global, AMS, ACC)
- Generate georeferenced bias maps using Cartopy
- Export figures and NetCDF files containing bias fields

Usage
-----
    python MONAN_Bias.py YEAR MONTH DAY HOUR FORECAST_HOURS

Arguments
---------
YEAR : int
    Year of model initialization (e.g., 2025)
MONTH : int
    Month of model initialization (e.g., 12)
DAY : int
    Day of model initialization (e.g., 01)
HOUR : int
    Hour of model initialization in UTC (e.g., 00)
FORECAST_HOURS : int
    Total forecast length in hours (e.g., 72, 120)

Example
-------
    python MONAN_Bias.py 2025 12 01 00 120

Reference
---------
The original script is treated as the *reference implementation*, and its numerical
results are preserved.
"""

#IMPORTS
import argparse
import datetime

#DATE PARSER
parser = argparse.ArgumentParser(
        description="Generates bias for each MONAN forecast cycle.", 
        formatter_class=argparse.RawTextHelpFormatter 
)
parser.add_argument('YEAR', type=str, help='Year of launch (e.g., 2025)')
parser.add_argument('MONTH', type=str, help='Month of launch (e.g., 12)')
parser.add_argument('DAY', type=str, help='Day of launch (e.g., 01)')
parser.add_argument('HOUR', type=str, help='Hour of launch (e.g., 00)')
parser.add_argument('LENGHT', type=int, help='Total forecast length in hours (e.g., 72, 120)')
args = parser.parse_args()

init_date = datetime.datetime(
    int(args.YEAR),
    int(args.MONTH),
    int(args.DAY),
    int(args.HOUR)
)

total_lenght = args.LENGHT
lead_times = range(24, total_lenght + 1, 24)

#MAIN LOOP
for lead in lead_times:
    print(f"\nProcessing lead {lead} of {total_lenght}")

    end_date = init_date + datetime.timedelta(hours=lead)

    ciclo_str = init_date.strftime("%Y%m%d%H")
    data_ini_mod_str = init_date.strftime("%Y%m%d%H")
    data_fim_mod_str = end_date.strftime("%Y%m%d%H")
    data_fim_obs_str = end_date.strftime("%Y%m%d%H")