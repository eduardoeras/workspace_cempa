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
# VARIABLES
debug = False

# IMPORTS
import argparse
import datetime
import os
import subprocess
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.feature as cfeature
from matplotlib.colors import ListedColormap, BoundaryNorm

import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

from pathlib import Path

# FILE PATHS
INPUT_PATH= "/home2/eduardo.eras/workspace/python/output/NetCDFs"
OUTPUT_PATH = "/home2/eduardo.eras/workspace/python/output"

# Color scale
colors_rgb = [
#   (90, 0, 0)        # 25
    (130, 0, 0),      # 20
    (192, 3, 0),      # 15
    (225, 18, 0),     # 12
    (255, 96, 2),     # 9
    (255, 193, 60),   # 6
    (255, 251, 190),  # 3
    (255, 255, 255),  # 0
    (179, 230, 249),  # -3
    (150, 200, 249),  # -6
    (75, 155, 244),   # -9
    (36, 116, 241),   # -12
    (25, 90, 234),    # -15
    (0,  45, 220),    # -20        
#   (0,  25, 200),    # -25
]

levels = [ -25, -20, -15, -12, -9, -6, -3, 3, 6, 9, 12, 15, 20, 25]
cmap = ListedColormap( [(r/255, g/255, b/255) for r, g, b in colors_rgb] )
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

# DATE PARSER
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

# MAIN LOOP
for lead in lead_times:
    print(f"\nProcessing lead {lead} of {total_lenght}")

    end_date = init_date + datetime.timedelta(hours=lead)

    ciclo_str = init_date.strftime("%Y%m%d%H")
    init_date_mod_str = init_date.strftime("%Y%m%d%H")
    end_date_mod_str = end_date.strftime("%Y%m%d%H")
    end_date_obs_str = end_date.strftime("%Y%m%d%H")

    # Input file paths
    monan_nc = (
        f"{INPUT_PATH}"
        f"/{ciclo_str}/"
        f"MONAN_Precipitation_24h_acum_"
        f"{init_date_mod_str}_{end_date_mod_str}_{lead:03d}h.nc"
    )

    gpm_nc = (
        f"{INPUT_PATH}"
        f"/GPM_IMERG/{end_date_obs_str}00/"
        f"GPM_IMERG_Precipitation_24h_accum_{end_date_obs_str}00.nc"
    )

    gsmap_nc = (
        f"{INPUT_PATH}"
        f"/GSMAP/{end_date_obs_str}00/"
        f"GSMAP_Precipitation_24h_accum_{end_date_obs_str}00.nc"
    )

    mswep_nc = (
        f"{INPUT_PATH}"
        f"/MSWEP/{end_date_obs_str}00/"
        f"MSWEP_Precipitation_24h_accum_{end_date_obs_str}00.nc"
    )

    # Input file name debug prints
    if debug:
        print("Model file:", monan_nc)
        print("GPM file:", gpm_nc)
        print("GSMAP file:", gsmap_nc)
        print("MSWEP file:", mswep_nc)

    # Output file path
    fig_dir = f"{OUTPUT_PATH}/Bias/{args.YEAR}{args.MONTH}/{ciclo_str}/"
    os.makedirs(fig_dir, exist_ok=True)

    # CDO Function (Climate Data Operators) remap observational data to model grid
    def remap_cdo(obs_nc, out_nc, ref_nc):
        if not os.path.exists(out_nc):
            cmd = [
                "cdo",
                "-f", "nc",
                f"-remapcon,{ref_nc}",
                obs_nc,
                out_nc
            ]
            subprocess.run(cmd, check=True)
            print(f"Remapping: {out_nc}")
        else:
            print(f"File already exists: {out_nc}")

    # Regrid
    gpm_remap = gpm_nc.replace(".nc", "_MONAN_grid.nc")
    gsmap_remap = gsmap_nc.replace(".nc", "_MONAN_grid.nc")
    mswep_remap = mswep_nc.replace(".nc", "_MONAN_grid.nc")

    remap_cdo(gpm_nc, gpm_remap, monan_nc)
    remap_cdo(gsmap_nc, gsmap_remap, monan_nc)
    remap_cdo(mswep_nc, mswep_remap, monan_nc)

    # NetCDF reading
    ds_monan = xr.open_dataset(monan_nc)
    ds_gpm = xr.open_dataset(gpm_remap)
    ds_gsmap = xr.open_dataset(gsmap_remap)
    ds_mswep = xr.open_dataset(mswep_remap)

    # Variable name adjustment
    var_monan = "prec"
    var_gpm = "prec"
    var_gsmap = "prec"
    var_mswep = "prec"

    monan = ds_monan[var_monan]
    gpm = ds_gpm[var_gpm]
    gsmap = ds_gsmap[var_gsmap]
    mswep = ds_mswep[var_mswep]

    lat = ds_monan["lat"]
    lon = ds_monan["lon"]

    # Differences
    diff_monan_gpm = monan - gpm
    diff_monan_gsmap = monan - gsmap
    diff_monan_mswep = monan - mswep

    abs_diff_monan_gpm = np.abs(monan - gpm)
    abs_diff_monan_gsmap = np.abs(monan - gsmap)
    abs_diff_monan_mswep = np.abs(monan - mswep)

    # Averages
    diff_med_gpm = np.nanmean(diff_monan_gpm)
    diff_med_gpm_AMS =  np.nanmean(diff_monan_gpm.sel(lat=slice(-55, 20),lon=slice(275, 340)))
    diff_med_gpm_ACC =  np.nanmean(diff_monan_gpm.sel(lat=slice(-10, 35),lon=slice(242, 335)))

    diff_med_gsmap = np.nanmean(diff_monan_gsmap)
    diff_med_gsmap_AMS =  np.nanmean(diff_monan_gsmap.sel(lat=slice(-55, 20),lon=slice(275, 340)))
    diff_med_gsmap_ACC =  np.nanmean(diff_monan_gsmap.sel(lat=slice(-10, 35),lon=slice(242, 335)))

    diff_med_mswep = np.nanmean(diff_monan_mswep)
    diff_med_mswep_AMS =  np.nanmean(diff_monan_mswep.sel(lat=slice(-55, 20),lon=slice(275, 340)))
    diff_med_mswep_ACC =  np.nanmean(diff_monan_mswep.sel(lat=slice(-10, 35),lon=slice(242, 335)))

    abs_diff_med_gpm = np.nanmean(abs_diff_monan_gpm)
    abs_diff_med_gpm_AMS =  np.nanmean(abs_diff_monan_gpm.sel(lat=slice(-55, 20),lon=slice(275, 340)))
    abs_diff_med_gpm_ACC =  np.nanmean(abs_diff_monan_gpm.sel(lat=slice(-10, 35),lon=slice(242, 335)))

    abs_diff_med_gsmap = np.nanmean(abs_diff_monan_gsmap)
    abs_diff_med_gsmap_AMS =  np.nanmean(abs_diff_monan_gsmap.sel(lat=slice(-55, 20),lon=slice(275, 340)))
    abs_diff_med_gsmap_ACC =  np.nanmean(abs_diff_monan_gsmap.sel(lat=slice(-10, 35),lon=slice(242, 335)))

    abs_diff_med_mswep = np.nanmean(abs_diff_monan_mswep)
    abs_diff_med_mswep_AMS =  np.nanmean(abs_diff_monan_mswep.sel(lat=slice(-55, 20),lon=slice(275, 340)))
    abs_diff_med_mswep_ACC =  np.nanmean(abs_diff_monan_mswep.sel(lat=slice(-10, 35),lon=slice(242, 335)))

    # Title values
    media_str_gpm = f"MAE: {abs_diff_med_gpm:.2f} mm, Bias: {diff_med_gpm:.2f} mm"
    media_str_gpm_AMS = f"MAE: {abs_diff_med_gpm_AMS:.2f} mm, Bias: {diff_med_gpm_AMS:.2f} mm"
    media_str_gpm_ACC = f"MAE: {abs_diff_med_gpm_ACC:.2f} mm, Bias: {diff_med_gpm_ACC:.2f} mm"

    media_str_gsmap = f"MAE: {abs_diff_med_gsmap:.2f} mm, Bias: {diff_med_gsmap:.2f} mm"
    media_str_gsmap_AMS = f"MAE: {abs_diff_med_gsmap_AMS:.2f} mm, Bias: {diff_med_gsmap_AMS:.2f} mm"
    media_str_gsmap_ACC = f"MAE: {abs_diff_med_gsmap_ACC:.2f} mm, Bias: {diff_med_gsmap_ACC:.2f} mm"

    media_str_mswep = f"MAE: {abs_diff_med_mswep:.2f} mm, Bias: {diff_med_mswep:.2f} mm"
    media_str_mswep_AMS = f"MAE: {abs_diff_med_mswep_AMS:.2f} mm, Bias: {diff_med_mswep_AMS:.2f} mm"
    media_str_mswep_ACC = f"MAE: {abs_diff_med_mswep_ACC:.2f} mm, Bias: {diff_med_mswep_ACC:.2f} mm"

    #############################################################################
    # Plotting Global Map                                                       #
    #############################################################################
    def plot_diff(data, titulo, nome_fig):
       fig = plt.figure(figsize=(10,5))
       ax = plt.axes(projection=ccrs.PlateCarree())

       im = data.plot(
       ax=ax,
       transform=ccrs.PlateCarree(),
       cmap=cmap,
       norm=norm,
       add_colorbar=False
       )

       ax.set_xticks(np.arange(-180, 180+1 , 60), crs=ccrs.PlateCarree())
       ax.set_yticks(np.arange(-60, 60, 30), crs=ccrs.PlateCarree())      
       ax.xaxis.set_major_formatter(LongitudeFormatter())
       ax.yaxis.set_major_formatter(LatitudeFormatter())
       ax.set_xlabel("")
       ax.set_ylabel("")
       ax.coastlines()
       ax.set_title(titulo) 

       cbar = plt.colorbar(im, orientation='vertical', pad=0.1, aspect=50, boundaries=levels, extend='both')
       cbar.set_label('(mm)')
       cbar.set_ticks( [-20, -15, -12, -9, -6, -3, 3, 6, 9, 12, 15, 20])                 

       plt.savefig(os.path.join(fig_dir, nome_fig), dpi=300, bbox_inches="tight")
       plt.close()
       print(f"Figura salva: {nome_fig}")

    plot_diff(
    diff_monan_gpm,
    f"Bias {lead:03d}h MONAN vs GPM IMERG \n Global - {media_str_gpm}",
    f"diff_MONAN_GPM_GLB_{lead:03d}h.png"
    )

    plot_diff(
    diff_monan_gsmap,
    f"Bias {lead:03d}h MONAN vs GSMAP \n Global - {media_str_gsmap}",
    f"diff_MONAN_GSMAP_GLB_{lead:03d}h.png"
    )

    plot_diff(
    diff_monan_mswep,
    f"Bias {lead:03d}h MONAN vs MSWEP \n Global - {media_str_mswep}",
    f"diff_MONAN_MSWEP_GLB_{lead:03d}h.png"
    )

    #############################################################################
    # Plotting AMS Map                                                          #
    #############################################################################
    def plot_diff_AMS(data, titulo, nome_fig):
       ax = plt.axes(projection=ccrs.PlateCarree())

       im = data.plot(
       ax=ax,
       transform=ccrs.PlateCarree(),
       cmap=cmap,
       norm=norm, 
       add_colorbar=False
       )

       ax.set_xticks(np.arange(-180, 180+1 , 10), crs=ccrs.PlateCarree())
       ax.set_yticks(np.arange(-60, 60, 10), crs=ccrs.PlateCarree())      
       ax.xaxis.set_major_formatter(LongitudeFormatter())
       ax.yaxis.set_major_formatter(LatitudeFormatter())
       ax.set_xlabel("")
       ax.set_ylabel("")
       ax.coastlines()
       ax.set_title(titulo) 
       ax.set_extent([-85, -20, -55, 20], crs=ccrs.PlateCarree())

       ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.6)
       ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.4)
       ax.add_feature(cfeature.LAND.with_scale('50m'), facecolor='none', edgecolor='black', linewidth=0.2)
       estados = cfeature.NaturalEarthFeature(
       category='cultural',
       name='admin_1_states_provinces_lines',
       scale='50m',
       facecolor='none'
       )
       ax.add_feature(estados, edgecolor='black', linewidth=0.4)
       
       cbar = plt.colorbar(im, orientation='vertical', pad=0.1, aspect=50, boundaries=levels, extend='both')
       cbar.set_label('(mm)')
       cbar.set_ticks( [-20, -15, -12, -9, -6, -3, 3, 6, 9, 12, 15, 20])           
       
       plt.savefig(os.path.join(fig_dir, nome_fig), dpi=300, bbox_inches="tight")
       plt.close()
       print(f"Figura salva: {nome_fig}")

    # Plotting AMS
    plot_diff_AMS(
    diff_monan_gpm,
    f"Bias {lead:03d}h MONAN {init_date_mod_str} vs GPM IMERG \n AMS - {media_str_gpm_AMS}",
    f"diff_MONAN_GPM_AMS_{lead:03d}h.png"
    )

    plot_diff_AMS(
    diff_monan_gsmap,
    f"Bias {lead:03d}h MONAN {init_date_mod_str} vs GSMAP \n AMS - {media_str_gsmap_AMS}",
    f"diff_MONAN_GSMAP_AMS_{lead:03d}h.png"
    )

    plot_diff_AMS(
    diff_monan_mswep,
    f"Bias {lead:03d}h MONAN {init_date_mod_str} vs MSWEP \n AMS - {media_str_mswep_AMS}",
    f"diff_MONAN_MSWEP_AMS_{lead:03d}h.png"
    )

    #############################################################################
    # Plotting ACC Map                                                          #
    #############################################################################
    def plot_diff_AMS(data, titulo, nome_fig):
       ax = plt.axes(projection=ccrs.PlateCarree())

       im = data.plot(
       ax=ax,
       transform=ccrs.PlateCarree(),
       cmap=cmap,
       norm=norm,
       add_colorbar=False
       )

       ax.set_xticks(np.arange(-180, 180+1 , 10), crs=ccrs.PlateCarree())
       ax.set_yticks(np.arange(-60, 60, 10), crs=ccrs.PlateCarree())      
       ax.xaxis.set_major_formatter(LongitudeFormatter())
       ax.yaxis.set_major_formatter(LatitudeFormatter())
       ax.set_xlabel("")
       ax.set_ylabel("")
       ax.coastlines()
       ax.set_title(titulo) 
       ax.set_extent([-118, -35, -10, 35], crs=ccrs.PlateCarree())


       ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.6)
       ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.4)
       ax.add_feature(cfeature.LAND.with_scale('50m'), facecolor='none', edgecolor='black', linewidth=0.2)
       estados = cfeature.NaturalEarthFeature(
       category='cultural',
       name='admin_1_states_provinces_lines',
       scale='50m',
       facecolor='none'
       )
       ax.add_feature(estados, edgecolor='black', linewidth=0.4)
   
       cbar = plt.colorbar(im, orientation='vertical', pad=0.1, aspect=50, boundaries=levels, extend='both')
       cbar.set_label('(mm)')
       cbar.set_ticks( [-20, -15, -12, -9, -6, -3, 3, 6, 9, 12, 15, 20])             
   
       plt.savefig(os.path.join(fig_dir, nome_fig), dpi=300, bbox_inches="tight")
       plt.close()
       print(f"Figura salva: {nome_fig}")

    # Plotting ACC
    plot_diff_AMS(
    diff_monan_gpm,
    f"Bias {lead:03d}h MONAN {init_date_mod_str} vs GPM IMERG \n ACC - {media_str_gpm_ACC}",
    f"diff_MONAN_GPM_ACC_{lead:03d}h.png"
    )

    plot_diff_AMS(
    diff_monan_gsmap,
    f"Bias {lead:03d}h MONAN {init_date_mod_str} vs GSMAP \n ACC - {media_str_gsmap_ACC}",
    f"diff_MONAN_GSMAP_ACC_{lead:03d}h.png"
    )

    plot_diff_AMS(
    diff_monan_mswep,
    f"Bias {lead:03d}h MONAN {init_date_mod_str} vs MSWEP \n ACC - {media_str_mswep_ACC}",
    f"diff_MONAN_MSWEP_ACC_{lead:03d}h.png"
    )

    #############################################################################
    # Saving NetCDF files with bias fields                                      #
    #############################################################################
    lat_da = xr.DataArray(
       lat,
       dims=['lat'],
       attrs={
        'units': 'degrees_north',
        'long_name': 'Latitude',
        '_CoordinateAxisType': 'Lat',
        'axis': 'Y'
       }
    )

    lon_da = xr.DataArray(
       lon,
       dims=['lon'],
       attrs={
        'units': 'degrees_east',
        'long_name': 'Longitude',
        '_CoordinateAxisType': 'Lon',
        'axis': 'X'
       }
    )

    # Output arrays
    diff_monan_gpm_out   = np.ascontiguousarray(diff_monan_gpm,   dtype='float32')
    diff_monan_gsmap_out = np.ascontiguousarray(diff_monan_gsmap, dtype='float32')
    diff_monan_mswep_out = np.ascontiguousarray(diff_monan_mswep, dtype='float32')

    #Dataset
    ds = xr.Dataset(
    {
        'bias_monangpm':   (['lat', 'lon'], diff_monan_gpm_out),
        'bias_monangsmap': (['lat', 'lon'], diff_monan_gsmap_out),
        'bias_monanmswep': (['lat', 'lon'], diff_monan_mswep_out),
    },
    coords={
        'lat': lat_da,
        'lon': lon_da
    }
    )

    # Encoding for all variables
    encoding = {
    'bias_monangpm': {
        'dtype': 'float32',
        '_FillValue': -999.0,
    },
    'bias_monangsmap': {
        'dtype': 'float32',
        '_FillValue': -999.0,
    },
    'bias_monanmswep': {
        'dtype': 'float32',
        '_FillValue': -999.0,
    }
    }

    # Output NetCDF
    file_name = f"Bias_MONAN_Prec_{init_date_mod_str}_{lead:03d}h.nc"
    netcdf_path = Path(OUTPUT_PATH) / "NetCDF" / init_date_mod_str / file_name
    netcdf_path.parent.mkdir(parents=True, exist_ok=True)
    ds.to_netcdf(netcdf_path, encoding=encoding, format='NETCDF4')  
    print(f"Netcdf: {netcdf_path}")

print("\nProcessing completed for all lead times.\n")

print("\n\nAll lead times processed. Calling figure cutting script...\n\n")

# Call figure cutting script after all lead times are processed
subprocess.run(
    ["bash", "MONAN_Bias.sh", f"{init_date_mod_str}", f"{total_lenght}", f"{OUTPUT_PATH}"],
    check=True
)  
