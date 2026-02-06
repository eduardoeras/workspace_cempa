'''
MONAN — 24-hour Accumulated Precipitation Analysis

**Purpose**

This script analyzes MONAN model diagnostic output to compute and visualize 24-hour
accumulated precipitation fields. The workflow is adapted from Dr. Andre Lyra existing
operational Python script.

**Scope of this script**
- Read MONAN diagnostic NetCDF files
- Compute 24-hour accumulated precipitation from cumulative fields
- Calculate global and regional statistics
- Generate georeferenced maps using Cartopy
- Export figures and NetCDF outputs

Usage:
    python MONAN_Prec.py YEAR MONTH DIA HOUR LENGHT

Example:
    python MONAN_Prec.py 2025 12 01 00 120

The original script is treated as the *reference implementation* and its results are preserved."
'''

#IMPORTS
import argparse
import datetime
import os
import xarray as xr
import numpy as np

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm

from pathlib import Path

#FILE PATHS
OUTPUT_PATH = "/home2/eduardo.eras/workspace/python/output"
BASE_PATH = "/p/projetos/monan_atm/eduardo.eras/dataout"
FILE_PREFIX = "MONAN_DIAG_G_POS_GFS_"
FILE_SUFIX = ".00.00.x655362L55.nc"

#COLOR PALETTE 
color_legend_rgb = [
    (255, 255, 255),  # Branco
    (220, 220, 220),  # Cinza Claro
    (180, 180, 180),  # Cinza
    (20, 0, 150),     # Azul Marinho/Roxo
    (0, 0, 255),      # Azul
    (0, 100, 100),    # Verde Escuro/Azul Petróleo
    (0, 200, 0),      # Verde
    (150, 255, 0),    # Verde Limão/Ciano
    (255, 255, 0),    # Amarelo Claro
    (255, 220, 0),    # Amarelo Escuro/Ouro
    (255, 130, 0),    # Laranja
    (230, 25, 25),    # Vermelho Claro
    (100, 0, 0),      # Vermelho Escuro/Borgonha
]

#PARSER
parser = argparse.ArgumentParser(
    description="Compute and visualize 24-hour accumulated precipitation from MONAN model output.",
    formatter_class=argparse.RawDescriptionHelpFormatter
)
parser.add_argument("YEAR", type=str, help="Year of the analysis (e.g., 2025)")
parser.add_argument("MONTH", type=str, help="Month of the analysis (1-12)")
parser.add_argument("DIA", type=str, help="Day of the analysis (1-31)")
parser.add_argument("HOUR", type=str, help="Hour of the analysis (0-23)")
parser.add_argument("LENGHT", type=int, help="Forecast lead time in hours (e.g., 120)")

args = parser.parse_args()

initial_base_date = datetime.datetime(
    int(args.YEAR),
    int(args.MONTH),
    int(args.DIA),
    int(args.HOUR)
)

initial_date = f"{initial_base_date:%Y%m%d%H}"
delta_hours = datetime.timedelta(hours=0)
delta_day = datetime.timedelta(hours=24)
days = args.LENGHT // 24

#MAIN LOOP
for i in range(days):
    print(f"\nProcessing day {i+1} of {days}")

    # Accumulation period
    accumulation_initial_date = initial_base_date + (delta_day * i) + (delta_hours * (i+1))
    accumulation_final_date = accumulation_initial_date + delta_day

    #Forecasting time
    Fct=f"{((i+1)*24):03d}"

    # Date strings
    date_format = '%Y%m%d%H'
    date_str_i = accumulation_initial_date.strftime(date_format)
    date_str_e = accumulation_final_date.strftime(date_format)

    #File paths
    DIR_DATE = f"{initial_date}/"
    file_name_i = f"{FILE_PREFIX}{initial_date}_{date_str_i}{FILE_SUFIX}" #Initial
    file_name_e = f"{FILE_PREFIX}{initial_date}_{date_str_e}{FILE_SUFIX}" #End

    full_path_i = os.path.join(BASE_PATH, DIR_DATE, file_name_i)
    full_path_e = os.path.join(BASE_PATH, DIR_DATE, file_name_e)

    fi = xr.open_dataset(full_path_i, engine="netcdf4")
    fe = xr.open_dataset(full_path_e, engine="netcdf4")

    # Sum of accumulated variables
    vari = (fi["rainnc"].isel(Time=0) + fi["rainc"].isel(Time=0)).values
    vare = (fe["rainnc"].isel(Time=0) + fe["rainc"].isel(Time=0)).values

    # Geographic coordinates (in degrees)
    lat = (fi["latitude"].values)
    lon = (fi["longitude"].values)
    
    # Difference (accumulated rainfall between the two periods)
    prec_acum = ((vare - vari) ).squeeze()

    # Calculating the Global Maximum
    prec_max = np.nanmax(prec_acum)
    mask_lat_AMS = (lat >= -55) & (lat <= 20)
    mask_lon_AMS = (lon >= 275) & (lon <= 340)
    prec_max_AMS = np.nanmax(prec_acum[np.ix_(mask_lat_AMS, mask_lon_AMS)])
    mask_lat_ACC = (lat >= -10) & (lat <= 35)
    mask_lon_ACC = (lon >= 242) & (lon <= 325)
    prec_max_ACC = np.nanmax(prec_acum[np.ix_(mask_lat_ACC, mask_lon_ACC)])

    # Calculating the Global Average
    weights = np.cos(np.deg2rad(lat))    
    w2d = weights[:, None]                 
    mask = np.isfinite(prec_acum)
    num = np.nansum(prec_acum * w2d * mask)
    den = np.nansum(w2d * mask)
    prec_media = num / den   
    
    prec_sub = prec_acum[np.ix_(mask_lat_AMS, mask_lon_AMS)]
    lat_sub = lat[mask_lat_AMS]
    weights = np.cos(np.deg2rad(lat_sub))[:, None]
    num = np.nansum(prec_sub * weights)
    den = np.nansum(weights * np.isfinite(prec_sub))
    prec_media_AMS = num / den

    prec_sub = prec_acum[np.ix_(mask_lat_ACC, mask_lon_ACC)]
    lat_sub = lat[mask_lat_ACC]
    weights = np.cos(np.deg2rad(lat_sub))[:, None]
    num = np.nansum(prec_sub * weights)
    den = np.nansum(weights * np.isfinite(prec_sub))
    prec_media_ACC = num / den

    # Format the values ​​for inclusion in the title.
    max_str = f"Max: {prec_max:.2f} mm"
    media_str = f"Mean: {prec_media:.2f} mm"
    max_str_AMS = f"Max: {prec_max_AMS:.2f} mm"
    media_str_AMS = f"Mean: {prec_media_AMS:.2f} mm"
    max_str_ACC = f"Max: {prec_max_ACC:.2f} mm"
    media_str_ACC = f"Mean: {prec_media_ACC:.2f} mm" 

    # Normalize RGB colors for Matplotlib
    normalized_colors_matplotlib = []
    for r, g, b in color_legend_rgb:
        normalized_colors_matplotlib.append((r / 255.0, g / 255.0, b / 255.0))

    # Creating the color map and the levels (clevs)
    cmap = ListedColormap(normalized_colors_matplotlib)
    clevs = [0, 1, 2, 4, 6, 10, 15, 25, 35, 50, 75, 100, 150]
    clevs_full = clevs + [40]
    norm = BoundaryNorm(clevs, ncolors=len(clevs), extend='max')

    # Cartopy projection
    proj = ccrs.PlateCarree()
    lon2d, lat2d = np.meshgrid(lon, lat)

    #############################################################################
    # Plotting Global Map                                                       #
    #############################################################################
    fig = plt.figure(figsize=(10,5))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())

    plt_triang = ax.pcolormesh(
        lon2d, lat2d, prec_acum,
        cmap=cmap,
        norm=norm,
        shading='auto',
        transform=proj
    )

    # Add coastlines and gridlines
    ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.6)
    ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.4)
    ax.add_feature(cfeature.LAND.with_scale('50m'), facecolor='none', edgecolor='black', linewidth=0.2)

    # gridlines
    gl = ax.gridlines(draw_labels=True, linewidth=0.3, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.bottom_labels = True
    gl.left_labels = True
    gl.xlines = False
    gl.ylines = False
    gl.xlocator = plt.FixedLocator(np.arange(-180, 181, 60))
    gl.ylocator = plt.FixedLocator(np.arange(-60, 61, 30))

    # Colorbar
    cbar = plt.colorbar(plt_triang, orientation='vertical', pad=0.1, aspect=50, boundaries=clevs, extend='max')
    cbar.set_label('(mm/day)')
    cbar.set_ticks(clevs)

    # Title
    plt.title(
        f'MONAN {initial_date}+{date_str_e} - 24h prec accum for {Fct}h \n' 
        f'Global - {max_str} | {media_str}', 
        fontsize=12
    )

    # Save figure
    nome_arquivo = f"MONAN_24precacum_{initial_date}_{date_str_e}_GLB.png"
    caminho_out = f"{OUTPUT_PATH}/MONAN/{args.YEAR}{args.MONTH}/{initial_date}/"
    os.makedirs(caminho_out, exist_ok=True)
    plt.savefig(os.path.join(caminho_out, nome_arquivo), dpi=300, bbox_inches="tight")
    plt.close(fig)

    #############################################################################
    # Plotting South America Map                                                #
    #############################################################################
    fig = plt.figure(figsize=(10,5))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent([-85, -20, -55, 20], crs=ccrs.PlateCarree())
    xticks = np.arange(-80, -19, 10)
    yticks = np.arange(-50, 21, 10)
    ax.set_xticks(xticks, crs=ccrs.PlateCarree())
    ax.set_yticks(yticks, crs=ccrs.PlateCarree())
    ax.xaxis.set_major_formatter(LongitudeFormatter(number_format='.0f'))
    ax.yaxis.set_major_formatter(LatitudeFormatter(number_format='.0f'))
    
    ax.tick_params(
        axis='both',
        which='major',
        direction='in',
        length=2,
        width=0.6,
        labelsize=8,
        top=True,
        right=True
    )
        
    plt_triang = ax.pcolormesh(
        lon2d, lat2d, prec_acum,
        cmap=cmap,
        norm=norm,
        shading='auto',
        transform=proj
    )

    # Add coastlines and gridlines
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

    #Gridlines
    gl = ax.gridlines(
        draw_labels=True,
        linewidth=0.3,
        color='gray',
        alpha=0.5,
        linestyle='--'
    )

    gl.top_labels = False
    gl.right_labels = False
    gl.bottom_labels = False
    gl.left_labels = False
    gl.xlines = False
    gl.ylines = False
    gl.xlocator = plt.FixedLocator(np.arange(-180, 181, 10))
    gl.ylocator = plt.FixedLocator(np.arange(-60, 61, 10))
    gl.xlabel_style = {'size': 9}
    gl.ylabel_style = {'size': 9}

    # Colorbar
    cbar = plt.colorbar(plt_triang, orientation='vertical', pad=0.1, aspect=50, boundaries=clevs, extend='max')
    cbar.set_label('(mm/day)')
    cbar.set_ticks(clevs)

    # Title
    plt.title(
        f'MONAN {initial_date}+{date_str_e} - prec24h {Fct}h\n' 
        f'AMS - {max_str_AMS} | {media_str_AMS}', 
        fontsize=12
    )

    # Save figure
    file_name = f"MONAN_24precacum_{initial_date}_{date_str_e}_AMS.png"
    out_path = f"{OUTPUT_PATH}/MONAN/{args.YEAR}{args.MONTH}/{initial_date}/"
    plt.savefig(os.path.join(out_path, file_name), dpi=300, bbox_inches="tight")
    plt.close(fig)

    #############################################################################
    # Plotting Central America and Caribbean                                    #
    #############################################################################
    fig = plt.figure(figsize=(10,5))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent([-118, -35, -10, 35], crs=ccrs.PlateCarree())
    xticks = np.arange(-110, -35, 10)
    yticks = np.arange(-10, 36, 10)
    ax.set_xticks(xticks, crs=ccrs.PlateCarree())
    ax.set_yticks(yticks, crs=ccrs.PlateCarree())
    ax.xaxis.set_major_formatter(LongitudeFormatter(number_format='.0f'))
    ax.yaxis.set_major_formatter(LatitudeFormatter(number_format='.0f'))
    
    ax.tick_params(
        axis='both',
        which='major',
        direction='in',
        length=2,
        width=0.6,
        labelsize=8,
        top=True,
        right=True
    )
        
    plt_triang = ax.pcolormesh(
    lon2d, lat2d, prec_acum,
        cmap=cmap,
        norm=norm,
        shading='auto',
        transform=proj
    )

    # Add coastlines and gridlines
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

    #Gridlines
    gl = ax.gridlines(
        draw_labels=True,
        linewidth=0.3,
        color='gray',
        alpha=0.5,
        linestyle='--'
    )

    gl.top_labels = False
    gl.right_labels = False
    gl.bottom_labels = False
    gl.left_labels = False
    gl.xlines = False
    gl.ylines = False
    gl.xlocator = plt.FixedLocator(np.arange(-180, 181, 10))
    gl.ylocator = plt.FixedLocator(np.arange(-60, 61, 10))
    gl.xlabel_style = {'size': 9}
    gl.ylabel_style = {'size': 9}

    # Colorbar
    cbar = plt.colorbar(plt_triang, orientation='vertical', pad=0.1, aspect=50, boundaries=clevs, extend='max')
    cbar.set_label('(mm/day)')
    cbar.set_ticks(clevs)

    # Title
    plt.title(
        f'MONAN {initial_date}+{date_str_e} - 24h prec accum for {Fct}h \n' 
        f'ACC - {max_str_ACC} | {media_str_ACC}', 
    fontsize=12
    )

    # Salvar figura
    nome_arquivo = f"MONAN_24precacum_{initial_date}_{date_str_e}_ACC.png"
    caminho_out = f"{OUTPUT_PATH}/MONAN/{args.YEAR}{args.MONTH}/{initial_date}/"
    plt.savefig(os.path.join(caminho_out, nome_arquivo), dpi=300, bbox_inches="tight")
    plt.close(fig)

    #############################################################################
    # Exporting NetCDF File                                                     #
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

    prec_out = np.ascontiguousarray(prec_acum, dtype='float32')

    ds = xr.Dataset(
       {
        'prec': (['lat', 'lon'], prec_out)
       },
       coords={
        'lat': lat_da,
        'lon': lon_da
       }
    )

    encoding = {
       'prec': {
               'dtype': 'float32',
               '_FillValue': -999.0,
       }
    }

    #Saving NetCDF
    file_name = f"MONAN_Precipitation_24h_acum_{initial_date}_{date_str_e}_{Fct}h.nc"
    netcdf_path = Path(OUTPUT_PATH) / "NetCDFs" / initial_date / file_name
    netcdf_path.parent.mkdir(parents=True, exist_ok=True)
    ds.to_netcdf(netcdf_path, encoding=encoding, format='NETCDF4')