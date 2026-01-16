from matplotlib.tri import Triangulation
import xarray as xr

monan_data = xr.open_dataset(monan_path, engine="h5netcdf")
monan_lons=monan_data['lon']
monan_lats =monan_data['lat']
# Create triangulation
tri = Triangulation(monan_lons, monan_lats)     
       
fig, ax = plt.subplots(figsize=(10, 8), subplot_kw=dict(projection=ccrs.PlateCarree()))

 # Plot
tpc = ax.tripcolor(
      tri, monan_data[var], 
      cmap="RdYlBu_r", 
      shading="flat"
                )
plt.colorbar(tpc, label= "var label", shrink= 0.5, aspect = 25, pad = 0.05)