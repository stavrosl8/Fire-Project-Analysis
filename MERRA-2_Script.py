# Info: http://hdfeos.org/zoo/LAADS_MOD_py.php
#Editor: Dr Logothetis Stavros-Andreas
#Affiliation: Laboratory of Atmospheric Physics, University of Patras
#Objective: MERRA-2 data files processing
#Merging, Transforming, Visualizing daily L2 data 
#The script starts here
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.io.shapereader as shpreader
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import xarray as xr
import matplotlib
import glob
import os

os.chdir('G:/Projects/Project_Forest_Fires_Evros/Forest_Fires')

fname = 'Europe.shp' 
adm1_shapes = list(shpreader.Reader(fname).geometries())
    
filenames = glob.glob('Data/MERRA-2/*.nc4')

colnames = ['DUEXTTAU', 'SSEXTTAU', 'TOTEXTTAU','BCEXTTAU','SUEXTTAU', 'OCEXTTAU']
coldir = {'DUEXTTAU':'Dust AOD', 'SSEXTTAU':'Sea Salt AOD', 'TOTEXTTAU':'Total AOD','BCEXTTAU': 'Black carbon AOD','SUEXTTAU':"Sulphate AOD", 'OCEXTTAU':'Organic Carbon AOD'}

plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
AOD = 'TOTEXTTAU'
for file in filenames:
        
    ds = xr.open_dataset(file)
    
    for hour in ds.time.values:
        ds_h = ds.sel(time=hour)
        hour = str(hour)
        title = hour.split('T')[0] + " " + hour.split('T')[1].split('.')[0]
        sname = hour.split('T')[0] + "_" + hour.split('T')[1].split('.')[0].split(':')[0]+ "_" + hour.split('T')[1].split('.')[0].split(':')[1]
       # for AOD in colnames:

        fig = plt.figure(figsize=(16, 10))
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
        ax.set_extent([19.5, 28.5, 34, 42], ccrs.PlateCarree())
        cmap = matplotlib.colormaps["jet"] # plt.cm.get_cmap("jet")
        cmap.set_over('darkred')
        img = ds_h[AOD].plot(cmap=cmap,vmax=0,vmin=1.2, add_colorbar=False)
        ax.add_geometries(adm1_shapes, ccrs.PlateCarree(), linewidth=0.4, edgecolor='black', facecolor='none', alpha=1)
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='k', alpha=0.5, linestyle='--')
        gl.xlabel_style = {'size': 40, 'color': 'gray'}
        gl.xlabel_style = {'color': 'black', 'weight': 'bold'}
        gl.ylabel_style = {'size': 40, 'color': 'gray'}
        gl.ylabel_style = {'color': 'black', 'weight': 'bold'}
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style={'size':16}
        gl.ylabel_style={'size':16}
        cbar = fig.colorbar(img, ax=ax, shrink=0.9, extend='max', orientation='vertical', fraction=0.04, pad=0.07)
        cbar.set_label('AOD', fontsize=20)
        cbar.ax.tick_params(labelsize=20)
        ax.set_title(title, fontsize=25, fontweight='bold', y=1.05)
        plt.savefig(f'Graphs/MERRA-2/{AOD}/{AOD}_{sname}.png', dpi=150, bbox_inches='tight')
        plt.close()