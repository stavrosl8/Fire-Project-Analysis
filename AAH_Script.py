#Editor: Dr Logothetis Stavros-Andreas
#Affiliation: Laboratory of Atmospheric Physics, University of Patras
#Objective: AC-SAF Metop-A/B/C GOME-2 Level 2 (L2) data files processing
#Merging, Transforming, Visualizing daily L2 data 
#The script starts here
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.io.shapereader as shpreader
from matplotlib import pyplot as plt
from netCDF4 import Dataset
import cartopy.crs as ccrs
import xarray as xr
import numpy as np
import tqdm
import glob
import os

def gen_xr_from_1D_vec(file, variable, parameter_name, longname, unit):
    latitude = file['GEOLOCATION/LatitudeCenter']
    longitude = file['GEOLOCATION/LongitudeCenter']
    param = variable 
    param_da = xr.DataArray(param[:],
        dims=["x","y"],
        coords={'latitude':(['x','y'], latitude[:]),
            'longitude':(['x','y'],longitude[:])},
        attrs={'long_name': longname, 'units': unit},
        name=parameter_name)        
    return (param_da)

def mask_data (xarray, mask):
    cloud_mask = xr.where(mask == 1, 1, 0)
    xarray_masked = xr.where(cloud_mask == 1, xarray, np.nan) #Apply mask onto the DataArray
    xarray_masked.attrs = xarray.attrs #Set DataArray attributes 
    xarray_masked.name = 'AAH_AbsorbingAerosolHeight'
    return(xarray_masked)

def visualize_scatter(xr_dataarray, conversion_factor, projection,point_size, color_scale, unit, 
                      title, shapefile, title_save):
  
    plt.rcParams["font.weight"] = "bold"
    plt.rcParams["axes.labelweight"] = "bold"
    fig = plt.figure(figsize=(16, 10))
    ax = fig.add_subplot(1, 1, 1, projection=projection)
   
    ax.set_extent([180, -180, 90, -90])
    ax.set_extent([19.5, 28.5, 34, 42], ccrs.PlateCarree())

    ax.add_geometries(shapefile, ccrs.PlateCarree(), linewidth=0.4,
                   edgecolor='black', facecolor='none', alpha=1)
    
    if (projection==ccrs.PlateCarree()):
        gl = ax.gridlines(draw_labels=True, linestyle='--')
        gl.top_labels=False
        gl.right_labels=False
        gl.xformatter=LONGITUDE_FORMATTER
        gl.yformatter=LATITUDE_FORMATTER
        gl.xlabel_style={'size':14}
        gl.ylabel_style={'size':14}

    from matplotlib.colors import BoundaryNorm
    
    lcmap=plt.cm.get_cmap('jet',12)
    bounds = np.arange(0,13,1)
    norm = BoundaryNorm(bounds, lcmap.N)

    # plot pixel positions
    img = ax.scatter(
        xr_dataarray.longitude.data,
        xr_dataarray.latitude.data,
        c=xr_dataarray.data*conversion_factor,
        cmap=lcmap,
        marker='o',
        s=point_size,norm=norm,
        transform=ccrs.PlateCarree())

    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel("Longitude", fontsize=16)
    plt.ylabel("Latitude", fontsize=16)
    cbar = fig.colorbar(img, ax=ax, orientation='horizontal', fraction=0.04, pad=0.1)
    cbar.set_label(unit, fontsize=16)
    cbar.ax.tick_params(labelsize=14)
    ax.set_title(title, fontsize=20, pad=10.0, fontweight="bold")
    plt.savefig(f'Graphs/GOME/{title_save}.png',dpi=150, bbox_inches='tight')
    plt.close()
    return()

os.chdir('G:/Projects/Project_Forest_Fires_Evros/Forest_Fires')
directory = os.getcwd()

fname = 'Europe.shp' 
adm1_shapes = list(shpreader.Reader(fname).geometries())

#Processing of each file individually
fileslist = glob.glob(directory+'/Data/GOME/*.hdf5')
for filename in tqdm.tqdm(fileslist):
    file = Dataset(filename)
    filename_ = filename.split("\\")[4][:-5]
    if 'AAH_AbsorbingAerosolHeight' in file['DATA'].variables.keys():
        var_aah = file['DATA/AAH_AbsorbingAerosolHeight']
        var_flag = file['DATA/AAH_RegimeFlag'] 
        time_ = file['GEOLOCATION/Time'][0][0].split('.')[0]
        satellite_id = filename_.split('_')[4]
        # to begin from here
        date_str = "Satellite ID: " + satellite_id + " | " + " ".join(time_.split('T')) +' UTC'
        ds_ahh = gen_xr_from_1D_vec(file = file, variable = var_aah, 
                                         parameter_name = 'AAH_AbsorbingAerosolHeight', 
                                         longname = 'Absorbing Aerosol Height', 
                                         unit = 'km')
                
        if (((ds_ahh.latitude>34) & (ds_ahh.latitude<42)) & ((ds_ahh.longitude>19.5) & (ds_ahh.longitude<28.5))).any():
            
            ds_ahh_masked = ds_ahh.where(ds_ahh >= -0.3)  
            ds_flag = gen_xr_from_1D_vec(file = file, variable = var_flag, 
                                             parameter_name = 'AAH_RegimeFlag', 
                                             longname = 'Flag Regime', unit = 'unitless')
            
            ds_flag_masked = ds_flag.where(ds_flag >=0)  
            ds_ahh_masked = mask_data (ds_ahh_masked, ds_flag_masked)
            
            #plot
            visualize_scatter(xr_dataarray=ds_ahh_masked, 
                         conversion_factor=1, 
                         projection=ccrs.PlateCarree(),
                         point_size=180,
                         color_scale='viridis',
                         unit=ds_ahh_masked.long_name + " ("+ds_ahh_masked.units+")", 
                         title=date_str,
                         shapefile=adm1_shapes,
                         title_save=filename_[:-5])