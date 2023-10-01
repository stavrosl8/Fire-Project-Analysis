# Info: http://hdfeos.org/zoo/LAADS_MOD_py.php
#Editor: Dr Logothetis Stavros-Andreas
#Affiliation: Laboratory of Atmospheric Physics, University of Patras
#Objective: MODIS Level 2 (L2) data files processing
#Merging, Transforming, Visualizing daily L2 data 
#The script starts here
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from pyresample.geometry import GridDefinition, SwathDefinition
from pyresample.kd_tree import resample_nearest
import cartopy.io.shapereader as shpreader
import matplotlib.pyplot as plt
from pyhdf.SD import SD, SDC
import cartopy.crs as ccrs
import xarray as xr
import numpy as np
import matplotlib
import datetime
import tqdm
import glob
import os

os.chdir('G:/Projects/Project_Forest_Fires_Evros/Forest_Fires')
directory = os.getcwd()
fileslist = glob.glob(directory+'/Data/MODIS/*.hdf')

def read_file (hdf_file, var_name):
    
    data2D = hdf_file.select(var_name)
    data = data2D[:,:].astype(np.double)

    # Retrieve attributes.
    attrs = data2D.attributes(full=1)
    aoa=attrs["add_offset"]
    add_offset = aoa[0]
    fva=attrs["_FillValue"]
    _FillValue = fva[0]
    sfa=attrs["scale_factor"]
    scale_factor = sfa[0]        
    data[data == _FillValue] = np.nan
    data = (data - add_offset) * scale_factor
    return(data)

fname = 'Europe.shp' 
adm1_shapes = list(shpreader.Reader(fname).geometries())
        
for file in tqdm.tqdm(fileslist):
   
    title_save = file.split('\\')[-1][:-4]
    date = datetime.datetime.strptime(file.split('\\')[-1].split('.')[1][1:], '%Y%j').date().strftime('%Y/%m/%d')
    hour = file.split('\\')[-1].split('.')[2][:2] + ":" + file.split('\\')[-1].split('.')[2][2:] + ":00 UTC"
    sat = file.split('\\')[-1].split('.')[0][:5]
    md = file.split('\\')[-1].split('.')[1][5:]
    title = 'Satellite ID: '+sat+" | "+date+" "+hour
   
    hdf = SD(file, SDC.READ)
    
    # Read geolocation dataset.
    lat = hdf.select('Latitude')
    latitude = lat[:,:]
    lon = hdf.select('Longitude')
    longitude = lon[:,:]
    
    AOD = read_file(hdf, 'AOD_550_Dark_Target_Deep_Blue_Combined')
    AOD_QC = read_file(hdf,  'AOD_550_Dark_Target_Deep_Blue_Combined_QA_Flag')
    AOD[AOD_QC != 3] = np.nan   

    # Regrid.
    # Define SwathDefinition.
    swathDef = SwathDefinition(lons=longitude, lats=latitude)    
    
    # Define GridDefinition.
    # 0.1 degree is about 10.11km, which is close enough to native resolution.
    cellSize = 0.1
    min_lon = np.min(longitude)
    max_lon = np.max(longitude)
    min_lat = np.min(latitude)
    max_lat = np.max(latitude)
    x0, xinc, y0, yinc = (min_lon, cellSize, max_lat, -cellSize)
    nx = int(np.floor((max_lon - min_lon) / cellSize))
    ny = int(np.floor((max_lat - min_lat) / cellSize))
    x = np.linspace(x0, x0 + xinc*nx, nx)
    y = np.linspace(y0, y0 + yinc*ny, ny)
    lon_g, lat_g = np.meshgrid(x, y)
    grid_def = GridDefinition(lons=lon_g, lats=lat_g)
    
    # Set radius_of_influence in meters.
    ri = 10000
    AOD_regrid = resample_nearest(swathDef, AOD, grid_def, 
                          radius_of_influence=ri, epsilon=0.5,
                          fill_value=np.nan)
    
    param_da = xr.DataArray(AOD_regrid[:],
        dims=["x","y"],
        coords={'latitude':(['x','y'], lat_g[:]),
            'longitude':(['x','y'],lon_g[:])},
        attrs={'long_name': 'AOD', 'units': 'Unitless'},
        name='AOD') 

    if (((param_da.latitude>34) & (param_da.latitude<42)) & ((param_da.longitude>19.5) & (param_da.longitude<28.5))).any():
    
        plt.rcParams["font.weight"] = "bold"
        plt.rcParams["axes.labelweight"] = "bold"
        fig = plt.figure(figsize=(16, 10))
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
        ax.set_extent([19.5, 28.5, 34, 42], ccrs.PlateCarree())
        cmap = matplotlib.colormaps["jet"] # plt.cm.get_cmap("jet")
        cmap.set_over('darkred')

        ax.add_geometries(adm1_shapes, ccrs.PlateCarree(), linewidth=0.4,
                   edgecolor='black', facecolor='none', alpha=1)  
        cbar_kwargs = {'shrink':0.95, 'pad':0.03,'extend':'max','orientation': 'vertical',}

        img = param_da.plot.pcolormesh(ax=ax, vmin=0, vmax=1.2, cbar_kwargs=cbar_kwargs,cmap=cmap, transform=ccrs.PlateCarree(), x="longitude", y="latitude")
              
        gl = ax.gridlines(draw_labels=True, linestyle='--')
        gl.top_labels=False
        gl.right_labels=False
        gl.xformatter=LONGITUDE_FORMATTER
        gl.yformatter=LATITUDE_FORMATTER
        gl.xlabel_style={'size':14}
        gl.ylabel_style={'size':14}
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)
        plt.xlabel("Longitude", fontsize=16)
        plt.ylabel("Latitude", fontsize=16)
        ax.set_title(title, fontsize=20, pad=10.0, fontweight="bold")
        plt.savefig(f'Graphs/MODIS/{title_save}.png',dpi=300)
        plt.close()    