#Editor: Dr Logothetis Stavros-Andreas
#Affiliation: Laboratory of Atmospheric Physics, University of Patras
#Objective: AC-SAF Metop-A/B/C GOME-2 Level 2 (L2) data files processing
#Merging, Transforming, Visualizing daily L2 data 
#The script starts here
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.io.shapereader as shpreader
from matplotlib import pyplot as plt
import matplotlib.colors as clr
import cartopy.crs as ccrs
import matplotlib as mpl
import pandas as pd
import tqdm
import os

os.chdir('G:/Projects\Project_Forest_Fires_Evros\Forest_Fires')
directory = os.getcwd()

df = pd.read_csv(directory+'/Data/PurpleAir/PM2.5_Hourly.csv',parse_dates=True,index_col=['Time'])
df_info = pd.read_csv(directory+'/Data/PurpleAir/Stations Info.csv')

fname = 'G:/Projects/Project_Forest_Fires_Evros/Forest_Fires/Europe.shp' 
adm1_shapes = list(shpreader.Reader(fname).geometries())

colors = ['#5ed5ff','#92d14f','#FFFF00','#fc3903']
bounds = [0,10,20,25,50]
cmap = mpl.colors.ListedColormap(colors)
cmap.set_over('#990100')
cbar_kwargs = {'shrink':0.8, 'extend':'max','ticks': bounds,'orientation': 'vertical','label': '$\mathbf{μg/m^{3}}$'}
norm = clr.BoundaryNorm(bounds, ncolors=4)

for idx in tqdm.tqdm(df.index):
    date = str(idx)
    sname = date.split(' ')[0] + "_" + date.split(' ')[1].split(':')[0] + "_" + date.split(' ')[1].split(':')[1]
    plt.rcParams["font.weight"] = "bold"
    plt.rcParams["axes.labelweight"] = "bold"
    fig = plt.figure(figsize=(16, 10))
        
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
       
    ax.set_extent([19.5, 28.5, 34, 42], ccrs.PlateCarree())
    
    ax.add_geometries(adm1_shapes, ccrs.PlateCarree(), linewidth=0.4,
                   edgecolor='black', facecolor='none', alpha=1)
    
    gl = ax.gridlines(draw_labels=True, linestyle='--')
    gl.top_labels=False
    gl.right_labels=False
    gl.xformatter=LONGITUDE_FORMATTER
    gl.yformatter=LATITUDE_FORMATTER
    gl.xlabel_style={'size':20}
    gl.ylabel_style={'size':20}
    
    # plot pixel positions
    for col in df.columns:
        lat = df_info['Lat'][df_info['Name'] == col].values[0]
        lon = df_info['Long'][df_info['Name'] == col].values[0]
        img = ax.scatter(lon,lat,
            c=df.loc[idx,col],
            cmap=cmap,
            marker='o',norm=norm,
            s=200,zorder=2)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel("Longitude", fontsize=16)
    plt.ylabel("Latitude", fontsize=16)
    cbar = fig.colorbar(img, ax=ax, shrink=0.8, extend='max', ticks=bounds, orientation='vertical', fraction=0.04, pad=0.03)
    cbar.set_label('$PM_{2.5} (\mathbf{μg/m^{3}})$', fontsize=20)
    cbar.ax.tick_params(labelsize=20)
    ax.set_title(str(idx), fontsize=25, pad=10.0, fontweight="bold")
    plt.savefig(f'Graphs/PurpleAir/{sname}.png',dpi=150, bbox_inches='tight')
    plt.close()
