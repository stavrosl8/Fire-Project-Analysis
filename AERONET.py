from datetime import datetime
import pandas as pd
import numpy as np
import tqdm
import glob

def read_file(station, product, ftype, columns, date_col):
    file = glob.glob(f'{product}/*{station}.{ftype}')[0]
    df = pd.read_csv(file,  parse_dates={'datetime':date_col},
                                  date_parser=dateparse, usecols=columns, skiprows=6, na_values=-999, sep=",")
    df.index = pd.to_datetime(df['datetime'])
    df = df.drop(columns=['datetime'])
    df = df.dropna(axis=0, how ='all')
    df = df.dropna(axis=1, how ='all')
    return(df)

#Read Files for Direct Sun Algorithm

AOD_columns = ['Date(dd:mm:yyyy)', 'Time(hh:mm:ss)','AOD_1020nm',
       'AOD_870nm', 'AOD_675nm', 'AOD_500nm', 'AOD_440nm', 'AOD_380nm',
       'AOD_340nm', 'Precipitable_Water(cm)','440-870_Angstrom_Exponent',
       'Site_Latitude(Degrees)', 'Site_Longitude(Degrees)',
       'Site_Elevation(m)', 'Solar_Zenith_Angle(Degrees)', 'Ozone(Dobson)',
       'NO2(Dobson)']

oneil_columns = ['Date_(dd:mm:yyyy)', 'Time_(hh:mm:ss)','Total_AOD_500nm[tau_a]',
                 'Fine_Mode_AOD_500nm[tau_f]','Coarse_Mode_AOD_500nm[tau_c]', 'FineModeFraction_500nm[eta]']

SSA_columns = ['Date(dd:mm:yyyy)', 'Time(hh:mm:ss)', 'Single_Scattering_Albedo[440nm]',
               'Single_Scattering_Albedo[675nm]','Single_Scattering_Albedo[870nm]',
                  'Single_Scattering_Albedo[1020nm]','Surface_Albedo[440m]',
                  'Surface_Albedo[675m]', 'Surface_Albedo[870m]','Surface_Albedo[1020m]']

VOL_columns = ['Date(dd:mm:yyyy)','Time(hh:mm:ss)', 'VolC-T','REff-T', 'VMR-T',
               'VolC-F', 'REff-F', 'VMR-F', 'VolC-C', 'REff-C', 'VMR-C']

SIZ_columns = ['Date(dd:mm:yyyy)','Time(hh:mm:ss)','0.050000','0.065604',
               '0.086077','0.112939','0.148184','0.194429','0.255105','0.334716',
               '0.439173','0.576227','0.756052','0.991996','1.301571','1.707757',
               '2.240702','2.939966','3.857452','5.061260','6.640745','8.713145',
               '11.432287','15.000000']

TAB_columns = ['Date(dd:mm:yyyy)','Time(hh:mm:ss)','Absorption_AOD[440nm]', 'Absorption_AOD[675nm]',
               'Absorption_AOD[870nm]', 'Absorption_AOD[1020nm]']

RIN_columns = ['Date(dd:mm:yyyy)','Time(hh:mm:ss)','Refractive_Index-Real_Part[440nm]',
'Refractive_Index-Real_Part[675nm]','Refractive_Index-Real_Part[870nm]','Refractive_Index-Real_Part[1020nm]',
'Refractive_Index-Imaginary_Part[440nm]','Refractive_Index-Imaginary_Part[675nm]','Refractive_Index-Imaginary_Part[870nm]',
'Refractive_Index-Imaginary_Part[1020nm]']

ASY_columns = ['Date(dd:mm:yyyy)','Time(hh:mm:ss)', 'Asymmetry_Factor-Total[440nm]',
               'Asymmetry_Factor-Total[675nm]','Asymmetry_Factor-Total[870nm]',
               'Asymmetry_Factor-Total[1020nm]','Asymmetry_Factor-Fine[440nm]',
               'Asymmetry_Factor-Fine[675nm]','Asymmetry_Factor-Fine[870nm]',
               'Asymmetry_Factor-Fine[1020nm]','Asymmetry_Factor-Coarse[440nm]',
               'Asymmetry_Factor-Coarse[675nm]','Asymmetry_Factor-Coarse[870nm]',
               'Asymmetry_Factor-Coarse[1020nm]']

oneil_date = ['Date(dd:mm:yyyy)', 'Time(hh:mm:ss)']
rest_date = ['Date_(dd:mm:yyyy)', 'Time_(hh:mm:ss)']

dateparse = lambda x: datetime.strptime(x, '%d:%m:%Y %H:%M:%S')

files_AOD = glob.glob('direct_product/*.lev15')

for station in tqdm.tqdm(files_AOD):
   
    file_name = station.split('_')[-1][:-6]
           
    df_AOD = read_file(file_name, 'direct_product', 'lev15', AOD_columns, rest_date)
    df_oneil = read_file(file_name, 'direct_product', 'ONEILL_lev15', oneil_columns, oneil_date)
    df_SSA = read_file(file_name, 'inversion_product', 'ssa', SSA_columns, rest_date)
    df_VOL = read_file(file_name, 'inversion_product', 'vol', VOL_columns, rest_date)
    df_SIZ = read_file(file_name, 'inversion_product', 'siz', SIZ_columns, rest_date)
    df_TAB = read_file(file_name, 'inversion_product', 'tab', TAB_columns, rest_date)
    df_RIN = read_file(file_name, 'inversion_product', 'rin', RIN_columns, rest_date)
    df_ASY = read_file(file_name, 'inversion_product', 'asy', ASY_columns, rest_date)

    df_st_DSA = pd.concat([df_AOD, df_oneil], axis=1)
    df_st_INV = pd.concat([df_SSA, df_VOL, df_SIZ, df_TAB, df_RIN, df_ASY], axis=1)

    df_merged = pd.merge_asof(df_st_DSA, df_st_INV, on="datetime", direction='nearest', tolerance=pd.Timedelta('3min'))
    df_merged.index = df_merged['datetime']
    df_merged = df_merged.resample('5min').mean().dropna(axis=0,how='all')

df_DSA = pd.concat(DSA,axis=1)


df_merged.to_csv(file_name+'.csv')
