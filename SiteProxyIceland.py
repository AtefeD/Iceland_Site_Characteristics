"""
This script creates a databse of site proxies and proxy-based site amplification values for stations in SISZ, Iceland from 

Local studies: 
# 1- our proposed site amplification units and the corresponding freq-dependent site amplification values (Darzi et al. 2024, SDEE)
# 2- empirical original amplification values from local studies (Rahpeyma et al. 2023, Rea23)
# 3- geology-based local Vs30 estimates [Vogfjörd et al. (2010), Vea10]

Large-scale, non-Icelandic, studies: 
# 1- site proxies e.g., slope-based Vs30 from USGS and ESRM20, Slope proxy, Geomorphological Sedimentray thickness (GST) proxy
# 2- frequency-dependent European proxy-based site amplification models, i.e., ESRM20 (Weatheril et al. 2023) and Loviknes23


first download "coefficient_table.csv" file from "https://zenodo.org/records/10686867" and put it in the root folder. 

*** References:
    1- Darzi 
    2- Rahpeyma, 
    
If you use this code, please cite Darzi et al. (2024). 
   

@author: Atefe Darzi (atefe@hi.is) 
Date: 06.2023- GFZ
Last update: 04.2024 at UICE-IMO

"""

import rasterio
import rasterio.features
from shapely.geometry import Point, shape
import geopandas as gpd
import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
import scipy.io as sio
import os
os.environ['USE_PYGEOS'] = '0'
import geopandas as gpd


import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter




# %% local site proxies for SISZ stations from: 
# 1- our proposed site amplification units and values
# 2- other local studies (Rea23)
# 3- geology-based local Vs30 estimates [Vea10]

# Define your directory path 
path = os.getcwd()

IceSPfile_path = os.path.join(path, '.\LocalSiteProxies_SISZstations.csv')
SP_df0 = pd.read_csv(IceSPfile_path)

SP_df0['Soil_type'] #--> traditional site classes in Iceland
SP_df0['SAFER_unit'] #--> geoological units associated with local Vs30 estimates, see fig3 of Darzi et al. 2024
SP_df0['SAFER_Vs30'] #-->  geological-based Vs30 estimate (Vogfjörd et al. (2010)), see fig3 of Darzi et al. 2024


SP_df0['GeoSiteAmp'] #--> the proposed geology-based site amplification units, see fig5 of Darzi et al. 2024
# 'HardRock': ['IS108', 'IS302', 'IS305', 'IS306', 'IS307', 'IS309'],
# 'Soil': ['IS103', 'IS105', 'IS109', 'IS610', 'IS304'],
# 'Lava': ['IS104', 'IS111', 'IS112', 'IS501', 'IS502', 'IS601', 'IS602', 'IS603', 'IS604', 'IS605', 'IS607', 'IS608', 'IS611'],
# 'Rock': ['IS101', 'IS102', 'IS106', 'IS107', 'IS113', 'IS301', 'IS303', 'IS401', 'IS402']

                        
# observed empirical δS2S (posterior median) at SISZ station locations in log10 units (Red star in Fig 11 of Darzi et al. 2024 corresponding to 2Hz)
# 1hz, 2hz, 5hz, 7hz, 10hz, 30hz, pga --> values correspond to periods of 1.07, 0.5187, 0.2157, 0.134, 0.1038, 0.0322, and 0 sec in Rea23
SP_df0['Ds2s_p50_1hz']
SP_df0['Ds2s_p50_2hz']
SP_df0['Ds2s_p50_5hz']
SP_df0['Ds2s_p50_7hz'] 
SP_df0['Ds2s_p50_10hz']
SP_df0['Ds2s_p50_30hz'] 
SP_df0['Ds2s_p50_pga']  


# average site amplification predictions corresponding to station location's key geological units 
# (black square in Fig 11 of Darzi et al. 2024, corresponding to 2Hz)
SP_df0['R23_p50_1hz']
SP_df0['R23_p50_2hz']
SP_df0['R23_p50_5hz']
SP_df0['R23_p50_7hz']
SP_df0['R23_p50_10hz']
SP_df0['R23_p50_30hz'] 
SP_df0['R23_p50_pga'] 

SP_df0['longitudes'] = SP_df0['st_Lon'].copy()
SP_df0['latitudes'] = SP_df0['st_Lat'].copy()


# %% read raster maps (Vs30 and site amplifictaion)
## Inferred VS30 from USGS-WA07
# edit the path according to your defined address
Vs30WA07tif = rasterio.open('D:\OneDrive\\OneDrive - Menntaský\\Atefe PC\\Iceland-GIS files\\usgs Vs30\\vs30\\vs30slope.tif') 
# if not works, unzip and use vs30.zip files: .\Iceland-GIS files\usgs Vs30\vs30\
print(Vs30WA07tif.bounds) # (minx, miny, maxx, maxy) tuple
dt_Vs30WA07 = Vs30WA07tif.read(1)

## slope-based inferredVs30 by Weatherill et al. 2023 
# edit the path according to your defined address
Vs30W23tif = rasterio.open('D:\OneDrive\\OneDrive - Menntaský\\Atefe PC\\Iceland-GIS files\W23_eSRm20\\European_Site_Model_Data\\layers\\Vs30 from slope (m_s).tif')
print(Vs30W23tif.bounds)
print(Vs30W23tif.nodata)
dt_Vs30W23 = Vs30W23tif.read(1)

## slope-based site amplification- Weatherill et al. 2023 - ref_slope=0.3m/m, ref_geo=Camberian
# edit the path according to your defined address
SAw23_path = 'D:\OneDrive\OneDrive - Menntaský\\Atefe PC\\Iceland-GIS files\\W23_eSRm20\\amplification_models\\'
file_SAw23_30hz = rasterio.open(SAw23_path + 'amplification_0.033.tif')
print(file_SAw23_30hz.bounds)  
file_SAw23_10hz = rasterio.open(SAw23_path + 'amplification_0.1.tif')
print(file_SAw23_10hz.bounds)
file_SAw23_7hz = rasterio.open(SAw23_path + 'amplification_0.14.tif')
print(file_SAw23_7hz.bounds)
file_SAw23_5hz = rasterio.open(SAw23_path + 'amplification_0.2.tif')
print(file_SAw23_5hz.bounds)
file_SAw23_2hz = rasterio.open(SAw23_path + 'amplification_0.5.tif')
print(file_SAw23_2hz.bounds)
file_SAw23_1hz = rasterio.open(SAw23_path + 'amplification_1.0.tif')
print(file_SAw23_1hz.bounds)
file_SAw23_pga = rasterio.open(SAw23_path + 'amplification_pga.tif')
print(file_SAw23_pga.bounds)

print(file_SAw23_pga.nodata)
print(file_SAw23_5hz.crs) #coord sys EPSG
print(file_SAw23_5hz.nodata) 
print(file_SAw23_5hz.shape) 

dt_SAw23_30hz = file_SAw23_30hz.read(1)
dt_SAw23_10hz = file_SAw23_10hz.read(1)
dt_SAw23_7hz = file_SAw23_7hz.read(1)
dt_SAw23_5hz = file_SAw23_5hz.read(1)
dt_SAw23_2hz = file_SAw23_2hz.read(1)
dt_SAw23_1hz = file_SAw23_1hz.read(1)
dt_SAw23_pga = file_SAw23_pga.read(1)

## slope-GEBCO14 in Weatherill et al. 2023 
# edit the path according to your defined address
slopeW23tif = rasterio.open('D:\OneDrive\OneDrive - Menntaský\Atefe PC\\Iceland-GIS files\\W23_eSRm20\\European_Site_Model_Data\\layers\\Slope (m_m).tif')
print(slopeW23tif.bounds)
print(slopeW23tif.nodata)
dt_W23slop = slopeW23tif.read(1)

## Sediment thickness
# edit the path according to your defined address
regolith_path = 'D:\OneDrive\\OneDrive - Menntaský\\Atefe PC\\Iceland-GIS files\\global geomorph\\Global_Soil_Regolith_Sediment_1304\\data\\'
file_sedt = rasterio.open(regolith_path + 'average_soil_and_sedimentary-deposit_thickness.tif')
print(file_sedt.bounds)  
dt_sedt = file_sedt.read(1)

# %% Extract the site proxies from large-scale studies at 34 stations in SISZ, Iceland
lats_st, longs_st = [], []
df1_SAw23_30hz, df1_SAw23_10hz, df1_SAw23_7hz ,df1_SAw23_5hz, df1_SAw23_2hz, df1_SAw23_1hz, df1_SAw23_pga= [], [], [], [], [], [], [] 
df1_Vs30WA07, df1_W23slop, df1_W23slopDeg, df1_Vs30W23, df1_sedt  = [], [], [], [], []
lat_st, lon_st =  SP_df0["st_Lat"].values, SP_df0["st_Lon"].values

for lat, long in zip(lat_st, lon_st): # for station locations 
    try:
        row7, col7 = Vs30WA07tif.index(long, lat)
        df1_Vs30WA07.append(dt_Vs30WA07[row7, col7])
    except:
        df1_Vs30WA07.append(np.nan)
    
    try:
        row8, col8 = slopeW23tif.index(long, lat)
        df1_W23slop.append(dt_W23slop[row8, col8])
        df1_W23slopDeg.append(math.degrees(dt_W23slop[row8, col8]))
    except:
        df1_W23slop.append(np.nan)
        df1_W23slopDeg.append(np.nan)
           
    try:
        row9, col9 = Vs30W23tif.index(long, lat)
        df1_Vs30W23.append(dt_Vs30W23[row9, col9])
    except:
        df1_Vs30W23.append(np.nan) 

    try:
        row10, col10 = file_sedt.index(long, lat)
        df1_sedt.append(dt_sedt[row10, col10])
    except:
        df1_sedt.append(np.nan)          

    row1, col1 = file_SAw23_30hz.index(long, lat)
    df1_SAw23_30hz.append(dt_SAw23_30hz[row1, col1])       
    row2, col2 = file_SAw23_10hz.index(long, lat)
    df1_SAw23_10hz.append(dt_SAw23_10hz[row2, col2])
    row3, col3 = file_SAw23_7hz.index(long, lat)
    df1_SAw23_7hz.append(dt_SAw23_7hz[row3, col3])
    row4, col4 = file_SAw23_5hz.index(long, lat)
    df1_SAw23_5hz.append(dt_SAw23_5hz[row4, col4])
    row4, col4 = file_SAw23_2hz.index(long, lat)
    df1_SAw23_2hz.append(dt_SAw23_2hz[row4, col4])
    row5, col5 = file_SAw23_1hz.index(long, lat)
    df1_SAw23_1hz.append(dt_SAw23_1hz[row5, col5])
    row6, col6 = file_SAw23_pga.index(long, lat)
    df1_SAw23_pga.append(dt_SAw23_pga[row6, col6])        

    lats_st.append(lat)
    longs_st.append(long)       
# len(df1_Vs30WA07) # == len(lat_st)*len(lon_st)           

file_SAw23_30hz.close()
file_SAw23_10hz.close()
file_SAw23_7hz.close()
file_SAw23_5hz.close()
file_SAw23_2hz.close()
file_SAw23_1hz.close()
file_SAw23_pga.close()
Vs30WA07tif.close()
Vs30W23tif.close()
slopeW23tif.close()
file_sedt.close()

# %% make a dataframe of available site proxies for Iceland

SP_df1 = pd.DataFrame({'longitudes':longs_st, 'latitudes':lats_st, 
                       'Vs30_usgs':df1_Vs30WA07,
                       'slope_W23':df1_W23slop,
                       'slope_W23_degrees':df1_W23slopDeg,
                       'Vs30_W23':df1_Vs30W23,
                       'sediment_thickness':df1_sedt,
                       
                       'site_amplification_W23_30hz':df1_SAw23_30hz,
                       'site_amplification_W23_10hz':df1_SAw23_10hz,
                       'site_amplification_W23_7hz':df1_SAw23_7hz, 
                       'site_amplification_W23_5hz':df1_SAw23_5hz, 
                       'site_amplification_W23_2hz':df1_SAw23_2hz,
                       'site_amplification_W23_1hz':df1_SAw23_1hz,
                       'site_amplification_W23_pga':df1_SAw23_pga })
SP_df= pd.merge(SP_df0, SP_df1, on=['longitudes', 'latitudes'], how='inner')
   
# modify the unvalid values to np.nan
SP_df.replace(-1, np.nan, inplace=True) 
SP_df.replace(-99.0, np.nan, inplace=True)
SP_df.replace('None', np.nan, inplace=True)
SP_df.replace('none', np.nan, inplace=True)
SP_df.replace('NAN', np.nan, inplace=True)
SP_df['sediment_thickness'].replace(file_sedt.nodata, np.nan, inplace=True)
SP_df['Vs30_W23'].replace(Vs30W23tif.nodata , np.nan, inplace=True)
SP_df['Vs30_usgs'].replace(Vs30WA07tif.nodata , np.nan, inplace=True)

print("count Null values of each column: ", SP_df.isnull().sum())

# check if any station is located on ocean as per Vs30-usgs data
ocean = np.full_like(SP_df['Vs30_usgs'], False, dtype=bool)
ocean[SP_df['Vs30_usgs'].astype(str)=='nan'] = True
SP_df['ocean_water'] = ocean



# %% European Geological era -Vilanovaet al 2018- add ERA0 attribute to our  dataframe
# edit the path according to your defined address
gdf_EuGeo_poly = gpd.read_file("D:\OneDrive\OneDrive - Menntaský\Atefe PC\Iceland-GIS files\Vilanova_Geo_W23\GEOL_V8_ERA2.shp")
gdf_EuGeo_poly = gdf_EuGeo_poly.to_crs(4326)  # transfer coordinate system

# converting longitude & latitude to geometry
SP_df['coordinates'] = list(zip(SP_df['longitudes'], SP_df['latitudes']))
SP_df.coordinates = SP_df.coordinates.apply(Point)
gdf_points = gpd.GeoDataFrame(SP_df, geometry='coordinates', crs=4326)

sjoin_EuGeo = gpd.sjoin(gdf_points, gdf_EuGeo_poly, how='left') #add EuGeo 
print(sjoin_EuGeo)
SP_df2 = pd.DataFrame(sjoin_EuGeo) # gdf to df

# remove unnecessary attributes 
SP_df = SP_df2.drop(['LTH_STR','Litho','Det_Strati','Litho_Simp',"Strati","System","Code_Syste","Syst_Litho","ERA1","ERA2"], axis='columns')
SP_df.keys().tolist()


# build EuGeo_ID : HOLOCENE= 0, CENOZOIC= 2, UNKNOWN= 7 for iceland
eraNum = np.full_like(SP_df.ERA0, np.nan, dtype=float)
for ie, era in enumerate(['HOLOCENE','PLEISTOCENE','CENOZOIC', 
                          'CRETACEOUS', 'JURASSIC-TRIASSIC','PRECAMBRIAN',
                          'PALEOZOIC', 'UNKNOWN']):
    eraNum[SP_df.ERA0 == era] = ie
SP_df['EuGeo_ID'] = eraNum


#  Vs-groups of Weatherill et al. 2023
print(SP_df[ ( SP_df['Vs30_W23']== 2)]['longitudes'].count()) # count values in each Vs30-class
# dict_Vs30_W23 = {2:[180, 239], 3:[240, 299], 4:[300, 359], 5:[360, 489], 6:[490, 619], 7:[620, 759] , 8: [760, 2500] }


    
        

# %% proxy-based predictions by Loviknes et al. 2024
extent = [-24.9, -13.1, 63.3, 66.6] # [llon, ulon, llat, ulat]
# Edit the path to the coefficient table file of Loviknes models 
coeff_table = pd.read_csv('.\coefficient_table.csv')
proxies = ['Vs30_WA', 'slope', 'GST1', 'Geology_slope']
proxies_label = ['Inferred $V_{\mathrm{S30}}$ (WA07)', 'Slope (Wea23)', 'GST', 'Geology and slope'] # subplot title

SP_df['Vs30_WA'] = SP_df['Vs30_usgs'] .copy() #Vs30_W23  
SP_df['slope'] = SP_df['slope_W23'] .copy() 
SP_df['GST1'] = SP_df['sediment_thickness'].copy() + 1
SP_df['slope_log'] = np.log(SP_df['slope_W23'].copy())

freqsP = ['1.062','1.991', '4.929','6.987','9.903'] # see coeff_table
freqsP_label = ['1','2', '5','7','10']

for t, tL in zip(freqsP, freqsP_label):
    fig, ax = plt.subplots(2,2, figsize=(35, 15),
              subplot_kw={'projection': ccrs.PlateCarree()})
    for i, j, pr, prL in zip([0, 0, 1, 1], [0, 1, 0, 1], proxies, proxies_label):
        b_lin = coeff_table[(coeff_table.Proxy==pr) & (coeff_table.Frequency.map('{:.3f}'.format)==t)]['b_lin'].values[0]
        b_lin_sd = coeff_table[(coeff_table.Proxy==pr) & (coeff_table.Frequency.map('{:.3f}'.format)==t)]['b_lin_std'].values[0]
           
        if pr=='Geology_slope':    
            pr1, pr2='slope_log' , 'ERA0'           
            xp = pd.get_dummies(data = SP_df[[pr1, pr2]])                    
            a_lin, a_lin_sd = [], []            
            for cla in xp.columns.values:
                if cla != 'slope_log':
                    clae = cla.replace('_','.').replace('-','.') #a_linslope_log vs a_linERA0.CENOZOIC
                else:
                    clae = cla
                a_linC = coeff_table[(coeff_table.Proxy==pr) & (coeff_table.Frequency.map('{:.3f}'.format)==t)]['a_lin'+clae].values[0]
                a_linCsd = coeff_table[(coeff_table.Proxy==pr) & (coeff_table.Frequency.map('{:.3f}'.format)==t)]['a_lin'+clae+'_std'].values[0]
                if np.isnan(a_linC):
                    a_linC = 0
                    a_linCsd = 0
                a_lin.append(a_linC)
                a_lin_sd.append(a_linCsd)
    
            amp_pred = np.dot(a_lin, xp.T.to_numpy()) + b_lin # xp --> slope_log
            
            a_lin_up= [x + y for x, y in zip(a_lin, a_lin_sd)]
            a_lin_lo= [x - y for x, y in zip(a_lin, a_lin_sd)]
            b_lin_up= b_lin + b_lin_sd
            b_lin_lo= b_lin - b_lin_sd
            
            amp_pred_up = np.dot(a_lin_up , xp.T.to_numpy()) + b_lin_up
            amp_pred_lo = np.dot(a_lin_lo , xp.T.to_numpy()) + b_lin_lo
            
            amp_pred[np.isnan(SP_df['EuGeo_ID'])]=np.nan  
            amp_pred[(SP_df['EuGeo_ID']==7)] = np.nan
            amp_pred_up[np.isnan(SP_df['EuGeo_ID'])]=np.nan  
            amp_pred_up[(SP_df['EuGeo_ID']==7)] = np.nan
            amp_pred_lo[np.isnan(SP_df['EuGeo_ID'])]=np.nan  
            amp_pred_lo[(SP_df['EuGeo_ID']==7)] = np.nan
            
            # np.min( epe_df['Lin_pred_Geology_slope_1.062']  )
        else:
            pr1=pr
            a_lin = coeff_table[(coeff_table.Proxy==pr) & (coeff_table.Frequency.map('{:.3f}'.format)==t)]['a_lin'].values[0]
            a_lin_sd = coeff_table[(coeff_table.Proxy==pr) & (coeff_table.Frequency.map('{:.3f}'.format)==t)]['a_lin_std'].values[0]
            
            xp = SP_df[pr]
            amp_pred = a_lin * np.log(xp) + b_lin #aln(x)+b
            a_lin_up= a_lin + a_lin_sd
            b_lin_up= b_lin + b_lin_sd
            a_lin_lo= a_lin - a_lin_sd
            b_lin_lo= b_lin - b_lin_sd
            
            amp_pred_up = a_lin_up * np.log(xp) + b_lin_up 
            amp_pred_lo = a_lin_lo * np.log(xp) + b_lin_lo          
            
        # have all ds2s in log10 unit
        SP_df['Lov23_' + pr + '_' + str(tL)] = np.log10(np.exp(amp_pred)) #:[0.5-3.5] 
        SP_df['Lov23up_' + pr + '_' + str(tL)] = np.log10(np.exp(amp_pred_up)) 
        SP_df['Lov23lo_' + pr + '_' + str(tL)] = np.log10(np.exp(amp_pred_lo)) 
        
        ax[i,j].coastlines(resolution='10m')
        ax[i,j].set_extent(extent, ccrs.PlateCarree()) 
        
        gl = ax[i,j].gridlines(draw_labels=True, linewidth=1, linestyle='--')
        gl.top_labels = False
        gl.right_labels = False
        gl.xformatter = LongitudeFormatter()
        gl.yformatter = LatitudeFormatter()
      
        sc = ax[i,j].scatter(SP_df['longitudes'], 
                        SP_df['latitudes'], 
                        c = np.exp(amp_pred), # SP_df['Lov23_' + pr + '_' + str(tL)]
                        cmap='seismic', s=100, zorder=2, 
                        linewidths=1, edgecolors='k', alpha=1,
                        marker='^', label='SISZ stations')
        sc.set_clim(0,3)
        cax1 = ax[i,j].inset_axes([1.01, 0.05, 0.05, 0.9])
        cbar = fig.colorbar(sc, cax=cax1, ticks=[0.5, 0.75, 1, 1.25,1.5, 1.75, 2, 2.25, 2.5, 2.75, 3])
        cbar.mappable.set_clim(0, 3)
        cbar.ax.tick_params(labelsize=12)

        # cbar.set_label('Predicted exp(' + r'$\delta \mathrm{S2S}_{s} $' + ')', fontsize = 15)
        ax[i,j].set_title('Predicted exp(' + r'$\delta \mathrm{S2S}_{s} $ ) using ' + prL, fontsize=30, loc='left')  
        
    ax[0,0].set_title('f = ' + t + ' Hz\n Predicted exp(' + r'$\delta \mathrm{S2S}_{s} $ ) using ' + proxies_label[0] , fontsize=30, loc='left')  
    fig.subplots_adjust(wspace=0.15)
    fig.subplots_adjust(hspace=0.05)
    plt.savefig(".\\Stations_Lea23_siteamp_f" + t + "Hz.png", orientation='landscape')


# SP_df = SP_df.dropna(subset=['Ds2s_p50_1hz'], how='any', axis=0)
       
SP_df.to_csv('SiteProxy_SiteAmp_DB.csv')


