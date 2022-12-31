#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 21 22:41:08 2022

@author: sambrandt
"""
# Import Libraries
# If you are stuck on how to use these, this may be able to help:
# https://docs.python.org/3/installing/index.html
# Note, there are sme libraries in here that you likely do not have already
import time
time0=time.time()
from datetime import datetime,timedelta
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np
from xarray import open_dataset
from xarray.backends import NetCDF4DataStore
from siphon.catalog import TDSCatalog
from netCDF4 import Dataset
import os
from scipy.interpolate import interp2d
from urllib.request import urlretrieve
from suntimes import SunTimes
# Inputs
central_lat=40 # In degrees
central_lon=-77 # In degrees
extent=3 # How far from the central coordinates the plot will go (in degrees)

script_directory='Your Directory Here' # Directory that you are
                                                    # running this script in
# Functions
def calc_latlon(data):
    # The math for this function was taken from 
    # https://makersportal.com/blog/2018/11/25/goes-r-satellite-latitude-and-
    # longitude-grid-projection-algorithm    
    x_1d = np.array(data['x'])*10**-6
    y_1d = np.array(data['y'])*10**-6
    x,y = np.meshgrid(x_1d,y_1d)
    goes_imager_projection=data.variables['fixedgrid_projection']
    r_eq=goes_imager_projection.semi_major_axis
    r_pol=goes_imager_projection.semi_minor_axis
    l_0=goes_imager_projection.longitude_of_projection_origin*(np.pi/180)
    h_sat=goes_imager_projection.perspective_point_height
    H=r_eq+h_sat
    a=np.sin(x)**2+(np.cos(x)**2*(np.cos(y)**2+(r_eq/r_pol)**2*np.sin(y)**2))
    b=-2*H*np.cos(x)*np.cos(y)
    c=H**2-r_eq**2
    r_s=(-b-(b**2-4*a*c)**0.5)/(2*a)
    print('^\nThis is expected behavior; the code is working as intended')
    s_x=r_s*np.cos(x)*np.cos(y)
    s_y=-r_s*np.sin(x)
    s_z=r_s*np.cos(x)*np.sin(y)
    lat=np.arctan((r_eq/r_pol)**2*(s_z/np.sqrt((H-s_x)**2+s_y**2)))*(180/np.pi)
    lon=(l_0-np.arctan(s_y/(H-s_x)))*(180/np.pi)
    return lon,lat,x_1d*h_sat,y_1d*h_sat
# Plot extent
north=central_lat+extent
south=central_lat-extent
east=central_lon+extent
west=central_lon-extent
# Determine whether to use visible or IR imagery
sun=SunTimes(longitude=central_lon,latitude=central_lat,altitude=0)
now=datetime.utcnow()
today=datetime(now.year,now.month,now.day)
sunrise=sun.riseutc(today)
sunset=sun.setutc(today)
if np.logical_and(now>sunrise,now<sunset)==True:
    channel='Channel02'
    wavelength='0.64'
    sat_type='Visible Red'
    lowerclip=0
    upperclip=1
    gamma=2
    cmap_sat='gray'
    print('\nDaytime at selected location, will use visible imagery'
          '\n(GOES-East ABI Channel 2, 0.64 μm Wavelength)')
elif np.logical_and(now>sunrise,now<sunset)==False:
    channel='Channel07'
    wavelength='3.9'
    sat_type='IR Shortwave'
    lowerclip=273.15-80
    upperclip=273.15+50
    gamma=1
    cmap_sat='gray_r'
    print('\nNighttime at selected location, will use infrared imagery'
          '\n(GOES-East ABI Channel 7, 3.9 μm Wavelength)')
time1=time.time()
# Download GOES data
sat_url='https://thredds.ucar.edu/thredds/catalog/satellite/goes/east/products/CloudAndMoistureImagery/CONUS/'+channel+'/current/catalog.xml'
sat_cat=TDSCatalog(sat_url)
sat_ds = sat_cat.datasets[-1]
sat_url=sat_ds.access_urls['HTTPServer']
tempfile=urlretrieve(sat_url,script_directory+'temp_sat.nc')
sat_data=Dataset(script_directory+'temp_sat.nc')
os.remove(script_directory+'temp_sat.nc')
sat_valid_datetime=datetime(int(sat_url[-17:-13]),1,1,int(sat_url[-10:-8]),int(sat_url[-8:-6]))+timedelta(int(sat_url[-13:-10])-1)
sat_valid=str(sat_valid_datetime)    
# Donwload RAP data
rap_url='https://thredds.ucar.edu/thredds/catalog/grib/NCEP/RAP/CONUS_13km/latest.xml'
best_rap=TDSCatalog(rap_url)
rap_ds=best_rap.datasets
ncss=rap_ds[0].subset()
query = ncss.query()
query.lonlat_box(north=north,south=south,east=east,west=west).time_range(sat_valid_datetime,sat_valid_datetime)
query.add_lonlat(value=True)
query.accept('netcdf4')
query.variables('Geopotential_height_isobaric','v-component_of_wind_height_above_ground','u-component_of_wind_height_above_ground','v-component_of_wind_isobaric','u-component_of_wind_isobaric','Pressure_surface','U-Component_Storm_Motion_height_above_ground_layer','V-Component_Storm_Motion_height_above_ground_layer')
rap_data=ncss.get_data(query)
rap_data=open_dataset(NetCDF4DataStore(rap_data))
rap_valid=ncss.metadata.time_span['begin'][0:10]+' '+str(sat_valid_datetime.hour).zfill(2)+':00:00'
hodo_lon=np.arange(west+0.1*extent,east,0.2*extent)
hodo_lat=np.arange(south+0.1*extent,north,0.2*extent)
# Time check #1
time2=time.time()
print('\nTime to download GOES and RAP from THREDDS:',time2-time1,'seconds\n')
# Retrieve, interpolate variables
lon=np.array(rap_data['lon'][0,:])
lat=np.array(rap_data['lat'][:,0])
p=np.arange(100,1025,25)*100
psfc=interp2d(lon,lat,np.array(rap_data['Pressure_surface'][0,:,:]))(hodo_lon,hodo_lat)
u_sfc=interp2d(lon,lat,np.array(rap_data['u-component_of_wind_height_above_ground'][0,0,:,:]))(hodo_lon,hodo_lat)
v_sfc=interp2d(lon,lat,np.array(rap_data['v-component_of_wind_height_above_ground'][0,0,:,:]))(hodo_lon,hodo_lat)
u_sm=interp2d(lon,lat,np.array(rap_data['U-Component_Storm_Motion_height_above_ground_layer'][0,0,:,:]))(hodo_lon,hodo_lat)
v_sm=interp2d(lon,lat,np.array(rap_data['V-Component_Storm_Motion_height_above_ground_layer'][0,0,:,:]))(hodo_lon,hodo_lat)
gpht=[]
u_iso=[]
v_iso=[]
for i in range(0,37):
    gpht.append(interp2d(lon,lat,np.array(rap_data['Geopotential_height_isobaric'][0,i,:,:]))(hodo_lon,hodo_lat))
    u_iso.append(interp2d(lon,lat,np.array(rap_data['u-component_of_wind_isobaric'][0,i,:,:]))(hodo_lon,hodo_lat))
    v_iso.append(interp2d(lon,lat,np.array(rap_data['v-component_of_wind_isobaric'][0,i,:,:]))(hodo_lon,hodo_lat))
gpht=np.array(gpht)
u_iso=np.array(u_iso)
v_iso=np.array(v_iso)
# Get surface height
gpht_sfc=np.zeros((10,10))
for i in range(0,10):
    for j in range(0,10):
        gpht_sfc[i,j]=np.interp([psfc[i,j]],p,gpht[:,i,j])[0]
        gpht[:,i,j]=gpht[:,i,j]-gpht_sfc[i,j]
# Gamma correction
sat=np.array(sat_data['Sectorized_CMI'])**(1/gamma)
sat=np.clip(sat,lowerclip,upperclip)
# Mask values not in the specified domain
# Speeds up plotting for later, especially with the higher-resolutino visible
# imagery
sat_lon,sat_lat,x_sat,y_sat=calc_latlon(sat_data)
sat_bool=np.logical_and(np.logical_and(sat_lat>=south-extent/2,sat_lat<=north+extent/2),np.logical_and(sat_lon>=west-extent/2,sat_lon<=east+extent/2))
sat=np.ma.masked_where(sat_bool==False,sat)
# Time check #2
time3=time.time()
print('\nTime to retrieve and perform calculations on variables:',time3-time2,'seconds')
# Create a tuple representing the bounds of the satellite data to be plotted
imextent_sat=(x_sat.min(),x_sat.max(),y_sat.min(),y_sat.max())
# Create a geostationary axis
globe_sat=ccrs.Globe(semimajor_axis=sat_data['fixedgrid_projection'].semi_major_axis,semiminor_axis=sat_data['fixedgrid_projection'].semi_minor_axis)
crs_sat=ccrs.Geostationary(central_longitude=sat_data['fixedgrid_projection'].longitude_of_projection_origin,globe=globe_sat)
# Create figure
# If you are running this through a terminal, you may need to add some code to
# save the figure to a specified directory
fig=plt.figure(dpi=500)
# Create cartopy axis
ax=plt.axes(projection=ccrs.PlateCarree())
# Make axis invisible
ax.axis('off')
# Set axis bounds
ax.set_extent([west,east,south,north])
# Add geographic features
ax.coastlines(linewidth=0.25,edgecolor='deepskyblue',zorder=1)
ax.add_feature(cfeature.STATES.with_scale('50m'),edgecolor='deepskyblue',linewidth=0.25,zorder=1)
# Display satellite imagery
im_sat=ax.imshow(sat,origin='upper',extent=imextent_sat,interpolation=None,transform=crs_sat,cmap=cmap_sat,vmin=lowerclip,vmax=upperclip,zorder=0)
# Create angles for ring plotting
angles=np.arange(0,2*np.pi+np.pi/32,np.pi/32)
# Loop through every hodograph axis
for i in range(0,10):
    for j in range(0,10):
        # Create hodograph axis
        axin = ax.inset_axes((j/10,i/10,0.1,0.1))
        # Make the axis invisible
        axin.axis('off')
        # Set axis limits (m/s)
        axin.set_xlim(-40,40)
        axin.set_ylim(-40,40)
        # Plot axes and rings
        axin.plot(10*np.cos(angles),10*np.sin(angles),color='black',lw=0.1,zorder=0)
        axin.plot(20*np.cos(angles),20*np.sin(angles),color='black',lw=0.1,zorder=0)
        axin.plot(30*np.cos(angles),30*np.sin(angles),color='black',lw=0.1,zorder=0)
        axin.plot([0,0],[-30,30],color='black',lw=0.1,zorder=0)
        axin.plot([-30,30],[0,0],color='black',lw=0.1,zorder=0)
        # Create the interpolated u,v,z arrays for plotting
        hgts=np.flip(gpht[:,i,j][gpht[:,i,j]>=0])
        u_plot=np.flip(u_iso[:,i,j][gpht[:,i,j]>=0])
        v_plot=np.flip(v_iso[:,i,j][gpht[:,i,j]>=0])
        hgts=np.insert(hgts,0,0)
        u_plot=np.insert(u_plot,0,u_sfc[i,j])
        v_plot=np.insert(v_plot,0,v_sfc[i,j])
        hgt_plot=np.concatenate((np.arange(0,1000,100),np.arange(1000,3000,500),np.arange(3000,13000,1000)),axis=0)
        u_plot=np.interp(hgt_plot,hgts,u_plot)
        v_plot=np.interp(hgt_plot,hgts,v_plot)
        # Plot the hodograph
        axin.plot(u_plot[np.logical_and(hgt_plot<=12000,hgt_plot>=9000)],v_plot[np.logical_and(hgt_plot<=12000,hgt_plot>=9000)],color='cyan',lw=0.5,clip_on=False,zorder=1)
        axin.plot(u_plot[np.logical_and(hgt_plot<=9000,hgt_plot>=6000)],v_plot[np.logical_and(hgt_plot<=9000,hgt_plot>=6000)],color='yellow',lw=0.5,clip_on=False,zorder=1)
        axin.plot(u_plot[np.logical_and(hgt_plot<=6000,hgt_plot>=3000)],v_plot[np.logical_and(hgt_plot<=6000,hgt_plot>=3000)],color='lime',lw=0.5,clip_on=False,zorder=1)
        axin.plot(u_plot[np.logical_and(hgt_plot<=3000,hgt_plot>=1000)],v_plot[np.logical_and(hgt_plot<=3000,hgt_plot>=1000)],color='red',lw=0.5,clip_on=False,zorder=1)
        axin.plot(u_plot[np.logical_and(hgt_plot<=1000,hgt_plot>=0)],v_plot[np.logical_and(hgt_plot<=1000,hgt_plot>=0)],color='magenta',lw=0.5,clip_on=False,zorder=1)
        # Plot the bunkers storm motion
        axin.scatter(u_sm[i,j],v_sm[i,j],s=1,c='white',linewidths=0,edgecolors=None,clip_on=False,zorder=2)
# Title text (needlessly complex because matplotlib doesn't think the ability
# to make multi-colored title objects is a worthwhile addition)
ax.set_title('GOES-16 '+wavelength+' µm ('+sat_type+'), Valid '+sat_valid+' UTC\nRAP Analysis Hodographs, Valid '+rap_valid+' UTC\n',fontsize=7)
ax.text(east-2*extent*0.04,north+extent/10,'km',fontsize=7,color='black',ha='right',va='top')
ax.text(east-5*extent*0.04,north+extent/10,'9-12',fontsize=7,color='cyan',ha='right',va='top')
ax.text(east-9*extent*0.04,north+extent/10,',',fontsize=7,color='black',ha='right',va='top')
ax.text(east-10*extent*0.04,north+extent/10,'6-9',fontsize=7,color='yellow',ha='right',va='top')
ax.text(east-13*extent*0.04,north+extent/10,',',fontsize=7,color='black',ha='right',va='top')
ax.text(east-14*extent*0.04,north+extent/10,'3-6',fontsize=7,color='lime',ha='right',va='top')
ax.text(east-17*extent*0.04,north+extent/10,',',fontsize=7,color='black',ha='right',va='top')
ax.text(east-18*extent*0.04,north+extent/10,'1-3',fontsize=7,color='red',ha='right',va='top')
ax.text(east-21*extent*0.04,north+extent/10,',',fontsize=7,color='black',ha='right',va='top')
ax.text(east-22*extent*0.04,north+extent/10,'0-1',fontsize=7,color='magenta',ha='right',va='top')
ax.text(west+2*extent*0.04,north+extent/10,'10 m/s Rings, BRM in White,',fontsize=7,color='black',ha='left',va='top')
ax.text(east,south,'Created By Sam Brandt (SamBrandtMeteo on GitHub)',fontsize=2,color='white',ha='right',va='bottom')
# Time check #3
time4=time.time()
print('\nTime to generate plot:',time4-time3,'seconds')
# Total runtime check
print('\nTotal script runtime:',time4-time0,'seconds')






