import pandas as pd
#import netCDF4 as nc
import numpy as np
import xarray as xr

#============== read data in ==================
mydata = pd.read_csv("./masks/CMIP6_NDVImax.csv",index_col=0)
#============== extract lonlat   ===================
def _unpackLatLon(index):
    lonST, latST = index.split(", ")
    return float(lonST[1:]), float(latST[:-1])
lonlat = np.array([_unpackLatLon(index) for index in mydata.index])
mydata["longitude"] = lonlat[:, 0]
mydata["latitude"]  = lonlat[:, 1]
mydata = mydata.reset_index(drop=True).set_index(["longitude", "latitude"])
lon    = lonlat[:, 0]
lat    = lonlat[:, 1]
lons   = lon[::60]
lats   = lat[::60]
#============== unstack var      ====================
ndvi = [[] for i in range(mydata.columns.shape[0])]
for i in range(mydata.columns.shape[0]):
  column = mydata.columns.values[i]
  ndvi[i] = mydata[column].unstack()
#============== 新建字典build dic for the data ======
nc_dict = {
    # nc文件的维度信息
    "dims": {"lat": lats.shape[0], "time": mydata.shape[1], "lon": lons.shape[0]},
    # nc文件的维度信息的坐标信息（lat,lon,time等）
    "coords": {
        "lat": {
            "dims": ("lat",),
            "attrs": {
                "standard_name": "latitude",
                "long_name": "Latitude",
                "units": "degrees_north",
                "axis": "Y",
            },
            "data":lats,
        },
        "lon": {
            "dims": ("lon",),
            "attrs": {
                "standard_name": "longitude",
                "long_name": "Longitude",
                "units": "degrees_east",
                "axis": "X",
            },
            "data":lons,
        },
        "time": {
            "dims": ("time",),
            "attrs": {"standard_name": "time", "long_name": "Time"},
            "data":pd.date_range(start='1982',periods=mydata.shape[1],freq='AS'),
        },
    },
    # nc文件中的变量
    "data_vars": {
        "ndvi": {
            "dims": ("time", "lat", "lon"),
            "attrs": {
                "long_name": "NDVImax",
                "var_desc": "max value of ndvi",
                "parent_stat": "Other",
            },
            "data":ndvi[0:33],
        }
    },
    # nc文件的全局属性信息
    "attrs": {
        "title": "NDVImax of CMIP6 historical scenario",
        "description": "Data calculated using parameters of TSSRESTREND",
    },
}

# 使用`from_dict`，将字典转化为xr.Dataset对象
ds = xr.Dataset.from_dict(nc_dict)
# 将xr.Dataset对象保存为nc文件
ds.to_netcdf("NDVImax_CMIP6.nc")









