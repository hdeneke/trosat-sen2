# Python standard library imports
import enum
import math
import os
import re
import zipfile

# SciPy importsR
import numpy as np
import xarray as xr

# OSGeo imports
from osgeo import gdal, ogr, osr
import pyproj

# misc. function or class imports
from functools import cached_property
from addict import Dict as adict
from lxml import etree as et

# package-local imports
from . import tilepar

# enable GDAL/OGR/OSR exceptions to silence warning
gdal.UseExceptions()
ogr.UseExceptions()
osr.UseExceptions()

# Mapping for strings to GDAL resampling algorithm IDs
resample_map = {
    'nearest'  : gdal.GRA_NearestNeighbour,
    'average'  : gdal.GRA_Average,
    'bilinear' : gdal.GRA_Bilinear,
    'cubic'    : gdal.GRA_Cubic,
    'cspline'  : gdal.GRA_CubicSpline,
    'lanczos'  : gdal.GRA_Lanczos,
    'min'      : gdal.GRA_Min,
    'max'      : gdal.GRA_Max,
    'mode'     : gdal.GRA_Mode,
    'median'   : gdal.GRA_Med,
    "Q1"       : gdal.GRA_Q1,
    "Q3"       : gdal.GRA_Q3
}


class _enum_meta(enum.EnumMeta):
    def __getitem__(self, item):
        if hasattr(self, "norm"): 
            item = self.norm(item)
        return super().__getitem__(item)


class l1c_band(enum.IntEnum, metaclass=_enum_meta):

    B01 =  0, 60.0, (2,1)
    B02 =  1, 10.0, (0,3)
    B03 =  2, 10.0, (0,2)
    B04 =  3, 10.0, (0,1)
    B05 =  4, 20.0, (1,1)
    B06 =  5, 20.0, (1,2)
    B07 =  6, 20.0, (1,3)
    B08 =  7, 10.0, (0,4)
    B8A =  8, 20.0, (1,4)    
    B09 =  9, 60.0, (2,2)    
    B10 = 10, 60.0, (2,3)    
    B11 = 11, 20.0, (1,5)    
    B12 = 12, 20.0, (1,6)    

    def __new__(cls, val, res, bidx):
        obj = int.__new__(cls, val)
        obj._value_ = val
        obj.res = res
        obj.bidx = bidx
        return obj

    def __str__(self):
        return self.name

    @property
    def band_id(self):
        return self._value_

    @classmethod
    def is_native_res(cls, b, res):
        b = b if isinstance(b, cls) else cls[b]
        return np.isclose(b.res, res)
        
    @classmethod
    def norm(self, s):
        if isinstance(s, str):
            m = re.match('^(B)?(\d+)(A)?$', s, re.IGNORECASE)
            if m:
                n = int(m[2])
                s = f'B{n}A'  if m[3] else f'B{n:0>2d}'
        elif isinstance(s, int):
            for k,v in l1c_band.__members__.items():
                if v._value_==s:
                    s = k
                    break
        return s

    @classmethod
    def parse_list(cls, s, sep=','):
        '''
        Parse a string into a list of bands
        '''
        return [cls[b] for b in s.split(sep)]

    @classmethod
    def match_bands(cls, regex, flags=re.IGNORECASE):
        '''
        
        '''
        p = re.compile(regex, flags)
        blist = [
            cls[b]
            for b in list(cls.__members__.keys())
            if p.match(b)
        ]
        return blist

# add bands to global namespace
globals().update(l1c_band.__members__)


def parse_data_obj(e, factory=dict):
    '''
    Parse data object elements in product metadata XML file
    '''
    d = factory()
    f = e.find('./byteStream')
    d['mimetype'] = f.attrib['mimeType']
    d['size']     = int(f.attrib['size'])
    d['relpath']  = f.find('fileLocation').attrib['href']
    d['checksum']   = f.find('checksum').text
    return d


def parse_zen_azi(e):
    '''
    Parse zenith/azimuth angles in tile metadata XML file
    '''
    zen  = np.array([ [float(v) for v in s.text.split()] for s in e.find('Zenith/Values_List') ])
    azi  = np.array([ [float(v) for v in s.text.split()] for s in e.find('Azimuth/Values_List')])
    return zen,azi


def slice_band(rb, sx, sy):
    '''
    Get a sub-region of a gdal.Band given by slices

    Parameters
    ----------
    rb : gdal.Band
       the GDAL raster band
    sx : slice or None
       the slice along the x-axis
    sy : slice or None
       the slice along the y-axis

    Returns
    -------
    
    '''

    # determine raster size
    ds = rb.GetDataset()
    nx,ny = ds.RasterXSize, ds.RasterYSize
    # get range of sub-window to read
    sx = sx or slice(None)
    sy = sy or slice(None)
    xl,xr,dx = sx.indices(nx)
    yu,yb,dy = sy.indices(ny)
    # sanity checks
    if dx<0 or dy<0: 
        raise ValueError("negative strides not supported")
    if xl>xr or yu>yb: 
        raise ValueError("start greater stop in slice not supported")            
    # read raster
    arr = rb.ReadAsArray(xoff=xl, yoff=yu, win_xsize=xr-xl,win_ysize=yb-yu)
    return arr[::dx,::dy]


def translate_band(rb, res, res_alg, sx=None, sy=None):
    '''
    Translate a gdal.Band to new resolution, optionally selecting a sub-region by slices
    '''        

    # get Dataset of raster and get raster shape / resolution
    ds = rb.GetDataset()
    xx,nat_res,_,yy,_,_ = ds.GetGeoTransform()
    nx,ny = ds.RasterXSize,ds.RasterYSize

    # get shape in new resolution
    nx = math.floor(nx*nat_res/res)
    ny = math.floor(ny*nat_res/res)

    # determine requested sub-window
    sx = sx or slice(None)
    sy = sy or slice(None)
    xl,xr,dx = sx.indices(nx)
    yu,yb,dy = sy.indices(ny)

    # sanity checks
    if dx<0 or dy<0: 
        raise ValueError("Negative strides not supported")
    if xl>xr or yu>yb: 
        raise ValueError("stop must be greater than start in slices")

    # prepare translate options
    topts = gdal.TranslateOptions(
        format="MEM",
        bandList=[rb.GetBand()],
        projWin=[xx+res*xl, yy-res*yu, xx+res*xr, yy-res*yb],
        xRes=float(res),
        yRes=float(res),
        resampleAlg=res_alg,
    )

    # do gdal.Translate
    mds = gdal.Translate('', ds, options=topts)

    # read translated array and return it
    arr = mds.GetRasterBand(1).ReadAsArray()
    return arr[::dx,::dy]


def read_band(rb, res, res_alg, sx, sy):
    # read raster band and return it
    if (res is None) or is_native_res(rb, res):
        # >>>> no resampling needed
        if sx is None and sy is None:
            # >>>> no slicing needed
            arr = rb.ReadAsArray()
        else:
            # >>>> slice raster band
            arr = slice_band(rb, sx, sy)
    else:
        # >>>> we must resample, use gdal.Translate
        arr = translate_band(rb, resolution, res_alg, sx, sy)        
    return arr


def is_native_res(rb, res):
    ds = rb.GetDataset()
    _,nat_res,_,_,_,_ = ds.GetGeoTransform()
    return np.isclose(nat_res, res)


class safe_reader(object):

    def read_file(self, name):
        '''
        Read a file contained in S2 SAFE archive

        Parameters
        ----------
        name : str
            The name of the file to read
        size : int,optional
            The maximum number of bytes to read
        '''

        if self.zipfile:
            p = os.path.normpath(os.path.join(self.safe_dir, name))
            s = self.zipfile.read(p)
        else:
            p = os.path.normpath(os.path.join(self.safe_dir, name))
            f = open(p, 'r')
            s = f.read()
            f.close()
        return s

    def parse_xml(self, fpath):
        '''
        Get XML file from SAFE archive/directory and get XML tree

        Parameters
        ----------
        fpath : str
            The path of the file to parse
        '''        

        s = self.read_file(fpath)
        return et.fromstring(s)

    @cached_property
    def inspire_mtd(self):
        return self.parse_xml(self.data_obj["INSPIRE_Metadata"].relpath)

    @cached_property
    def manifest(self):
        return self.parse_xml('manifest.safe')

    @cached_property
    def gds(self):
        return gdal.Open(self.fpath, gdal.GA_ReadOnly)

    @cached_property
    def sds(self):
        return [
            gdal.Open(_sds, gdal.GA_ReadOnly)
            for _sds,_ in self.gds.GetSubDatasets()
            
        ]

    @cached_property
    def data_obj(self):
        xpath_data_obj = './dataObjectSection/dataObject'
        return adict({
            e.attrib['ID']: parse_data_obj(e, factory=adict)
            for e in self.manifest.findall(xpath_data_obj) 
        })


class l1c_reader(safe_reader):

    def __init__(self, fpath):
        # init fpath/zipfile/safe_dir attributes
        if zipfile.is_zipfile(fpath):
            self.zipfile  = zipfile.ZipFile(fpath, "r")
            self.fpath    = os.path.abspath(fpath)
            self.safe_dir = os.path.splitext(os.path.basename(fpath))[0]+".SAFE"
        elif os.path.isdir(fpath):
            self.zipfile  = False
            self.safe_dir = os.path.abspath(fpath)
            self.fpath    = os.path.join(self.safe_dir, "MTD_MSIL1C.xml")
        elif os.path.isfile(fpath) and os.path.basename(fpath)=="MTD_MSIL1C.xml":
            self.zipfile  = False
            self.safe_dir = os.path.dirname(self.fpath)
            self.fpath    = os.path.abspath(fpath)
        else:
            raise ValueError(f"Unrecognized SAFE file: {fpath}")
        return
    
    @cached_property
    def l1c_mtd(self):
        return self.parse_xml(self.data_obj["S2_Level-1C_Product_Metadata"].relpath)

    @cached_property
    def tile_mtd(self):
        return self.parse_xml(self.data_obj["S2_Level-1C_Tile1_Metadata"].relpath)

    @cached_property
    def datastrip_mtd(self):
        return self.parse_xml(self.data_obj["S2_Level-1C_Datastrip1_Metadata"].relpath)

    @cached_property
    def detfoo_ds(self):
        return {
            m[1]:self.open_gdal(v.relpath)
            for v in rdr.data_obj.values()
            if (m:=re.match('^MSK_DETFOO_(B\d[\dA]).jp2$', os.path.basename(v.relpath)))
        }
        return 

    @cached_property
    def qamask_ds(self):
        return {
            m[1]:self.open_gdal(v.relpath)
            for v in rdr.data_obj.values()
            if (m:=re.match('^MSK_QUALIT_(B\d[\dA]).jp2$', os.path.basename(v.relpath)))
        }

    def get_band_ref(self, b):
        b = l1c_band[b] if not isinstance(b, l1c_band) else b
        i,j = b.bidx
        return self.sds[i].GetRasterBand(j)


    def read_band(self, b, *, resolution=None, resample_alg=None, sx=None, sy=None):
        '''
        Read Sentinel-2 band

        Parameters
        ----------
        bnd : str, int, or l1c_band
            The band, e.g. 'B01'...'B11' and 'B8A'
        resolution: None, int or str
            The resolution as None, int or string(e.g. '20m','native')
        resampling: str, default= 'average' if resolution>'native' else 'cubic'
            The resampling algorithm, valid choices are 'nearest','bilinear','cubic','cubic_spline',
            'average','lanczos','min','max','mode','median','Q1','Q3'
        sx : None or slice
            Slice along the x-axis, applied after resampling
        sy : None or slice
            Slice along the y-axis, applied after resampling
        '''
        # get gdal raster band ref.
        rb = self.get_band_ref(b)
        res_alg = resample_map[resample_alg] if isinstance(resample_alg, str) else resample_alg
        return read_band(rb, resolution, res_alg, sx, sy)

    
    def open_gdal(self, relpath):
        if self.zipfile:
            fpath = "/vsizip/" + os.path.normpath(self.fpath) \
                + "/"+os.path.normpath(
                    os.path.join(self.safe_dir, relpath)
                )
        else:
            fpath = os.path.normpath(os.path.join(rdr.safe_dir, relpath))
        ds = gdal.Open(fpath)
        return ds


    def read_detfoo(self, b, resolution=None, sx=None, sy=None):
        '''
        Read Sentinel-2 detector footprint mask

        Parameters
        ----------
        bnd : str, int, or l1c_band
            The band, e.g. 'B01'...'B11' and 'B8A'
        resolution: None or int
            The resolution as None or int
        sx : None or slice
            Slice along the x-axis, applied after resampling
        sy : None or slice
            Slice along the y-axis, applied after resampling
        '''
        b   = l1c_band[b] if not isinstance(b, l1c_band) else b
        rb  = self.detfoo_ds[b.name].GetRasterBand(1)
        arr = read_band(rb, resolution, gdal.GRA_NearestNeighbour, sx, sy)
        return arr


    def read_qamask(self, b, resolution=None, sx=None, sy=None):
        '''
        Read Sentinel-2 band quality mask

        Parameters
        ----------
        bnd : str, int, or l1c_band
            The band, e.g. 'B01'...'B11' and 'B8A'
        resolution: None or int
            The resolution as None or int
        sx : None or slice
            Slice along the x-axis, applied after resampling
        sy : None or slice
            Slice along the y-axis, applied after resampling
        '''
        b   = l1c_band[b] if not isinstance(b, l1c_band) else b
        rb  = self.qamask_ds[b.name].GetRasterBand(1)
        arr = read_band(rb, resolution, gdal.GRA_NearestNeighbour, sx, sy)
        return arr

    def get_tile_id(self):
        tm = self.tile_mtd
        return tm.find('n1:General_Info/TILE_ID',tm.nsmap).text

    def get_sensing_time(self):
        tm = self.tile_mtd
        tstr = tm.find('n1:General_Info/SENSING_TIME',tm.nsmap).text
        return np.datetime64(tstr[:-1])
    
    def get_epsg_code(self):
        '''
        Get EPSG code
        '''
        tm = self.tile_mtd
        s = tm.find('n1:Geometric_Info/Tile_Geocoding/HORIZONTAL_CS_CODE',tm.nsmap).text
        return int(s.split(':')[1])
    
    def get_proj_str(self):
        sr = osr.SpatialReference() 
        sr.ImportFromWkt(self.sds[0].GetProjection())
        return sr.ExportToProj4()

    def get_proj(self):
        '''
        Get a pyproj Projection object
        '''
        return pyproj.Proj(self.get_proj_str())

    def get_xy_ul(self):
        '''
        Get the xy coordinates of the upper left corner
        '''
        x0,_,_,y0,_,_ = self.sds[0].GetGeoTransform()
        return x0,y0

    def get_xy_axes(self, res=10.0, bounds=False):
        if np.isclose(res,10.0):
            i=0
        elif np.isclose(res,20.0):
            i=1
        elif np.isclose(res, 60.0):
            i=2
        x0,dx,_,y0,_,dy = self.sds[i].GetGeoTransform()
        nx = self.sds[i].RasterXSize
        ny = self.sds[i].RasterYSize
        if bounds:
            xax = x0+dx*np.arange(nx+1)
            yax = y0+dy*np.arange(ny+1)
        else:
            xax = x0+dx*(0.5+np.arange(nx))
            yax = y0+dy*(0.5+np.arange(ny))
        return xax,yax
    
    def get_lonlat(self, res=10.0, bounds=False):
        xax,yax = self.get_xy_axes(res,bounds)
        xg,yg = np.meshgrid(xax,yax)
        proj = self.get_proj()
        return proj(xg,yg,inverse=True)

    def get_sun_angles(self, axes=True):
        tm = self.tile_mtd
        e  = tm.find('n1:Geometric_Info/Tile_Angles/Sun_Angles_Grid', tm.nsmap)
        szen,sazi = parse_zen_azi(e)
        retval = (szen, sazi)
        if axes:
            x,y = self.get_xy_ul()
            dx,dy = ( float(e.find('Zenith/'+s).text) for s in ('COL_STEP','ROW_STEP'))
            xax = x+dx*np.arange(szen.shape[1],dtype=float)
            yax = y-dy*np.arange(szen.shape[0],dtype=float)
            return szen,sazi,(xax,yax)
        else:
            return szen,sazi

    def get_view_angles(self, b, axes=True):
        '''
        Get band-specific viewing angles
        '''
        b   = l1c_band[b] if not isinstance(b, l1c_band) else b
        tm = self.tile_mtd
        vzen = {}
        vazi = {}
        xp_fmt = "n1:Geometric_Info/Tile_Angles/Viewing_Incidence_Angles_Grids[@bandId='{0}']"
        for e in tm.iterfind(xp_fmt.format(b.band_id), tm.nsmap):
            det_id = int(e.attrib['detectorId'])
            zn, az = parse_zen_azi( e )
            vzen[det_id] = zn
            vazi[det_id] = az
            s = zn.shape
            dx,dy = ( float(e.find('Zenith/'+s).text) for s in ('COL_STEP','ROW_STEP'))
        if axes:
            x,y = self.get_xy_ul()
            xax = x+dx*np.arange(s[1],dtype=float)
            yax = y-dy*np.arange(s[0],dtype=float)
            return vzen,vazi,(xax,yax)
        else:
            return vzen,vazi

    def get_detector_starttime(self):
        dm = self.datastrip_mtd
        # create empty dictionary as return value
        det_stime = dict()
        # iterate over the Band_Time_Stamp elements
        xp = 'n1:Image_Data_Info/Sensor_Configuration/Time_Stamp/Band_Time_Stamp'
        for bnd_tstamp in dm.iterfind(xp, dm.nsmap):
            bnd_id = int(bnd_tstamp.attrib['bandId'])
            # add empty dict
            det_stime[bnd_id] = dict()
            # iterate over detectors
            for det in bnd_tstamp.iterfind('Detector'):
                # parse time and detector ID
                ts     = np.datetime64(det.find('GPS_TIME').text)
                det_id = int(det.attrib['detectorId'])
                det_stime[bnd_id][det_id] = ts
        return det_stime

    def get_gps_ephem(self):

        def parse_gps_elem(e):
            time     = np.datetime64(e.find('GPS_TIME').text)
            x,y,z    = np.fromstring(e.find('POSITION_VALUES').text, float, sep=' ')/1e3
            vx,vy,vz = np.fromstring(e.find('VELOCITY_VALUES').text, float, sep=' ')/1e3
            return time,x,y,z,vx,vy,vz

        gps_ephem_dtype = [
            ('time','M8[ns]'), ('x','f8'),('y','f8'),('z','f8'), 
            ('vx','f8'),('vy','f8'),('vz','f8')
        ]
    
        dsm = self.datastrip_mtd
        xp =  'n1:Satellite_Ancillary_Data_Info/Ephemeris/GPS_Points_List/GPS_Point'
        pts = np.array(
            [ 
                parse_gps_elem(p)
                for p in dsm.iterfind(xp, dsm.nsmap)
            ],
            dtype=ephem_dtype,
        ).view(np.recarray)
        return pts

    
