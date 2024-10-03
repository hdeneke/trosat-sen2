# Python standard library imports
import enum
import os
import re
import zipfile

# SciPy imports
import numpy as np
import xarray as xr

# OSGeo imports
from osgeo import gdal, ogr, osr
import pyproj

# misc. function or class imports
from functools import cached_property
from addict import Dict as adict

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
        raise ValueError("start greater stop in slice not supported")            
    # prepare translate options
    topts = gdal.TranslateOptions(
        format="MEM",
        bandList=[rb.GetBand()],
        projWin=[xx+res*xl, yy-res*yu, xx+res*xr, yy-res*yb],
        xRes=float(res),
        yRes=float(res),
        resampleAlg=res_alg,
    )
    # apply translation
    mds = gdal.Translate('', ds, options=topts)
    # read translated array and return it
    arr = mds.GetRasterBand(1).ReadAsArray()
    return arr[::dx,::dy]


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
            p = os.path.normpath(os.path.join(self.safe_id+'.SAFE', name))
            s = self.zipfile.read(p)
        else:
            p = os.path.normpath(os.path.join(self.path, name))
            f = open(os.path.normpath(p), 'r')
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

