# Python standard library imports
import os
import re
import zipfile

# SciPy imports
import numpy as np
import xarray as xr
import pyproj

from osgeo import gdal, ogr, osr

from functools import cached_property
from addict import Dict as adict

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


class safe_reader(object):

    # list of band names
    bands = [ 'B01', 'B02', 'B03', 'B04', 'B05', 'B06', 'B07', \
              'B08', 'B8A', 'B09', 'B10', 'B11', 'B12'         ]

    def __init__(self, name, mode='r'):

        # init path/name attributes
        self.path    = os.path.abspath(name)
        self.name    = os.path.basename(name)
        self.safe_id = os.path.splitext(self.name)[0]

        # check if 1. argument is a zipfile or directory
        if zipfile.is_zipfile(self.path):
            self.zipfile = zipfile.ZipFile(name,mode)
        elif os.path.isdir(self.path):
            self.zipfile = None
        else:
            raise ValueError('name is neither zipfile nor directory')

        # read manifest file
        self.manifest = self.read_xml('manifest.safe')

        # get data objects
        xpath_data_obj = './dataObjectSection/dataObject'
        self.data_obj = { e.attrib['ID']:adict(parse_data_obj(e)) for e in  self.manifest.findall(xpath_data_obj) }
        if 'S2_Level-1C_Product_Metadata' in self.data_obj:
            self.proc_level = 'Level-1C'
            self.init_l1c()
        elif 'S2_Level-2A_Product_Metadata' in self.data_obj:
            self.proc_level = 'Level-2A'
            self.init_l2a()
        else:
            self.proc_level = 'Unknown'
        return

    def init_l1c(self):
        
        # parse tile/product xml metadata
        self.tile_meta = self.read_xml(self.data_obj['S2_Level-1C_Tile1_Metadata'].relpath)
        self.prod_meta = self.read_xml(self.data_obj['S2_Level-1C_Product_Metadata'].relpath)

        # open gdal dataset
        p = self.path if self.zipfile else os.path.join(self.path,'MTD_MSIL1C.xml')
        self.gdal_ds = gdal.Open(p, gdal.GA_ReadOnly)
        # open gdal sub-datasets
        self.gdal_sds = {}
        for sds,_ in self.gdal_ds.GetSubDatasets():
            res = sds.split(':')[-2]
            self.gdal_sds[res] = gdal.Open(sds, gdal.GA_ReadOnly)
        # set band meta data
        self.band_meta = adict()
        for res in ('10m','20m','60m'):
            for i in range(self.gdal_sds[res].RasterCount):
                raster = self.gdal_sds[res].GetRasterBand(i+1)
                bmeta =  {k.lower():v for k,v in raster.GetMetadata_Dict().items()}
                bmeta['resolution'] = res
                bmeta['index'] = i+1
                bmeta['bandname'] = re.sub('B(\d)$',r'B0\1',bmeta['bandname'])
                self.band_meta[bmeta['bandname']] = adict(bmeta)
        return

    def init_l2a(self):
        return

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

    def read_xml(self, name):
        '''
        Get XML file from SAFE archive and return parsed XML tree
        '''
        s = self.read_file(name)
        return et.fromstring(s)

    def read_band(self, bnd, resolution=None, resampling=None, sx=None, sy=None):
        '''
        Read Sentinel-2 band

        Parameters
        ----------
        bnd : str
            The band name, 'B01'...'B11' and 'B8A'
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

        # get band metadata
        bm = self.band_meta[bnd]
        res = bm['resolution']
        i = bm['index']

        kwargs = {}

        # if we do not need to resample...
        if resolution in {None,'native'} or resolution==res:
            # ... we can return a virtual array view
            # get kwargs dict from slices
            if sx: kwargs.update({'xoff':sx.start,'xsize':sx.stop-sx.start, 'bufxsize':sx.stop-sx.start})
            if sy: kwargs.update({'yoff':sy.start,'ysize':sy.stop-sy.start, 'bufysize':sy.stop-sy.start})
            return self.gdal_sds[res].GetRasterBand(i).GetVirtualMemArray(**kwargs)
        # otherwise do resampling
        # get resolution
        r = int(resolution[:-1])
        # set default resampling
        if resampling is None:
            rorig = int(res[:-1])
            resampling = 'average' if r>rorig else 'nearest'
        # resample
        d1 = gdal.Translate('', self.gdal_sds[res], format='VRT', bandList=[i])
        d2 = gdal.Warp('', d1, format='VRT', xRes=r, yRes=r, resampleAlg=resample_map[resampling])
        if sx: kwargs.update({'xoff':sx.start,'win_xsize':sx.stop-sx.start})
        if sy: kwargs.update({'yoff':sy.start,'win_ysize':sy.stop-sy.start})
        a = d2.GetRasterBand(1).ReadAsArray(**kwargs)
        del d1,d2
        return a

    def get_tile_id(self):
        return self.tile_meta.find('n1:General_Info/TILE_ID',self.tile_meta.nsmap).text


    def get_sensing_time(self):
        tstr = self.tile_meta.find('n1:General_Info/SENSING_TIME', self.tile_meta.nsmap).text
        return np.datetime64(tstr[:-1])
    

    def get_epsg_code(self):
        '''
        Get EPSG code
        '''
        s = self.tile_meta.find('n1:Geometric_Info/Tile_Geocoding/HORIZONTAL_CS_CODE', self.tile_meta.nsmap).text
        return int(s.split(':')[1])


    def get_ul_xy(self):
        '''
        Get the xy coordinates for the upper left corner
        '''
        x,_,_,y,_,_ = self.gdal_sds["10m"].GetGeoTransform()
        return x,y


    def get_sun_angles(self, axes=False):
        e = self.tile_meta.find('n1:Geometric_Info/Tile_Angles/Sun_Angles_Grid', self.tile_meta.nsmap)
        szen,sazi = parse_zen_azi(e)
        if not axes:
            return szen,sazi
        x,y = self.get_ul_xy()
        dx,dy = ( float(e.find('Zenith/'+s).text) for s in ('COL_STEP','ROW_STEP'))
        xax = x+dx*np.arange(szen.shape[1],dtype=np.float)
        yax = y-dy*np.arange(szen.shape[0],dtype=np.float)
        return szen,sazi,(xax,yax)

    def get_axes_xy(self, resolution="10m", bounds=False):
        x0,dx,_,y0,_,dy = self.gdal_sds[resolution].GetGeoTransform()
        nx = self.gdal_sds[resolution].RasterXSize
        ny = self.gdal_sds[resolution].RasterYSize
        if bounds:
            xax = x0+dx*np.arange(nx+1)
            yax = y0+dy*np.arange(ny+1)
        else:
            xax = x0+dx*(0.5+np.arange(nx))
            yax = y0+dy*(0.5+np.arange(ny))
        return xax,yax

    def get_proj_str(self):
        sr = osr.SpatialReference() 
        sr.ImportFromWkt(self.gdal_sds['10m'].GetProjection())
        return sr.ExportToProj4()

    def get_proj(self):
        sr = osr.SpatialReference() 
        sr.ImportFromWkt(self.gdal_sds['10m'].GetProjection())
        return pyproj.Proj(sr.ExportToProj4())
        
    def get_lonlat(self, resolution="10m", bounds=False):
        xax,yax = self.get_axes_xy(resolution,bounds)
        xg,yg = np.meshgrid(xax,yax)
        proj = self.get_proj()
        return proj(xg,yg,inverse=True)
