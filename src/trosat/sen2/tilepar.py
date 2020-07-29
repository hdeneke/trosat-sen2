import pkg_resources as pkg_res
import numpy as np
import pyproj
import ogr

# Latitude bands for UTM projection, generated with the following code:
utm_lat_bands = 'CDEFGHJKLMNPQRSTUVWX'
# Generated as follows:
# > utm_lat_bands = [chr(i) for i in range(67,89)]
# > utm_lat_bands.remove('I'); utm_lat_bands.remove('O')
# > utm_lat_bands = ''.join(utm_lat_bands)

# Format for Proj string for UTM projection with WGS84 ellipsoid (GRS80 is Proj default)
utm_proj_str_fmt  = '+proj=utm +zone={0} +ellps=WGS84 +units=m +no_defs'

def numpy2poly(pts):
    ring = ogr.Geometry(ogr.wkbLinearRing)
    for x,y in pts[:,:2].astype(np.float):
        ring.AddPoint(float(x),float(y),0.0)
    poly = ogr.Geometry(ogr.wkbPolygon)
    poly.AddGeometry(ring)
    return poly


def get_bbox_dict(roi=None, return_type="ndarray"):
    '''
    Get dictionary containing bounding box of Sentinel2 tiles
    
    Parameters
    ----------
    roi : None, str, or ndarray
        the region of interest, either None,or given as WKT-formatted polygon or numpy array
    return_type: str
        accepted arguments are "ndarray", "wkt", "json", "gml"

    Returns
    -------
    bbox_dict : dict
        a dict containing the tile string as keys, and the bounding box as values. The Representation
        of the bounding box depends on return_type, either as wkt/JSON/GML-formatted string, or as numpy
        array.
    '''

    # get path to TILPAR file provided as package resource
    fn = pkg_res.resource_filename('trosat.sen2','share/S2A_OPER_GIP_TILPAR_MPC.kmz')
    # open file with LIBKML driver
    drv = ogr.GetDriverByName('LIBKML')
    ds = drv.Open(fn,0)
    # get feature layer
    lyr = ds.GetLayer()
    # define return function, based on return_fmnt
    if return_type=="wkt":
        feat2ret = lambda f: f.geometry().GetGeometryRef(0).ExportToWkt()
    elif return_type=="gml":
        feat2ret = lambda f: f.geometry().GetGeometryRef(0).ExportToGML()
    elif return_type=="json":
        feat2ret = lambda f: f.geometry().GetGeometryRef(0).ExportToJson()
    elif return_type=='ndarray':
        def feat2ret(f):
            g = f.geometry().GetGeometryRef(0)
            assert g.GetGeometryName()=='POLYGON'
            return  np.array(g.GetGeometryRef(0).GetPoints())[:,:2]
    else:
        raise ValueError('Unkown argument specified for return_type: ', )
    # extract boundary boxes from KML file feature layer
    if roi is None:
        # no region of interest given, return bounding box for all tiles
        bbox_dict = { f.GetField(0):feat2ret(f) for f in lyr }
    else:
        # handle different specifications of ROI
        if isinstance(roi, str):
            roi = ogr.CreateGeometryFromWkt(roi)
        elif isinstance(roi,ogr.Geometry):
            pass
        elif isinstance(roi,np.ndarray):
            roi = numpy2poly(roi)
        else:
            raise ValueError('Unkown argument specified for roi: ', roi)
        # return bounding box for tiles intersecting ROI
        bbox_dict = { f.GetField(0):feat2ret(f) for f in lyr if f.geometry().Intersect(roi)}
    # cleanup and return
    del lyr, ds
    return bbox_dict

def get_proj(tile):
    '''
    Get pyproj Proj object for given tile
    '''
    return pyproj.Proj(get_proj_str(tile))


def get_proj_str(tile):
    '''
    Get Proj projection string for tile

    Parameters
    ----------
    tile : str
        the tile string

    Returns
    -------
    proj_str : str
        a string that describes the projection of the tile
    '''
    # extract numeric/East-West part of UTM zone
    zone=tile[:2]
    # test if we are in the Southern hemisphere
    if utm_lat_bands.index(tile[2])<10:
        zone+= " +south"
    return utm_proj_str_fmt.format(zone)


def is_north_hemisphere(c):
    if ord(c.upper())>=78:
        return True
    return False

def get_epsg_code(tile):
    '''
    Get EPSG code for tile

    Parameters
    ----------
    tile : str
        the tile string

    Returns
    -------
    code : int
        The EPSG code of the tile projection
    '''

    zone = int(tile[:2])
    return 32600+zone if is_north_hemisphere(tile[2]) else 32700+zone
