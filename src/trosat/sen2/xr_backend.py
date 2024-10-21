import os
import re

import numpy as np
import xarray as xr

from xarray.core.utils import Frozen
from xarray.backends.common import (
    BackendEntrypoint,
    BackendArray,
    AbstractDataStore
)

from osgeo import gdal

from .core import l1c_reader, l1c_band

class l1c_store(AbstractDataStore, l1c_reader):

    def __init__(
            self, 
            fpath,
            *, 
            res = 60.0,
            detfoo = False, 
            qmask = False,
            bands = ["B02","B03","B04"],
    ):
        l1c_reader.__init__(self, fpath)
        self.res = res
        #self.rs_up = gdal.GRA_Lanczos
        #self.rs_dn = gdal.GRA_Average
        self.detfoo = detfoo
        self.qmask = qmask
        #self.lazy = lazy
        self.bands = bands
        return

    @classmethod
    def open(cls, fpath, **kwargs):
        return cls(fpath, **kwargs)

    def get_attrs(self):
        return Frozen(self.get_global_attrs())

    def get_band_var(self, b):
        b = l1c_band[b] if not isinstance(b, l1c_band) else b
        arr  = self.read_band(
            b, res=self.res, rs_algo=gdal.GRA_Bilinear
        )
        atts = self.get_band_attrs(b)
        return xr.Variable(
            dims=("y","x"),
            data=arr,
            attrs=atts,
        )

    def get_qmask_var(self, b):
        arr  = self.read_qmask(b, res=self.res)
        atts = {"comment": "band quality mask: "+b.lower()}
        return xr.Variable(
            dims=("y","x"),
            data=arr,
            attrs=atts,
        )

    def get_detfoo_var(self, b):
        arr = self.read_detfoo(b, res=self.res)
        atts = {"comment": "band quality mask: "+b.lower()}
        return xr.Variable(
            dims=("y","x"),
            data=arr,
            attrs=atts,
        )

    def get_variables(self):
        vars = {}
        for b in self.bands:
            k = "refl_"+b.lower()
            vars[k] = self.get_band_var(b)
            if self.detfoo:
                k = "detfoo_"+b.lower()
                vars[k] = self.get_detfoo_var(b)
            if self.qmask:
                k = "qmask_"+b.lower()
                vars[k] = self.get_qmask_var(b)
        xax,yax = self.get_xy_axes(res=self.res, bounds=False)
        vars["x"] = xr.Variable(
            dims="x",
            data=xax/1e3,
            attrs={"standard_name":"projection_x_coordinate", "units":"km"},
        )
        vars["y"] = xr.Variable(
            dims="y",
            data=yax/1e3,
            attrs={"standard_name":"projection_y_coordinate", "units":"km"},
        )
        return xr.core.utils.Frozen(vars)

class sen2_backend(BackendEntrypoint):

    description = "Open Sentinel-2 L1C SAFE files in Xarray"
    url = "https://github.com/hdeneke/trosat-sen2"

    fpat =  '^(?P<sat>S2B)_(?P<prod>MSI(L1C|L2A))_'               \
            '(?P<time>\d{8}T\d{6})_(N\w+)_(R\w+)_(?P<tile>T\w+)_' \
            '(?P<ptime>\d{8}T\d{6}).zip'

    def guess_can_open(self, fpath):
        m = re.match(self, fpat, os.path.basename(fpath))
        retval = True if m else False
        return retval

    def open_dataset(
        self,
        fpath,
        *,
        #mask_and_scale=False,
        drop_variables=None,
        #decode_times=True,
        #decode_coords=True,
        res=60.0,
        bands=["B02","B03","B04","B08","B11"],
        detfoo=False,
        qmask=False,
    ):

        # build dataset using data store
        with l1c_store(
                fpath,
                res=res,
                bands=bands,
                detfoo=detfoo,
                qmask=qmask
        ) as store:
            # get attributes
            vars, atts = store.load()
            ds = xr.Dataset(vars, attrs=atts)
            ds.set_close(store.close)
        return ds
