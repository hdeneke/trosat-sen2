import os
import re

import numpy as np
import xarray as xr

from xarray.backends.common import (
    BackendEntrypoint,
    BackendArray,
    AbstractDataStore
)

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
        mask_and_scale=False,
        drop_variables=None,
        decode_times=True,
        decode_coords=True,
        res=60.0,
        bands=["B02","B03","B04","B08","B11"],
        lazy=False,
    ):

        # build dataset using data store
        with l1c_store(fpath, res=res, bands=bands) as store:
            # get attributes
            vars, atts = store.load()
            ds = xr.Dataset(vars, attrs=atts)
            ds.set_close(store.close)
        return ds
