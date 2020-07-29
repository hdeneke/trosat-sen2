#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# This file is part of the TROPOS Satellite Suite (TROOSAT) developed 
# by the satellite group at TROPOS.
#
# TROSAT is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#
# Authors/Copyright(2012-2020):
# *Hartwig Deneke (deneke@tropos.de)
#
# With various contributions by:
# * Frank Werner (werner@tropos.de)
# * Fabian Senf (senf@tropos.de)
# * Sebastian Bley (bley@tropos.de)
# * Jonas Witthuhn (witthuhn@tropos.de)

import os
import zipfile
import numpy as np
import xml.etree.ElementTree as et
import pyproj
from osgeo import osr

class xml_meta(object):
    
    def __init__(self, s):
        self._root = et.fromstring( s )
        self._ns  = { 'n1': self._root.tag[1:].split('}')[0] }
        return

    def __getattr__(self, att):
        return getattr(self._root, att)

    @classmethod
    def from_zipfile(cls, fn, name):
        zip = zipfile.ZipFile(fn,'r')
        for f in zip.filelist:
            if os.path.basename(f.filename)==name:
                s = zip.read(f.filename)
                return cls(s)
        return None

class mtd_tl( xml_meta ):

    name = 'MTD_TL.xml'

    @classmethod
    def from_zipfile(cls, fn):
        return super().from_zipfile(cls, fn, cls.name)

    def get_tile_id( self ):
        '''
        Get TILE_ID element from general info sub-tree
        '''
        e = self.find( 'n1:General_Info/TILE_ID', self._ns )
        return e.text

    def get_sensing_time(self):
        '''
        Get SENSING_TIME element
        '''
        e = self.find( 'n1:General_Info/SENSING_TIME', self._ns )
        return np.datetime64( e.text[:-1] )

    def get_ul_xy(self):
        '''
        Get the xy coordinates for the upper left corner
        '''
        e = self.find( 'n1:Geometric_Info/Tile_Geocoding/Geoposition[1]', self._ns)
        x = float( e.find('ULX').text )
        y = float( e.find('ULY').text )
        return x,y

    def get_epsg_code(self):
        s = self.find( 'n1:Geometric_Info/Tile_Geocoding/HORIZONTAL_CS_CODE', self._ns).text
        return int(s.split(':')[1])

    def get_qi_masks(self):
        return [e.text for e in self.iterfind( 'n1:Quality_Indicators_Info/Pixel_Level_QI/'   \
                                               'MASK_FILENAME[@type="MSK_DETFOO"]', self._ns) ]

    @staticmethod
    def get_zen_azi( elem ):
        '''
        Return zenith/azimuth angle fields from  Values_List elements
        '''
        zen  = np.array([ [float(v) for v in e.text.split()] for e in elem.find('Zenith/Values_List') ])
        azi  = np.array([ [float(v) for v in e.text.split()] for e in elem.find('Azimuth/Values_List')])
        return zen, azi

    def get_view_angles( self ):
        '''
        Get the view angles from the MTD root element
        '''
        vzen = {}
        vazi = {}
        for e in self.iterfind( 'n1:Geometric_Info/Tile_Angles/Viewing_Incidence_Angles_Grids', self._ns ):
            b = int(e.attrib['bandId'])
            d = int(e.attrib['detectorId'])
            a1, a2 = get_zen_azi( e )
            if b not in vzen:
                vzen[b] = dict()
                vazi[b] = dict()
                vzen[b][d] = a1
                vazi[b][d] = a2
        return vzen, vazi

    def get_sun_angles(self, axes=False):
        '''
        Get the sun angles from the MTD root element
        '''
        e = self.find( 'n1:Geometric_Info/Tile_Angles/Sun_Angles_Grid', self._ns )
        szen,sazi = self.get_zen_azi(e)
        if axes:
            x,y = self.get_ul_xy()
            dx = float(e.find('Zenith/COL_STEP').text)
            dy = float(e.find('Zenith/ROW_STEP').text)
            xax = x+np.arange(szen.shape[1])*dx
            yax = y-np.arange(szen.shape[0])*dy
            return szen,sazi,(xax,yax)
        else:
            return szen,sazi

    
