#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# This file is part of the TROPOS Satellite utilities (TROSAT) developed 
# within the remote sensing department at the Leibniz Institute for
# Tropospheric research (TROPOS).
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
# Authors/Copyright(2012-2024):
# * Hartwig Deneke (deneke@tropos.de)

import glob
from setuptools import setup
import versioneer

setup(
    name         = 'trosat-sen2',
    version      = versioneer.get_version(),
    cmdclass     = versioneer.get_cmdclass(),
    author       = 'Hartwig Deneke', 
    author_email = 'deneke@tropos.de',
    license      = 'GPLv3',
    description  = 'Package for working with Sentinel-2 data',
    #scripts     = glob.glob('bin/*.py'),
    packages     = [ 'trosat', 'trosat.sen2' ],
    package_dir  = { '': 'src'},
    package_data = { 'trosat.sen2': ['share/*']},
    install_requires = [
        'gdal',
    ],
    keywords = ['Sentinel-2'],
    classifiers = [
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 3',
    ]
)
