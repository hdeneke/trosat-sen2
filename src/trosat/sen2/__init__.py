#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# This file is part of the TROPOS Satellite Suite (TROSAT) developed 
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
# Authors/Copyright(2011-2020):
# *Hartwig Deneke (deneke@tropos.de)
#
# With contributions by:
# * Sebastian Bley (bley@tropos.de)
# * Fabian Senf (senf@tropos.de)
# * Frank Werner (werner@tropos.de)
# * Jonas Witthuhn (witthuhn@tropos.de)

import pkg_resources as pkg_res
import ogr

from .core import *
from . import meta
from . import tilepar


from . import _version
__version__ = _version.get_versions()['version']
