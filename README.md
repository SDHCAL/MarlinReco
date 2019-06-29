# Note for this branch

The software on this branch seems to be the version to deal with the previous SDHCAL prototype simulation where
the step information needed for the SDHCAL digitizer was saved in two different collections :

- the SimCalorimeterHit collection with the steps center stored
- A GenericObject collection with the steps begin and end points stored.

The LCIO format have been evolved to include all the steps information in the SimCalorimeterHit objects.
This branch should be considered only if one wants to reprocess old simulation files from SDHCAL prototype simulation



# MarlinReco
[![Build Status](https://travis-ci.org/iLCSoft/MarlinReco.svg?branch=master)](https://travis-ci.org/iLCSoft/MarlinReco)
[![Coverity Scan Build Status](https://scan.coverity.com/projects/12363/badge.svg)](https://scan.coverity.com/projects/ilcsoft-marlinreco)

Assembly of various Marlin processor for reconstruction.

MarlinReco is distributed under the [GPLv3 License](http://www.gnu.org/licenses/gpl-3.0.en.html)

[![License](https://www.gnu.org/graphics/gplv3-127x51.png)](https://www.gnu.org/licenses/gpl-3.0.en.html)


## License and Copyright
Copyright (C), MarlinReco Authors

MarlinReco is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License long with this program.  If not, see <http://www.gnu.org/licenses/>.
