# MarlinKinfitTutorial
[![linux](https://github.com/ILDAnaSoft/MarlinKinfitTutorial/actions/workflows/linux.yml/badge.svg)](https://github.com/ILDAnaSoft/MarlinKinfitTutorial/actions/workflows/linux.yml)
[![Coverity Scan Build Status](https://scan.coverity.com/projects/12362/badge.svg)](https://scan.coverity.com/projects/ilcsoft-marlinkinfitprocessors)

Tutorial to MarlinKinfit

MarlinKinfitTutorial is distributed under the [GPLv3 License](http://www.gnu.org/licenses/gpl-3.0.en.html)

[![License](https://www.gnu.org/graphics/gplv3-127x51.png)](https://www.gnu.org/licenses/gpl-3.0.en.html)

## Content

* This tutorial performs a fit with four hard (px, py, pz, E) and and one soft (Z mass) constraint on ZH->mumubb events.  

* It comprises a workflow based on a DST file from a recent ILD production:
  - MC level selection of H->bb with "user" processor contained here
  - a fix to the covariance matrices of neutral PFOs (c.f. PhD Thesis Yasser Radkhorrami, in prep) with "user" processor contained here 
  - IsolatedLeptonTagging with standard iLCSoft/MarlinRecoAnalysis
  - Jet clustering with FastJet with standard iLCSoft/MarlinRecoAnalysis
  - ErrorFlow with standard iLCSoft/MarlinRecoAnalysis
  - the actual kinematic fit processor ZHllqq5CFit from MarlinKinfitProcessors

* more examples can be found in MarlinKinfitProcessors

## License and Copyright
Copyright (C), MarlinKinfitTutorial Authors

MarlinKinfitProcessors is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License long with this program.  If not, see <http://www.gnu.org/licenses/>.
