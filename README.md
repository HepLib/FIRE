## About
`FIRE` was developed by [A. Smirnov](https://gitlab.com/feynmanintegrals/fire). All credit for the original code goes to him and we would like to take this opportunity to thank [A. Smirnov](https://gitlab.com/feynmanintegrals/fire) for his great work!

More information, please visit [Wiki Section](https://github.com/HepLib/FIRE/wiki).

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.


## Installation
```
git clone https://github.com/HepLib/FIRE.git
cd FIRE
make -j 8 dep
make -j 8
make test
```
## 3-versions
- `M/FIRE` for polynormial version with ration polynomial coefficients in IBP equations.
- `Q/FIRE` for ration version with rational integer coefficients in IBP equations.
- `F/FIRE` for float version with fix-precision float coefficients in IBP equations.
