% pimrc2011_polynomial_svd
% Copyright (C) 2011, Rasmus Brandt

% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.

% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.

% You should have received a copy of the GNU General Public License along
% with this program; if not, write to the Free Software Foundation, Inc.,
% 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
if ispc
    blaslib = fullfile(matlabroot, ...
      'extern', 'lib', 'win32', 'lcc', 'libmwblas.lib');

    mex('-v', '-largeArrayDims', 'palomar_am_loop.c', blaslib);
elseif isunix
	% https://stackoverflow.com/questions/22367516/mex-compile-error-unknown-type-name-char16-t
    mex -v -Dchar16_t=uint16_T -largeArrayDims -lmwblas palomar_am_loop.c
end
