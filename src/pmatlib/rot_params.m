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
function [ alpha, theta, phi ] = rot_params(a_null, a_rot)
    % ROT_PARAMS Calculate rotation parameters for Givens rotation

    if a_rot == 0
        theta = pi/2;
    else
        theta = atan(abs(a_null)/abs(a_rot));
    end
        
    alpha = -angle(a_rot);
    phi = -angle(a_null);
end