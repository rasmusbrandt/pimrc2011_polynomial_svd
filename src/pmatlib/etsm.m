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
function r = etsm(dims, row, t)
    % ETSM Returns an elementary time shift matrix
    
    if t > 0
        coefs = zeros(dims, dims, abs(t)+1);
        coefs(:,:,end) = eye(dims);
        coefs(row,row,end) = 0;
        coefs(row,row,1) = 1;

        r = PMatrix(coefs, t+1);
    elseif t < 0
        coefs = zeros(dims, dims, abs(t)+1);
        coefs(:,:,1) = eye(dims);
        coefs(row,row,1) = 0;
        coefs(row,row,end) = 1;

        r = PMatrix(coefs, 1);
    else
        r = PMatrix(eye(dims));
    end
end