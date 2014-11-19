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
function [ em, e2 ] = epgrm(dims, row, col, t, alpha, theta, phi)
    % EPGRM Elementary polynomial Givens rotation
    
    % Get epgr
    e2 = epgr(t, alpha, theta, phi);
    
    % Extend into larger matrix
    if ndims(e2.coefs) < 3
        ss = 1;
    else
        s2 = size(e2.coefs);
        ss = s2(3);
    end
    
    coefs = zeros(dims, dims, ss);
    const_ind = e2.const_ind;
    em = PMatrix(coefs, const_ind);
    
    em.coefs(:,:,em.lead2ind(0)) = eye(dims);
    em.coefs(col,col,:) = e2.coefs(1,1,:);
    em.coefs(col,row,:) = e2.coefs(1,2,:);
    em.coefs(row,col,:) = e2.coefs(2,1,:);
    em.coefs(row,row,:) = e2.coefs(2,2,:);
end