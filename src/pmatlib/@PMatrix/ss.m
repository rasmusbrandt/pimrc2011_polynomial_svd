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
function r = ss(obj)
    % FIXME: optimize the code
    
    [p,q,m] = size(obj.coefs);
    n = m - 1;
    
    if obj.const_ind ~= 1
        error('Can only deal with causal matrices, starting at lag 0!');
    end
    
    % Form A
    A = [];
    for i = 1:(n-1)
        A = blkdiag(eye(p),A);
    end
    A = [zeros(p,(n-1)*p);A];
    A = [A zeros(n*p,p)];
    
    % Form B
    B = [];
    for i = 1:n
        B = [obj.coefs(:,:,i+1);B]; % +1 to not take the direct feedthrough matrix
    end
    
    % Form C
    C = [zeros(p,(n-1)*p) eye(p)];
    
    % Form D
    D = obj.coefs(:,:,1);
    
    r = ss(A,B,C,D,[]);
end