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
function [A,B,C,D] = ss_ctrb(obj)
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
    B = [eye(p);zeros((n-1)*p,p)];
    
    % Form C
    C = [];
    for i = 2:m
        C = [C obj.coefs(:,:,i)];
    end
    
    % Form D
    D = obj.coefs(:,:,1);

end