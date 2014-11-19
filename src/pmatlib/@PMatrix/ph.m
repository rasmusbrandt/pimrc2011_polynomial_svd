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
function r = ph(obj)
    % PH Parahermitian conjugate

    % Hermitian conjugate
    r = ctranspose(obj);
    
    % Flip z -> 1/z
    s = size(r.coefs);
    
    if length(s) > 2
        % Not scalar matrix
        
        r.coefs = flipdim(r.coefs, 3);
        r.const_ind = s(3) - r.const_ind + 1;
    end
end

