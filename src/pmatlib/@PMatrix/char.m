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
function c = char(obj, approx)
    % CHAR Creates a string representation of the polynomial matrix

    s = size(obj.coefs);
    
    % Preallocate storage
    c = cell(s(1), s(2));

    for i = 1:s(1)
        for j = 1:s(2)
            % Convert coefficients to cells containing text
            tmp = PMatrix.polynomial_as_cells( ...
                squeeze(obj.coefs(i,j,:)), obj.const_ind, approx);
            
            % Store as array, in a cell
            c(i,j) = {[tmp{:}]};
        end
    end
end