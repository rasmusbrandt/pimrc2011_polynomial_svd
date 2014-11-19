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
function r = get(obj, rows, cols, leads)
    % GET Returns a select part of the matrix

    if leads ~= ':'
        r = PMatrix(obj.coefs(rows,cols,obj.lead2ind(leads)),obj.const_ind + leads(1));
    else
        r = PMatrix(obj.coefs(rows,cols,:),obj.const_ind);
    end
end