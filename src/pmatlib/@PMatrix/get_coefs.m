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
function r = get_coefs(obj, row, col, lead)
    % GET_COEFS Return coefficient(s) for selected indices
    
    if nargin == 1
        r = obj.coefs;
    elseif nargin == 2
        r = obj.coefs(row,:,:);
    elseif nargin == 3
        r = obj.coefs(row,col,:);
    elseif nargin == 4
        if obj.lead_exists(lead)
            r = obj.coefs(row,col,obj.lead2ind(lead));
        else
            r = 0;
        end
    else
        ME = MException('PMatrix:bad_arguments', ...
                        'Incorrect number of arguments to function call.');
        throw(ME);
    end
end