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
function [ row, col, lead, mag ] = flodc_upper_all(obj)
    % FLODC_UPPER_ALL Finds largest coefficient over the main diagonal
    
    % Set up variables
    s = size(obj);
    col_store = zeros(3,s(2));
    
    % Loop over all columns
    for col = 2:s(2)
        [ row, lead, mag ] = obj.flodc_upper(col);
        col_store(:,col) = [ row, lead, mag ].';
    end
    
    % Find max, of all column maxes
    [ mag, col ] = max(col_store(3, :));
    
    % Return it
    row = col_store(1, col);
    lead = col_store(2, col);
end