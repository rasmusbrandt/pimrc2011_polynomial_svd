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
function [ row, col, lead, mag ] = flodc_all(obj)
    % FLODC_ALL Finds largest off-diagonal coefficient
    
    % Set up variables
    s = size(obj);
    col_store_upper = zeros(3,s(2)-1);
    col_store_lower = zeros(3,min(s(1) - 1, s(2)));
    
    % Loop over all columns
    for col = 2:s(2)
        [ row, lead, mag ] = obj.flodc_upper(col);
        col_store_upper(:,col) = [ row, lead, mag ].';
    end
    for col = 1:min(s(1) - 1, s(2))    
        [ row, lead, mag ] = obj.flodc_lower(col);
        col_store_lower(:,col) = [ row, lead, mag ].';
    end
    
    % Find max, of all column maxes
    [ mag_upper, col_upper ] = max(col_store_upper(3, :));
    [ mag_lower, col_lower ] = max(col_store_lower(3, :));
    
    % Return it
    if isempty(mag_lower)
        mag = mag_upper;
        col = col_upper;
        row = col_store_upper(1, col);
        lead = col_store_upper(2, col);
    elseif isempty(mag_upper)
        mag = mag_lower;
        col = col_lower;
        row = col_store_lower(1, col);
        lead = col_store_lower(2, col);
    elseif mag_upper > mag_lower
        mag = mag_upper;
        col = col_upper;
        row = col_store_upper(1, col);
        lead = col_store_upper(2, col);
    elseif mag_upper < mag_lower
        mag = mag_lower;
        col = col_lower;
        row = col_store_lower(1, col);
        lead = col_store_lower(2, col);   
    end
end