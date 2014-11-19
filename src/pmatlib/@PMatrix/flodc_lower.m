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
function [ row, lead, mag ] = flodc_lower(obj, col)
    % FLODC_LOWER Find Largest Off-Diagonal Coefficient (magnitude, under
    % main diagonal)
    
    % Select coefficients to check
    cc = abs(obj.coefs(col+1:end,col,:));
    s = size(cc);
    
    % Put them in matrix, so we can use MAX
    ccs = squeeze(cc);
    
    % Find max
    if s(1) == 0
        ME = MException('PMatrix:bad_column', ...
                        'Incorrectly selected column');
        throw(ME);
    elseif s(1) > 1
        % Find max, in rows
        [Y, rows2] = max(ccs);

        % Find max, in leads
        [mag, ind] = max(Y);

        % What (true) row was that?
        row = col + rows2(ind);

        % Convert to lead
        lead = obj.ind2lead(ind);
    else
        % Find lead with largest coefficient
        [mag, ind] = max(ccs);
        
        row = col + 1;
        lead = obj.ind2lead(ind);
    end
end