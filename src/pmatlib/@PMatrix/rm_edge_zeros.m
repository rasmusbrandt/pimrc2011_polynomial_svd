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
function r = rm_edge_zeros(obj)
    % RM_EDGE_ZEROS Removes edging zeros. Needed for proper division to
    % work.
    
    % Is this the empty polynomial?
    if norm(sum(abs(obj.coefs),3)) == 0
        r = PMatrix(zeros(1,1,1));
        return;
    end
    
    r = PMatrix(obj.coefs, obj.const_ind);
    ss = size(obj);
    
    % Find leading zeros
    for lead = 1:ss(3)
        if norm(r.coefs(:,:,lead))
           % This is not a zero coefficient, break here
           break;
        end
    end
    
    % Find trailing zeros
    for last = ss(3):-1:1
        if norm(r.coefs(:,:,last))
           % This is not a zero coefficient, break here
           break;
        end
    end
    
    % Return trimmed polynomial
    r.coefs = r.coefs(:,:,lead:last);
    r.const_ind = r.const_ind - (lead-1);
end