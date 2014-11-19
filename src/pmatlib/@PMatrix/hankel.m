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
function hank = hankel(obj)
    % HANK Returns Hankel matrix of filter
    
    sH = size(obj.coefs);
    
    hank = zeros(sH(1)*(sH(3)-1), sH(2)*(sH(3)-1));
    for j = 1:(sH(3)-1)
        for i = 1:(sH(3) - j)
            hank( ((i-1)*sH(1)+1):i*sH(1), ((j-1)*sH(2)+1):j*sH(2) ) = ...
                obj.coefs(:,:, i + j);
        end
    end 
end