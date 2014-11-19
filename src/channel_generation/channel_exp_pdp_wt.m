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
function H = channel_exp_pdp_wt(rows, cols, lags, WT)
    % CHANNEL_EXP_PDP Returns a channel with exponential power delay
    % profile, given a certain bandwidth-coherence time product.

    b = exp(-1/WT); a = 1 - b;
    
    coefs = zeros(rows,cols,lags);
    
    for i = 1:lags
        coefs(:,:,i) = (1/sqrt(2))*sqrt(a*b^i)*(randn(rows,cols) + ...
                          1i*randn(rows,cols));
    end

    % PMatrix
    H = PMatrix(coefs);
    
    % Normalize channel
    H = (1/H.fnorm)*H;
end
