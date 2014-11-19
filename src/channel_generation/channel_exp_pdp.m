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
function H = channel_exp_pdp(rows, cols, lags, psi)
    % CHANNEL_EXP_PDP Returns a channel with exponential power delay
    % profile
    
    % Get Gaussian channel
    H = PMatrix((1/sqrt(2))*(randn(rows,cols,lags) + ...
                          1i*randn(rows,cols,lags)));
    
    % Apply exponential PDP
    for n = 1:lags
        H.coefs(:,:,n) = exp(-psi*(n-1))*H.coefs(:,:,n);
    end
    
    % Normalize channel
    H = (1/H.fnorm)*H;
end
