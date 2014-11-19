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
function r = coefstep(obj, row, col, t, alpha, theta, phi)
    % COLSTEP Applies a coefficient step
    
    sC = size(obj.coefs);
    if numel(sC) == 2
        sC(3) = 1;
    end
    
    % Set up 2x2 scalar Givens rotation
    c = cos(theta); s = sin(theta);
    G = [c*exp(1i*alpha),s*exp(1i*phi);-s*exp(-1i*phi),c*exp(-1i*alpha)];
    
    % Get affected elements
    B = obj.coefs([col row],:,:);
    
    % Extend with zeros, so we can circularly shift the lower row
    z = zeros(2,sC(2),abs(t));
    B = cat(3, z, B); B = cat(3, B, z); sB = size(B);
    
    % Shift lower row
    B(2,:,:) = circshift(B(2,:,:),[0 0 -t]);
    
    % Apply Givens rotation
    if exist('convmat', 'file')
        % Use mex implementation
        B = convmat(G, 1, B, obj.const_ind + abs(t));
    else
        % Loop over leads
        if length(sB) == 3
            for k = 1:sB(3)
                B(:,:,k) = G*B(:,:,k);
            end
        else
            B = G*B;
        end
    end

    % Unshift lower row
    B(2,:,:) = circshift(B(2,:,:),[0 0 t]);
    
    % Set up new coefficient matrix
    unaffected = 1:sC(1); unaffected([col row]) = [];
    new_coefs = zeros(sC(1),sC(2),sC(3) + 2*abs(t));
    new_coefs(unaffected,:,1+abs(t):end-abs(t)) = obj.coefs(unaffected,:,:);
    new_coefs([col row],:,:) = B;
    
    % Return
    r = PMatrix(new_coefs, obj.const_ind + abs(t));
end