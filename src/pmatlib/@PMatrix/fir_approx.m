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
function r = fir_approx(obj, common_den)
    % FIR_APPROX Return FIR approximation of matrix
    
    % PARAMETERS
    K = 1000;
    mu = 1e-6;
    
    % Sizes
    s = size(obj); p = s(1); q = s(2); r = s(3);
    
    % Get trimmed common numerator
    ap = common_den.rm_edge_zeros();
    a = squeeze(ap.coefs);
    na1lead = common_den.ind2lead(1);
    
    % Loop over all matrix elements
    r = PMatrix(obj);
    for i = 1:p
        for j = 1:q
            % Get denominator
            bp = r.get(i,j,:);
            
            % Break out so that first coefficient is constant
            nb1lead = bp.ind2lead(1);
            
            % Set in correct format
            b = squeeze(bp.coefs);
            
            % Get impulse response
            tmp = filter(b,a,[1 zeros(1,K-1)]);
            
            % Set up PMatrix for impulse response
            ptmp = PMatrix(reshape(tmp,1,1,K),na1lead-nb1lead + 1);
            
            % Truncate away the non-important taps
            %ptmp2 = ptmp.truncate(mu);
            ptmp2 = ptmp;
            
            % Set it
            r = r.set_element(i,j,ptmp2);
        end
    end
end