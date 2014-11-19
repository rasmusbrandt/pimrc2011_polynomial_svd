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
function g = gcd(A, B)
    % GCD Returns Greatest Common Denominator of two polynomials
    %     Uses the Euclidean algorithm
    
    % Check ordering
    sA = size(A); sB = size(B);
    if sB(3) > sA(3)
        C = A; A = B; B = C; clear C;
    end
    
    % Perform first long division
    [~,r2] = A/B;
    
    % Loop it
    r1 = B; iter = 0;
    while (any(abs(r2.coefs) > 0) && iter < 10000)
        [~,r3] = r1/r2;
        r1 = r2; r2 = r3;
    end
    
    g = r1;
end