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
function [ U, D, V, c, d ] = ...
    psvd_lcm(A, MaxPSVDIter)
    % PSVD_LCM PSVD through LCM method
    
    % Set up initial conditions
    [p,q,r] = size(A.coefs);
    
    U = PMatrix(eye(p));
    V = PMatrix(eye(q));
    
    c = PMatrix(ones(1,1,1));
    d = PMatrix(ones(1,1,1));

    iter = 0;
   
    % Iterate
    while (iter < MaxPSVDIter)
        iter = iter + 1;
        
        % PQRD 1
        [ U1, R1, c1 ] = pqrd_lcm(A);

        % Truncate 1
        %U1 = U1.truncate(mu);
        %R1 = R1.truncate(mu);
        %c1 = c1.truncate(mu);

        % Orthogonal 1
        U = U1*U;

        % Numerator 1
        c = c1*c;

        % Flip 1
        Aprim = R1.ph();

        % PQRD 2
        [ V1, R2, d1 ] = pqrd_lcm(Aprim);

        % Truncate 2
        %V1 = V1.truncate(mu);
        %R2 = R2.truncate(mu);
        %d1 = d1.truncate(mu);

        % Orthogonal 2
        V = V1*V;

        % Numerator 2
        d = d1*d;

        % Flip 2
        A = R2.ph();
    end
    
    % Singular values
    D = A;
end