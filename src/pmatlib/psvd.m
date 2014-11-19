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
function [ U, D, V, cpxdata, cnvdata ] = ...
    psvd(A, eps, mu, MaxIter, MaxSweeps, MaxPSVDIter)
    % PSVD Performs a polynomial QR decomposition
    
    % Should we measure complexity and convergence?
    nao = nargout;
    
    % Complexity data structure
    if nao > 2
        cpxdata = struct('NumPSVDIters' , 0,  ...
                         'NumPQRDSweeps', zeros(2,MaxPSVDIter), ...
                         'NumPQRDIters' , zeros(MaxSweeps,MaxPSVDIter,2));
    end
    
    % Convergence data structure
    if nao > 3
        cnvdata = struct('OD_fnorms',[]); 
    end
               
    % Set up initial conditions
    s = size(A);

    U = PMatrix(eye(s(1)));
    V = PMatrix(eye(s(2)));
    g = 1 + eps;
    iter = 0;
    
    % Iterate
    while (iter < MaxPSVDIter) && (g > eps)
        [ ~, ~, ~, g ] = A.flodc_all();

        if g > eps
            iter = iter + 1;

            if nao > 3
                [ U1, R1, cpx1, cnv1 ] = ...
                    pqrd(A, eps, mu, MaxIter, MaxSweeps);
            elseif nao > 2
                [ U1, R1, cpx1 ] = ...
                    pqrd(A, eps, mu, MaxIter, MaxSweeps);
            else
                [ U1, R1 ] = ...
                    pqrd(A, eps, mu, MaxIter, MaxSweeps);
            end
            
            Aprim = R1.ph();
            U = U1*U;
            
            if nao > 3
                [ V1, R2, cpx2, cnv2 ] = ...
                    pqrd(Aprim, eps, mu, MaxIter, MaxSweeps);
            elseif nao > 2
                [ V1, R2, cpx2 ] = ...
                    pqrd(Aprim, eps, mu, MaxIter, MaxSweeps);
            else
                [ V1, R2 ] = ...
                    pqrd(Aprim, eps, mu, MaxIter, MaxSweeps);
            end
            
            V = V1*V;
            A = R2.ph();
            
            A.truncate(mu);
            U.truncate(mu);
            V.truncate(mu);
            
            % Complexity data
            if nao > 2
                cpxdata.NumPQRDSweeps(1,iter) = cpx1.NumSweeps;
                cpxdata.NumPQRDSweeps(2,iter) = cpx2.NumSweeps;
                cpxdata.NumPQRDIters(:,iter,1) = ...
                    cpx1.NumIters;
                cpxdata.NumPQRDIters(:,iter,2) = ...
                    cpx2.NumIters;
            end
            
            % Convergence data
            if nao > 3
                cnvdata.OD_fnorms = ...
                    [cnvdata.OD_fnorms;cnv1.OD_fnorms;cnv2.OD_fnorms];
            end
        end
    end
    
    % Complexity data
    if nao > 2
        cpxdata.NumPSVDIters = iter;
    end
    
    D = A;
end