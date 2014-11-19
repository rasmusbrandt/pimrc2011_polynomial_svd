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
function [ Q, R, cpxdata, cnvdata ] = pqrd(A, eps, mu, MaxIter, MaxSweeps)
    % PQRD Performs a polynomial QR decomposition
    
    % Should we measure complexity and convergence?
    nao = nargout;
    
    % Complexity data structure
    if nao > 2
        cpxdata = struct('NumSweeps',0, ...
                         'NumIters', zeros(MaxSweeps,1));
    end
    
    % Convergence data structure
    if nao > 3
        cnvdata = struct('BD_fnorms',[], ...
                         'OD_fnorms',[]); 
        BD_tmp = zeros(MaxIter,1); OD_tmp = zeros(MaxIter,1);
    end
    
    % Set up initial conditions
    s = size(A); 
    Q = PMatrix(eye(s(1)));
    g1 = 1 + eps;
    n = 0;
    
    % Algorithm sweeps
    while (n < MaxSweeps) && (g1 > eps)
        n = n + 1;
        
        % Loop over columns
        for col = 1:min(s(1) - 1, s(2))
            iter = 0;
            g2 = 1 + eps;
            
            % Loop over leads
            while (iter < MaxIter) && (g2 > eps)
                [ row, lead, g2 ] = A.flodc_lower(col); t = -lead;
                
                if g2 > eps
                    iter = iter + 1;
                    
                    % Calculate rotation parameters
                    a_null = A.get_coefs(row,col,lead);
                    a_rot = A.get_coefs(col,col,0);
                    [ alpha, theta, phi ] = rot_params(a_null, a_rot);
                    
                    % Coefficient step
                    A = A.coefstep(row, col, t, alpha, theta, phi);
                    Q = Q.coefstep(row, col, t, alpha, theta, phi);
                    
                    % Truncate
                    A = A.truncate(mu);
                    Q = Q.truncate(mu);
                    
                    % Convergence data
                    if nao > 3
                        BD_tmp(iter) = A.fnorm_bd;
                        OD_tmp(iter) = A.fnorm_od;
                    end
                end
            end
            
            % Complexity data
            if nao > 2
                cpxdata.NumIters(n) = cpxdata.NumIters(n) + iter;
            end
            
            % Convergence data
            if nao > 3
                cnvdata.BD_fnorms = [cnvdata.BD_fnorms;BD_tmp(1:iter)];
                cnvdata.OD_fnorms = [cnvdata.OD_fnorms;OD_tmp(1:iter)];
            end
        end

        % Check max coefficient under main diagonal
        [ ~, ~, ~, g1 ] = A.flodc_lower_all();
    end

    % Complexity data
    if nao > 2
        cpxdata.NumSweeps = n;
    end
    
    R = A;
end