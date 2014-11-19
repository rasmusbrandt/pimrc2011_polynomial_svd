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
function [ Q, R, cpxdata, cnvdata ] = pqrd_mod(A, eps_r, mu, rho, MaxSweeps)
    % PQRD_MOD Modified Polynomial QRD By Columns
    
    % Should we measure complexity and convergence?
    nao = nargout;

    % Set up initial conditions
    [p,q,r] = size(A.coefs);
    Nsub = p*min(p-1,q) - sum(1:min(p-1,q));
    MaxIter = round(rho*abs(log10(eps_r))*Nsub*r);
    Q = PMatrix(eye(p));
    h1 = 1 + eps_r;
    n = 0;
    Aorig = A;
    
    % Complexity data structure
    if nao > 2
        cpxdata = struct('NumSweeps',0, ...
                         'NumIters', zeros(MaxSweeps,1), ...
                         'SumKflogKf', 0, ...
                         'SumKf', 0);
    end
    
    % Convergence data structure
    if nao > 3
        cnvdata = struct('BD_fnorms',[], ...
                         'OD_fnorms',[]); 
        BD_tmp = zeros(MaxIter,1); OD_tmp = zeros(MaxIter,1);
    end
    
    % Algorithm sweeps
    while (n < MaxSweeps) && (h1 > eps_r)
        n = n + 1;
        
        % Loop over columns
        for col = 1:min(p-1,q)
            iter = 0;
            h2 = 1 + eps_r;
            
            % Loop over leads
            while (iter < MaxIter) && (h2 > eps_r)
                % Find coefficient with largest magnitude
                [ row, lead, g2 ] = A.flodc_lower(col); t = -lead;
                
                % Estimate proportion of energy below main diagonal,
                % if all coefficients had this magnitude
                h2 = r*Nsub*(g2/Aorig.fnorm)^2;
                
                % Is it too large?
                if h2 > eps_r
                    iter = iter + 1;
                    
                    % Calculate rotation parameters
                    a_null = A.get_coefs(row,col,lead);
                    a_rot = A.get_coefs(col,col,0);
                    [ alpha, theta, phi ] = rot_params(a_null, a_rot);
                    
                    % Coefficient step
                    A = A.coefstep(row, col, t, alpha, theta, phi);
                    Q = Q.coefstep(row, col, t, alpha, theta, phi);
                    
                    % Complexity data
                    if nao > 2
                        Kfi = size(A.coefs,3);
                        
                        cpxdata.SumKf = cpxdata.SumKf + Kfi;
                        cpxdata.SumKflogKf = cpxdata.SumKflogKf + ...
                                                       Kfi*log2(Kfi);
                    end
                    
                    % Truncate
                    A = A.truncate(mu);
                    Q = Q.truncate(mu);
                    
                    % Get new number of lags
                    [~,~,r] = size(A.coefs);
                    
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

        % Check proportion of energy under main diagonal
        h1 = (A.fnorm_bd/Aorig.fnorm)^2;
    end

    % Complexity data
    if nao > 2
        cpxdata.NumSweeps = n;
    end
    
    R = A;
end