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
function [ Q, R, c, cpxdata, cnvdata ] = pqrd_lcm(A)
    % PQRD_LCM Performs a polynomial QR decomposition through the LCM
    % method
    
    % Should we measure complexity and convergence?
    nao = nargout;
    
    % Set up initial conditions
    [p,q,r] = size(A.coefs); 
    Nsub = p*min(p-1,q) - sum(1:min(p-1,q));
    Q = PMatrix(eye(p));
    cacc = PMatrix(1);
    
    % Complexity data structure
    if nao > 2
        cpxdata = struct('NumIters', 0);
    end
    
    % Convergence data structure
    if nao > 3
        cnvdata = struct('BD_fnorms',zeros(Nsub,1), ...
                         'OD_fnorms',zeros(Nsub,1)); 
    end

    % Loop over columns
    iter = 0;
    for col = 1:min(p-1,q)
        % Loop over rows
        for row = col+1:p
            iter = iter + 1;
            
            % Get diagonal element
            d = A.get(col,col,:);
            
            % Get element to null
            f = A.get(row,col,:);
            
            % Find gcd
            %g = gcd(d,f);
            
            %Is gcd a monic polynomial?
            %if sum(g.coefs ~= 0) == 1
            %   g = PMatrix(1);
            %else
            %    disp('GCD found!');
            %    g
            %end
            
            g = PMatrix(1);
            
            % Calculate temporary polynomials
            alpha = f/g;
            beta  = d/g;
            cc = alpha*alpha.ph + beta*beta.ph;
            c = cc.spf('kalman');
            
            % Set up Given's rotation
            G = PMatrix(eye(p,p));
            for n = 1:p
                G = G.set_element(n,n,c);
            end
            G = G.set_element(col,col,-beta.ph);
            G = G.set_element(col,row,-alpha.ph);
            G = G.set_element(row,col,alpha);
            G = G.set_element(row,row,-beta);
            
            % Normalize
            nmp = G.fnorm;
            G.coefs = G.coefs/nmp;
            c.coefs = c.coefs/nmp;
            
            % Apply it
            A = G*A;
            Q = G*Q;
            cacc = c*cacc;
            
            % Convergence data
            if nao > 3
                cnvdata.BD_fnorms(iter) = A.fnorm_bd/A.fnorm;
                cnvdata.OD_fnorms(iter) = A.fnorm_od/A.fnorm;
            end
        end
    end
    
    % Complexity data
    if nao > 2
        cpxdata.NumIters = iter;
    end
    
    % Set triangular matrix
    R = A;
    
    % Get accumulated spectral factor. DO NOT TRUNCATE THE FACTOR!!!
    % SPF algorithm will not converge then.
%     QQph = Q*Q.ph;
%     QQph11 = QQph.get(1,1,:);
%     c = QQph11.spf('kalman');
    c = cacc;
end