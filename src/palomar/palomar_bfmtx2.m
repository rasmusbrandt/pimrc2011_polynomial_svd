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
function [Uopt, Vopt] = palomar_bfmtx2(Uf, Vf, Lu, Lv, M)
    % Joint phase optimization of two equidimensional matrices
    
    % Variables
    [r1, c1, Nu] = size(Uf);
    [r2, c2, Nv] = size(Vf);
    
    if (r1 ~= r2) || (c1 ~= c2)
        error('Only equidimensional matrices allowed (yet)!');
    else
        rows = r1; cols = c1;
    end
    
    N = max(Nu,Nv);
    X = zeros(cols,cols,N);
    
    % Build projection matrices
    F = (1/sqrt(N))*fft(eye(N));
    
    FLu = F(:,1:Lu); 
    PF1 = FLu*FLu';
    PFT1 = eye(size(PF1)) - PF1;
    
    FLv = F(:,1:Lv); 
    PF2 = FLv*FLv';
    PFT2 = eye(size(PF2)) - PF2;
    
    % Loop over columns and optimize each separately
    for j = 1:cols
        
        % Build Q
        Q = zeros(N,N);
        for i = 1:rows
            Dm1 = diag(squeeze(Uf(i,j,:)));
            Q = Q + Dm1'*PFT1*Dm1;
            
            Dm2 = diag(squeeze(Vf(i,j,:)));
            Q = Q + Dm2'*PFT2*Dm2;
        end
        
        % Optimize
        if nargin == 5
            % Aimed initial guess
            x0 = exp(-1i*2*pi*M/N*(0:N-1)');
        else
            % Uniform initial guess
            x0 = exp(-1i*2*pi*rand(N,1));
        end
        x_am = palomar_am(Q, x0);
        
        % Store the results
        X(j,j,:) = reshape(x_am,1,1,N);
    end
    
    % Multiply all frequencies with obtained phases
    for n = 1:N
        Uf(:,:,n) = Uf(:,:,n)*X(:,:,n);
        Vf(:,:,n) = Vf(:,:,n)*X(:,:,n);
    end
    
    % Build optimized impulse response matrices
    coefs = ifft(Uf,[],3); Uopt = PMatrix(coefs(:,:,1:Lu));
    coefs = ifft(Vf,[],3); Vopt = PMatrix(coefs(:,:,1:Lv));
end
