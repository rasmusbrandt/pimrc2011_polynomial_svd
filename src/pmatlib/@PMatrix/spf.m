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
function r = spf(obj, algorithm)
    % SPF Spectral factorization
    
    % Based on 5. SPECTRAL FACTORIZATION VIA A RICCATTI EQUATION in
    % Sayed, Kailath: A survey of spectral factorization methods
    
    % PARAMETER
    CONV_P = 1e-5;
    
    if nargin == 1
        algorithm = 'dare';
    end
    
    s = size(obj);
    if ~obj.is1x1()
        error('Only works for 1x1 polynomials.');
    end
    
    if obj.isconstant()
        r = PMatrix(sqrt(obj.coefs));
        return;
    end

    % Degree
    m = max(abs(obj.ind2lead(1)), abs(obj.ind2lead(s(3))));

    % Set up matrices
    F = diag(ones(m-1,1),-1);
    h = [zeros(1,m-1), 1];
    p0 = obj.coefs(1,1,obj.lead2ind(0));

    % Get coefficients (in a way so we DON'T get the small coefficients 
    % on the side).
    if abs(abs(obj.ind2lead(1))) > abs(obj.ind2lead(s(3)))
        Nbar = squeeze(obj.coefs(1,1,1:obj.lead2ind(1)));
        Nbar = conj(Nbar);
    else
        Nbar = squeeze(obj.coefs(1,1,obj.lead2ind(-1):s(3)));
        Nbar = flipud(Nbar);
    end
    
    if strcmp(algorithm, 'dare')
        % DARE ALGORITHM
        
        sigmabar = dare(F', h', zeros(size(F')), -p0, -Nbar);
        re = p0 - h*sigmabar*h';
        kp = (Nbar - F*sigmabar*h')*(1/re);

        % Check that the correct solution was obtained
        if any(eig(F - (F*sigmabar*h' - Nbar)/(h*sigmabar*h' - p0)*h) > 1)
            warning('Non-stable spectral factorization! DARE.');
        end
    elseif strcmp(algorithm, 'kalman')
        % KALMAN ALGORITHM
        
        L = Nbar; K = Nbar;
        Re = p0; Rr = p0;
        
        reldiff = 1;
        while reldiff > CONV_P
%             ORIGINAL VERSION
%             K_new  = K   - F*L*inv(Rr)*L'*h';
%             L_new  = F*L - K*inv(Re)*h*L;
%             Re_new = Re  - h*L*inv(Rr)*L'*h';
%             Rr_new = Rr  - L'*h'*inv(Re)*h*L;

            K_new  = K   - F*L/Rr*L'*h';
            L_new  = F*L - K/Re*h*L;
            Re_new = Re  - h*L/Rr*L'*h';
            Rr_new = Rr  - L'*h'/Re*h*L;

            reldiff = norm(K_new - K)/norm(K);
            
            K  = K_new;  L  = L_new;
            Re = Re_new; Rr = Rr_new;
        end
        
        re = Re;
        kp = K*(1/re);
    elseif strcmp(algorithm, 'array')
        % ARRAY KALMAN
        
        Lbar = Nbar/sqrt(p0); K = Nbar; Re = p0;
        
        A = [ sqrt(Re) h*Lbar ; K*1/sqrt(conj(Re)) F*Lbar ];
        
        reldiff = 1;
        
        while reldiff > CONV_P
            [~, Upper] = qr(A'); Lower = Upper';
            a = Lower(1,1); b = Lower(2:end,1); c = Lower(2:end,2);

            A = [ a h*c ; b F*c ];
            
            K_new = conj(a)*b;
            
            reldiff = norm(K_new - K)/norm(K);
            
            K = K_new;
            Re = a^2;
        end
        
        re = Re;
        kp = K*(1/re);
    end
    
    % Return factorization
    r = sqrt(re)*PMatrix(reshape([1;flipud(kp)],1,1,1+length(kp)),1);
    
    % Check that the filter is in fact stable
    if any(abs(roots(r.coefs)) > 1)
        warning('Non-stable spectral factorization! General');
    end
end