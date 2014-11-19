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
% Order comparison

% Testing or real?
[~,hn] = system('hostname');
if strcmp(hn(1:3), 'sim')
    TESTING = false;

    L = 30;
    N = 100;
    rows = 3;
    cols = 3;
    WTau1 = 5;
    WTau2 = 10;

    eps_r = 1e-3;
    mu = 1e-6;
    mu1 = 1e-6;
    mu2 = 1e-3;
    rho = 100;
    MaxSweeps = 1000;
    MaxPSVDIter = 1000;

    nw = matlabpool('size'); if nw == 0;nw = 1;end
    fprintf('Starting production run with %d workers.\n', nw);
else
    TESTING = true;
    
    L = 45;
    N = 2;
    rows = 3;
    cols = 3;
    WTau1 = 5;
    WTau2 = 10;

    eps_r = 1e-3;
    mu = 1e-6;
    mu1 = 1e-6;
    mu2 = 1e-3;
    rho = 100;
    MaxSweeps = 1000;
    MaxPSVDIter = 1000;

    
    nw = matlabpool('size'); if nw == 0;nw = 1;end
    fprintf('Starting testing run with %d workers.\n', nw);
end



results1 = zeros(N,L);
results2 = zeros(N,L);
results_trunc16 = zeros(N,L);
results_trunc13 = zeros(N,L);
results_trunc26 = zeros(N,L);
results_trunc23 = zeros(N,L);

for l = 2:L
    fprintf('Starting L %d out of %d\n', l, L);
    parfor n = 1:N
        H = channel_exp_pdp_wt(rows,cols,l,WTau1);
        
         % Get PSVD
        [U_ph,Dp,V_ph] = psvd_mod(H, eps_r, mu, rho, MaxSweeps, MaxPSVDIter);
        
        Ut6 = U_ph.truncate(mu1);
        Ut3 = U_ph.truncate(mu2);
        
        results1(n,l) = size(U_ph.coefs,3);
        results_trunc16(n,l) = size(Ut6.coefs,3);
        results_trunc13(n,l) = size(Ut3.coefs,3);
    end
    
%     parfor n = 1:N
%         H = channel_exp_pdp_wt(rows,cols,l,WTau2);
%         
%          % Get PSVD
%         [U_ph,Dp,V_ph] = psvd_mod(H, eps_r, mu, rho, MaxSweeps, MaxPSVDIter);
%         
%         Ut6 = U_ph.truncate(mu1);
%         Ut3 = U_ph.truncate(mu2);
%         
%         results2(n,l) = size(U_ph.coefs,3);
%         results_trunc26(n,l) = size(Ut6.coefs,3);
%         results_trunc23(n,l) = size(Ut3.coefs,3);
%     end
end

if ~TESTING
    mkdir('/nobackup/rabr5411/poly/order_results');
    cl = clock; 
    HH = num2str(cl(4)); if length(HH) == 1;HH = ['0' HH]; end
    MM = num2str(cl(5)); if length(MM) == 1;MM = ['0' MM]; end
    save(['/nobackup/rabr5411/poly/order_results/' date '-' HH '_' MM '.mat']);
end