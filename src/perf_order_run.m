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
clear;
addpath('pmatlib');

% Testing or real?
[~,hn] = system('hostname');
if strcmp(hn(1:3), 'sim')
    TESTING = false;
    
    Nsim = 100;
    rows = 3;
    cols = 3;
    L1 = 7;  WTau1 = 2.5;
    L2 = 15; WTau2 = 5;

    eps_r = 1e-3;
    mu = 1e-6;
    rho = 100;
    MaxSweeps = 1000;
    MaxPSVDIter = 1000;

    nw = matlabpool('size'); if nw == 0;nw = 1;end
    fprintf('Starting production run with %d workers.\n', nw);
else
    TESTING = true;
    
    Nsim = 2;
    rows = 3;
    cols = 3;
    L1 = 7;  WTau1 = 2.5;
    L2 = 15; WTau2 = 5;

    eps_r = 1e-3;
    mu = 1e-6;
    rho = 100;
    MaxSweeps = 1000;
    MaxPSVDIter = 1000;

    
    nw = matlabpool('size'); if nw == 0;nw = 1;end
    fprintf('Starting testing run with %d workers.\n', nw);
end

results1U = zeros(Nsim,1);
results2U = zeros(Nsim,1);
results1V = zeros(Nsim,1);
results2V = zeros(Nsim,1);

parfor ii_Nsim = 1:Nsim
    fprintf('Starting sim %d out of %d\n', ii_Nsim, Nsim);
    
    % L1, WTau2
    H = channel_exp_pdp_wt(rows,cols,L1,WTau1);
    [U_ph,Dp,V_ph] = psvd_mod(H, eps_r, mu, rho, MaxSweeps, MaxPSVDIter);

    Ut6 = U_ph.truncate(mu); 
    Vt6 = V_ph.truncate(mu);

    results1U(ii_Nsim) = size(Ut6.coefs,3);
    results1V(ii_Nsim) = size(Vt6.coefs,3);

    % L2, WTau2
    H = channel_exp_pdp_wt(rows,cols,L2,WTau2);
    [U_ph,Dp,V_ph] = psvd_mod(H, eps_r, mu, rho, MaxSweeps, MaxPSVDIter);

    Ut6 = U_ph.truncate(mu);
    Vt6 = V_ph.truncate(mu);

    results2U(ii_Nsim) = size(Ut6.coefs,3);
    results2V(ii_Nsim) = size(Vt6.coefs,3);
end

if ~TESTING
    mkdir('/nobackup/rabr5411/poly/perf_order_results');
    cl = clock; 
    HH = num2str(cl(4)); if length(HH) == 1;HH = ['0' HH]; end
    MM = num2str(cl(5)); if length(MM) == 1;MM = ['0' MM]; end
    save(['/nobackup/rabr5411/poly/perf_order_results/' date '-' HH '_' MM '.mat']);
end