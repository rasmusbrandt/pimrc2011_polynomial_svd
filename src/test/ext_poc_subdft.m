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
% PROOF OF CONCEPT

clear;

% Parameters
rows = 3;
cols = 3;


lags = 15;
WT = 5;

eps_r = 1e-2;
mu = 1e-6;
rho = 5;
MaxSweeps = 10;
MaxPSVDIter = 10;

Ls = round(linspace(1,lags,10));
ERR2 = zeros(2,length(Ls));
ERR3 = zeros(2,length(Ls));


% Generate channel
H = channel_exp_pdp_wt(rows, cols, lags, WT);

% Normalize channel
H = (1/H.fnorm)*H;

% Select L
N = lags;

% Get PSVD
[U_ph,Dp,V_ph] = psvd_mod(H, eps_r, mu, rho, MaxSweeps, MaxPSVDIter);
Up = U_ph.ph; Vp = V_ph.ph;
Up = Up.truncate(mu); Dp = Dp.truncate(mu); Vp = Vp.truncate(mu);

[~,~,Nu] = size(Up.coefs); [~,~,Nd] = size(Dp.coefs); [~,~,Nv] = size(Vp.coefs); 
N2 = Nu + Nd + Nv - min(min(Nu,Nd),Nv) - 1;
N = N2;

Upf = Up.to_freq(N2);    Vpf = Vp.to_freq(N2); Dpf = Dp.to_freq(N2);

% Get frequency domain SVDs
Hf = H.to_freq(N);       Df = zeros(rows,cols,N);
Uf = zeros(rows,rows,N); Vf = zeros(cols,cols,N);
for n = 1:N
    [Uf(:,:,n), Df(:,:,n), Vf(:,:,n)] = svd(Hf(:,:,n));
end

for ii_l = 1:length(Ls)
    ii_l
    L = Ls(ii_l); Lu = L; Lv = L;
    
    

    % Optimize phases
    [Ua,Da,Va] = palomar_bfmtx2adv(Uf,Df,Vf,Lu,Lv);
    [Uo, Vo] = palomar_bfmtx2(Uf,Vf,Lu,Lv);
		[Ut, Vt] = palomar_bfmtx2(Upf,Vpf,Lu,Lv);
		[Ub,Db,Vb] = palomar_bfmtx2adv(Upf,Dpf,Vpf,Lu,Lv);

    % Get frequency domain representations
    Uof = Uo.to_freq(N); Vof = Vo.to_freq(N);
    Utf = Ut.to_freq(N); Vtf = Vt.to_freq(N);
		Uaf = Ua.to_freq(N); Vaf = Va.to_freq(N); Daf = Da.to_freq(N);
		Ubf = Ub.to_freq(N); Vbf = Vb.to_freq(N); Dbf = Db.to_freq(N);

    % Unitarity errors

    % Decomposition errors
    ERR2_p = 0;
    ERR2_o = 0;
    ERR2_t = 0;
    ERR2_a = 0;
    ERR2_b = 0;
    ERR2_2o = 0;
    k = 0;
    for n = 1:N
        tru = Hf(:,:,n);

        k = k + norm(tru, 'fro')^2;

        %ERR2_p = ERR2_p + norm(Upf(:,:,n) * Dpf(:,:,n) * Vpf(:,:,n)' - tru, 'fro')^2;
        ERR2_o = ERR2_o + norm(Uof(:,:,n) * Df(:,:,n)  * Vof(:,:,n)' - tru, 'fro')^2;
        ERR2_t = ERR2_t + norm(Utf(:,:,n) * Dpf(:,:,n) * Vtf(:,:,n)' - tru, 'fro')^2;
        ERR2_a = ERR2_a + norm(Uaf(:,:,n) * Daf(:,:,n) * Vaf(:,:,n)' - tru, 'fro')^2;
        ERR2_b = ERR2_b + norm(Ubf(:,:,n) * Dbf(:,:,n) * Vbf(:,:,n)' - tru, 'fro')^2;

        %ERR2_2o = ERR2_2o + norm(U2of(:,:,n) * Df(:,:,n)  * V2of(:,:,n)' - tru, 'fro')^2;
    end
    %ERR2_p = (1/k)*ERR2_p;
    ERR2_o = (1/k)*ERR2_o;
    ERR2_t = (1/k)*ERR2_t;
    ERR2_a = (1/k)*ERR2_a;
    ERR2_b = (1/k)*ERR2_b;
    %ERR2_2o = (1/k)*ERR2_2o;


    % Effective channel errors
    ERR3_p = 0;
    ERR3_o = 0;
    ERR3_t = 0;
    ERR3_a = 0;
    ERR3_b = 0;
    k = 0;
    for n = 1:N
                tru = Df(:,:,n);

                k = k + norm(tru, 'fro')^2;

                %M = abs(Upf(:,:,n)' * Hf(:,:,n) * Vpf(:,:,n));  [~, Mind] = sort(diag(M), 1, 'descend');
                %ERR3_p = ERR3_p + norm(M(Mind,Mind) - Df(:,:,n), 'fro')^2;
                M = abs(Uof(:,:,n)' * Hf(:,:,n) * Vof(:,:,n));  [~, Mind] = sort(diag(M), 1, 'descend');
                ERR3_o = ERR3_o + norm(M(Mind,Mind) - Df(:,:,n), 'fro')^2;
                M = abs(Utf(:,:,n)' * Hf(:,:,n) * Vtf(:,:,n));  [~, Mind] = sort(diag(M), 1, 'descend');
                ERR3_t = ERR3_t + norm(M(Mind,Mind) - Df(:,:,n), 'fro')^2;
                M = abs(Uaf(:,:,n)' * Hf(:,:,n) * Vaf(:,:,n));  [~, Mind] = sort(diag(M), 1, 'descend');
                ERR3_a = ERR3_a + norm(M(Mind,Mind) - Df(:,:,n), 'fro')^2;
                M = abs(Ubf(:,:,n)' * Hf(:,:,n) * Vbf(:,:,n));  [~, Mind] = sort(diag(M), 1, 'descend');
                ERR3_b = ERR3_b + norm(M(Mind,Mind) - Df(:,:,n), 'fro')^2;
    end
    %ERR3_p = (1/k)*ERR3_p;
    ERR3_o = (1/k)*ERR3_o;
    ERR3_t = (1/k)*ERR3_t;
    ERR3_a = (1/k)*ERR3_a;
    ERR3_b = (1/k)*ERR3_b;
    
    ERR2(1,ii_l) = ERR2_a;
    %ERR2(2,ii_l) = ERR2_o;
    ERR2(2,ii_l) = ERR2_b;
    %ERR2(4,ii_l) = ERR2_t;
    
    ERR3(1,ii_l) = ERR3_a;
    %ERR3(2,ii_l) = ERR3_o;
    ERR3(2,ii_l) = ERR3_b;
    %ERR3(4,ii_l) = ERR3_t;
end

figure;
plot(Ls,ERR2);
legend('DFT-SVD', 'PSVD');
title(['Decomposition error, WT = ' num2str(WT)]);

figure;
plot(Ls,ERR3);
legend('DFT-SVD', 'PSVD');
title(['Magnitude eff.chan. error, WT = ' num2str(WT)]);

