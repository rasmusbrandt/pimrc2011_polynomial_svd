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
rows = 2;
cols = 2;
lags = 50;
psi = 0.05;

eps_r = 1e-2;
mu = 1e-6;
rho = 5;
MaxSweeps = 10;
MaxPSVDIter = 10;

% Generate channel
H = channel_exp_pdp(rows, cols, lags, psi);

% Normalize channel
H = (1/H.fnorm)*H;

% Get PSVD
[U_ph,Dp,V_ph] = psvd_mod(H, eps_r, mu, rho, MaxSweeps, MaxPSVDIter);
Up = U_ph.ph; Vp = V_ph.ph;

% Select L
[~,~,Nu] = size(Up.coefs); [~,~,Nv] = size(Vp.coefs); N = max(Nu,Nv);
% Lu = ceil(alpha*Nu);
% Lv = ceil(alpha*Nv);
N = 512;
Lu = 45;
Lv = 45;

% Get frequency domain SVDs
Upf = Up.to_freq(N);    Vpf = Vp.to_freq(N); Dpf = Dp.to_freq(N);
Hf = H.to_freq(N);       Df = zeros(rows,cols,N);
Uf = zeros(rows,rows,N); Vf = zeros(cols,cols,N);
for n = 1:N
    [Uf(:,:,n), Df(:,:,n), Vf(:,:,n)] = svd(Hf(:,:,n));
end

% Optimize phases
[Ua,Da,Va] = palomar_bfmtx2adv(Uf,Df,Vf,Lu,Lv);
[Uo, Vo] = palomar_bfmtx2(Uf,Vf,Lu,Lv);
[Ut, Vt] = palomar_bfmtx2(Upf,Vpf,Lu,Lv);
[Ub,Db,Vb] = palomar_bfmtx2adv(Upf,Dpf,Vpf,Lu,Lv);

N2 = lags;
H2f = H.to_freq(N2);       D2f = zeros(rows,cols,N2);
U2f = zeros(rows,rows,N2); V2f = zeros(cols,cols,N2);
for n = 1:N2
    [U2f(:,:,n), D2f(:,:,n), V2f(:,:,n)] = svd(H2f(:,:,n));
end

[U2o, V2o] = palomar_bfmtx2(U2f,V2f,Lu,Lv);

% Stupid truncation
Us = stupid_truncate(Up, Lu);
Vs = stupid_truncate(Vp, Lv);

% Get frequency domain representations
Uof = Uo.to_freq(N); Vof = Vo.to_freq(N);
U2of = U2o.to_freq(N); V2of = V2o.to_freq(N);
Utf = Ut.to_freq(N); Vtf = Vt.to_freq(N);
Usf = Us.to_freq(N); Vsf = Vs.to_freq(N);
Uaf = Ua.to_freq(N); Vaf = Va.to_freq(N); Daf = Da.to_freq(N);
Ubf = Ub.to_freq(N); Vbf = Vb.to_freq(N); Dbf = Db.to_freq(N);

% Unitarity errors
ERR1_Up = 0;
ERR1_Uo = 0;
ERR1_Ut = 0;
ERR1_Vp = 0;
ERR1_Vo = 0;
ERR1_Vt = 0;
k = 0;
for n = 1:N
    tru = eye(rows,rows);
    
    k = k + norm(tru, 'fro')^2;
    
    ERR1_Up = ERR1_Up + norm(Upf(:,:,n)'*Upf(:,:,n) - tru,'fro')^2;
    ERR1_Uo = ERR1_Uo + norm(Uof(:,:,n)'*Uof(:,:,n) - tru,'fro')^2;
    ERR1_Ut = ERR1_Ut + norm(Utf(:,:,n)'*Utf(:,:,n) - tru,'fro')^2;
    ERR1_Vp = ERR1_Vp + norm(Vpf(:,:,n)'*Vpf(:,:,n) - tru,'fro')^2;
    ERR1_Vo = ERR1_Vo + norm(Vof(:,:,n)'*Vof(:,:,n) - tru,'fro')^2;
    ERR1_Vt = ERR1_Vt + norm(Vtf(:,:,n)'*Vtf(:,:,n) - tru,'fro')^2;
end
ERR1_Up = (1/k)*ERR1_Up;
ERR1_Uo = (1/k)*ERR1_Uo;
ERR1_Ut = (1/k)*ERR1_Ut;
ERR1_Vp = (1/k)*ERR1_Vp;
ERR1_Vo = (1/k)*ERR1_Vo;
ERR1_Vt = (1/k)*ERR1_Vt;

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
    
    ERR2_p = ERR2_p + norm(Upf(:,:,n) * Dpf(:,:,n) * Vpf(:,:,n)' - tru, 'fro')^2;
    ERR2_o = ERR2_o + norm(Uof(:,:,n) * Df(:,:,n)  * Vof(:,:,n)' - tru, 'fro')^2;
    ERR2_t = ERR2_t + norm(Utf(:,:,n) * Dpf(:,:,n) * Vtf(:,:,n)' - tru, 'fro')^2;
    ERR2_a = ERR2_a + norm(Uaf(:,:,n) * Daf(:,:,n) * Vaf(:,:,n)' - tru, 'fro')^2;
    ERR2_b = ERR2_b + norm(Ubf(:,:,n) * Dbf(:,:,n) * Vbf(:,:,n)' - tru, 'fro')^2;
    
    ERR2_2o = ERR2_2o + norm(U2of(:,:,n) * Df(:,:,n)  * V2of(:,:,n)' - tru, 'fro')^2;
end
ERR2_p = (1/k)*ERR2_p;
ERR2_o = (1/k)*ERR2_o;
ERR2_t = (1/k)*ERR2_t;
ERR2_a = (1/k)*ERR2_a;
ERR2_b = (1/k)*ERR2_b;
ERR2_2o = (1/k)*ERR2_2o;


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
            
            M = abs(Upf(:,:,n)' * Hf(:,:,n) * Vpf(:,:,n));  [~, Mind] = sort(diag(M), 1, 'descend');
            ERR3_p = ERR3_p + norm(M(Mind,Mind) - Df(:,:,n), 'fro')^2;
            M = abs(Uof(:,:,n)' * Hf(:,:,n) * Vof(:,:,n));  [~, Mind] = sort(diag(M), 1, 'descend');
            ERR3_o = ERR3_o + norm(M(Mind,Mind) - Df(:,:,n), 'fro')^2;
            M = abs(Utf(:,:,n)' * Hf(:,:,n) * Vtf(:,:,n));  [~, Mind] = sort(diag(M), 1, 'descend');
            ERR3_t = ERR3_t + norm(M(Mind,Mind) - Df(:,:,n), 'fro')^2;
            M = abs(Uaf(:,:,n)' * Hf(:,:,n) * Vaf(:,:,n));  [~, Mind] = sort(diag(M), 1, 'descend');
            ERR3_a = ERR3_a + norm(M(Mind,Mind) - Df(:,:,n), 'fro')^2;
            M = abs(Ubf(:,:,n)' * Hf(:,:,n) * Vbf(:,:,n));  [~, Mind] = sort(diag(M), 1, 'descend');
            ERR3_b = ERR3_b + norm(M(Mind,Mind) - Df(:,:,n), 'fro')^2;
end
ERR3_p = (1/k)*ERR3_p;
ERR3_o = (1/k)*ERR3_o;
ERR3_t = (1/k)*ERR3_t;
ERR3_a = (1/k)*ERR3_a;
ERR3_b = (1/k)*ERR3_b;

ERR3_2o = 0;
k = 0;
for n = 1:N
            tru = Df(:,:,n);

            k = k + norm(tru, 'fro')^2;
            
            M = abs(U2of(:,:,n)' * Hf(:,:,n) * V2of(:,:,n));  [~, Mind] = sort(diag(M), 1, 'descend');
            ERR3_2o = ERR3_2o + norm(M(Mind,Mind) - Df(:,:,n), 'fro')^2;
end
ERR3_2o = (1/k)*ERR3_2o;