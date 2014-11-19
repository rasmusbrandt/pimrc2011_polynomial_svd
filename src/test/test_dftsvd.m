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

rows = 3;
cols = 3;
lags = 15;
WT = 5;
N = 2048;

H = channel_exp_pdp_wt(rows, cols, lags, WT);
H = (1/H.fnorm)*H;

Hf = H.to_freq(N);       Df = zeros(rows,cols,N);
Uf = zeros(rows,rows,N); Vf = zeros(cols,cols,N);
for n = 1:N
    [Uf(:,:,n), Df(:,:,n), Vf(:,:,n)] = svd(Hf(:,:,n));
end

U = PMatrix(ifftshift(ifft(Uf,[],3),3), N/2 + 1);
D = PMatrix(ifftshift(ifft(Df,[],3),3), N/2 + 1);
V = PMatrix(ifftshift(ifft(Vf,[],3),3), N/2 + 1);

Lu = floor(N/3); Lv = floor(N/3);
[Ua, Da, Va, ~, ~, ~, Xua, Xva] = palomar_bfmtx2adv(Uf,Df,Vf,Lu,Lv);

xua = PMatrix(ifft(Xua,[],3)); xva = PMatrix(ifft(Xva,[],3));
xua_shift = xua;
xua_shift.coefs = ifftshift(xua_shift.coefs,3); xua_shift.const_ind = N/2 + 1;
xva_shift = xva;
xva_shift.coefs = ifftshift(xva_shift.coefs,3); xva_shift.const_ind = N/2 + 1;
truA = xua_shift.ph*D*xva_shift;
    
%fnorm(H - U*D*V.ph)
%fnorm_od(U.ph*H*V)
%fnorm(U.ph*U - eye(rows,cols))/fnorm(PMatrix(eye(rows,cols)))
fnorm(truA - Ua.ph*H*Va)
