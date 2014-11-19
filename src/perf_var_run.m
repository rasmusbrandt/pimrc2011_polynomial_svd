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
% Performance VARIABLE
% Shows error performance for N = 2*L

clear;
addpath('pmatlib');

sim_id = datestr(now,'MMSSFFF');

fprintf('Performance VARIABLE.\n');

% Testing or real?
[~,hn] = system('hostname');
if strcmp(hn(1:3), 'sim')
    TESTING = false;
    
    fprintf('Choose parameter set:\n');
    fprintf('0. Testing\n');
    fprintf('21.  H ~ 3x3x7, exp_pdp_wt(2.5), taps(5:10:115), Nsim(100), x h \n');
    fprintf('51.  H ~ 3x3x15, exp_pdp_wt(5), taps(10:20:250), Nsim(100), 0.5 h \n');
    fprintf('3051.  H ~ 3x3x30, exp_pdp_wt(10), taps(10:35:500), Nsim(100), 3 h \n');
    
    pset = input('Select: ');
    
    nw = matlabpool('size'); if nw == 0;nw = 1;end
    fprintf('Starting production run with parameter set %d and %d workers.\n', pset, nw);
else
    TESTING = true;
    
    pset = 0;
    
    nw = matlabpool('size'); if nw == 0;nw = 1;end
    fprintf('Starting testing run with %d workers.\n', nw);
end

% Parameters
switch pset
    case 0
        % TESTING
        
        % General
        Nsim = 2;
        taps = 5:10:115;
        
        % Channel
        rows = 3;
        cols = 3;
        ch_type = 'exp_pdp_wt';
        lags = 7;
        WTau = 2.5;
       
        
    case 21
        % General
        Nsim = 100;
        taps = 5:10:115;
        
        % Channel
        rows = 3;
        cols = 3;
        ch_type = 'exp_pdp_wt';
        lags = 7;
        WTau = 2.5;
        
    case 51
        % General
        Nsim = 100;
        taps = 10:20:250;
        
        % Channel
        rows = 3;
        cols = 3;
        ch_type = 'exp_pdp_wt';
        lags = 15;
        WTau = 5;
        
    case 3051
        % General
        Nsim = 100;
        taps = 10:35:500;
        
        % Channel
        rows = 3;
        cols = 3;
        ch_type = 'exp_pdp_wt';
        lags = 30;
        WTau = 10;
end

% Storage variables
varRES_ERR1_U0 = zeros(Nsim,length(taps));
varRES_ERR1_Ua = zeros(Nsim,length(taps));
varRES_ERR1_V0 = zeros(Nsim,length(taps));
varRES_ERR1_Va = zeros(Nsim,length(taps));

varRES_ERR2_0 = zeros(Nsim,length(taps));
varRES_ERR2_a = zeros(Nsim,length(taps));

varRES_ERR3_0 = zeros(Nsim,length(taps));
varRES_ERR3_a = zeros(Nsim,length(taps));

tic;
for ii_Nsim = 1:Nsim
    fprintf('Starting sim run %d of %d.\n', ii_Nsim, Nsim);
    
    % Generate channel
    if strcmp(ch_type, 'exp_pdp_wt')
        H = channel_exp_pdp_wt(rows, cols, lags, WTau);
    elseif strcmp(ch_type, 'wssus')
        H = channel_wssus(rows, cols, lags);
    else
        error('Incorrect channel type!');
    end

    % Normalize channel
    H = (1/H.fnorm)*H;

    % Calculate approximations
    parfor ii_tap = 1:length(taps)
        % Select L
        Lu = taps(ii_tap); Lv = taps(ii_tap);
        
        % Get frequency domain SVDs for DFT-SVD
        N = taps(ii_tap);
        Hf = H.to_freq(N);       Df = zeros(rows,cols,N);
        Uf = zeros(rows,rows,N); Vf = zeros(cols,cols,N);
        for n = 1:N
            [Uf(:,:,n), Df(:,:,n), Vf(:,:,n)] = svd(Hf(:,:,n));
        end

        % Get time domain SVD 
        U = PMatrix(ifftshift(ifft(Uf,[],3),3), floor(N/2));
        D = PMatrix(ifftshift(ifft(Df,[],3),3), floor(N/2));
        V = PMatrix(ifftshift(ifft(Vf,[],3),3), floor(N/2));
        
        % Optimize DFT-SVD solution
        N = 2*taps(ii_tap);
        
        % Get frequency domain SVDs for phopt DFT-SVD
        Hf = H.to_freq(N);       Df = zeros(rows,cols,N);
        Uf = zeros(rows,rows,N); Vf = zeros(cols,cols,N);
        for n = 1:N
            [Uf(:,:,n), Df(:,:,n), Vf(:,:,n)] = svd(Hf(:,:,n));
        end
        [Ua, Da, Va] = palomar_bfmtx2adv(Uf,Df,Vf,Lu,Lv);
        Da.coefs = ifftshift(Da.coefs,3); Da.const_ind = N/2 + 1;
        
        % Unitarity errors
        varRES_ERR1_U0(ii_Nsim,ii_tap) = fnorm(U.ph*U - eye(rows,rows))/fnorm(PMatrix(eye(rows,rows)));
        varRES_ERR1_Ua(ii_Nsim,ii_tap) = fnorm(Ua.ph*Ua - eye(rows,rows))/fnorm(PMatrix(eye(rows,rows)));
        varRES_ERR1_V0(ii_Nsim,ii_tap) = fnorm(V.ph*V - eye(rows,rows))/fnorm(PMatrix(eye(rows,rows)));
        varRES_ERR1_Va(ii_Nsim,ii_tap) = fnorm(Va.ph*Va - eye(cols,cols))/fnorm(PMatrix(eye(cols,cols)));
        
        % Decomposition errors
        varRES_ERR2_0(ii_Nsim,ii_tap) = fnorm(H - U*D*V.ph);
        varRES_ERR2_a(ii_Nsim,ii_tap) = fnorm(H - Ua*Da*Va.ph);
        
        % Effective channel error
        varRES_ERR3_0(ii_Nsim,ii_tap) = fnorm_od(U.ph*H*V);
        varRES_ERR3_a(ii_Nsim,ii_tap) = fnorm_od(Ua.ph*H*Va);
    end
end
toc

if ~TESTING
    mkdir('/nobackup/rabr5411/poly/perf_var_results');
    cl = clock; 
    HH = num2str(cl(4)); if length(HH) == 1;HH = ['0' HH]; end
    MM = num2str(cl(5)); if length(MM) == 1;MM = ['0' MM]; end
    save(['/nobackup/rabr5411/poly/perf_var_results/' date '-' HH '_' MM '-pset' num2str(pset) '-id' num2str(sim_id) '.mat']);
end