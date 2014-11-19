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
% Performance FIXED
% Shows error performance for fixed algorithm parameters

clear;
addpath('pmatlib');

sim_id = datestr(now,'MMSSFFF');

fprintf('Performance FIXED.\n');

% Testing or real?
[~,hn] = system('hostname');
if strcmp(hn(1:3), 'sim')
    TESTING = false;
    
    fprintf('Choose parameter set:\n');
    fprintf('0. Testing\n');
    fprintf('21.  H ~ 3x3x7, exp_pdp_wt(2.5), eps = 1e-3, taps(5:10:115), N = 125,  Nsim(100), x  h \n');
    fprintf('51.  H ~ 3x3x15, exp_pdp_wt(5), eps = 1e-3, taps(10:20:250), N = 512,  Nsim(100), 5  h \n');
    fprintf('101. H ~ 3x3x15, exp_pdp_wt(5), eps = 1e-3, taps(10:20:250), N = 1024, Nsim(50),  8 h \n');
    fprintf('102. H ~ 3x3x15, exp_pdp_wt(5), eps = 1e-3, taps(10:20:250), N = 1024, Nsim(100), 16 h \n');
    fprintf('3051.  H ~ 3x3x30, exp_pdp_wt(10), eps = 1e-3, taps(10:35:500), N = 512,  Nsim(100), 15  h \n');
    fprintf('30101. H ~ 3x3x30, exp_pdp_wt(10), eps = 1e-3, taps(10:35:500), N = 1024, Nsim(50),  15 h \n');
    fprintf('30102. H ~ 3x3x30, exp_pdp_wt(10), eps = 1e-3, taps(10:35:500), N = 1024, Nsim(100), 30 h \n');
    
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
        N = 115;
        taps = 5:10:115;
        
        % Channel
        rows = 3;
        cols = 3;
        ch_type = 'exp_pdp_wt';
        lags = 15;
        WTau = 5;

        % PSVD
        eps_r = 1e-3;
        mu = 1e-6;
        rho = 100;
        MaxSweeps = 100;
        MaxPSVDIter = 100;
       
        
    case 21
        % General
        Nsim = 100;
        N = 125;
        taps = 5:10:115;
        
        % Channel
        rows = 3;
        cols = 3;
        ch_type = 'exp_pdp_wt';
        lags = 7;
        WTau = 2.5;

        % PSVD
        eps_r = 1e-3;
        mu = 1e-6;
        rho = 100;
        MaxSweeps = 100;
        MaxPSVDIter = 100;
        
    case 51
        % General
        Nsim = 100;
        N = 512;
        taps = 10:20:250;
        
        % Channel
        rows = 3;
        cols = 3;
        ch_type = 'exp_pdp_wt';
        lags = 15;
        WTau = 5;

        % PSVD
        eps_r = 1e-3;
        mu = 1e-6;
        rho = 100;
        MaxSweeps = 100;
        MaxPSVDIter = 100;
        
        
    case 101
        % General
        Nsim = 50;
        N = 1024;
        taps = 10:20:250;
        
        % Channel
        rows = 3;
        cols = 3;
        ch_type = 'exp_pdp_wt';
        lags = 15;
        WTau = 5;

        % PSVD
        eps_r = 1e-3;
        mu = 1e-6;
        rho = 100;
        MaxSweeps = 100;
        MaxPSVDIter = 100;
        
    case 102
        % General
        Nsim = 100;
        N = 1024;
        taps = 10:20:250;
        
        % Channel
        rows = 3;
        cols = 3;
        ch_type = 'exp_pdp_wt';
        lags = 15;
        WTau = 5;

        % PSVD
        eps_r = 1e-3;
        mu = 1e-6;
        rho = 100;
        MaxSweeps = 100;
        MaxPSVDIter = 100;
        
        
    case 3051
        % General
        Nsim = 100;
        N = 512;
        taps = 10:35:500;
        
        % Channel
        rows = 3;
        cols = 3;
        ch_type = 'exp_pdp_wt';
        lags = 30;
        WTau = 10;

        % PSVD
        eps_r = 1e-3;
        mu = 1e-6;
        rho = 100;
        MaxSweeps = 100;
        MaxPSVDIter = 100;
        
        
    case 30101
        % General
        Nsim = 50;
        N = 1024;
        taps = 10:35:500;
        
        % Channel
        rows = 3;
        cols = 3;
        ch_type = 'exp_pdp_wt';
        lags = 30;
        WTau = 10;

        % PSVD
        eps_r = 1e-3;
        mu = 1e-6;
        rho = 100;
        MaxSweeps = 100;
        MaxPSVDIter = 100;
        
    case 30102
        % General
        Nsim = 100;
        N = 1024;
        taps = 10:35:500;
        
        % Channel
        rows = 3;
        cols = 3;
        ch_type = 'exp_pdp_wt';
        lags = 30;
        WTau = 10;

        % PSVD
        eps_r = 1e-3;
        mu = 1e-6;
        rho = 100;
        MaxSweeps = 100;
        MaxPSVDIter = 100;
        
end

% Storage variables
RES_ERR1_U0 = zeros(Nsim,length(taps));
RES_ERR1_Up = zeros(Nsim,length(taps));
RES_ERR1_Ua = zeros(Nsim,length(taps));
RES_ERR1_Ub = zeros(Nsim,length(taps));
RES_ERR1_Us = zeros(Nsim,length(taps));
RES_ERR1_V0 = zeros(Nsim,length(taps));
RES_ERR1_Vp = zeros(Nsim,length(taps));
RES_ERR1_Va = zeros(Nsim,length(taps));
RES_ERR1_Vb = zeros(Nsim,length(taps));
RES_ERR1_Vs = zeros(Nsim,length(taps));

RES_ERR2_0 = zeros(Nsim,length(taps));
RES_ERR2_p = zeros(Nsim,length(taps));
RES_ERR2_a = zeros(Nsim,length(taps));
RES_ERR2_b = zeros(Nsim,length(taps));
RES_ERR2_s = zeros(Nsim,length(taps));

RES_ERR3_0 = zeros(Nsim,length(taps));
RES_ERR3_p = zeros(Nsim,length(taps));
RES_ERR3_a = zeros(Nsim,length(taps));
RES_ERR3_b = zeros(Nsim,length(taps));
RES_ERR3_s = zeros(Nsim,length(taps));

RES_ORDER_Up = zeros(Nsim,1);
RES_ORDER_Dp = zeros(Nsim,1);
RES_ORDER_Vp = zeros(Nsim,1);

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

    % Get PSVD
    [U_ph,Dp,V_ph] = psvd_mod(H, eps_r, mu, rho, MaxSweeps, MaxPSVDIter);
    Up = U_ph.ph; Vp = V_ph.ph;

    % Truncate with mu, to remove some zeros (speed up simulations)
    Up = Up.truncate(mu); Dp = Dp.truncate(mu); Vp = Vp.truncate(mu);
    
    % Get frequency domain SVDs
    Hf = H.to_freq(N);       Df = zeros(rows,cols,N);
    Uf = zeros(rows,rows,N); Vf = zeros(cols,cols,N);
    for n = 1:N
        [Uf(:,:,n), Df(:,:,n), Vf(:,:,n)] = svd(Hf(:,:,n));
    end
    
    % Get time domain SVD 
    U = PMatrix(ifftshift(ifft(Uf,[],3),3), N/2 + 1);
    D = PMatrix(ifftshift(ifft(Df,[],3),3), N/2 + 1);
    V = PMatrix(ifftshift(ifft(Vf,[],3),3), N/2 + 1);
    
    % Calculate Npsvd
    [~,~,Nu] = size(Up.coefs); [~,~,Nd] = size(Dp.coefs); [~,~,Nv] = size(Vp.coefs); 
    Npsvd = Nu + Nd + Nv - min(min(Nu,Nd),Nv) - 1;
    
    % Consistency check, so that we can truncate without a problem!
    if Npsvd <= taps(end)
        Npsvd = taps(end);
    end
    
    % Store orders
    RES_ORDER_Up(ii_Nsim) = Nu;
    RES_ORDER_Dp(ii_Nsim) = Nd;
    RES_ORDER_Vp(ii_Nsim) = Nv;
    
    % Get non-modified PSVD frequency domain
    Upf = Up.to_freq(Npsvd); Vpf = Vp.to_freq(Npsvd); Dpf = Dp.to_freq(Npsvd);

    % Calculate approximations
    parfor ii_tap = 1:length(taps)
        % Select L
        Lu = taps(ii_tap); Lv = taps(ii_tap);
        
        % Optimize DFT-SVD solution
        [Ua, Da, Va] = palomar_bfmtx2adv(Uf,Df,Vf,Lu,Lv);
        Da.coefs = ifftshift(Da.coefs,3); Da.const_ind = N/2 + 1;
        
        % Optimize PSVD solution
        [Ub, Db, Vb] = palomar_bfmtx2adv(Upf,Dpf,Vpf,Lu,Lv,ceil(Lu/2),ceil(Lv/2)); 
        Db.coefs = ifftshift(Db.coefs,3); Db.const_ind = N/2 + 1;

        % Stupid truncation
        Us = stupid_truncate(Up, Lu);
        Ds = Dp;
        Vs = stupid_truncate(Vp, Lv);
        
        % Unitarity errors
        RES_ERR1_U0(ii_Nsim,ii_tap) = fnorm(U.ph*U - eye(rows,rows))/fnorm(PMatrix(eye(rows,rows)));
        RES_ERR1_Up(ii_Nsim,ii_tap) = fnorm(Up.ph*Up - eye(rows,rows))/fnorm(PMatrix(eye(rows,rows)));
        RES_ERR1_Ua(ii_Nsim,ii_tap) = fnorm(Ua.ph*Ua - eye(rows,rows))/fnorm(PMatrix(eye(rows,rows)));
        RES_ERR1_Ub(ii_Nsim,ii_tap) = fnorm(Ub.ph*Ub - eye(rows,rows))/fnorm(PMatrix(eye(rows,rows)));
        RES_ERR1_Us(ii_Nsim,ii_tap) = fnorm(Us.ph*Us - eye(rows,rows))/fnorm(PMatrix(eye(rows,rows)));
        RES_ERR1_V0(ii_Nsim,ii_tap) = fnorm(V.ph*V - eye(rows,rows))/fnorm(PMatrix(eye(rows,rows)));
        RES_ERR1_Vp(ii_Nsim,ii_tap) = fnorm(Vp.ph*Vp - eye(cols,cols))/fnorm(PMatrix(eye(cols,cols)));
        RES_ERR1_Va(ii_Nsim,ii_tap) = fnorm(Va.ph*Va - eye(cols,cols))/fnorm(PMatrix(eye(cols,cols)));
        RES_ERR1_Vb(ii_Nsim,ii_tap) = fnorm(Vb.ph*Vb - eye(cols,cols))/fnorm(PMatrix(eye(cols,cols)));
        RES_ERR1_Vs(ii_Nsim,ii_tap) = fnorm(Vs.ph*Vs - eye(cols,cols))/fnorm(PMatrix(eye(cols,cols)));
        
        % Decomposition errors
        RES_ERR2_0(ii_Nsim,ii_tap) = fnorm(H - U*D*V.ph);
        RES_ERR2_p(ii_Nsim,ii_tap) = fnorm(H - Up*Dp*Vp.ph);
        RES_ERR2_a(ii_Nsim,ii_tap) = fnorm(H - Ua*Da*Va.ph);
        RES_ERR2_b(ii_Nsim,ii_tap) = fnorm(H - Ub*Db*Vb.ph);
        RES_ERR2_s(ii_Nsim,ii_tap) = fnorm(H - Us*Ds*Vs.ph);
        
        % Effective channel error
        RES_ERR3_0(ii_Nsim,ii_tap) = fnorm_od(U.ph*H*V);
        RES_ERR3_p(ii_Nsim,ii_tap) = fnorm_od(Up.ph*H*Vp);
        RES_ERR3_a(ii_Nsim,ii_tap) = fnorm_od(Ua.ph*H*Va);
        RES_ERR3_b(ii_Nsim,ii_tap) = fnorm_od(Ub.ph*H*Vb);
        RES_ERR3_s(ii_Nsim,ii_tap) = fnorm_od(Us.ph*H*Vs);

    end
end
toc

if ~TESTING
    mkdir('/nobackup/rabr5411/poly/perf_fixed_results');
    cl = clock; 
    HH = num2str(cl(4)); if length(HH) == 1;HH = ['0' HH]; end
    MM = num2str(cl(5)); if length(MM) == 1;MM = ['0' MM]; end
    save(['/nobackup/rabr5411/poly/perf_fixed_results/' date '-' HH '_' MM '-pset' num2str(pset) '-id' num2str(sim_id) '.mat']);
end