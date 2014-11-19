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
% Polynomial extension project, simulation CPXBOUND

clear;
addpath('pmatlib');

sim_id = datestr(now,'MMSSFFF');

fprintf('Poly extension CPXBOUND.\n');

% Testing or real?
[~,hn] = system('hostname');
if strcmp(hn(1:3), 'sim')
    TESTING = false;
    
    fprintf('Choose parameter set:\n');
    fprintf('0. Testing\n');
    fprintf('11. Prod, 2:5:32 3.3h\n');
    fprintf('21. Prod, 2:2:32 7.6h\n');
    fprintf('31. Prod, 2:5:32 Wtau = [2.5,5]\n');
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
        rows = 3;
        cols = 3;
        ch_type = 'exp_pdp_wt';
        lags = 2:5:32;
        WTaus = [5 10];

        eps_r = 1e-3;
        mu = 1e-6;
        
        rho = 100;
        MaxSweeps = 1000;
        MaxPSVDIter = 1000;

        Nsim = 2;
        
        
    case 11
        rows = 3;
        cols = 3;
        ch_type = 'exp_pdp_wt'; 
        lags = 2:5:32;
        WTaus = [5 10];

        eps_r = 1e-3;
        mu = 1e-6;
        
        rho = 100;
        MaxSweeps = 1000;
        MaxPSVDIter = 1000;
        
        Nsim = 100;   
      
        
    case 21
            
        rows = 3;
        cols = 3;
        ch_type = 'exp_pdp_wt'; 
        lags = 2:2:32;
        WTaus = [5 10];

        eps_r = 1e-3;
        mu = 1e-6;
        
        rho = 100;
        MaxSweeps = 1000;
        MaxPSVDIter = 1000;
        
        Nsim = 100;
        
        
    case 31
        rows = 3;
        cols = 3;
        ch_type = 'exp_pdp_wt'; 
        lags = 2:5:32;
        WTaus = [2.5 5];

        eps_r = 1e-3;
        mu = 1e-6;
        
        rho = 100;
        MaxSweeps = 1000;
        MaxPSVDIter = 1000;
        
        Nsim = 100;   
      
end


% Storage variables
psvd = zeros(Nsim,length(lags),length(WTaus));
dftsvd_L = zeros(Nsim,length(lags),length(WTaus));
dftsvd_Cn = zeros(Nsim,length(lags),length(WTaus));

tic;
for ii_lags = 1:length(lags)
    fprintf('Starting lag %d of %d.\n', ii_lags, length(lags));
    
    for ii_wtau = 1:length(WTaus)

        % Calculate approximantions
        parfor ii_Nsim = 1:Nsim
            H = PMatrix(0);

            % Generate channel
            if strcmp(ch_type, 'exp_pdp_wt')
                H = channel_exp_pdp_wt(rows, cols, lags(ii_lags), WTaus(ii_wtau));
            elseif strcmp(ch_type, 'wssus')
                H = channel_wssus(rows, cols, lags);
            else
                error('Incorrect channel type!');
            end

            % Normalize channel
            H = (1/H.fnorm)*H;

            % Get PSVD
            [U_ph,Dp,V_ph,cpxdata] = ...
                psvd_mod(H, eps_r, mu, rho, MaxSweeps, MaxPSVDIter);

            % Complexity lower bound
            psvd(ii_Nsim,ii_lags,ii_wtau) = 4*cols*cpxdata.SumKf;
            
            % DFT-SVD optimization complexity
            L = lags(ii_lags); % Starting L
            while true
                N = 2*L;
                
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

                % Optimize DFT-SVD solution
                [Ua, Da, Va, ~, ~, ~, ~, ~, Cn] = palomar_bfmtx2adv(Uf,Df,Vf,L,L);

                % Effective channel error
                eff_chan = fnorm_od(Ua.ph*H*Va);
                
                if eff_chan^2 < eps_r
                    dftsvd_L(ii_Nsim,ii_lags,ii_wtau) = L;
                    dftsvd_Cn(ii_Nsim,ii_lags,ii_wtau) = Cn;
                    break;
                else
                    L = L + 16;
                end
            end
        end
    end
end
toc

if ~TESTING
    mkdir('/nobackup/rabr5411/poly/cpxbound_results');
    cl = clock; 
    HH = num2str(cl(4)); if length(HH) == 1;HH = ['0' HH]; end
    MM = num2str(cl(5)); if length(MM) == 1;MM = ['0' MM]; end
    save(['/nobackup/rabr5411/poly/cpxbound_results/' date '-' HH '_' MM '-pset' num2str(pset) '-id' num2str(sim_id) '.mat']);
end