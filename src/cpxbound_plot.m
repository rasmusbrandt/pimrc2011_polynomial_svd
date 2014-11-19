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
text_size = 8;
axis_size = 6;

% PSVD
m_psvd = 2*mean(psvd);

% FIXED DFT
dftlengths = [2^8 2^9 2^10 2^11];
dftsvd_x = linspace(0,max(lags));
dftsvd_y = zeros(100,length(dftlengths));

for ii_dft = 1:length(dftlengths)
    L = dftlengths(ii_dft);
    tmp_svd = L*rows*rows*cols + L*rows*cols*cols + L*cols*cols*cols;
    tmp_fft_chan = L*rows*cols*log2(L);

    dftsvd_y(:,ii_dft) = (tmp_svd + tmp_fft_chan)*ones(100,1);
end

% DFT-SVD with phase optimization
dftsvd_N = 2*dftsvd_L;
Cinre = dftsvd_N;

% Channel 1
tmp_fft_chan = dftsvd_N*rows*cols.*log2(dftsvd_N); 
tmp_svd = dftsvd_N*(4*rows*rows*cols + 8*rows*cols*cols + 9*cols*cols*cols);
tmp_phopt = Cinre.*dftsvd_Cn.*dftsvd_N;
m_dftsvd = mean(tmp_fft_chan + tmp_svd + tmp_phopt);

%dftsvd_phopt1 = mean(dftsvd_Cn(:,:,1).*dftsvd_L(:,:,1).^2 + dftsvd_L(:,:,1)*rows*cols.*log2(dftsvd_L(:,:,1)) + dftsvd_L(:,:,1)*(rows*rows*cols + rows*cols*cols + cols*cols*cols));

% Channel 2
%dftsvd_phopt2 = mean(dftsvd_Cn(:,:,2).*dftsvd_L(:,:,2).^2 + dftsvd_L(:,:,1)*rows*cols.*log2(dftsvd_L(:,:,1)) + dftsvd_L(:,:,2)*(rows*rows*cols + rows*cols*cols + cols*cols*cols));

% Plot it
figure;
semilogy(lags,m_psvd(:,:,1),'-o',lags,m_psvd(:,:,2),'-v', ...
         lags,m_dftsvd(:,:,1),'-s',lags,m_dftsvd(:,:,2),'-^'); %,dftsvd_x,dftsvd_y);
h = xlabel('Number of taps in H(n)'); set(h,'FontSize',text_size);
h = ylabel('Number of flops'); set(h,'FontSize',text_size);
h = title(sprintf('$%d \\times %d$ channel. $\\epsilon_r = 10^{%d}$, $\\mu = 10^{%d}$', ...
    rows,cols,log10(eps_r),log10(mu))); set(h,'FontSize',text_size);
h = gca; set(h,'FontSize',axis_size);
%xlim([10 lags(end)]); ylim([1e6 0.4e8]);
legend('PSVD lower bound, $W\tau = 5$', 'PSVD lower bound, $W\tau = 10$', ... 
       'Ph. Opt. DFT-SVD, $W\tau = 5$', 'Ph. Opt. DFT-SVD, $W\tau = 10$',  ...
       'Location', [0.55 0.5 0.1 0.1]);
%text(5,2.7e4,['DFT-SVD ' num2str(dftlengths(1))],'FontSize',text_size);
%text(10,6e4,['DFT-SVD ' num2str(dftlengths(2))],'FontSize',text_size);
%text(15,1.2e5,['DFT-SVD ' num2str(dftlengths(3))],'FontSize',text_size);
%text(20,2.5e5,['DFT-SVD ' num2str(dftlengths(4))],'FontSize',text_size);
print('-depsc','cpxbound.eps')