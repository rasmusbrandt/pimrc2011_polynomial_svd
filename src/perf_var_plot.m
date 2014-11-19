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
xlims = [taps(1) taps(end)];
xticks = [10 50:50:taps(end)];
%ylims = [1e-3 1];

% Process the errors
m_ERR1_U0 = mean(varRES_ERR1_U0);
m_ERR1_Ua = mean(varRES_ERR1_Ua);
m_ERR1_V0 = mean(varRES_ERR1_V0);
m_ERR1_Va = mean(varRES_ERR1_Va);

m_ERR2_0 = mean(varRES_ERR2_0);
m_ERR2_a = mean(varRES_ERR2_a);

m_ERR3_0 = mean(varRES_ERR3_0);
m_ERR3_a = mean(varRES_ERR3_a);

% Unitarity errors
figure;
semilogy(taps, m_ERR1_U0, '.g-', ...
         taps, m_ERR1_Ua, 'om-');         
h = xlabel('Model order $L$'); set(h,'FontSize',text_size);
h = ylabel('Unitarity error'); set(h,'FontSize',text_size);
h = title(sprintf('$%d \\times %d$ channel, $%d$ taps, $W\\tau = %d$.', ...
    rows,cols,lags,WTau)); set(h,'FontSize',text_size);
h = gca; set(h,'FontSize',axis_size);
set(h, 'XTick', xticks);
xlim(xlims); ylim([1e-2 1]);
legend('DFT-SVD', 'Ph. Opt. DFT-SVD, $N = 2*L$',  ...
       'Location', 'NorthEast');
print('-depsc','unitarity.eps');

% Effective channel errors
figure;
semilogy(taps, m_ERR3_0, '.g-', ...
         taps, m_ERR3_a, 'om-');
h = xlabel('Model order $L$'); set(h,'FontSize',text_size);
h = ylabel('Effective channel off-diagonal error'); set(h,'FontSize',text_size);
h = title(sprintf('$%d \\times %d$ channel, $%d$ taps, $W\\tau = %d$.', ...
    rows,cols,lags,WTau)); set(h,'FontSize',text_size);
h = gca; set(h,'FontSize',axis_size);
set(h, 'XTick', xticks);
xlim(xlims); ylim([0.4e-2 1]);
legend('DFT-SVD', 'Ph. Opt. DFT-SVD, $N = 2*L$',  ...
       'Location', 'NorthEast');
print('-depsc','effective_channel.eps');