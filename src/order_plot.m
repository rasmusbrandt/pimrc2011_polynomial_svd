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

figure;
plot(1:L,mean(results1), '-b', ...
     1:L,mean(results_trunc16), '--g', ...
     1:L,mean(results_trunc13), '-.r', ...
     1:L,mean(results2), '-c', ...
     1:L,mean(results_trunc26), '--m', ...
     1:L,mean(results_trunc23), '-.y');
h = xlabel('Number of taps in $\pmat{H}$'); set(h,'FontSize',text_size);
h = ylabel('Number of taps in $\pmat{U}$'); set(h,'FontSize',text_size);
h = title(['$' num2str(rows) '\times' num2str(cols) '$ channel. $\epsilon_r = 10^{-3}$']); set(h,'FontSize',text_size);
h = gca; set(h,'FontSize',axis_size);
legend('PSVD', 'Trunc. PSVD ($\mu = 10^{-6}$)', 'Trunc. PSVD ($\mu = 10^{-3}$)', ...
       'PSVD', 'Trunc. PSVD ($\mu = 10^{-6}$)', 'Trunc. PSVD ($\mu = 10^{-3}$)', ...
       'Location', 'NorthWest');
print('-depsc', 'order.eps');