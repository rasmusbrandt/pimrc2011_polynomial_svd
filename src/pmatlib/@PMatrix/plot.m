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
function plot(obj, matlab_axis, xlbl, ylbl, ttl)
    % PLOT Displays a stemplot over the matrix FIR filters

    if nargin == 1
        matlab_axis = false;
    elseif nargin == 2
        % Do nothing
    elseif nargin == 5
        % Do nothing
    else
        % Incorrect number of arguments!
        error('Incorrect number of arguments!');
    end
    
    % Variables
    ss = size(obj);
    mid_row = floor(ss(1)/2) + 1;
    mid_col = floor(ss(2)/2) + 1;
    
    s = size(obj.coefs); k = 1;
    for i = 1:s(1)
        for j = 1:s(2)
            if length(s) > 2
                % Polynomial matrix
                lags = -obj.ind2lead(1:s(3));
            else
                % Scalar matrix
                lags = 0;
            end
            stems = fliplr(abs(squeeze(obj.coefs(i,j,:))));

            subplot(s(1), s(2), k);
            stem(lags, stems, 'k');
            if ~matlab_axis
                axis([lags(1)-1 lags(end)+1 0 1.1*max(max(max((abs(obj.coefs)))))]);
            end
            
            % FontSize
            set(gca, 'FontSize', 6);

            % Display text
            if nargin == 5
                if j == mid_col
                    switch i
                        case 1
                            title(ttl);
                        case ss(1)
                            hx = xlabel(xlbl);
                    end
                end

                if (i == mid_row) && (j == 1)
                    ylabel(ylbl);
                end
            end
            
            k = k + 1;
        end
    end
end