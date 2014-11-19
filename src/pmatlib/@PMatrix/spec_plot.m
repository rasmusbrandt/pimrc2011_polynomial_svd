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
function spec_plot(obj, den, ax_val, xlbl, ylbl, ttl)
    % SPEC_PLOT Frequency plot
    
    % Consistency check
    N = 2^15;
    
    % Variables
    ss = size(obj);
    mid_row = floor(ss(1)/2) + 1;
    mid_col = floor(ss(2)/2) + 1;
    
    % Get frequency representation
    Hf = obj.freq(N);
    
    % Denominator (if given)
    if nargin > 1 && ~isempty(den)
        denf = squeeze(den.freq(N));
    end
    
    % Plot power gains
    figure; ind = 0;
    for i = 1:ss(1)
        for j = 1:ss(2)
            % Increment index
            ind = ind + 1;
            
            % Reshape vector
            h = squeeze(Hf(i,j,:));
            
            % Divide by denominator if given
            if nargin > 1 && ~isempty(den)
                h = h./denf;
            end
            
            % Plot it
            subplot(ss(1),ss(2),ind); 
            plot(linspace(0,1,N),20*log10(abs(h)), 'k');
            ax = axis;
            if nargin > 2
                axis([0 1 ax_val]);
            else
                axis([0 1 ax(3:4)]);
            end
            
            % Display text
            if j == mid_col
                switch i
                    case 1
                        if nargin <= 5 || isempty(ttl)
                            title('Frequency response');
                        else
                            title(ttl);
                        end
                    case ss(1)
                        if nargin <= 3 || isempty(xlbl)
                            xlabel('Normalized frequency');
                        else
                            xlabel(xlbl);
                        end
                end
            end

            if (i == mid_row) && (j == 1)
                if nargin <= 4 || isempty(ylbl)
                    ylabel('Power gain (dB)');
                else
                    ylabel(ylbl);
                end
            end
        end
    end
    
    % Don't plot phase response
    return;
    
    % Plot phases
    figure; ind = 0;
    for i = 1:ss(1)
        for j = 1:ss(2)
            % Increment index
            ind = ind + 1;
            
            % Reshape vector
            h = squeeze(Hf(i,j,:));
            
            % Divide by denominator if given
            if nargin == 2
                h = h./denf;
            end
            
            % Plot it
            subplot(ss(1),ss(2),ind); 
            plot(linspace(0,1,N),unwrap(angle(h)));
            ax = axis; axis([0 1 ax(3:4)]);
            
            % Display text
            if j == mid_col
                switch i
                    case 1
                        title('Frequency response');
                    case ss(1)
                        xlabel('Normalized frequency');
                end
            end

            if (i == mid_row) && (j == 1)
                ylabel('Phase (rad))');
            end
        end
    end
end