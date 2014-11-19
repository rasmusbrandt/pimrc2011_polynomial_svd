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
function [x,n] = palomar_am(Q, x0)
    % Variables
    N = length(x0);
    MAXITER = 4000;
    CRIT = 1e-3;
    
    % Initial guess
    x = x0; x_old = zeros(size(x0));
    
    % Choose implementation
    if exist('palomar_am_loop', 'file');
        % Iterate
        n = 0;
        while (norm(x - x_old)/norm(x)) > CRIT && (n < MAXITER)
            n = n + 1;

            % Keep old x around
            x_old = x;

            x = palomar_am_loop(Q,x);
        end
    else
        warning('Matlab implementation is being used. There is a MEX implementation, but that crashes on the authors new computer at the time of this writing (Nov 2014).');
        
        % Iterate
        n = 0;
        while (norm(x - x_old)/norm(x)) > CRIT && (n < MAXITER)
            n = n + 1;

            % Keep old x around
            x_old = x;

            % Loop over all components
            for k = 1:N
                % Optimize this particular phase

    %             Qtkk = Q(k,:)*x - Q(k,k)*x(k);
    %             phi_k = angle(Qtkk) + pi;
    %             x(k) = exp(1i*phi_k);
                x(k) = exp(1i*(angle(Q(k,:)*x - Q(k,k)*x(k)) + pi));
            end
        end
    end
    
    if n == MAXITER
        warning('Max iterations reached.');
    end
end
