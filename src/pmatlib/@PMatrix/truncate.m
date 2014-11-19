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
function r = truncate(obj, mu, Tstart, Tstop)
    % TRUNCATE Removes outer coefficient matrices according to the appendix
    % of Foster, McWhirter, Davies and Chambers
    
    % Consistency check
    if mu < 0
        ME = MException('PMatrix:bad_mu', 'Please provide positive mu');
        throw(ME);
    elseif mu > 1
        ME = MException('PMatrix:bad_mu', ...
                        'Please provide mu less than one');
        throw(ME);
    end
    
    % Size
    s = size(obj.coefs);
    
    if ndims(obj.coefs) < 3
        % Scalar matrix
        r = obj;
        return;
    end
    
    % Total F-norm squared of the matrix
    F = (obj.fnorm())^2;
    
    if nargin == 2
        % Find T1
        for T1 = s(3):-1:1
            step1 = abs(obj.coefs(:,:,T1:end)).^2;
            step2 = sum(sum(sum(step1)));

            if step2/F > mu/2
                break;
            end
        end

        % Find T2
        for T2 = 1:s(3)
            step1 = abs(obj.coefs(:,:,1:T2)).^2;
            step2 = sum(sum(sum(step1)));

            if step2/F > mu/2
                break;
            end
        end
    else
        T1 = Tstop;
        T2 = Tstart;

        if T2 < 1
            T2 = 1;
        end
        
        if T1 > s(3)
            T1 = s(3);
        end
    end
    
    % Consistency check
    if T1 < T2
        ME = MException('PMatrix:bad_truncation', ...
                        'mu was specified too big');
        throw(ME);
    end
    
    % Truncate!
    obj.coefs = obj.coefs(:,:,T2:T1);
    obj.const_ind = obj.const_ind - T2 + 1;
    
    r = obj;
end