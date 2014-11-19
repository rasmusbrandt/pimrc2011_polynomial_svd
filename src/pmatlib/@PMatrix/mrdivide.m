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
function [q,r] = mrdivide(obj1, obj2)
    % MRDIVIDE Polynomial division. 
    %          ONLY WORKS FOR POLYNOMIALS WITH SCALAR COEFFICIENTS!!!
    
    if obj1.is1x1() && obj2.is1x1()
        % Trim the polynomials
        p1 = obj1.rm_edge_zeros(); p2 = obj2.rm_edge_zeros();
        sp1 = size(p1);            sp2 = size(p2);
        
        % Get coefficients in correct format
        c1 = squeeze(p1.coefs); c2 = squeeze(p2.coefs);
        
        % Perform long division
        [qc,rc] = deconv(c1,c2);
        
        % What fraction was pulled out? Only keep negative exponents!
        k_num = p1.ind2lead(sp1(3)); if k_num > 0, k_num = 0; end
        k_den = p2.ind2lead(sp2(3)); if k_den > 0, k_den = 0; end
        
        % Get quotient
        qci = length(qc) + k_num - k_den;
        q = PMatrix(reshape(qc, 1, 1, length(qc)), qci);
        
        % Get remainder
        rci = length(rc) + k_num - k_den;
        r = PMatrix(reshape(rc, 1, 1, length(rc)), rci);
        
        % Trim quotient and remainder
        r = r.rm_edge_zeros();
        q = q.rm_edge_zeros();
    else
        error('Not 1x1 polynomials!');
    end
end