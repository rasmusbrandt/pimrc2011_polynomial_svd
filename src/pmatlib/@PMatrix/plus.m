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
function r = plus(obj1, obj2)
    % PLUS Addition of PMatrices.

    % Make sure both arguments are PMatrices
    obj1 = PMatrix(obj1);
    obj2 = PMatrix(obj2);
    
    % Check dimensions
    s1 = size(obj1); s2 = size(obj2);
    if(s1(1:2) ~= s2(1:2))
        ME = MException('PMatrix:bad_dimensions', ...
                        'Dimension mismatch!');
        throw(ME);
    end

    % Get alignment parameters
    nb1 = obj1.const_ind - 1; nb2 = obj2.const_ind - 1;
    na1 = s1(3) - nb1 - 1;    na2 = s2(3) - nb2 - 1;

    % Align coefficient lists
    c1b = zeros(s1(1), s1(2), nb2-nb1);
    c1a = zeros(s1(1), s1(2), na2-na1);
    c2b = zeros(s1(1), s1(2), nb1-nb2);
    c2a = zeros(s1(1), s1(2), na1-na2);

    aligned_coefs1 = cat(3, c1b, obj1.coefs);
    aligned_coefs1 = cat(3, aligned_coefs1, c1a);
    aligned_coefs2 = cat(3, c2b, obj2.coefs);
    aligned_coefs2 = cat(3, aligned_coefs2, c2a);

    aligned_const_ind = max(obj1.const_ind, obj2.const_ind);

    % Return new PMatrix
    r = PMatrix(aligned_coefs1 + aligned_coefs2, aligned_const_ind);

    % Trim it
    r = r.rm_edge_zeros();
end % plus
