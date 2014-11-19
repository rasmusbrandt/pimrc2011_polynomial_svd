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
function r = circonv(obj1, obj2)
    % MTIMES Multiplication of PMatrices.

    % Make sure both arguments are PMatrices
    obj1 = PMatrix(obj1);
    obj2 = PMatrix(obj2);
   
    % Get coefficients
    c1 = obj1.coefs; c2 = obj2.coefs;
    s1 = size(c1); s2 = size(c2);
    N = max(s1(3),s2(3));

    % Enter frequency domain
    ff1 = fft(c1,N,3); 
    ff2 = fft(c2,N,3);
    ff3 = zeros(s1(1),s2(2),N);

    % Perform component-wise multiplication (.* doesn't work)
    for i = 1:N
        ff3(:,:,i) = ff1(:,:,i)*ff2(:,:,i);
    end

    % Leave frequency domain
    c3 = ifft(ff3,N,3);

    % Return it
    r = PMatrix(c3, obj1.const_ind);

end % mtimes