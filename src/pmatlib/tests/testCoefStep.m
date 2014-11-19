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
function test_suite = testCoefStep
  initTestSuite;
end

function testCoefStep1
    % Generate matrix
    A = PMatrix(randn(3,3,3)); s = size(A);
    
    % Select coefficient to null
    row = 2; col = 1; lead = -1; t = -lead;
    
    % Calculate rotation parameters
    a_null = A.get_coefs(row,col,lead);
    a_rot = A.get_coefs(col,col,0);
    [ alpha, theta, phi ] = rot_params(a_null, a_rot);
    
    % Colstep method
    Acs = A.coefstep(row, col, t, alpha, theta, phi);
    
    % Straight-forward method
    B = etsm(s(1), row, -t);
    G = epgrm(s(1), row, col, t, alpha, theta, phi);

    % Rotate
    Asf = B*G*A;
    
    % Check consistency
    assertElementsAlmostEqual(Acs.coefs, Asf.coefs, 'absolute', 10*eps);
end