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
function test_suite = testMult
  initTestSuite;
end

function testMult1
    coefs = ones(2,2,2);
    p1 = PMatrix(coefs);
    p2_conv = p1.mtimes(p1,'conv');
    p2_fft = p1.mtimes(p1,'fft');
    
    coefs2 = 2*ones(2,2,3);
    coefs2(:,:,2) = 4*ones(2,2);
    
    assertElementsAlmostEqual(p2_conv.coefs, coefs2, 'absolute', 2*eps);
    assertElementsAlmostEqual(p2_fft.coefs, coefs2, 'absolute', 2*eps);
    assertElementsAlmostEqual(p2_conv.coefs, p2_fft.coefs, 'absolute', 2*eps);
end

function testMult2
    c1 = zeros(2,2,3);
    c1(:,:,1) = diag([1 2]);
    c1(1,2,2) = 1;
    c1(:,1,3) = ones(2,1);
    p1 = PMatrix(c1);
    
    c2 = zeros(2,2,4);
    c2(:,1,1) = ones(2,1);
    c2(:,:,2) = diag([5 2]);
    c2(1,2,4) = 1;
    p2 = PMatrix(c2, 2);
    
    p3_conv = p1.mtimes(p2, 'conv');
    p3_fft = p1.mtimes(p2, 'fft');
    
    c3 = zeros(2,2,6);
    c3(:,1,1) = [1 2]';
    c3(:,:,2) = diag([6 4]);
    c3(:,:,3) = [1 2;1 0];
    c3(:,:,4) = [5 1;5 0];
    c3(:,2,6) = ones(2,1);
    
    assertElementsAlmostEqual(p3_conv.coefs, c3, 'absolute', 2*eps);
    assertElementsAlmostEqual(p3_fft.coefs, c3, 'absolute', 2*eps);
    assertElementsAlmostEqual(p3_conv.coefs, p3_fft.coefs, 'absolute', 2*eps);
end

function testMultScalar
    dims = randi(100,1,2);
    c1 = randn(dims) + j*randn(dims);
    c2 = randn(fliplr(dims)) + j*randn(fliplr(dims));
    p1 = PMatrix(c1);
    p2 = PMatrix(c2);
    sum_conv = p1.mtimes(p2,'conv');
    sum_fft = p1.mtimes(p2,'fft');
    
    assertElementsAlmostEqual(sum_conv.coefs, c1*c2, 'absolute', 2*eps);
    assertElementsAlmostEqual(sum_fft.coefs, c1*c2, 'absolute', 2*eps);
    assertElementsAlmostEqual(sum_conv.coefs, sum_fft.coefs, 'absolute', 2*eps);
end

function testExample1
    % Matrices
    coefs = zeros(3,3,3);
    coefs(2,1,1) = 1;
    coefs(1,3,1) = 2;
    coefs(:,:,2) = diag([2 1 1]);
    coefs(3,2,3) = 1;
    A = PMatrix(coefs, 2);

    coefs = zeros(3,3,3);
    coefs(:,:,1) = [0 0 0;-0.2981 0 0.7454;0.3333 0 0.6667];
    coefs(:,:,2) = [0.8944 0 0;0 0.5963 0;0 -0.6667 0];
    coefs(1,2,3) = 0.4472;
    Q = PMatrix(coefs, 2);

    coefs = zeros(3,3,4);
    coefs([2 3],3,1) = [-0.5963;0.6667];
    coefs(:,3,2) = [1.7889;0.7454;0.6667];
    coefs(1,1,3) = 2.2361;
    coefs(2,2,3) = 1.3416;
    coefs(1,2,4) = 0.4472;
    R = PMatrix(coefs, 3);
    
    % Compute multiplications
    Qph = Q.ph;
    diff_conv = A - Qph.mtimes(R,'conv');
    diff_fft = A - Qph.mtimes(R,'fft');
    diff_diff = diff_conv - diff_fft;
    
    assertElementsAlmostEqual(diff_conv.coefs, zeros(size(diff_conv.coefs)), 'absolute', 1e-3);
    assertElementsAlmostEqual(diff_fft.coefs, zeros(size(diff_fft.coefs)), 'absolute', 1e-3);
    assertElementsAlmostEqual(diff_diff.coefs, zeros(size(diff_diff.coefs)), 'absolute', 2*eps);
end