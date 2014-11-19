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
function test_suite = testAddSub
  initTestSuite;
end

function testAddCausal
    dims = randi(100,1,3);
    c1 = randn(dims) + j*randn(dims);
    c2 = randn(dims) + j*randn(dims);
    p1 = PMatrix(c1);
    p2 = PMatrix(c2);
    sum = p1 + p2;
    
    assertEqual(sum.coefs, c1 + c2);
end

function testSubCausal
    dims = randi(100,1,3);
    c1 = randn(dims) + j*randn(dims);
    c2 = randn(dims) + j*randn(dims);
    p1 = PMatrix(c1);
    p2 = PMatrix(c2);
    sum = p1 - p2;
    
    assertEqual(sum.coefs, c1 - c2);
end

function testAddNonCausal
    c1 = zeros(2,2,6);
    c1(2,1,1) = 1;
    c1(:,:,3) = eye(2);
    c1(:,:,4) = diag([1 5]);
    c1(:,2,6) = ones(2,1);
    p1 = PMatrix(c1, 4);
    
    c2 = zeros(2,2,8);
    c2(:,:,1) = [0 1;0 1];
    c2(2,2,2) = 1;
    c2(2,:,3) = [5 1];
    c2(2,:,4) = [1 1];
    c2(1,1,8) = 1;
    p2 = PMatrix(c2, 3);
    
    sum = p1 + p2;
    
    cs = zeros(2,2,9);
    cs(2,1,1) = 1;
    cs(:,2,2) = ones(2,1);
    cs(:,:,3) = diag([1 2]);
    cs(:,:,4) = [1 0;5 6];
    cs(2,:,5) = [1 1];
    cs(:,2,6) = ones(2,1);
    cs(1,1,9) = 1;
    
    assertEqual(sum.coefs, cs);
end

function testSubNonCausal
    c1 = zeros(2,2,6);
    c1(2,1,1) = 1;
    c1(:,:,3) = eye(2);
    c1(:,:,4) = diag([1 5]);
    c1(:,2,6) = ones(2,1);
    p1 = PMatrix(c1, 4);
    
    c2 = zeros(2,2,8);
    c2(:,:,1) = [0 1;0 1];
    c2(2,2,2) = 1;
    c2(2,:,3) = [5 1];
    c2(2,:,4) = [1 1];
    c2(1,1,8) = 1;
    p2 = PMatrix(c2, 3);
    
    sum = p1 - p2;
    
    cs = zeros(2,2,9);
    cs(2,1,1) = 1;
    cs(:,2,2) = [-1 -1]';
    cs(:,:,3) = diag([1 0]);
    cs(:,:,4) = [1 0;-5 4];
    cs(2,:,5) = -[1 1];
    cs(:,2,6) = [1 1]';
    cs(1,1,9) = -1;
    
    assertEqual(sum.coefs, cs);
end

function testAddScalar
    dims = randi(100,1,2);
    c1 = randn(dims) + j*randn(dims);
    c2 = randn(dims) + j*randn(dims);
    p1 = PMatrix(c1);
    p2 = PMatrix(c2);
    sum = p1 + p2;
    
    assertEqual(sum.coefs, c1 + c2);
end

function testSubScalar
    dims = randi(100,1,2);
    c1 = randn(dims) + j*randn(dims);
    c2 = randn(dims) + j*randn(dims);
    p1 = PMatrix(c1);
    p2 = PMatrix(c2);
    sum = p1 - p2;
    
    assertEqual(sum.coefs, c1 - c2);
end