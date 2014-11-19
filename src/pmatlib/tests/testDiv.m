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
function test_suite = testDiv
  initTestSuite;
end

function testDiv1
    len1 = ceil(10*rand);
    len2 = ceil(10*rand);
    pos1 = ceil(10*rand);
    pos2 = ceil(10*rand);
    
    p1 = PMatrix(randn(1,1,len1),pos1);
    p2 = PMatrix(randn(1,1,len2),pos2);
    
    p1 = p1.rm_edge_zeros();
    p2 = p2.rm_edge_zeros();
    
    p1p2 = p1*p2;
    
    [q1,r1] = p1p2/p1;
    assertElementsAlmostEqual(squeeze(q1.coefs)', squeeze(p2.coefs)', 'absolute', 1e-8);
    assertElementsAlmostEqual(squeeze(r1.coefs)', squeeze(zeros(size(r1)))', 'absolute', 1e-8);
    
    [q2,r2] = p1p2/p2;
    assertElementsAlmostEqual(squeeze(q2.coefs)', squeeze(p1.coefs)', 'absolute', 1e-8);
    assertElementsAlmostEqual(squeeze(r2.coefs)', squeeze(zeros(size(r2)))', 'absolute', 1e-8);
end