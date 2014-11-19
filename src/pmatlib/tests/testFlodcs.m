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
function test_suite = testFlodcs
  initTestSuite;
end

function testFlodcLower
    for i = 1:5
    for j = 1:5
    for k = 1:5
            
    % Generate matrix
    dims = [i j k];
    c1 = randn(dims) + j*randn(dims);
    p1 = PMatrix(c1);
    
    % Iterate over all elements
    for row = 2:dims(1)
        for col = 1:min(row-1, dims(2))
            for ind = 1:dims(3)
                lead = p1.ind2lead(ind);

                % Increase the selected coefficient
                old_coef = p1.coefs(row,col,ind);
                new_coef = 1e3 + sum(sum(sum(c1))); mag = abs(new_coef);
                p1.coefs(row,col,ind) = new_coef;

                % Find largest coefficient under main diagonal, 
                % in the correct column
                [ row2, lead2, mag2 ] = p1.flodc_lower(col);
                assertEqual(row2, row);
                assertEqual(lead2, lead);
                assertEqual(mag2, mag);

                % Find largest coefficient under main diagonal
                [ row2, col2, lead2, mag2 ] = p1.flodc_lower_all();
                assertEqual(row2, row);
                assertEqual(col2, col);
                assertEqual(lead2, lead);
                assertEqual(mag2, mag);
                
                % Return coefficient to old value
                p1.coefs(row,col,ind) = old_coef;
            end
        end
    end
    
    end
    end
    end
end

function testFlodcUpper
    for i = 1:5
    for j = 1:5
    for k = 1:5
            
    % Generate matrix
    dims = [i j k];
    c1 = randn(dims) + j*randn(dims);
    p1 = PMatrix(c1);
    
    % Iterate over all elements
    for row = 1:min(dims(1),dims(2)-1)
        for col = row+1:dims(2)
            for ind = 1:dims(3)
                lead = p1.ind2lead(ind);

                % Increase the selected coefficient
                old_coef = p1.coefs(row,col,ind);
                new_coef = 1e3 + sum(sum(sum(c1))); mag = abs(new_coef);
                p1.coefs(row,col,ind) = new_coef;

                % Find largest coefficient over main diagonal,
                % in the correct column
                [ row2, lead2, mag2 ] = p1.flodc_upper(col);
                assertEqual(row2, row);
                assertEqual(lead2, lead);
                assertEqual(mag2, mag);

                % Find largest coefficient over main diagonal
                [ row2, col2, lead2, mag2 ] = p1.flodc_upper_all();
                assertEqual(row2, row);
                assertEqual(col2, col);
                assertEqual(lead2, lead);
                assertEqual(mag2, mag);

                % Return coefficient to old value
                p1.coefs(row,col,ind) = old_coef;
            end
        end
    end
    
    end
    end
    end
end

function testFlodcAll
    for i = 1:5
    for j = 1:5
    for k = 1:5
            
    % Generate matrix
    dims = [i j k];
    c1 = randn(dims) + j*randn(dims);
    p1 = PMatrix(c1);
    
    % Iterate over all elements
    for row = 1:dims(1)
        for col = 1:dims(2)
            % Don't look at main diagonal
            if row == col
                continue;
            end
            
            for ind = 1:dims(3)
                lead = p1.ind2lead(ind);
    
                % Increase the selected coefficient
                old_coef = p1.coefs(row,col,ind);
                new_coef = 1e3 + sum(sum(sum(c1))); mag = abs(new_coef);
                p1.coefs(row,col,ind) = new_coef;

                % Find largest coefficient under main diagonal
                [ row2, col2, lead2, mag2 ] = p1.flodc_all();
                assertEqual(row2, row);
                assertEqual(col2, col);
                assertEqual(lead2, lead);
                assertEqual(mag2, mag);
                
                % Return coefficient to old value
                p1.coefs(row,col,ind) = old_coef;
            end
        end
    end
    
    end
    end
    end
end