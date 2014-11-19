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
classdef PMatrix
    % PMatrix
    % A class which represents a polynomial matrix.
    %
    % Rasmus Brandt, <rabr5411@kth.se>
    % 2010-05-03
    
    properties
        coefs           % Numerator coefficient matrices
        const_ind       % Index of constant term
    end
    
    methods
        function obj = PMatrix(arg1, arg2)
            % PMATRIX Constructor
            
            if nargin == 0
                % The empty matrix
                obj.coefs = 0;
                obj.const_ind = 1;
            elseif nargin == 1
                if(isa(arg1, 'PMatrix'))
                    % Copy other Pmatrix
                    obj.coefs = arg1.coefs;
                    obj.const_ind = arg1.const_ind;
                else
                    % Assume coefficient matrices given
                    obj.coefs = arg1;
                    obj.const_ind = 1; 
                end
            elseif nargin == 2
                % Assume the given arguments are ok
                obj.coefs = arg1;
                obj.const_ind = arg2;
            else
                ME = MException('PMatrix:bad_constructor_call', ...
                    'Incorrect number of constructor arguments!');
                throw ME;
            end
        end
    end
    
    methods (Static)
        s = polynomial_as_cells(coef, const_ind, approx);
    end
end