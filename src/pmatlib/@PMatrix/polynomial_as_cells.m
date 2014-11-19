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
function s = polynomial_as_cells(coef, const_ind, approx)
   % POLYNOMIAL_AS_STRING Returns a string formatted version of a
   % polynomial
   
   % Make sure the coefs are in a row vector
   coef = coef(:).';
   
   % Preallocate some space (otherwise MATLAB might freak out)
   s = cell(1,8);
   
   % Iterate over coefficients, and display them accordingly
   if all(coef == 0)
      s = {'0'};
   else
      d = const_ind - 1; % First exponent
      ind = 1; % String component index
      
      for k = 1:length(coef)
         % For each coefficient
         a = coef(k);
         
         % Remove small terms, for the sake of clarity
         if approx
             if abs(a) < 1e-9
                 a = 0;
             end
             if abs(real(a)) < 1e-9
                 a = imag(a);
             end
             if abs(imag(a)) < 1e-9
                 a = real(a);
             end
         end
         
         if a ~= 0
            % Print sign
            if real(a) >= 0
               if ind ~= 1
                 % Except for first coefficient, if positive
                 s(ind) = {' + '};
                 ind = ind + 1;
               end
            else
               if ind ~= 1
                   s(ind) = {' - '};
                   a = -a;
                   ind = ind + 1;
               else
                   % No leading space
                   s(ind) = {'-'};
                   a = -a;
                   ind = ind + 1;
               end
            end
            
            % Print coefficient, unless it is a 1 which is not the constant
            % coefficient
            if a ~= 1 || d == 0
              if(isreal(a))
                s(ind) = {num2str(a)};
              else
                s(ind) = {['(' num2str(a) ')']}; 
              end
              ind = ind + 1;

              % Print * for non-constant terms
              if d ~= 0
                s(ind) = {'*'};
                ind = ind + 1;
              end
            end
            
            % Print the indeterminate variable z, but not for the constant
            % term
            if d < 0
               s(ind) = {['z^(' int2str(d) ')']};
               ind = ind + 1; 
            elseif d > 1
               s(ind) = {['z^' int2str(d)]};
               ind = ind + 1;
            elseif d == 1
               s(ind) = {'z'};
               ind = ind + 1;
            end
         end
         
         % Go to next exponent
         d = d - 1;
      end
   end
   
   % If nothing, display a zero
   if isempty(s)
       s = {'0'};
   end
end