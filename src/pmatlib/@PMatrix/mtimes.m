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
function r = mtimes(obj1, obj2, algorithm)
    % MTIMES Multiplication of PMatrices.

    % Make sure both arguments are PMatrices
    obj1 = PMatrix(obj1);
    obj2 = PMatrix(obj2);
    
    % Default arguments
    if nargin == 2
        if exist('convmat', 'file')
            algorithm = 'conv';
        else
            algorithm = 'fft';
        end
        
        if (obj1.is1x1() || obj2.is1x1())
            algorithm = 'fft';
        end
        
        if (obj1.isconstant()) || (obj2.isconstant())
            algorithm = 'loop';
        end
    end
    
    % Select appropriate algorithm
    if strcmp(algorithm, 'conv')
        % MEX implementation (convolution)
        
        [ new_coefs, new_const_ind ] = ...
           convmat(obj1.coefs, obj1.const_ind, obj2.coefs, obj2.const_ind);
        
        r = PMatrix(new_coefs, new_const_ind);
        
    elseif strcmp(algorithm, 'fft')
        % FFT implementation
        
        % Get coefficients
        c1 = obj1.coefs; c2 = obj2.coefs;
        s1 = obj1.size; s2 = obj2.size;
        
        % Degree of FFT
        d = obj1.max_degree + obj2.max_degree + 1;
        
        % Zero-pad to the right size
        cc1 = cat(3,c1,zeros(s1(1),s1(2),d-s1(3)));
        cc2 = cat(3,c2,zeros(s2(1),s2(2),d-s2(3)));
        
        % Enter frequency domain
        ff1 = fft(cc1,[],3); 
        ff2 = fft(cc2,[],3);
        ff3 = zeros(s1(1),s2(2),d);
        
        % Perform component-wise multiplication (.* doesn't work)
        for i = 1:d
            ff3(:,:,i) = ff1(:,:,i)*ff2(:,:,i);
        end
        
        % Leave frequency domain
        cc3 = ifft(ff3,[],3);
        
        % Calculate new constant index
        new_ind_one_lead = obj1.ind2lead(1) + obj2.ind2lead(1);
        new_const_ind = new_ind_one_lead + 1;
        
        % Return it
        r = PMatrix(cc3, new_const_ind);
        
    elseif strcmp(algorithm, 'loop')
        % Multiply using Matlab loops. SLOW.

        s1 = obj1.size; s2 = obj2.size;

        % Allocate new PMatrix
        new_ind_one_lead = obj1.ind2lead(1) + obj2.ind2lead(1);
        new_const_ind = new_ind_one_lead + 1;

        N_leads = s1(3) + s2(3) - 1;

        if obj1.is1x1()
            new_coefs = zeros(s2(1), s2(2), N_leads);
        elseif obj2.is1x1()
            new_coefs = zeros(s1(1), s1(2), N_leads);
        else
            new_coefs = zeros(s1(1), s2(2), N_leads);
        end

        r = PMatrix(new_coefs, new_const_ind);

        % Multiply!
        for i = 1:s1(3)
            for j = 1:s2(3)
                % Method calls take a long time apparently!

                % ORIGINAL
                % lead = obj1.ind2lead(i) + obj2.ind2lead(j);
                % ind = r.lead2ind(lead);

                % OPTIMIZED
                ind = r.const_ind - ...
                      (obj1.const_ind - i + obj2.const_ind - j);

                r.coefs(:,:,ind) = r.coefs(:,:,ind) + ...
                    obj1.coefs(:,:,i) * obj2.coefs(:,:,j);
            end
        end
    else
        error('No algorithm selected.');
    end
end % mtimes