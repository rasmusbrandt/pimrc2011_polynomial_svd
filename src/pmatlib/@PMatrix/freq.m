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
function out = freq(obj, Nfft, c)
    % FREQ Returns frequency domain representation of filter
    % If c is present, the spectrum is divided by c's spectrum
    
    % Get coefficients
    coefs = obj.coefs; ss = size(obj);
    
    % Check arguments
    if nargin == 1
        Nfft = ss(3);
    end
    
    % Zero-pad for alignment
    nb = obj.const_ind - 1; na = ss(3) - nb - 1;
    ca = zeros(ss(1), ss(2), nb-na);
    cb = zeros(ss(1), ss(2), na-nb);
    coefs = cat(3, cb, coefs);
    coefs = cat(3, coefs, ca);
    
    % Zero-pad for FFT
    diff = Nfft - length(coefs);
    if diff > 0
        ca = zeros(ss(1), ss(2), floor(diff/2));
        cb = zeros(ss(1), ss(2), ceil(diff/2));
        coefs = cat(3, cb, coefs);
        coefs = cat(3, coefs, ca);
    end
    
    % Align for FFT
    coefs = ifftshift(coefs, 3);
    
    % Do it
    out = fft(coefs,Nfft,3);
    
    % Divide by c spectrum?
    if nargin == 3
        cf = c.freq(Nfft);
        
        for i = 1:ss(1)
            for j = 1:ss(2)
                out(i,j,:) = out(i,j,:)./cf;
            end
        end
    end
end