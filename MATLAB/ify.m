function a = ify( b )
%IFY Summary of this function goes here
%   Detailed explanation goes here
a = fftshift(fft(fftshift(b.'))).';

end

