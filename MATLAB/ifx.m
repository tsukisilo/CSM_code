function a = ifx( b )
%IFX Summary of this function goes here
%   Detailed explanation goes here
c = fftshift(fft(fftshift(b)));
a=c;
end

