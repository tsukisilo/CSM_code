function a = iifx( b )
%IIFX Summary of this function goes here
%   Detailed explanation goes here
a = fftshift(ifft(fftshift(b)));

end

