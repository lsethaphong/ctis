function Fm = filtration(Im)
% filtration(Im, RES)
% Im = the Radon transform data
% ===================================================
RES = length(Im);
Fm = zeros(RES,RES);
mins = -1.0;
maxs =  1.0;
step = (maxs-mins)/RES;
theta = [2*pi/RES:2*pi/RES:2*pi];
pp  = sqrt(2.0)*[mins+step/2:step:maxs-step/2];
for w=1:RES
    Fm(w,:) = ifft(ifftshift(abs(pp).*(fftshift(fft(Im(w,:))))));
end
Fm = real(Fm);
    
