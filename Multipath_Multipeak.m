%% Multipath_multi-peak detection for delay testing

clear;clc;close all;
fb=1e3;
bit_period=1/fb; %bit period is 0.001 second 
fc=1e4; %Carrier Frequency
fs=2e4; %Sampling Rate
bits=[1,0,1,0,0,1];
t=0:1/(fs*fb):1/fb;
for i=1:length(bits)
    analog_sig=(bits(i)*2-1)*cos(2*pi*fc*t);
    if i==1
        analog=analog_sig;
    else
        analog=[analog,analog_sig];
    end
end
signal=analog;
path0=signal;
path1=3*[zeros(1,5600),signal];
path2=0.5*[zeros(1,6600),signal];
path3=2*[zeros(1,7200),signal];
path4=0.3*[zeros(1,9400),signal];
signal_r=path0+path1(1:120006)+path2(1:120006)+path3(1:120006)+path4(1:120006);

cross=fftshift(ifft(fft(signal).*conj(fft(signal))));
peak_index=find(cross==max(cross));
cross_r=fftshift(ifft(fft(signal_r).*conj(fft(signal))));
peak_index_r=find(cross_r==max(cross_r));
cross_r_PHAT=1e3*abs(fftshift(ifft((fft(signal_r).*conj(fft(signal)))./(abs(fft(signal_r)).*abs(fft(signal))))));
peak_index_PHAT=find(cross_r_PHAT==max(cross_r_PHAT));

signal_r=path0+path2(1:120006)+path3(1:120006)+path4(1:120006);
cross_r=fftshift(ifft(fft(signal_r).*conj(fft(signal))));
peak_index_r2=find(cross_r==max(cross_r));
cross_r_PHAT=1e3*abs(fftshift(ifft((fft(signal_r).*conj(fft(signal)))./(abs(fft(signal_r)).*abs(fft(signal))))));
peak_index_PHAT2=find(cross_r_PHAT==max(cross_r_PHAT));



