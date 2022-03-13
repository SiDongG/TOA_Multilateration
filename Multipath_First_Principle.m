%% Multipath without channel and wave objects

%%
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
subplot (3,3,1)
plot(1:length(signal),signal)
cross=fftshift(ifft(fft(signal).*conj(fft(signal))));
subplot (3,3,2)
plot(1:length(cross),cross);

xcor_PHAT=1e3*abs(fftshift(ifft((fft(signal).*conj(fft(signal)))./(abs(fft(signal)).*abs(fft(signal))))));
hold on
subplot (3,3,3)
plot(1:length(xcor_PHAT),xcor_PHAT);

% multipath_gain(Non-dB) [1,3,0.5,2,0.3]
% multipath_delay(misec) [0,280us,330us,360us,470us]
% delay_samples [0,5600,6600,7200,9400]
path0=signal;
path1=3*[zeros(1,5600),signal];
path2=0.5*[zeros(1,6600),signal];
path3=2*[zeros(1,7200),signal];
path4=0.3*[zeros(1,9400),signal];
signal_r=path0+path1(1:120006)+path2(1:120006)+path3(1:120006)+path4(1:120006);

subplot(3,3,4)
plot(1:length(signal_r),signal_r);

cross2=fftshift(ifft(fft(signal_r).*conj(fft(signal))));
subplot (3,3,5)
plot(1:length(cross2),cross2);

xcor_PHAT2=1e3*abs(fftshift(ifft((fft(signal_r).*conj(fft(signal)))./(abs(fft(signal_r)).*abs(fft(signal))))));

subplot (3,3,6)
plot(1:length(xcor_PHAT2),xcor_PHAT2);

