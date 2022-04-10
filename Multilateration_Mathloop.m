clear;clc;close all;
SNR_List=[1e-5,1e-4,1e-3,1e-2,1e-1,1e0,1e1,1e2];
Error=zeros(5,length(SNR_List));
loop=20000;
for Mode=1:4
    
    disp(Mode)
    for v=1:length(SNR_List)
        SNR=SNR_List(v);
        sum=0;
        Measurement=zeros(1,loop);
        for i=1:loop
            rx=randi([100,900]);
            ry=randi([100,900]);
            [X,Y]=Multilateration_Math(Mode,rx,ry,SNR);
            sum=sum+(X-rx)^2+(Y-ry)^2;
            Measurement(i)=X;
        end
        Error(Mode,v)=var(Measurement);
    end

end  
for v=1:length(SNR_List)
    SNR=SNR_List(v);
    [X,Y]=Multilateration_Math(5,rx,ry,SNR);
    Error(5,v)=X;
end

figure;
hold on; box on;
plot(1:length(SNR_List),log(Error(1,:)),'k-','LineWidth',1);
plot(1:length(SNR_List),log(Error(2,:)),'g-','LineWidth',1);
plot(1:length(SNR_List),log(Error(3,:)),'r-','LineWidth',1);
plot(1:length(SNR_List),log(Error(4,:)),'b-','LineWidth',1);
plot(1:length(SNR_List),log(Error(5,:)),'LineWidth',1);


legend('LLS1','LLS2','WLLS','2SWLLS','CRLB')