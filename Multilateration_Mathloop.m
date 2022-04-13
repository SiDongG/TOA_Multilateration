%clear;clc;close all;
SNR_List=[1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1e0,1e1,1e2,1e3];
rx=300;
ry=600;
Error=zeros(6,length(SNR_List));
loop=10;
for Mode=1:6
    
    disp(Mode)
    for v=1:length(SNR_List)
        SNR=SNR_List(v);
        sum=0;
%         Measurement=zeros(1,loop);
        for i=1:loop
            [X,Y]=Multilateration_Math(Mode,rx,ry,SNR);
            sum=sum+(X-rx)^2+(Y-ry)^2;
%             Measurement(i)=X-rx;
        end
%         Error(Mode,v)=var(Measurement);
        Error(Mode,v)=sum/loop;
    end

end  
% for v=1:length(SNR_List)
%     SNR=SNR_List(v);
%     [X,Y]=Multilateration_Math(10,rx,ry,SNR);
%     Error(6,v)=X;
% end

figure;
hold on; box on;
semilogy(-60:10:30,log10(Error(1,:)),'k-','LineWidth',1);
semilogy(-60:10:30,log10(Error(2,:)),'g-','LineWidth',1);
semilogy(-60:10:30,log10(Error(3,:)),'r-','LineWidth',1);
semilogy(-60:10:30,log10(Error(4,:)),'b-','LineWidth',1);
semilogy(-60:10:30,log10(Error(5,:)),'LineWidth',1);
semilogy(-60:10:30,log10(Error(6,:)),'LineWidth',1);
% plot(-40:10:50,log10(Error(1,:)),'k-','LineWidth',1);
% plot(-40:10:50,log10(Error(2,:)),'g-','LineWidth',1);
% plot(-40:10:50,log10(Error(3,:)),'r-','LineWidth',1);
% plot(-40:10:50,log10(Error(4,:)),'b-','LineWidth',1);
% plot(-40:10:50,log10(Error(6,:)),'LineWidth',1);


legend('LLS1','LLS2','WLLS','2SWLLS','CWLLS','CRLB')