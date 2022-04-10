Error=zeros(5,6);
rx=100;ry=800;
Variance=[1,10,100,1000,10000,100000,1e6,1e7];
loop=20000;
for Mode=1:4
    
    disp(Mode)
    for v=1:length(Variance)
        Var=Variance(v);
        sum=0;
        for i=1:loop
            [X,Y]=Multilateration_Math(Mode,rx,ry,Var);
            sum=sum+(X-rx)^2+(Y-ry)^2;
        end
        Error(Mode,v)=sum/loop;
    end

end  
for v=1:length(Variance)
    Var=Variance(v);
    Error(5,v)=Multilateration_Math(5,rx,ry,Var);
end

figure;
hold on; box on;
plot(1:length(Variance),log(Error(1,:)),'k-','LineWidth',1);
plot(1:length(Variance),log(Error(2,:)),'g-','LineWidth',1);
plot(1:length(Variance),log(Error(3,:)),'r-','LineWidth',1);
plot(1:length(Variance),log(Error(4,:)),'b-','LineWidth',1);
plot(1:length(Variance),log(Error(5,:)),'LineWidth',1);


legend('LLS1','LLS2','WLLS','2SWLLS','CRLB')