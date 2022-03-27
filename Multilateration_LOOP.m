Mode=1;
Error=zeros(1,9);
index=1;
for r=5:5:45
    sum=0;
    for x=-r:r/10:r
        y=sqrt(r^2-x^2);
        [X,Y]=Multilateration_LLS(Mode,x,y);
        sum=abs(X-x)+abs(Y-y);
    end
    Error(index)=sum;
    index=index+1;
end
        
        