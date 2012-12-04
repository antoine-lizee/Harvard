%This script enables one to test and see that the convolution with a gaussian 
%of std Dt selects the frequencies wich are smaller than 1/(2*pi*Dt), and
%therefore periods wich are greater than 2*pi*Dt.


n=1000
m=n+1;
x=-n:n;
X=x/n*10

w=0.02; %std of the gaussian
c1=sqrt(2*pi)*w;
c2=2*w.^2;
z=zeros(1,length(x));
z(m+floor(-4*w*100):m+floor(4*w*100))=1/c1*exp(-(X(m+floor(-4*w*100):m+floor(4*w*100))).^2/c2);

figure(6); clf
plot(X,z,'Linewidth',2)

p=20
s=cos((1:p)'*2*pi*X); s(p+1,:)=sin(2*pi*X);
f=sum(s);
figure(7); clf;
plot(X,f)
hold all
plot(X,z*w*p)
plot([-10, -w*pi-eps, -w*pi, w*pi, w*pi+eps, 10],[0,0,1,1,0,0])

C=zeros(2*n+1);
C(:,m)=z;
for i=1:n
    C(:,i)=(circshift(z,[0,-m+i]))';
    C(:,end+1-i)=circshift(z,[0,m-i])';
end

F=f*C*1/100;
plot(X,F,'Linewidth',2)

Xf=X*length(X)/range(X).^2;
FF=abs(ft(F));
Ff=abs(ft(f));
figure(8); clf;
plot(Xf,Ff); hold all
plot(Xf,FF,'Linewidth',2);

Fz=abs(ft(z))
FP=Ff.*Fz;
%figure(9);
%plot(Xf,[FP;Fz;Ff]);