
x0=100;
t=1;
for i=1:x0
    t1=[0,t(end,:),0];
    t2=t1+circshift(t1,[0,1])+circshift(t1,[0,-1]);
    t=[zeros(i,1),t,zeros(i,1);t2];
end
T=diag(1./sum(t,2))*t;


x=1:1:2*x0+1;
m=x0+1;

X2=T.*(((ones(x0+1,1)*x)-m).^2);
V=sqrt(sum(X2,2));
%V2=sqrt(var(x,T(end,:))); %Just a verification of the V computation


figure(5); clf
plot(x,T);

%display(V(1:end))

w=V(end);

c1=sqrt(2*pi)*w;
c2=2*w.^2;
z=1/c1*exp(-(x-m).^2/c2);
hold all;
plot(x,z,'Linewidth',2)

X=[2:10,15,22,30,39,61,95];
disp(sprintf('%4d',X'))
disp(sprintf('%.1f|',V(X+1)'*pi))