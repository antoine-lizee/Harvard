%%
tissue28_2=tissue28;
tissue28_2.e=circshift(tissue28_2.e,[-2,0]);
tissue28_2.v(1,:)=V(1,:);
tissue28_2.e(3,3)=-A(2);
tissue28_2.e(4,2:3)=[3,A(3)];
tissue28_2.c{1}=[1,2,-3,4];
tissue28_2=clean(tissue28_2);
figure(103);
clf
plot(tissue28_2);
axis equal

%%
%M=findcarac(tissue36,1);
fnplt(rsmak('circle',M.R_e(2),M.C(2,:)));
%plot(M.C(1,1),M.C(1,2),'or');

%%
%P=findP([V(2,:);M.v1(2,:)],A(1),0:0.01:1);
P=findP([M.v1(2,:);V(1,:)],A(2),0:0.01:1);
p=plot(P(:,1),P(:,2),'-k','Linewidth',1*2);

%%Test
d2=[-10 -7 -5 -3 0 3 5 7 10 15];
disp(d2);
d2=d2+(d2<-2*pi)*2*pi+(d2<0)*2*pi-(d2>2*pi)*2*pi-(d2>4*pi)*2*pi;
disp(d2);


%%
x=0:0.01:pi;
figure(100);
clf
hold all
plot(x,x-sin(x));
plot(x,sin(x/2).*sin(x).^2/2);

%%
gamma=0:0.02:2*pi;
figure(100);
cla
a=(pi-gamma/2)/pi;
y=asin(a);
% y2=atan(sqrt(a.^2./(1-a.^2)));

plot(gamma*180/pi,y*180/pi);
% hold all
% plot(gamma*180/pi,y2*180/pi);