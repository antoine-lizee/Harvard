
x=0:0.01:5;
T=2;
y=cos(2*pi*x/T)+y;
Fy=abs(fft(y));
figure(3);
f=x*length(x)/range(x).^2;
plot(f(1:floor(end/2)),Fy(1:floor(end/2)));
