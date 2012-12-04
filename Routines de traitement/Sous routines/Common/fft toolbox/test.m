[nu t] = ftAxis(256,200);

w=0.05
c1=sqrt(2*pi)*w;
c2=2*w^2;
y = 1/c1*exp(-t.^2/c2);
figure(1);
plot(t,y);
xlabel('Temps (s)');
figure(2);
Fy=abs(ft(y));
plot(nu,Fy);
xlabel('Fréquence (Hz)');