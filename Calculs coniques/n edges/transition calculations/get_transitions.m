figure(1)
B2=findtransitions(2, [0,200;200,300], 0.01);
B21=findtransitionsfor2(2, [0,20; 200, 320], 0.01);
B2(2,:)=B21(2,:);
figure(2)
B3=findtransitions(3, [0,150;150,266], 0.01);
figure(3)
B4=findtransitions(4, [1,100;100,186], 0.01);
figure(4)
B5=findtransitions(5, [0,37;37,90], 0.01);
figure(5)
B52=findtransitions_2(5, [0,20; 20, 120], 0.01);

string={'the boundaries of strict self similarity for the solutions are : \n'};
index=[1:5, 52];
for ii=2:6
    eval(['B=B' num2str(index(ii)) ';']);
    string{ii}=[num2str(B([3,2]), '%.1f and %.1f for ') num2str(index(ii),'%g edges') '\n'];
end

s=[string{:}]
s=sprintf(s)
disp(s)
    
    
    