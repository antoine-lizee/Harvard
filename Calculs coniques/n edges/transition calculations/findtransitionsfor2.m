function [ bounds ] = findtransitionsfor2( n, bounds0 , p )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

a=1;

bounds=bounds0;
[ Vt At ~ ] = findfixedpoint( n, bounds0(2,1), a );
cell0=getcellfromsol(Vt,At);
[~,~,~, index]=findwall(cell0,n,bounds0(2,1));
walls0=index;

k=0;
while k==0
    bnext=(bounds(2,2)+bounds(2,1))/2;
    [ Vt At ~ ] = findfixedpoint( n, bnext, a );
    cell=getcellfromsol(Vt,At);
    [V A, ~, index]=findwall(cell,n,bnext);
    if index==walls0,
        bounds(2,1)=bnext;
    else
        bounds(2,2)=bnext;
    end
    
    if (bounds(2,2)-bounds(2,1))<p
        k=1;
    end
    disp(bounds(2,:));
    disp(index);
end

%% Plot the solutions

gcf;
bounds_bis=bounds';

for ii=3:4
    subplot(1,4,ii); 
    cla;
    [ Vt At ~ ] = findfixedpoint( n, bounds_bis(ii), a );
    cell=getcellfromsol(Vt,At);
    [V A, ~, ~]=findwall(cell,n,bounds_bis(ii));   
    plot(cell,'-g*','Linewidth',2);
    axis equal
    P=findP(V,A(3),0:0.01:1);
    plot(P(:,1),P(:,2),'-r','Linewidth',2);

    title(num2str(bounds_bis(ii),'gamma = %.1f'));
end

end




function [V A length index]=findwall(cell,n, gamma)

V=[];
A=[];
length=100;

M=findcarac(cell,1);

for i=2:n+1
    for j=i:n+2
        if i==j, continue, end
        [ ~,~, Vi, Ai, ~, lengthi] = findplane( M, i, j );
        if lengthi<length
            index=[i j];
            V=Vi;
            A=Ai;
            length=lengthi;
        end
    end
end

gamma2=2*pi-gamma/180*pi;
[Vix Viy]=pol2cart([-gamma2/2; gamma2/2],ones(2,1)*1/sqrt(gamma2));
Vi=[Vix Viy];
Ai=[1,1,gamma2/2];
lengthi=sqrt(gamma2);
disp([length,lengthi]);
if lengthi<length
    index=[0,0];
    V=Vi;
    A=Ai;
    length=lengthi;
end

end

