function [ bounds] = findtransitions_2( n, bounds0 , p)
%FINDTRANSITIONS Summary of this function goes here
%   This function compute the values of gamma that define the boundaries of
%   the existence space of a solution.
%   gamma0 is a starting point, p is the precision.

a=1;

bounds=bounds0;
[ Vt At ~ ] = findfixedpoint( n+0.2, bounds0(2,1), a );
cell0=getcellfromsol(Vt,At);
[~,~,~, index]=findwall(cell0,n);
walls0=index;

k=0;
while k==0
    bnext=(bounds(1,2)+bounds(1,1))/2;
    [ Vt At ~ ] = findfixedpoint( n+0.2, bnext, a );
    cell=getcellfromsol(Vt,At);
    [V A, ~, index]=findwall(cell,n);
    if index==walls0,
        bounds(1,2)=bnext;
    else
        bounds(1,1)=bnext;
    end
    
    if (bounds(1,2)-bounds(1,1))<p
        k=1;
	end
    disp(bounds(1,:));
    disp(index);
end

k=0;
while k==0
    bnext=(bounds(2,2)+bounds(2,1))/2;
    [ Vt At ~ ] = findfixedpoint( n+0.2, bnext, a );
    cell=getcellfromsol(Vt,At);
    [V A, ~, index]=findwall(cell,n);
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
clf
bounds_bis=bounds';

for ii=1:4
    subplot(1,4,ii); 
    [ Vt At ~ ] = findfixedpoint( n+0.2, bounds_bis(ii), a );
    cell=getcellfromsol(Vt,At);
    [V A, ~, ~]=findwall(cell,n);   
    plot(cell,'-g*','Linewidth',2);
    axis equal
    P=findP(V,A(3),0:0.01:1);
    plot(P(:,1),P(:,2),'-r','Linewidth',2);

    title(num2str(bounds_bis(ii),'gamma = %.1f'));
end

end

function [V A length index]=findwall(cell,n)

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

end

