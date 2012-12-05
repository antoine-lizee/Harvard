% This simple script compute the dividing length of the solution cells for
% different values of the number of edges. This dividing length is in fact
% the length of the first edge, which is the dividing wall to be, divided by the
% square of two ratio.
%
%Copyright Antoine Lizée @ Dumais Laboratory, Harvard University
%July 2011

a=1;
n=2:5;
gamma0=[0,0,0,0,0,-50,-150];
gamma=0:0.1:360;
dl=zeros(length(gamma),numel(n)+1);
dl(:,1)=sqrt(2*(360-gamma)*pi/180);
for ii=n
    k1=0;
    k2=0;
    
    j=find(gamma==gamma0(ii));    
    while k1==0 && j<=length(gamma)
        try
            [V, beta_final, ~]=findfixedpoint(ii,gamma(j),a);
            j=j+1;
            cell=getcellfromsol(V,beta_final);
            V=[cell.v(3,:)-cell.v(2,:);cell.v(ii+3,:)-cell.v(ii+2,:)];
            A=cell.e([2,ii+2],3);
            dl(j,ii)=sum(sqrt(sum(V.^2,2)).*A./sin(A));
        catch err
            if (strcmp(err.identifier,'FIND:endofsol'))
                disp(['for ' num2str(ii) ' edges, calculation has stopped at gamma=' num2str(gamma(j-1))]);
                k1=1;
                gamma_range(ii,2)=gamma(j-1);
            else
                rethrow(err);
            end
        end
    end
    
    j=find(gamma==gamma0(ii));    
    while k2==0 && j>=1
        try
            [V, beta_final, ~]=findfixedpoint(ii,gamma(j),a);
            j=j-1;
            cell=getcellfromsol(V,beta_final);
            V=[cell.v(3,:)-cell.v(2,:);cell.v(ii+3,:)-cell.v(ii+2,:)];
            A=cell.e([2,ii+2],3);
            dl(j,ii)=sum(sqrt(sum(V.^2,2)).*A./sin(A));
        catch err
            disp(['for ' num2str(ii) ' edges, calculation has stopped at gamma=' num2str(gamma(j+1))]);
            k2=1;
            gamma_range(ii,1)=gamma(j+1);
        end
    end
    
end
dl=dl/sqrt(2);
figure(2);
plot(gamma,dl,'Linewidth',1.5);
L=num2str([1,n]','%d edges');
legend(L)
    