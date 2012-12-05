% This simple script compute the perimeter of the solution cells for
% different values of the number of edges. It is very useful to determine
% the maximum gamma angle possible for each solution.
%
%Copyright Antoine Lizée @ Dumais Laboratory, Harvard University
%July 2011

a=1;
n=5;
gamma0=[0,0,0,0,0,-50,-150];
gamma=-500:0.1:360;
perim=zeros(length(gamma),n);
perim(:,1)=sqrt(2*(360-gamma)*pi/180);
for ii=n
    k1=0;
    k2=0;
    
    j=find(gamma==gamma0(ii));    
    while k1==0 && j<=length(gamma)
        try
            [V, beta_final, perim(j,ii)]=findfixedpoint(ii+0.2,gamma(j),a);
            j=j+1;
        catch err
            if (strcmp(err.identifier,'FIND:endofsol'))
                disp(['for ' num2str(ii) ' edges, calculation has stopped at gamma=' num2str(gamma(j-1))]);
                k1=1;
                gamma_range_2(ii,2)=gamma(j-1);
            else
                rethrow(err);
            end
        end
    end
    
    j=find(gamma==gamma0(ii));    
    while k2==0 && j>=1
        try
            [V, beta_final, perim(j,ii)]=findfixedpoint(ii+0.2,gamma(j),a);
            j=j-1;
        catch err
            disp(['for ' num2str(ii) ' edges, calculation has stopped at gamma=' num2str(gamma(j+1))]);
            k2=1;
            gamma_range_2(ii,1)=gamma(j+1);
        end
    end
    
end

figure(3);
plot(gamma,perim,'Linewidth',1.5);
L=num2str(n','%d edges');
legend(L)
    