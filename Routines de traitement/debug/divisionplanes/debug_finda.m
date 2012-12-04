for h=21:24
    figure(h);
    clf
    for j=1:3:10
        name=['originaltissue' num2str(h) '_' num2str(j) ];
        eval(['tissue=' name ';']);
        M=findcarac(tissue,1); 
        subplot(2,2,ceil(j/3));  
        %figure(j)
        %clf
        plot(tissue)
        axis equal

        % disp([M.I_e M.I_v M.v1 M.v2 M.a M.xy M.l_e M.R_e M.r_e M.T M.N M.CV M.C])
        for i=0:0.2:1
            [ s0 s1 a V A]=finda(i,M,2,4);
            P=findP(V,A(3),0:0.01:1);
            plot(P(:,1),P(:,2),'-k');
            text(V(1:2,1),V(1:2,2),{num2str(s1,'\\color{red}  %.1f'); num2str(s0,'\\color{red}  %.1f')});
            [ s0 s1 a V A]=finda(i,M,3,4);
            P=findP(V,A(3),0:0.01:1);
            plot(P(:,1),P(:,2),'-k');
            text(V(1:2,1),V(1:2,2),{num2str(s1,'  %.1f'); num2str(s0)});
            
        end
        %plot([M.C(2,1),M.C(4,1)],[M.C(2,2),M.C(4,2)],'-r*','MarkerSize',10)
    end
end