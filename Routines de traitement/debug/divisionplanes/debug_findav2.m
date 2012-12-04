load('tissue_test')
n=12;
new=1; %choice between creating new debugging tissues or concentrating on existing ones to understand issues with the algorithm.
sub=0;
for h=1:n
    if sub
        f=ceil(h/4);figure(f);subplot(2,2,h-(f-1)*4);
    else
        figure(h);clf;
    end
    if new
        tissue=originaltissue3;
        tissue.e([1,2,5],3)=rand(3,1)*2-1;
        tissue.v([1,2],1)=rand(2,1)+[-0.5; 0.5];
        eval(['tissue' num2str(h) '=tissue']);
    end
    name=['tissue' num2str(h) ];
    eval(['tissue=' name ';']);
    M=findcarac(tissue,1); 
    plot(tissue,'-r*','Linewidth',2)
    axis equal

    % disp([M.I_e M.I_v M.v1 M.v2 M.a M.xy M.l_e M.R_e M.r_e M.T M.N M.CV M.C])
    for i=0:0.1:1
        [ s0 s1 a V A]=finda(i,M,1,2);
        P=findP(V,A(3),0:0.01:1);
        plot(P(:,1),P(:,2),'-k','Linewidth',1);
%             text(V(1:2,1),V(1:2,2),{num2str(s1,'\\color{red}  % .1f'); num2str(s0,'\\color{red}  % .1f')});
        text(V(2,1),V(2,2),num2str(a,'\\color{black}  % .2f'));
        [ s0 s1 a V A]=finda(i,M,2,5);
        P=findP(V,A(3),0:0.01:1);
        plot(P(:,1),P(:,2),'-g','Linewidth',1);
%             text(V(1:2,1),V(1:2,2),{num2str(s1,'  % .1f'); num2str(s0,'  % .1f')});
        text(V(2,1),V(2,2),num2str(a,'\\color{green}  % .2f'));
    end
%         plot([M.C(2,1),M.C(1,1)],[M.C(2,2),M.C(1,2)],'-r*','MarkerSize',8)
end
maximize(1:ceil(n/4));