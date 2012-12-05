load('tissue_test')
n=20;
for h=[20]
    %f=ceil(h/4);
    figure(h);
    %subplot(2,2,h-(f-1)*4);
    clf;
    %You can comment these four following lines and uncomment the two next
    %in order to focus on one tissue that seems to have pb (see the debug
    %part of 'finda')
%     tissue=originaltissue3;
%     tissue.e([1,2,5],3)=rand(3,1)*2-1;
%     tissue.v([1,2],1)=rand(2,1)+[-0.5; 0.5];
%     eval(['tissue' num2str(h) '=tissue']);
    for j=1
        name=['tissue' num2str(h) ];
        eval(['tissue=' name ';']);
        M=findcarac(tissue,1); 
        plot(tissue,'-r*','Linewidth',2)
        axis equal
        s1=zeros(1000,1);
        s2=s1;
        szero=s;
        % disp([M.I_e M.I_v M.v1 M.v2 M.a M.xy M.l_e M.R_e M.r_e M.T M.N M.CV M.C])
        for i=0:0.01:1
            index=int8(i*100)+1;
            szero(index)=i;
%             [ s0 s1(index) a V A]=finda(i,M,1,2);
%             P=findP(V,A(3),0:0.01:1);
%             plot(P(:,1),P(:,2),'-k','Linewidth',1);
% %             text(V(1:2,1),V(1:2,2),{num2str(s1,'\\color{black}  % .1f'); num2str(s0,'\\color{black}  % .1f')});
%             text(V(2,1),V(2,2),num2str(a,'\\color{black} \\bf % .2f  \\rm'),'VerticalAlignment','Top');
            [ s0 s2(index) a V A]=finda(i,M,2,5);
            P=findP(V,A(3),0:0.01:1);
            plot(P(:,1),P(:,2),'-g','Linewidth',1);
%             text(V(1:2,1),V(1:2,2),{num2str(s1,'\\color{green}  % .1f'); num2str(s0,'\\color{green}  % .1f')});
            text(V(2,1),V(2,2),num2str(a,'\\color{green} \\bf % .2f  \\rm'),'VerticalAlignment','Top');
        end
        figure(100+h);
        plot(szero,[s1,s2]);
%         plot([M.C(2,1),M.C(1,1)],[M.C(2,2),M.C(1,2)],'-r*','MarkerSize',8)
    end
end
maximize(1:n);
maximize(101:100+n);
%clear('-regexp','^tissue\w*');