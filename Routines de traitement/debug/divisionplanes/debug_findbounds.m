load('tissue_test')
n=12;
k=2;
l=5;
new=1; %choice between creating new debugging tissues or concentrating on existing ones to understand issues with the algorithm.
for h=1:n
    f=ceil(h/4);
    figure(f);
    subplot(2,2,h-(f-1)*4);
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
%     for i=0:0.05:1
%         [ s0 ~, a V A]=finda(i,M,k,l);
%         P=findP(V,A(3),0:0.01:1);
%         plot(P(:,1),P(:,2),'-g','Linewidth',1);
% %         text(V(1:2,1),V(1:2,2),{num2str(s1,'\\color{black}  % .1f'); num2str(s0,'\\color{black}  % .1f')});
% %         text(V(2,1),V(2,2),num2str(a,'\\color{black} \\bf % .2f  \\rm'),'VerticalAlignment','Top');
%     end
    s1=findbounds(M,k,l);        
    for i=linspace(s1(1),s1(2),20)
        [ s0 s1 a V A]=finda(i,M,k,l);
        P=findP(V,A(3),0:0.01:1);
        plot(P(:,1),P(:,2),'-b','Linewidth',2);
        text(V(1:2,1),V(1:2,2),{num2str(s1,'\\color{black}  % .1f'); num2str(s0,'\\color{black}  % .1f')});
%         text(V(2,1),V(2,2),num2str(a,'\\color{black} \\bf % .2f  \\rm'),'VerticalAlignment','Top');
    end
%     plot([M.C(2,1),M.C(1,1)],[M.C(2,2),M.C(1,2)],'-r*','MarkerSize',8)

end
maximize(1:ceil(n/4));
%% Cleaning lines
% delete(findobj('Type','Figure'));
%clear('-regexp','^tissue\w*');