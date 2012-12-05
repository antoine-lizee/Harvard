load('tissue_test')
n=40;
for h=1:n
    f=ceil(h/4);
    figure(f);
    subplot(2,2,h-(f-1)*4);
    %clf;
    %You can comment these four following lines and uncomment the two next
    %in order to focus on one tissue that seems to have pb (see the debug
    %part of 'finda')
    tissue=originaltissue3;
    tissue.e([1,2,5],3)=rand(3,1)*2-1;
    tissue.v([1,2],1)=rand(2,1)+[-0.5; 0.5];
    eval(['tissue' num2str(h) '=tissue']);
    for j=1
%         name=['tissue' num2str(h) ];
%         eval(['tissue=' name ';']);
        M=findcarac(tissue,1); 
        plot(tissue,'-r*','Linewidth',2)
        axis equal
        s=zeros(1,2);
        s2=s;
        s1=s;
        [ s0 s2(1) ]=finds2(0,M,2,5);
        [ s0 s1(1) ]=finds2(1,M,5,2);
        [ s0 s2(2) ]=finds2(1,M,2,5);
        [ s0 s1(2) ]=finds2(0,M,5,2);
        if  s2(2)>s2(1), s1=circshift(s1, [0,1]), end
        s1=[0,1]+(s1-[0,1]).*(s2>1 | s2<0);        
        for i=linspace(s1(1),s1(2),20)
%             [ s0 s1 a V A]=finda(i,M,1,2);
%             P=findP(V,A(3),0:0.01:1);
%             plot(P(:,1),P(:,2),'-k','Linewidth',1);
%             text(V(1:2,1),V(1:2,2),{num2str(s1,'\\color{black}  % .1f'); num2str(s0,'\\color{black}  % .1f')});
%             text(V(2,1),V(2,2),num2str(a,'\\color{black} \\bf % .2f  \\rm'),'VerticalAlignment','Top');
            [ s0 s1 a V A]=finda(i,M,2,5);
            P=findP(V,A(3),0:0.01:1);
            plot(P(:,1),P(:,2),'-g','Linewidth',1);
%             text(V(1:2,1),V(1:2,2),{num2str(s1,'\\color{green}  % .1f'); num2str(s0,'\\color{green}  % .1f')});
            text(V(2,1),V(2,2),num2str(a,'\\color{green} \\bf % .2f  \\rm'),'VerticalAlignment','Top');
        end
%         plot([M.C(2,1),M.C(1,1)],[M.C(2,2),M.C(1,2)],'-r*','MarkerSize',8)
    end
end
maximize(1:ceil(n/4));
%clear('-regexp','^tissue\w*');