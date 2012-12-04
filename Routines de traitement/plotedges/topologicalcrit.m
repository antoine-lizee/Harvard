% load
% ne_M1=nedges
% nedges=cellfun(@(X) length(X),originaltissue.c)
% load
% ne_M2=nedges
% nedges=cellfun(@(X) length(X),originaltissue.c)
% load
% ne_M3=nedges
% load
% nedges=cellfun(@(X) length(X),originaltissue.c)
% ne_M4=nedges
% load
% nedges=cellfun(@(X) length(X),originaltissue.c)
% ne_M5=nedges
% load
% nedges=cellfun(@(X) length(X),originaltissue.c)
% ne_Z=nedges
% ne_M={ne_M1, ne_M2, ne_M3, ne_M4, ne_M5}
% ne_M=[ne_M{:}]
%%
name={'M1' 'M2' 'M3' 'M4' 'M5' 'Z'};
n=[1 2 3 4 5 6 7 8 9 10];

for i=1:6
    eval(['t=sort(ne_' name{i} ');']);
    tb=find(t==1,1,'last');
    if isempty(tb)
        n(1+i,1)=0;
    else
        n(1+i,1)=tb;
    end
    for j=2:length(n)
        tb=find(t==j,1,'last');
        if isempty(tb)
            n(1+i,j)=0;
        else
            n(1+i,j)=find(t==j,1,'last');
        end
    end
    nef(1+i,2:10)=n(1+i,2:10)-n(1+i,1:9);    
    nef(nef<0)=0;
    n(1+i,:)=nef(1+i,:)/sum(nef(1+i,:));
end
figure(1);
bar(n(1,4:9),n(2:end,4:9)'*100);
legend(name);
n2=sum(nef(2:6,:))
n2=[n2/sum(n2); n(7,:)];
figure(2);
bar(n(1,4:9),n2(:,4:9)'*100);
legend('Microsorum', 'Zinnia')

    