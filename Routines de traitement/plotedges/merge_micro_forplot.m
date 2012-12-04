ot=cellNetwork();
le=zeros(1,5);
ls=zeros(1,5);
Ne=zeros(1,5);
Nv=zeros(1,5);
for i=1:5
    filename=['Microsorum' int2str(i) '_1-treated'];
    load(filename);
    le(i)=length(originaltissue.e);
    lv(i)=length(originaltissue.v);
    Ne(i)=length(ot.e);
    Nv(i)=length(ot.v);
    ot.e(Ne(i)+1:Ne(i)+le(i),:)=originaltissue.e;
    ot.e(Ne(i)+1:Ne(i)+le(i),1:2)=ot.e(Ne(i)+1:Ne(i)+le(i),1:2)+Nv(i);
    ot.v(Nv(i)+1:Nv(i)+lv(i),:)=originaltissue.v;
end
Ne(i+1)=length(ot.e);
Nv(i+1)=length(ot.v);
save('Micro_1-5_forplotedges.mat', 'ot');