%%
for i=1:10
    tissue=originaltissue20;
    tissue.e([2,4],3)=[0.5-2*i*0.041 -0.4];
    name=['originaltissue22_' num2str(i) ];
    eval([name '=tissue;']);
    save('tissue_test', name, '-append');
end

%%
originaltissue20=originaltissue3;
originaltissue20.e(:,3)=[0.2; 0.5; 0.4; 0.1; 0.5];