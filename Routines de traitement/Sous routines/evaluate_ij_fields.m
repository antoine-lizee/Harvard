input_var=planes;
i='P';
j='Q';
iv1=size(input_var,1);
iv2=size(input_var,2);
h=cell(iv1,iv2);

for m=1:iv1
    for n=1:iv2
       % h(m,n)={[mat2str(input_var(m,n).(i)) ' , ' mat2str(input_var(m,n).(j))]};
        h{m,n}= [input_var(m,n).(i),input_var(m,n).(i)];
    end
end

h