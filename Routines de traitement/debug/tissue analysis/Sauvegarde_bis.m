%% Saves all the variables with the _bis to compare them
Who=whos; %% create the struct
for w=1:length(Who)
    eval([Who(w).name '_bis = ' Who(w).name ';']);
    %clear((Who(w).name));
end


