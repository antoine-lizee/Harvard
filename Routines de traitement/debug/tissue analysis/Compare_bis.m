%% Debugging : comparing the similar variables (with or without "-bis") and showing them when they are different

deb_var=who; % Fetching the name of the variables in an array of string
nbr_var=length(deb_var);
% d=[strcmp(deb_var(1:nbr_var-1), [deb_var(2:nbr_var) '_bis']); false];
% good try, but it seems impossible to get this | properly...
d1= false(nbr_var,1); % Initializing the two logical vectors
d2=d1;
for z=1:nbr_var-1 % Feeding the two logical vectors
    d1(z)=strcmp(deb_var{z+1}, [deb_var{z} '_bis']) && ~(eval(['isequal(' deb_var{z} ',' deb_var{z+1} ')']));
    d2(z)=strcmp(deb_var{z+1}, [deb_var{z} '_bis']);
end
common_variables=deb_var(d2); length(common_variables) % Displaying the number of common variables
different_variables=deb_var(d1) % Displaying the variables which have been backed up and which are different
