function duplicate = check_duplicates(y,all_y)

duplicate = 0;

ny = size(all_y,2)-1;

%res = zeros(ny,1);

for y_i = 1:ny
    %if all(sort(table2array(y)) == sort(table2array(all_y(:,y_i))))
    if all(sort(y) == sort(table2array(all_y(:,y_i))))
        duplicate = 1;
        break
    end
    
    %n_same = sum(sort(table2array(y)) == sort(table2array(all_y(:,y_i))));
    n_same = sum(sort(y) == sort(table2array(all_y(:,y_i))));
    
    if n_same > (size(all_y,1) - 2)
        duplicate = 1;
        break
    end    
end

if duplicate == 1
    fprintf("duplicate, skipping...");
end

end