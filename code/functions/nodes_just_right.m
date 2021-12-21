function [sampled_id,xyz_coords] = nodes_just_right(nodes,node_id,n_to_sample,empirical_net)

dist_emp = coord2dist(table2array(empirical_net(:,1:3)));
dist_emp = dist_emp(triu(true(size(dist_emp)),+1));

while true
    sampled_id = table2array(datasample(node_id,n_to_sample,'Replace',false));
   
    nodes_i = nodes(sampled_id,:);
    
    dist = coord2dist(table2array(nodes_i(:,1:3)));
    dist = dist(triu(true(size(dist)),+1));
    
    if min(dist) < min(dist_emp)
        %fprintf('\nskipping, some nodes are too close\n')
        continue
    elseif (mean(dist) < (mean(dist_emp) - std(dist_emp)) || mean(dist) > (mean(dist_emp) + std(dist_emp)))
        %fprintf('skipping, mean node distance too large or small\n')
        continue
    elseif (max(dist) < (max(dist_emp) - std(dist_emp)) || max(dist) > (max(dist_emp) + std(dist_emp)))
        %fprintf('skipping, some nodes are too far\n')
        continue
    else
        xyz_coords = nodes_i;
        break
    end
end

% resulting net properties
%fprintf('\n random net min dist = %i\n mean dist = %i\n max dist = %i\n',min(dist),mean(dist),max(dist))

end