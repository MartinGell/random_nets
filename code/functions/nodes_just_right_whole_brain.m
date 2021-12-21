function [sampled_id,xyz_coords] = nodes_just_right_whole_brain(n_to_sample,empirical_net)

dist_emp = coord2dist(table2array(empirical_net(:,1:3)));
dist_emp = dist_emp(triu(true(size(dist_emp)),+1));

while true
    %y = datasample(node_id,n_to_sample,'Replace',false);
    %nodes_i = nodes(table2array(y),:);
    
    GM_mask = spm_vol('/home/mgell/Matlab_toolboxes/Juspace_v1/PETatlas/mask/PET_GM_mask_bin.nii');
    GM_data = spm_read_vols(GM_mask);
    GM_id = find(GM_data);

    sampled_id = datasample(GM_id,n_to_sample,'Replace',false);
    [i,j,k] = ind2sub(GM_mask(1).dim,sampled_id);
    phys_space = [i j k]';
    %phys_space = phys_space';
    vox_space = GM_mask(1).mat * [phys_space; ones(1,size(phys_space,2))];
    vox_space = vox_space(1:3,:);
    % the above is equivalent to:
    %x = srow_x(1) * i + srow_x(2) * j + srow_x(3) * k + srow_x(4);
    %y = srow_y(1) * i + srow_y(2) * j + srow_y(3) * k + srow_y(4);
    %z = srow_z(1) * i + srow_z(2) * j + srow_z(3) * k + srow_z(4);
    
    dist = coord2dist(vox_space');
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
        xyz_coords = vox_space';
        break
    end
end

% resulting net properties
%fprintf('\n random net min dist = %i\n mean dist = %i\n max dist = %i\n',min(dist),mean(dist),max(dist))

end