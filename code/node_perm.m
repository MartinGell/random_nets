
clear

addpath('/home/mgell/Matlab_toolboxes/spm12');

% directory with project
wd = '/home/mgell/Work/random_nets';
addpath(genpath(wd));
cd(wd);
%addpath('/home/mgell/Matlab_toolboxes/NIFTI_toolbox-master');

% n of random nets
node_size = 5 % sphere size
use_mask = 0
nperm = 10
%receptor = 'NAT' % NAT or Sert5HT
all_masks = 0 % print all networks as nii?
community = ''; % whole network or specific community? Need to provide community assignment. 
              % Use same value as in community assignment. e.g. 1 if only
              % want community called 1.

% Template network
imp_net = readtable(strcat(wd, '/networks/All_localmax.txt'));
imp_net([2,24],:) = [];

% Atlas to sample from - needs to be xyz coordiantes + name/id in last column
nodes = readtable(strcat(wd, '/networks/power_nodes.node'), 'FileType','text');
node_id = nodes(:,4);

% Output dir
outputdir = strcat(wd, '/res/');

% switch receptor
%     case 'NAT'
%         PET_map = '/home/mgell/Matlab_toolboxes/Juspace_v1/PETatlas/NAT_MRB_c11_2x2x2.nii.gz';
%         outputdir = '/home/mgell/Work/FC/PET/permutation/NAT/';
% 
%     case 'Sert5HT'
%         PET_map = '/home/mgell/Matlab_toolboxes/Juspace_v1/PETatlas/5HT1a_WAY_HC36_2x2x2.nii.gz';
%         outputdir = '/home/mgell/Work/FC/PET/permutation/Sert5HT/';
% end
% 
% % Load PET map
% V = spm_vol(PET_map);

% subsample network to community if requested
if ~isempty(community)
    communities = readtable('/home/mgell/Work/FC/hcp/text_files/res/replic_5mm_net_communities.csv','ReadVariableNames',false);
    imp_net(:,5) = communities;
    rows = imp_net.Var5 == community;
    imp_net = imp_net(rows,:);
end

node_n = size(imp_net,1);
all_y = table(zeros(node_n,1));

% Do the thing
for perm_i = 1:nperm

    fprintf('permutation: %i\n', perm_i) % permutation number
    
    duplicate = 1;
    while duplicate == 1 
        % create random net
        [y, nodes_i] = nodes_just_right(nodes,node_id,node_n,imp_net);
        %[y, nodes_i] = nodes_just_right_whole_brain(node_n,imp_net);
        %nodes_i = array2table([nodes_i, y],'VariableNames',{'x','y','z','id'});
        
        % check if random network already exists
        duplicate = check_duplicates(y,all_y); 
    end
    
    all_y(:,perm_i) = array2table(y); % save all ys so we can avoid duplicates of random samples
    
    % save random network
    coords_i = [outputdir 'coords' num2str(perm_i) '.txt'];
    if all_masks == 1
        net_i = [outputdir 'net' num2str(perm_i) '.nii'];
    else
        net_i = [outputdir 'net' '.nii'];
    end
    writetable(nodes_i,coords_i,'WriteVariableNames',false,'Delimiter','tab')
    %create_net(node_size,use_mask,coords_i,V,net_i);
    %parcellate(PET_map,mask_i,fullfile(outputdir, [receptor '_' num2str(perm_i) '.txt'])); 
end

