function [ nodes_vec ] = nodes( tree_struct )
% Generates an alternative representation for a given partition tree. 
%
% Input:
%     tree_struct - cell of struct arrays representing a tree 
% 
% Output:
%     nodes_vec - vector of parents pointers, in the sense that tree_vec(i)
%         represents the index of the i-th node parent, and tree_vec(1)=0 
%         for the root.
%--------------------------------------------------------------------------

n_levels = length(tree_struct);
nodes_vec = [0]; %root
pointer = 0;

for ii = (n_levels-1):(-1):1,
    temp = tree_struct{ii}.super_folders + pointer;
    pointer = length(nodes_vec) ;
    nodes_vec = [nodes_vec, temp];
end

end