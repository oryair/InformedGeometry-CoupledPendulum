function [ ] = PlotEmbedding( embedding, true_perm, title_str )
% Displays the 3-D diffusion embedding, i.e. the 3 eigenvectors corresponding
% to the 3 largest eigenvalues, colored according to the ground truth.
%
% Inputs:
%     embedding - the embedding of the data set into Euclidean space
%     true_perm - permutation executed on the original data
%     title_str - title string (optional)
%--------------------------------------------------------------------------

% color_mat = jet(length(embedding(:,1)));
% color_mat = color_mat(true_perm,:);

switch size(embedding, 2)
    case 1
        plot(true_perm, embedding, 'x'); %'MarkerEdgeColor', 'k'
        xlabel(title_str), ylabel('\psi_1');
    case 2
        scatter(embedding(:,1), embedding(:,2), 30, true_perm, 'filled'); %'MarkerEdgeColor', 'k'
        xlabel('\psi_1'), ylabel('\psi_2');

    otherwise
        scatter3(embedding(:,1), embedding(:,2), embedding(:,3), 30, true_perm, 'filled'); %'MarkerEdgeColor', 'k'
        xlabel('\psi_1'), ylabel('\psi_2'), zlabel('\psi_3'), hold on

end

if nargin < 3,
    title('Diffusion Embedding'), hold off
else
    title([title_str ' Embedding']), hold off
end
colorbar;
% [~, indsort] = sort(true_perm);
% if size(embedding, 2) > 2
%     figure;subplot(3,1,1);plot(embedding(indsort,1));
%     title('Diffusion Embedding');ylabel('\psi_1')
%     subplot(3,1,2);plot(embedding(indsort,2));
%     ylabel('\psi_2')
%     subplot(3,1,3);plot(embedding(indsort,3));
%     ylabel('\psi_3')
% else
%     figure;subplot(2,1,1);plot(embedding(indsort,1));
%     ylabel('\psi_1');
%     subplot(2,1,2);plot(embedding(indsort,2));
%     ylabel('\psi_2');
% end
end
