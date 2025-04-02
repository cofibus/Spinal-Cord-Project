%% Connectivity matrix WITHOUT genetic bias
clear
load('N_wo.mat');
load('N_w.mat');
load('N_r.mat');
%% 
function ct_mat = connectivity_plot(N, ptitle, filename)
    unq = unique(N.Types);
    
    types = length(unq);
    ct_mat = NaN(types);
    
    for pre = 1:types
        for post = 1:types
            pre_ind = N.Types == unq(pre);
            post_ind = N.Types == unq(post);
    
            %(post_ind&pre_ind');
    
            % Type-pair sub-matrix (pre w/ post)
            submatrix = abs(N.ConnMat(post_ind, pre_ind)); 
    
            % mean of cell-pairs
            ct_mat(post, pre) = mean(submatrix,'all');
        end
    end
    
    % Normalize
    ct_mat = ct_mat ./ max(ct_mat,[],1); 
  
    figure;
    imagesc(ct_mat); 
    colormap('gray');
    colorbar;
    xticks(1:types);
    yticks(1:types);
    xticklabels(unq);
    yticklabels(unq);
    xlabel('Presynaptic');
    ylabel('Postsynaptic');
    title(ptitle);
    print(filename, '-dpng', '-r300');
end
%% 
function diff_mat_plot(ct_mat_wo, ct_mat_w, ct_mat_r)
    diffMat1 = ct_mat_wo - ct_mat_w;
    % diffMat2 = ct_mat_wo - ct_mat_r;
    % diffMat3 = ct_mat_w - ct_mat_r;
    
    figure();
    % subplot(1, 3, 1);
    imagesc(abs(diffMat1));
    colormap('cool');
    colorbar;
    xlabel('Presynaptic');
    ylabel('Postsynaptic');
    title('Without/With');
    print('diff_conn_w_wo', '-dpng', '-r300')

    % subplot(1, 3, 2);
    % imagesc(abs(diffMat2));
    % colormap('cool');
    % colorbar;
    % xlabel('Presynaptic');
    % ylabel('Postsynaptic');
    % title('Without/Random');
    % 
    % subplot(1, 3, 3);
    % imagesc(abs(diffMat3));
    % colormap('cool');
    % colorbar;
    % xlabel('Presynaptic');
    % ylabel('Postsynaptic');
    % title('With/Random');
    % sgtitle('Difference in connectivity');
    % print('diff_conn_3', '-dpng', '-r300');
 end
%% 
ct_mat_wo = connectivity_plot(N_wo, 'Without genetic bias', 'CT_conn_wo');
ct_mat_w = connectivity_plot(N_w, 'With genetic bias', 'CT_conn_w');
ct_mat_r = connectivity_plot(N_r, 'Randomized projectome', 'CT_conn_r');
diff_mat_plot(ct_mat_wo, ct_mat_w, ct_mat_r);