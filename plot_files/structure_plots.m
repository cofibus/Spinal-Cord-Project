% MN pool without:
clear
load('N_wo.mat');
load('N_w.mat');
load('N_r.mat');

wo_col = [0.996 0.380 0];
w_col =  [0.47 0.368 0.941];
r_col = [0.862 0.149 0.498];
%% 
% ---- Tibialis anterior pool ----
function MN_distribution(N_wo, N_w, N_r, wo_col, w_col, r_col)

    figure('Position', [100 100 850 1000]);
    set(gcf, 'InvertHardcopy', 'off');
    set(gcf, 'Color', 'w');

    MN_ind_wo = N_wo.MnID == 'Tibialis Anterior' & N_wo.Latera == 1; % Keeping not only tib but righ-side tib
    % Monosinaptically connected to the MN pool
    IN_ind_wo = any(N_wo.ConnMat(MN_ind_wo, :), 1); % selecting all columns & 
    % tibialis rows and checking along the 1st dim
    
    % Positions
    MN_pos_wo = N_wo.Position(MN_ind_wo, :); % Positions for the right-side MN pool
    IN_pos_wo = N_wo.Position(IN_ind_wo, :); % Positions for the Interneurons 
    % (all pre no matter the latera)
    
    % Scatter plot
    sc1 = axes('Position', [0.08 0.11 0.18 0.85]);
    scatter(sc1, IN_pos_wo(:, 1), IN_pos_wo(:, 2), 12, [0.529 0.529 0.529], 'filled');
    hold(sc1, 'on');
    scatter(sc1, MN_pos_wo(:, 1), MN_pos_wo(:, 2), 20, wo_col ,'o', 'filled');
    set(sc1, 'YDir', 'reverse', 'YTickLabel',[], 'XTickLabel',[], ...
        'YTick',[], 'XTick',[], 'Color', [0.9 0.9 0.9]);
    xline(sc1, 0, ':', 'LineWidth', 1);
    title(sc1, 'Without genetic bias');
    
    % Segment lines
    seg = unique(N_wo.Segment, 'stable');
    x = N_wo.Position(:, 1);
    y = N_wo.Position(:, 2);
    [seg_min, seg_max] = arrayfun(@(s) bounds(y(N_wo.Segment == s)), seg, 'UniformOutput', false);
    for i = 1:length(seg)
        lim = seg_max{i};
        yline(sc1, lim, ':', 'LineWidth', 1);
        meanY = (seg_min{i} + seg_max{i}) / 2;
        offset = min(x) * 1.35;
        text(sc1, offset, meanY, char(seg(i)), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
    end

    hold (sc1, 'off');
   
    sections = 30;
    y_lim = linspace(min(y), max(y), sections + 1); % Limit between each section
    IN_x_sec = histcounts(y(IN_ind_wo), y_lim); % INs in each section
    norm_wo = IN_x_sec / sum(IN_x_sec);
    meanY = (y_lim(1:end-1) + y_lim(2:end)) / 2; % better than in the limit of a section
    offset = max(x) * 1.2;
    
    % Distribution line plot
    line1 = axes('Position', [0.27 0.11 0.04 0.85]);
    plot(offset + norm_wo * 8000, meanY, 'k-', 'LineWidth', 2, 'Color', wo_col);
    set(line1, 'YTickLabel',[], 'XTickLabel',[], 'YTick',[], 'XTick',[], ...
         'Box', 'off', 'XColor', 'none', 'YColor', 'none', 'Color', 'none'); % Invisible axes
    
    linkaxes([sc1, line1], 'y');
    set([sc1, line1], 'YDir', 'reverse');
    
    % With model
    MN_ind_w = N_w.MnID == 'Tibialis Anterior' & N_w.Latera == 1;
    IN_ind_w = any(N_w.ConnMat(MN_ind_w, :), 1);
    MN_pos_w = N_w.Position(MN_ind_w, :); 
    IN_pos_w = N_w.Position(IN_ind_w, :); 
    
    sc2 = axes('Position', [0.37 0.11 0.18 0.85]);
    scatter(sc2, IN_pos_w(:, 1), IN_pos_w(:, 2), 12, [0.529 0.529 0.529], 'filled');
    hold(sc2, 'on');
    scatter(sc2, MN_pos_w(:, 1), MN_pos_w(:, 2), 20, w_col , 'o', 'filled');
    hold(sc2, 'off');
    set(sc2, 'YDir', 'reverse', 'YTickLabel',[], 'XTickLabel',[], ...
        'YTick',[], 'XTick',[], 'Color', [0.9 0.9 0.9]);
    xline(sc2, 0, ':', 'LineWidth', 1);
    title(sc2, 'With genetic bias');
    
    seg = unique(N_w.Segment);
    x = N_w.Position(:, 1);
    y = N_w.Position(:, 2);
    [seg_min, seg_max] = arrayfun(@(s) bounds(y(N_w.Segment == s)), seg, 'UniformOutput', false);
    
    for i = 1:length(seg)
        lim = seg_max{i};
        yline(sc2, lim, ':', 'LineWidth', 1);
        meanY = (seg_min{i} + seg_max{i}) / 2;
        offset = min(x) * 1.35;
        text(sc2, offset, meanY, char(seg(i)), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
    end
    
    y_lim = linspace(min(y), max(y), sections + 1);
    IN_x_sec = histcounts(y(IN_ind_w), y_lim);
    norm_w = IN_x_sec / sum(IN_x_sec);
    meanY = (y_lim(1:end-1) + y_lim(2:end)) / 2;
    offset = max(x) * 1.2;
    
    line2 = axes('Position', [0.56 0.11 0.04 0.85]);
    plot(offset + norm_wo * 8000, meanY, '-', 'LineWidth', 2, 'Color', wo_col);
    hold(line2, 'on');
    plot(offset + norm_w * 8000, meanY, 'k-', 'LineWidth', 2, 'Color', w_col);
    hold(line2, 'off');
    set(line2, 'YTickLabel',[], 'XTickLabel',[], 'YTick',[], 'XTick',[], ...
         'Box', 'off', 'XColor', 'none', 'YColor', 'none', 'Color', 'none');
    
    linkaxes([sc2, line2], 'y');
    set([sc2, line2], 'YDir', 'reverse');
    
    % Random model
    MN_ind_r = N_r.MnID == 'Tibialis Anterior' & N_r.Latera == 1;
    IN_ind_r = any(N_r.ConnMat(MN_ind_r, :), 1);
    MN_pos_r = N_r.Position(MN_ind_r, :);
    IN_pos_r = N_r.Position(IN_ind_r, :);
    
    sc3 = axes('Position', [0.66, 0.11, 0.18, 0.85]);
    scatter(sc3, IN_pos_r(:, 1), IN_pos_r(:, 2), 12, [0.529 0.529 0.529], 'filled');
    hold(sc3, 'on');
    scatter(sc3, MN_pos_r(:, 1), MN_pos_r(:, 2), 20, r_col , 'o', 'filled');
    hold(sc3, 'off');
    set(sc3, 'YDir', 'reverse', 'YTickLabel',[], 'XTickLabel',[], ...
        'YTick',[], 'XTick',[], 'Color', [0.9 0.9 0.9]);
    xline(sc3, 0, ':', 'LineWidth', 1);
    title(sc3, 'Random bias');
    
    seg = unique(N_r.Segment);
    x = N_r.Position(:, 1);
    y = N_r.Position(:, 2);
    [seg_min, seg_max] = arrayfun(@(s) bounds(y(N_r.Segment == s)), seg, 'UniformOutput', false);
    for i = 1:length(seg)
        lim = seg_max{i};
        yline(sc3, lim, ':', 'LineWidth', 1);
        meanY = (seg_min{i} + seg_max{i}) / 2;
        offset = min(x) * 1.35;
        text(sc3, offset, meanY, char(seg(i)), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
    end
    
    y_lim = linspace(min(y), max(y), sections + 1);
    IN_x_sec = histcounts(y(IN_ind_r), y_lim);
    norm_r = IN_x_sec / sum(IN_x_sec);
    meanY = (y_lim(1:end-1) + y_lim(2:end)) / 2;
    offset = max(x) * 1.2;
    
    line3 = axes('Position', [0.85, 0.11, 0.04, 0.85]);
    plot(offset + norm_wo * 8000, meanY, '-', 'LineWidth', 2, 'Color', wo_col);
    hold(line3, 'on');
    plot(offset + norm_r * 8000, meanY, '-', 'LineWidth', 2, 'Color', r_col);
    plot(offset + norm_w * 8000, meanY, '-', 'LineWidth', 2, 'Color', w_col);
    hold(line3, 'off');
    set(line3, 'YTickLabel',[], 'XTickLabel',[], 'YTick',[], 'XTick',[], ...
         'Box', 'off', 'XColor', 'none', 'YColor', 'none', 'Color', 'none');
    
    linkaxes([sc3, line3], 'y');
    set([sc3, line3], 'YDir', 'reverse');

    print('MN_distribution', '-dpng', '-r300');
end

% ---- FRACTIONS ----
function fractions(N_wo,  N_w, N_r, wo_col, w_col, r_col)
    figure('Position', [100 500 600 1000]);
    subplot(3, 1, 1);
    MN_ind_wo = N_wo.MnID == 'Tibialis Anterior' & N_wo.Latera == 1;
    MN_ind_w = N_w.MnID == 'Tibialis Anterior' & N_w.Latera == 1;
    MN_ind_r = N_r.MnID == 'Tibialis Anterior' & N_r.Latera == 1;
    IN_ind_wo = any(N_wo.ConnMat(MN_ind_wo, :), 1);
    IN_ind_w = any(N_w.ConnMat(MN_ind_w, :), 1);
    IN_ind_r = any(N_r.ConnMat(MN_ind_r, :), 1);
    
    % -------- Celltype fractions --------
    % without
    %unqT = unique(N_wo.Types); % Written by hand so that I can specify the order
    unqT = categorical(cellstr(["V2a-2", "V2a-1", "V2b", "V0d", "V0v", "DI6", "V3", "MN", "V1-Foxp2", "V1-Pou6f2", "V1-Rensh", "V1-Sp8"]));
    IN_f_wo = zeros(size(unqT));
    for i = 1:length(unqT)
        IN_f_wo(i) = sum(IN_ind_wo & (N_wo.Types' == unqT(i)));
    end
    
    IN_f_wo = IN_f_wo / sum(IN_ind_wo);
    
    % with
    IN_f_w = zeros(size(unqT));
    for i = 1:length(unqT)
        IN_f_w(i) = sum(IN_ind_w & (N_w.Types' == unqT(i)));
    end
    
    IN_f_w = IN_f_w / sum(IN_ind_w);
    
    % random
    IN_f_r = zeros(size(unqT));
    for i = 1:length(unqT)
        IN_f_r(i) = sum(IN_ind_r & (N_r.Types' == unqT(i)));
    end
    
    IN_f_r = IN_f_r / sum(IN_ind_r);
    
    fractions = [IN_f_wo' * 100, IN_f_w' * 100, IN_f_r' * 100]; 
    % column 1 is without, column 2 is with, column 3 is random
    
    b = bar(fractions, 0.9);
    b(1).FaceColor = wo_col;
    b(2).FaceColor = w_col;
    b(3).FaceColor = r_col;
    
    xlabel('Celltypes');
    ylabel('Fraction of premotor (%)');
    set(gca,'TickLength',[0 0])
    ylim padded;
    
    xticks(1:length(unqT));
    xticklabels(unqT);
    
    % -------- Laminae fractions --------
    subplot(3, 1, 2);
    % without
    first_char_num = str2double(extractBefore(string(N_wo.Layers), 2));
    
    % Create layers matrix, it assigns 1 to 9 checking first character, and 10
    % checking the 2 first characters. 
    layers_mat = nan(size(N_wo.Layers));
    layers_mat(first_char_num >= 1 & first_char_num <= 9) = first_char_num(first_char_num >= 1 & first_char_num <= 9);
    layers_mat(strncmp(string(N_wo.Layers), '10', 2)) = 10;
    
    unqT = 1:10;
  
    IN_f_wo = zeros(size(unqT));
    for i = 1:length(unqT)
        IN_f_wo(i) = sum(IN_ind_wo & (layers_mat' == unqT(i)));
    end
    IN_f_wo = IN_f_wo / sum(IN_ind_wo);
    
    % with
    first_char_num = str2double(extractBefore(string(N_w.Layers), 2));
    
    layers_mat = nan(size(N_w.Layers));
    layers_mat(first_char_num >= 1 & first_char_num <= 9) = first_char_num(first_char_num >= 1 & first_char_num <= 9);
    layers_mat(strncmp(string(N_w.Layers), '10', 2)) = 10;
    
    IN_f_w = zeros(size(unqT));
    for i = 1:length(unqT)
        IN_f_w(i) = sum(IN_ind_w & (layers_mat' == unqT(i)));
    end
    IN_f_w = IN_f_w / sum(IN_ind_w);

    % random
    first_char_num = str2double(extractBefore(string(N_r.Layers), 2));
    
    layers_mat = nan(size(N_r.Layers));
    layers_mat(first_char_num >= 1 & first_char_num <= 9) = first_char_num(first_char_num >= 1 & first_char_num <= 9);
    layers_mat(strncmp(string(N_r.Layers), '10', 2)) = 10;
    
    IN_f_r = zeros(size(unqT));
    for i = 1:length(unqT)
        IN_f_r(i) = sum(IN_ind_r & (layers_mat' == unqT(i)));
    end
    IN_f_r = IN_f_r / sum(IN_ind_r);
    
    fractions = [IN_f_wo' * 100, IN_f_w' * 100, IN_f_r' * 100]; 
    
    b = bar(fractions, 0.9);
    b(1).FaceColor = wo_col;
    b(2).FaceColor = w_col;
    b(3).FaceColor = r_col;
    
    xlabel('Laminae');
    ylabel('Fraction of premotor (%)');
    set(gca,'TickLength',[0 0])
    ylim padded;
    
    xticks(1:length(unqT));
    xticklabels({'I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX', 'X'});
    
    % -------- Segment fractions --------
    subplot(3, 1, 3);
    seg = unique(N_wo.Segment, 'stable');
    
    % without
    IN_f_wo = zeros(size(seg));
    for i = 1:length(seg)
        IN_f_wo(i) = sum(IN_ind_wo & (N_wo.Segment' == seg(i)));
    end
    IN_f_wo = IN_f_wo / sum(IN_ind_wo);
    
    % with
    IN_f_w = zeros(size(seg));
    for i = 1:length(seg)
        IN_f_w(i) = sum(IN_ind_w & (N_w.Segment' == seg(i)));
    end
    
    IN_f_w = IN_f_w / sum(IN_ind_w);

    % random
    IN_f_r = zeros(size(seg));
    for i = 1:length(seg)
        IN_f_r(i) = sum(IN_ind_r & (N_r.Segment' == seg(i)));
    end
    IN_f_r = IN_f_r / sum(IN_ind_r);
    
    fractions = [IN_f_wo * 100, IN_f_w * 100, IN_f_r * 100]; 

    b = bar(fractions, 0.9);
    b(1).FaceColor = wo_col; 
    b(2).FaceColor = w_col; 
    b(3).FaceColor = r_col;

    xlabel('Segments');
    ylabel('Fraction of premotor (%)');
    set(gca,'TickLength',[0 0])
    ylim padded;
    
    xticks(1:length(seg));
    xticklabels(seg);
    legend({'Without genetic bias', 'With genetic bias', 'Randomized projectome'}, 'FontSize',12);
    print('MN_fractions', '-dpng', '-r300');
end
%% 
MN_distribution(N_wo, N_w, N_r, wo_col, w_col, r_col);

fractions(N_wo, N_w, N_r, wo_col, w_col, r_col);
