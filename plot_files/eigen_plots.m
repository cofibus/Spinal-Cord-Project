clear
load('N_wo.mat');
load('N_w.mat');
load('N_r.mat');

r_col = [0.862 0.149 0.498];
wo_col = [0.996 0.380 0];
w_col =  [0.47 0.368 0.941];

[evec_wo, eval_wo] = eig(N_wo.ConnMat);
[evec_w, eval_w] = eig(N_w.ConnMat);
[evec_r, eval_r] = eig(N_r.ConnMat);
eval_wo = diag(eval_wo);
eval_w = diag(eval_w);
eval_r = diag(eval_r);

[~, sorted_idx_wo] = sort(real(eval_wo), 'descend');
[~, sorted_idx_w] = sort(real(eval_w), 'descend');
[~, sorted_idx_r] = sort(real(eval_r), 'descend');

function [d_evec_wo, d_evec_w, d_evec_r] = dominant_evec(mode, sorted_idx_wo, sorted_idx_w, sorted_idx_r, evec_wo, evec_w, evec_r)
d_evec_wo = evec_wo(:, sorted_idx_wo(mode));
d_evec_w = evec_w(:, sorted_idx_w(mode));
d_evec_r = evec_r(:, sorted_idx_r(mode));
end

[d_evec_wo, d_evec_w, d_evec_r] = dominant_evec(1, sorted_idx_wo, sorted_idx_w, sorted_idx_r, evec_wo, evec_w, evec_r);
[d_evec_wo2, d_evec_w2, d_evec_r2] = dominant_evec(3, sorted_idx_wo, sorted_idx_w, sorted_idx_r, evec_wo, evec_w, evec_r);

d_e_wo = eval_wo(sorted_idx_wo(1));
d_e_w = eval_w(sorted_idx_w(1));
d_e_r = eval_r(sorted_idx_r(1));
%%
function eval_plots(eval_wo, eval_w, eval_r, d_e_wo, d_e_w, d_e_r, wo_col, w_col, r_col)
    % Without
    figure;
    hold on;
    xline(1, 'r:', 'LineWidth', 2);
    xline(0, 'k:', 'LineWidth', 1);
    yline(0, 'k:', 'LineWidth', 1);
    scatter(real(eval_wo), imag(eval_wo), 40, 'o', 'MarkerEdgeColor', wo_col);
    
    xlabel('Real');
    ylabel('Imaginary');
    title('Eigen spectrum');
    title('Without genetic bias');
    hold off;
    print('wo_d_eval', '-dpng', '-r300');
    
    % With
    figure;
    hold on;
    xline(1, 'r:', 'LineWidth', 2);
    xline(0, 'k:', 'LineWidth', 1);
    yline(0, 'k:', 'LineWidth', 1);
    scatter(real(eval_w), imag(eval_w), 40, 'o',  'MarkerEdgeColor', w_col);
    
    xlabel('Real');
    ylabel('Imaginary');
    title('Eigen spectrum');
    title('With genetic bias');
    hold off;
    print('w_d_eval', '-dpng', '-r300');
    
    % Random
    figure;
    hold on;
    xline(1, 'r:', 'LineWidth', 2);
    xline(0, 'k:', 'LineWidth', 1);
    yline(0, 'k:', 'LineWidth', 1);
    scatter(real(eval_r), imag(eval_r), 40, 'o',  'MarkerEdgeColor', r_col);
    
    xlabel('Real');
    ylabel('Imaginary');
    title('Eigen spectrum');
    title('Randomised projectome');
    hold off;
    print('r_d_eval', '-dpng', '-r300');
    
    % --- Eigenspectrum ---
    figure;
    hold on;
    xline(1, 'r--', 'LineWidth', 1);
    xline(0, 'k:', 'LineWidth', 1);
    yline(0, 'k:', 'LineWidth', 1);
    l1 = scatter(real(eval_wo), imag(eval_wo), 40, 'x', 'LineWidth', 2.3, 'DisplayName', 'Without genetic bias', 'MarkerEdgeColor', wo_col);
    l2 = scatter(real(eval_w), imag(eval_w), 40, 'o', 'filled', 'LineWidth', 1.2, 'DisplayName', 'With genetic bias', 'MarkerFaceColor', w_col);
    
    xlabel('Real');
    ylabel('Imaginary');
    title('Eigen spectrum');
    
    scatter(real(d_e_wo), imag(d_e_wo), 60, 'kx', 'LineWidth', 2, 'DisplayName', 'Dominant (Without)');
    scatter(real(d_e_w), imag(d_e_w), 60, 'ko', 'filled', 'LineWidth', 2, 'DisplayName', 'Dominant (With)');
    legend([l1, l2], 'Location', 'northwest', 'FontSize', 12);
    hold off;
    print('wo_w_d_eval', '-dpng', '-r300');
    
    % with random too
    figure;
    hold on;
    xline(1, 'r--', 'LineWidth', 1);
    xline(0, 'k:', 'LineWidth', 1);
    yline(0, 'k:', 'LineWidth', 1);
    l3 = scatter(real(eval_r), imag(eval_r), 40, 'd', 'LineWidth', 1.2, 'DisplayName', 'Random bias', 'MarkerEdgeColor', r_col, 'MarkerFaceColor', r_col);
    l2 = scatter(real(eval_w), imag(eval_w), 40, 'o', 'filled', 'LineWidth', 1.2, 'DisplayName', 'With genetic bias', 'MarkerFaceColor', w_col);
    l1 = scatter(real(eval_wo), imag(eval_wo), 40, 'x', 'LineWidth', 2.3, 'DisplayName', 'Without genetic bias', 'MarkerEdgeColor', wo_col);
    
    xlabel('Real');
    ylabel('Imaginary');
    title('Eigen spectrum');
    
    scatter(real(d_e_wo), imag(d_e_wo), 60, 'kx', 'LineWidth', 2, 'DisplayName', 'Dominant (Without)');
    scatter(real(d_e_w), imag(d_e_w), 60, 'ko', 'filled', 'LineWidth', 2, 'DisplayName', 'Dominant (With)');
    scatter(real(d_e_r), imag(d_e_r), 60, 'kd', 'MarkerEdgeColor', r_col, 'MarkerFaceColor', r_col,'LineWidth', 2, 'DisplayName', 'Dominant (Random)');
    legend([l1, l2, l3], 'Location', 'northwest', 'FontSize', 12);
    hold off;
    print('3_d_eval', '-dpng', '-r300');
    
    % --- Margin plot ---
    figure('Position', [100 100 800 600]);
    
    hold on;
    xline(1, 'r--', 'LineWidth', 1);
    xline(0, 'k:', 'LineWidth', 1);
    yline(0, 'k:', 'LineWidth', 1);
    scatter(real(eval_wo), imag(eval_wo), 40, 'kx', 'LineWidth', 1.2, 'DisplayName', 'Without genetic bias');
    scatter(real(eval_w), imag(eval_w), 40, 'ko', 'LineWidth', 1.2, 'DisplayName', 'With genetic bias');
    scatter(real(eval_r), imag(eval_r), 40, 'kd', 'LineWidth', 1.2, 'DisplayName', 'Randomised projectome');
    
    xlabel('Real');
    ylabel('Imaginary');
    title('Eigen spectrum');
    
    diff_magnitude1 = abs(d_e_wo - d_e_w);
    diff_magnitude2 = abs(d_e_wo - d_e_r);
    diff_magnitude3 = abs(d_e_w - d_e_r);
    
    text((real(d_e_wo) + real(d_e_w))/2.1, (imag(d_e_wo) + imag(d_e_w))/2 + 0.1, sprintf('%.2f', diff_magnitude1), 'HorizontalAlignment', 'center');
    text((real(d_e_wo) + real(d_e_r))/2.1, (imag(d_e_wo) + imag(d_e_r))/2 + 0.1, sprintf('%.2f', diff_magnitude2), 'HorizontalAlignment', 'center');
    text((real(d_e_w) + real(d_e_r))/2.1, (imag(d_e_w) + imag(d_e_r))/2 + 0.1, sprintf('%.2f', diff_magnitude3), 'HorizontalAlignment', 'center');
    
    distance_wo = abs(real(d_e_wo)-1);
    distance_w = abs(real(d_e_w)-1);
    distance_r = abs(real(d_e_r)-1);
    
    % text labels for distances above the arrows
    text(real(d_e_wo) - distance_wo/2, imag(d_e_wo) + 0.2, sprintf('%.2f', distance_wo), 'HorizontalAlignment', 'center', 'Color', wo_col);
    text(real(d_e_w) - distance_w/2, imag(d_e_w) + 0.2, sprintf('%.2f', distance_w), 'HorizontalAlignment', 'center', 'Color', w_col);
    text(real(d_e_r) - distance_r/2, imag(d_e_r) + 0.2, sprintf('%.2f', distance_r), 'HorizontalAlignment', 'center', 'Color', r_col);
    
    l1 = scatter(real(d_e_wo), imag(d_e_wo), 100, 'x', 'LineWidth', 3, 'DisplayName', 'Without genetic bias', 'MarkerEdgeColor', wo_col);
    l2 = scatter(real(d_e_w), imag(d_e_w), 80, 'o', 'filled', 'LineWidth', 2, 'DisplayName', 'With genetic bias', 'MarkerFaceColor', w_col);
    l3 = scatter(real(d_e_r), imag(d_e_r), 80, 'd', 'filled', 'LineWidth', 2, 'DisplayName', 'Randomised projectome', 'MarkerFaceColor', r_col, 'MarkerEdgeColor', r_col);
    leg = legend([l1, l2, l3], 'Location', 'northwest', 'FontSize', 12);
    title(leg, 'Dominant Eigenvalue (Mode 1)')
    hold off;
    
    print('margin_d_eval', '-dpng', '-r300');
end
%% 
function plot_both_mode(d_evec_wo, d_evec_w, d_evec_r, d_evec_wo2, d_evec_w2, d_evec_r2, N_wo, N_w, N_r)

% ---- First Mode ----
phase_wo = angle(d_evec_wo);
dotSizes_wo = rescale(phase_wo, 1, 10);

% --- Without 1 ---
figure('Position', [100 100 1300 1000]);
set(gcf, 'InvertHardcopy', 'off');
set(gcf, 'Color', 'w');

subplot(1, 6, 1);
scatter(N_wo.Position(:, 1), N_wo.Position(:, 2), dotSizes_wo, 'k', 'filled');
hold on;
seg_wo = unique(N_wo.Segment, 'stable');
x_wo = N_wo.Position(:, 1);
y_wo = N_wo.Position(:, 2);
[seg_min_wo, seg_max_wo] = arrayfun(@(s) bounds(y_wo(N_wo.Segment == s)), seg_wo, 'UniformOutput', false);
for i = 1:length(seg_wo)
    lim = seg_max_wo{i};
    yline(lim, ':', 'LineWidth', 1);
    meanY = (seg_min_wo{i} + seg_max_wo{i}) / 2;
    offset = min(x_wo) * 1.3;
    text(offset, meanY, char(seg_wo(i)), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
end
set(gca, 'YTickLabel', [], 'Color', [0.98 0.98 0.98]);
hold off;

% ---- Second Mode ----
phase_wo = angle(d_evec_wo2);
dotSizes_wo = rescale(phase_wo, 1, 10);

% --- Without 2 ---
subplot(1, 6, 2);
scatter(N_wo.Position(:, 1), N_wo.Position(:, 2), dotSizes_wo, 'k', 'filled');
hold on;
seg_wo = unique(N_wo.Segment, 'stable');
x_wo = N_wo.Position(:, 1);
y_wo = N_wo.Position(:, 2);
[seg_min_wo, seg_max_wo] = arrayfun(@(s) bounds(y_wo(N_wo.Segment == s)), seg_wo, 'UniformOutput', false);
for i = 1:length(seg_wo)
    lim = seg_max_wo{i};
    yline(lim, ':', 'LineWidth', 1);
    meanY = (seg_min_wo{i} + seg_max_wo{i}) / 2;
    offset = max(x_wo) * 1.1;
    %text(offset, meanY, char(seg_wo(i)), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
end
set(gca, 'Color', [0.98 0.98 0.98]);
title({'Without genetic bias';' Mode 1                                           Mode2'}, ...
    'Position', [-1400 14000]);
hold off;

% --- With 1 ---
phase_w = angle(d_evec_w);
dotSizes_w = rescale(phase_w, 1, 10);

subplot(1, 6, 3);
scatter(N_w.Position(:, 1), N_w.Position(:, 2), dotSizes_w, 'k', 'filled');
hold on;
seg_w = unique(N_w.Segment, 'stable');
x_w = N_w.Position(:, 1);
y_w = N_w.Position(:, 2);
[seg_min_w, seg_max_w] = arrayfun(@(s) bounds(y_w(N_w.Segment == s)), seg_w, 'UniformOutput', false);
for i = 1:length(seg_w)
    lim = seg_max_w{i};
    yline(lim, ':', 'LineWidth', 1);
    meanY = (seg_min_w{i} + seg_max_w{i}) / 2;
    offset = max(x_w) * 1.1;
    %text(offset, meanY, char(seg_w(i)), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
end
set(gca, 'Color', [0.98 0.98 0.98]);
hold off;

% --- With 2 ---
phase_w = angle(d_evec_w2);
dotSizes_w = rescale(phase_w, 1, 10);

subplot(1, 6, 4);
scatter(N_w.Position(:, 1), N_w.Position(:, 2), dotSizes_w, 'k', 'filled');
hold on;
seg_w = unique(N_w.Segment, 'stable');
x_w = N_w.Position(:, 1);
y_w = N_w.Position(:, 2);
[seg_min_w, seg_max_w] = arrayfun(@(s) bounds(y_w(N_w.Segment == s)), seg_w, 'UniformOutput', false);
for i = 1:length(seg_w)
    lim = seg_max_w{i};
    yline(lim, ':', 'LineWidth', 1);
    meanY = (seg_min_w{i} + seg_max_w{i}) / 2;
    offset = max(x_w) * 1.1;
    %text(offset, meanY, char(seg_w(i)), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
end
set(gca, 'Color', [0.98 0.98 0.98]);
title({'With genetic bias';' Mode 1                                           Mode2'}, ...
    'Position', [-1400 14000]);
hold off;

% ---- Random -----
phase_r = angle(d_evec_r);
dotSizes_r = rescale(phase_r, 1, 10);

% Mode 1
subplot(1, 6, 5);
scatter(N_r.Position(:, 1), N_r.Position(:, 2), dotSizes_r, 'k', 'filled');
hold on;
seg_r = unique(N_r.Segment, 'stable');
x_r = N_r.Position(:, 1);
y_r = N_r.Position(:, 2);
[seg_min_r, seg_max_r] = arrayfun(@(s) bounds(y_r(N_r.Segment == s)), seg_r, 'UniformOutput', false);
for i = 1:length(seg_r)
    lim = seg_max_r{i};
    yline(lim, ':', 'LineWidth', 1);
    meanY = (seg_min_r{i} + seg_max_r{i}) / 2;
    offset = max(x_r) * 1.1;
    %text(offset, meanY, char(seg_r(i)), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
end
set(gca, 'Color', [0.98 0.98 0.98]);
hold off;

phase_r = angle(d_evec_r2);
dotSizes_r = rescale(phase_r, 1, 10);

% Mode 2
subplot(1, 6, 6);
scatter(N_r.Position(:, 1), N_r.Position(:, 2), dotSizes_r, 'k', 'filled');
hold on;
seg_r = unique(N_r.Segment, 'stable');
x_r = N_r.Position(:, 1);
y_r = N_r.Position(:, 2);
[seg_min_r, seg_max_r] = arrayfun(@(s) bounds(y_r(N_r.Segment == s)), seg_r, 'UniformOutput', false);
for i = 1:length(seg_r)
    lim = seg_max_r{i};
    yline(lim, ':', 'LineWidth', 1);
    meanY = (seg_min_r{i} + seg_max_r{i}) / 2;
    offset = max(x_r) * 1.3;
    text(offset, meanY, char(seg_r(i)), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
end
set(gca, 'Color', [0.98 0.98 0.98]);
hold off;
title({'Randomised projectome';' Mode 1                                           Mode2'}, ...
    'Position', [-1400 14000]);
print('all_d_modes', '-dpng', '-r300');
end
%% 
function d_evec_distribution(mode, d_evec_wo, d_evec_w, d_evec_r, N_wo, N_w, N_r)
    % ---- Second Mode ----
    path = '/Users/sabrina/Library/CloudStorage/OneDrive-UniversityofCopenhagen/SpinalCord_Project/Final figures';
    r_col = [0.862 0.149 0.498];
    wo_col = [0.996 0.380 0];
    w_col =  [0.47 0.368 0.941];
    phase_wo = angle(d_evec_wo);
    phase_w = angle(d_evec_w);
    phase_r = angle(d_evec_r);
    dotSizes_wo = rescale(phase_wo, 1, 10);
    dotSizes_w = rescale(phase_w, 1, 10);
    dotSizes_r = rescale(phase_r, 1, 10);

    % --- Calculations for other figures --- WITHOUT
    seg_wo = unique(N_wo.Segment, 'stable');
    x_wo = N_wo.Position(:, 1);
    y_wo = N_wo.Position(:, 2);
    [seg_min_wo, seg_max_wo] = arrayfun(@(s) bounds(y_wo(N_wo.Segment == s)), seg_wo, 'UniformOutput', false);

    left_indices = find(N_wo.Latera == -1);
    right_indices = find(N_wo.Latera == 1);
    
    sections = 50;
    y_lim = linspace(min(y_wo), max(y_wo), sections + 1);
    meanY = (y_lim(1:end-1) + y_lim(2:end)) / 2;
    devec_means_left_wo = arrayfun(@(i) mean(phase_wo(left_indices(y_wo(left_indices) >= y_lim(i) & y_wo(left_indices) < y_lim(i+1)))), 1:sections);
    devec_means_right_wo = arrayfun(@(i) mean(phase_wo(right_indices(y_wo(right_indices) >= y_lim(i) & y_wo(right_indices) < y_lim(i+1)))), 1:sections);

    % --- Figure 1: Scatter Plots & distribution line plot ---
    figure;
    fig = gcf;
    set(fig, 'Units', 'pixels');
    set(fig, 'Position', [100 100 900 1000]);

    % Without
    subplot(1, 3, 1);
    scatter(N_wo.Position(:, 1), N_wo.Position(:, 2), dotSizes_wo, 'k', 'filled');
    xlim("padded");
    title(['Mode ', num2str(mode), ': Without genetic bias']);
    hold on;
    for i = 1:length(seg_wo)
        lim = seg_max_wo{i};
        yline(lim, ':', 'LineWidth', 1);
        meanY2 = (seg_min_wo{i} + seg_max_wo{i}) / 2;
        offset = max(x_wo) * 1.35;
        %text(offset, meanY2, char(seg_wo(i)), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
    end
    plot(max(x_wo) * 1.2 + devec_means_right_wo * 50, meanY, '-', 'LineWidth', 2, 'Color', wo_col);
    plot(min(x_wo) * 1.4 + devec_means_left_wo * 50, meanY, '-', 'LineWidth', 2, 'Color', wo_col);
    hold off;

    % With
    % --- Calculations for other figures --- WITH
    seg_w = unique(N_w.Segment, 'stable');
    x_w = N_w.Position(:, 1);
    y_w = N_w.Position(:, 2);
    [seg_min_w, seg_max_w] = arrayfun(@(s) bounds(y_w(N_w.Segment == s)), seg_w, 'UniformOutput', false);

    left_indices = find(N_w.Latera == -1);
    right_indices = find(N_w.Latera == 1);
    
    y_lim = linspace(min(y_w), max(y_w), sections + 1);
    meanY = (y_lim(1:end-1) + y_lim(2:end)) / 2;
    devec_means_left_w = arrayfun(@(i) mean(phase_w(left_indices(y_w(left_indices) >= y_lim(i) & y_w(left_indices) < y_lim(i+1)))), 1:sections);
    devec_means_right_w = arrayfun(@(i) mean(phase_w(right_indices(y_w(right_indices) >= y_lim(i) & y_w(right_indices) < y_lim(i+1)))), 1:sections);
    
    subplot(1, 3, 2);
    scatter(N_w.Position(:, 1), N_w.Position(:, 2), dotSizes_w, 'k', 'filled');
    xlim("padded");
    title(['Mode ', num2str(mode), ': With genetic bias']);
    hold on;
    for i = 1:length(seg_w)
        lim = seg_max_w{i};
        yline(lim, ':', 'LineWidth', 1);
        meanY2 = (seg_min_w{i} + seg_max_w{i}) / 2;
        offset = max(x_w) * 1.35;
        %text(offset, meanY2, char(seg_w(i)), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
    end
    plot(max(x_w) * 1.2 + devec_means_right_w * 50, meanY, '-', 'LineWidth', 2, 'Color', w_col);
    plot(min(x_w) * 1.4 + devec_means_left_w * 50, meanY, '-', 'LineWidth', 2, 'Color', w_col);
    hold off;

    % --- Calculations for other figures --- RANDOM
    seg_r = unique(N_r.Segment, 'stable');
    x_r = N_r.Position(:, 1);
    y_r = N_r.Position(:, 2);
    [seg_min_r, seg_max_r] = arrayfun(@(s) bounds(y_r(N_r.Segment == s)), seg_r, 'UniformOutput', false);

    left_indices = find(N_r.Latera == -1);
    right_indices = find(N_r.Latera == 1);
    
    y_lim = linspace(min(y_r), max(y_r), sections + 1);
    meanY = (y_lim(1:end-1) + y_lim(2:end)) / 2;
    devec_means_left_r = arrayfun(@(i) mean(phase_r(left_indices(y_r(left_indices) >= y_lim(i) & y_r(left_indices) < y_lim(i+1)))), 1:sections);
    devec_means_right_r = arrayfun(@(i) mean(phase_r(right_indices(y_r(right_indices) >= y_lim(i) & y_r(right_indices) < y_lim(i+1)))), 1:sections);
    
    subplot(1, 3, 3);
    scatter(N_r.Position(:, 1), N_r.Position(:, 2), dotSizes_r, 'k', 'filled');
    xlim("padded");
    title(['Mode ', num2str(mode), ': Random bias']);
    hold on;
    for i = 1:length(seg_r)
        lim = seg_max_r{i};
        yline(lim, ':', 'LineWidth', 1);
        meanY2 = (seg_min_r{i} + seg_max_r{i}) / 2;
        offset = max(x_r) * 1.65;
        text(offset, meanY2, char(seg_r(i)), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
    end
    plot(max(x_r) * 1.2 + devec_means_right_r * 50, meanY, '-', 'LineWidth', 2, 'Color', r_col);
    plot(min(x_r) * 1.4 + devec_means_left_r * 50, meanY, '-', 'LineWidth', 2, 'Color', r_col);
    hold off;
    filename = [path, '/dist_d_evec', num2str(mode)];
    print(filename, '-dpng', '-r300');

    % --- Figure 2: Power spectrum ---
    figure;
    fig = gcf;
    set(fig, 'Units', 'pixels');
    set(fig, 'Position', [100 100 900 700]);
    power_data = struct();
    % Without Genetic Bias
    fft_r = fft(devec_means_right_wo);
    fft_l = fft(devec_means_left_wo);
    power_data.wo.r = abs(fft_r).^2;
    power_data.wo.l = abs(fft_l).^2;
    subplot(1, 3, 1);
    plot(power_data.wo.l, 'Color', [0.392 0.560 1], 'LineWidth', 1.4, 'LineStyle', '-.');
    hold on;
    plot(power_data.wo.r, 'Color', [1 0.690 0], 'LineWidth', 1.4, 'LineStyle', '--');
    xlabel('Frequency'); 
    ylabel('Power');
    title('Without genetic bias');
    hold off;

    % With Genetic Bias
    fft_r = fft(devec_means_right_w);
    fft_l = fft(devec_means_left_w);
    power_data.w.r = abs(fft_r).^2;
    power_data.w.l = abs(fft_l).^2;
    subplot(1, 3, 2);
    plot(power_data.w.l, 'Color', [0.392 0.560 1], 'LineWidth', 1.4, 'LineStyle', '-.');
    hold on;
    plot(power_data.w.r, 'Color', [1 0.690 0], 'LineWidth', 1.4, 'LineStyle', '--');
    xlabel('Frequency'); 
    ylabel('Power');
    title('With genetic bias');
    hold off;

    % Random
    fft_r = fft(devec_means_right_r);
    fft_l = fft(devec_means_left_r);
    power_data.r.r = abs(fft_r).^2;
    power_data.r.l = abs(fft_l).^2;
    subplot(1, 3, 3);
    plot(power_data.r.l, 'Color', [0.392 0.560 1], 'LineWidth', 1.4, 'LineStyle', '-.');
    hold on;
    plot(power_data.r.r, 'Color', [1 0.690 0], 'LineWidth', 1.4, 'LineStyle', '--');
    title('Randomised projectome');
    xlabel('Frequency'); 
    ylabel('Power');
    legend('Left hemimode', 'Right hemimode', 'FontSize', 14);
    legend('boxoff');
    sgt=sgtitle(['Mode ', num2str(mode)]); 
    sgt.FontWeight = 'bold';
    hold off;
    filename = [path, '/power_spectrum_hemi', num2str(mode)];
    print(filename, '-dpng', '-r300');

    % Model Comparison
    figure;
    fig = gcf;
    set(fig, 'Units', 'pixels');
    set(fig, 'Position', [100 100 900 700]);
    subplot(2, 2, 1);
    plot(power_data.wo.l/max(power_data.wo.l), 'Color', wo_col, 'LineWidth', 1.4, 'LineStyle', '-.');
    hold on;
    plot(power_data.w.l/max(power_data.w.l), 'Color', w_col, 'LineWidth', 1.4, 'LineStyle', '--');
    plot(power_data.r.l/max(power_data.r.l), 'Color', r_col, 'LineWidth', 1.4, 'LineStyle', ':');
    xlabel('Frequency'); 
    ylabel('Normalized power');
    title('Left hemimode power spectrum');
    legend('Without genetic bias', 'With genetic bias', 'Randomised projectome', 'Location','north', 'FontSize', 14);
    legend('boxoff');
    hold off;

    subplot(2, 2, 2);
    plot(power_data.wo.r/max(power_data.wo.r), 'Color', wo_col, 'LineWidth', 1.4, 'LineStyle', '-.');
    hold on;
    plot(power_data.w.r/max(power_data.w.r), 'Color', w_col, 'LineWidth', 1.4, 'LineStyle', '--');
    plot(power_data.r.r/max(power_data.r.r), 'Color', r_col, 'LineWidth', 1.4, 'LineStyle', ':');
    xlabel('Frequency'); 
    ylabel('Normalized power');
    title('Right hemimode power spectrum'); %divided by the max
    hold off;
    %print('power_spectrum_models', '-dpng', '-r300');

    % Model Comparison -- Log scale
    subplot(2, 2, 3);
    plot(power_data.wo.l, 'Color', wo_col, 'LineWidth', 1.4, 'LineStyle', '-.');
    hold on;
    plot(power_data.w.l, 'Color', w_col, 'LineWidth', 1.4, 'LineStyle', '--');
    plot(power_data.r.l, 'Color', r_col, 'LineWidth', 1.4, 'LineStyle', ':');
    xlabel('Frequency'); 
    ylabel('Log(Power)');
    title('Left hemimode power spectrum (log scale)');
    set(gca, 'YScale', 'log')
    hold off;

    subplot(2, 2, 4);
    plot(power_data.wo.r, 'Color', wo_col, 'LineWidth', 1.4, 'LineStyle', '-.');
    hold on;
    plot(power_data.w.r, 'Color', w_col, 'LineWidth', 1.4, 'LineStyle', '--');
    plot(power_data.r.r, 'Color', r_col, 'LineWidth', 1.4, 'LineStyle', ':');
    title('Right hemimode power spectrum (log scale)');
    set(gca, 'YScale', 'log');
    xlabel('Frequency'); 
    ylabel('Log(Power)');
    sgt=sgtitle(['Mode ', num2str(mode)]); 
    sgt.FontWeight = 'bold';
    hold off;
    %print('power_spectrum_models_log', '-dpng', '-r300');
    filename = [path, '/power_spectrum_models_both', num2str(mode)];
    print(filename, '-dpng', '-r300');
end
%% 
eval_plots(eval_wo, eval_w, eval_r, d_e_wo, d_e_w, d_e_r, wo_col, w_col, r_col);
plot_both_mode(d_evec_wo, d_evec_w, d_evec_r, d_evec_wo2, d_evec_w2, d_evec_r2, N_wo, N_w, N_r);
d_evec_distribution(1, d_evec_wo, d_evec_w, d_evec_r, N_wo, N_w, N_r);
d_evec_distribution(2, d_evec_wo2, d_evec_w2, d_evec_r2, N_wo, N_w, N_r);
%% 
%figure;
%scatter3(N_w.Position(:,1), N_w.Position(:,2), N_w.Position(:,3));