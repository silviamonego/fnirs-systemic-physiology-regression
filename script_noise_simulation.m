%% SCRIPT: PINK NOISE (1/f) SPECTRAL ANALYSIS & SPURIOUS CORRELATION SIMULATION
% aim: to demonstrate the 1/f power-law spectrum in low-frequency physiological
% noise and prove the necessity of block-by-block linear detrending

clearvars;
close all;
clc;

% =========================================================================
% Academic Context: Master's Thesis in Bioengineering for Neuroscience
% University of Padova (DEI) | NIRx Medical Technologies LLC
% Reference: Chapter 3 (Signal Processing and Data Analysis Pipeline, Sec. 3.4)
% =========================================================================

% Add external toolboxes if present in local or parent directories
external_dirs = {'cvxEDA', '../cvxEDA', 'Homer3_analysis', '../Homer3_analysis', ...
                 'homer2', '../homer2', fullfile('Homer3_analysis', 'cvxEDA-main', 'matlab')};
for d = 1:numel(external_dirs)
    if exist(external_dirs{d}, 'dir')
        addpath(genpath(external_dirs{d}));
    end
end

% Force the MATLAB Signal Processing Toolbox path to the beginning of the search path
% to resolve the shadowing of the built-in findpeaks function by Homer2's custom version.
spt_path = fullfile(matlabroot, 'toolbox', 'signal', 'signal');
if exist(spt_path, 'dir')
    addpath(spt_path, '-begin');
end

% Simulation parameters
noise_exponent = 1.7; % Exponent beta for 1/f^beta noise (1.0 = pink, 2.0 = red/brownian)

% =========================================================================
% DYNAMIC FILE PATH CONFIGURATION 
% =========================================================================
fprintf('Waiting for subject folder selection...\n');
default_data_dir = fullfile(pwd, 'DATASET');
if ~exist(default_data_dir, 'dir'), default_data_dir = pwd; end
subj_dir = uigetdir(default_data_dir, 'Select the Subject Folder');
if isequal(subj_dir, 0)
    disp('Operation canceled by the user');
    return;
end

[~, subject_name] = fileparts(subj_dir);

file_snirf     = fullfile(subj_dir, [subject_name, '_Satori2.snirf']);
file_snirf_raw = fullfile(subj_dir, [subject_name, '.snirf']);
file_csv       = fullfile(subj_dir, [subject_name, '_ProcessedPhysio.csv']);

if ~exist(file_snirf, 'file')
    error('Satori SNIRF file not found: %s', file_snirf);
end
if ~exist(file_snirf_raw, 'file')
    error('Raw SNIRF file not found: %s', file_snirf_raw);
end
if ~exist(file_csv, 'file')
    error('Processed physiological CSV file not found: %s', file_csv);
end

fprintf('\nSubject Directory: %s\n', subj_dir);
fprintf('Found Satori SNIRF File: %s\n', [subject_name, '_Satori2.snirf']);
fprintf('Found Raw SNIRF File: %s\n', [subject_name, '.snirf']);
fprintf('Found Processed Physio CSV File: %s\n', [subject_name, '_ProcessedPhysio.csv']);

% =========================================================================
%% PART 1: LOAD REAL DATA FOR TIMELINE AND REAL PSD ANALYSIS
% =========================================================================
fprintf('\n1. Loading fNIRS and physiological time series...\n');

% Load fNIRS timeline
time_vector = h5read(file_snirf, '/nirs/data1/time');
num_time_points = length(time_vector);
fs_fnirs = 1 / mean(diff(time_vector));

% Load Stimuli onsets
stim1 = h5read(file_snirf, '/nirs/stim1/data');
stim2 = h5read(file_snirf, '/nirs/stim2/data');
stim3 = h5read(file_snirf, '/nirs/stim3/data');
onsets_stim1 = stim1(1,:)'; % Fingertapping (FT)
onsets_stim2 = stim2(1,:)'; % Breath Holding (BH)
onsets_stim3 = stim3(1,:)'; % Combined (FT+BH)

% Read Labels and Indices to map HbO and HbR
data_time_series = h5read(file_snirf, '/nirs/data1/dataTimeSeries')'; 
info_snirf = h5info(file_snirf, '/nirs/data1');
num_elements = length(info_snirf.Groups);
data_labels = cell(num_elements, 1);
for i = 1:num_elements
    path_ml = info_snirf.Groups(i).Name; 
    raw_label = h5read(file_snirf, [path_ml, '/dataTypeLabel']);
    data_labels{i} = strtrim(char(raw_label));
end
idx_hbo = find(strcmp(data_labels, 'HbO')); 
num_channels = length(idx_hbo);

% Initialize HbO data matrix
y_hbo = zeros(num_time_points, num_channels);
for ch = 1:num_channels
    y_hbo(:, ch) = data_time_series(:, idx_hbo(ch));
end

% Dynamically identify Short Channels (< 15 mm) from 3D probe geometry
% In the 54-channel motor montage used in this study, these correspond to channels [3, 14, 18, 28, 31, 42, 45, 54]
src_pos = h5read(file_snirf_raw, '/nirs/probe/sourcePos3D')';
det_pos = h5read(file_snirf_raw, '/nirs/probe/detectorPos3D')';
ml_src = zeros(num_channels, 1);
ml_det = zeros(num_channels, 1);
sd_dist = zeros(num_channels, 1);

for ch = 1:num_channels
    path_ch = sprintf('/nirs/data1/measurementList%d', ch);
    ml_src(ch) = double(h5read(file_snirf, [path_ch, '/sourceIndex']));
    ml_det(ch) = double(h5read(file_snirf, [path_ch, '/detectorIndex']));
    sd_dist(ch) = norm(src_pos(ml_src(ch), :) - det_pos(ml_det(ch), :));
end

short_channels = find(sd_dist < 15)';

if isempty(short_channels)
    error('No short channels (< 15 mm) identified in the probe geometry. Please verify SNIRF 3D optode coordinates.');
end

fprintf('   -> Identified %d Short Channels (< 15 mm): [%s]\n', length(short_channels), num2str(short_channels));

sc_data_hbo = y_hbo(:, short_channels);
R_sc_hbo = corr(zscore(sc_data_hbo), 'Rows', 'complete');
mean_corr_hbo = zeros(length(short_channels), 1);
for ch = 1:length(short_channels)
    other_idx = setdiff(1:length(short_channels), ch);
    mean_corr_hbo(ch) = mean(R_sc_hbo(ch, other_idx), 'omitnan');
end
bad_sc_hbo = mean_corr_hbo < 0.1;
short_channels_clean_hbo = short_channels(~bad_sc_hbo);
sc_data_clean_hbo = y_hbo(:, short_channels_clean_hbo);
[~, score_sc_hbo] = pca(zscore(sc_data_clean_hbo));
sc_pca1_hbo = zscore(score_sc_hbo(:, 1));

% Load Physio
physio_table = readtable(file_csv); 
time_physio = physio_table.Time_s;
hr_ecg_raw = physio_table.HR_ECG;
hr_ecg_aligned = interp1(time_physio, hr_ecg_raw, time_vector, 'linear', 'extrap');

% =========================================================================
%% PART 2: DEFINE EXPERIMENTAL SEGMENTS
% =========================================================================
fprintf('2. Segmenting timeline based on experimental paradigm...\n');
% to create the temporal masks for the 4 different tasks (Baseline,
% Fingertapping, Breath holding, and Combined)
% these masks will be applied to the synthetic noise to replicate the exact
% concatenation process

% Baseline: 30s before first stimulus onset
all_onsets = [onsets_stim1; onsets_stim2; onsets_stim3];
first_onset = min(all_onsets);
idx_baseline = time_vector >= (first_onset - 30) & time_vector < first_onset;

% Tasks (21s blocks)
task_duration = 21;

% Fingertapping (FT)
idx_ft = false(num_time_points, 1);
for i = 1:length(onsets_stim1)
    onset_t = onsets_stim1(i);
    idx_block = time_vector >= onset_t & time_vector <= (onset_t + task_duration);
    idx_ft = idx_ft | idx_block;
end

% Breath Holding (BH)
idx_bh = false(num_time_points, 1);
for i = 1:length(onsets_stim2)
    onset_t = onsets_stim2(i);
    idx_block = time_vector >= onset_t & time_vector <= (onset_t + task_duration);
    idx_bh = idx_bh | idx_block;
end

% Combined (COMB)
idx_comb = false(num_time_points, 1);
for i = 1:length(onsets_stim3)
    onset_t = onsets_stim3(i);
    idx_block = time_vector >= onset_t & time_vector <= (onset_t + task_duration);
    idx_comb = idx_comb | idx_block;
end

% =========================================================================
%% PART 3: POWER SPECTRAL DENSITY (PSD) ANALYSIS
% =========================================================================
fprintf('3. Performing Spectral Analysis (PSD) and replication of 1/f slope...\n');
% to visually prove that SC blood flow follows the exact same 1/f decay
% slope as the synthetic noise, justifying why we use a B=1.7 model to
% simulate fNIRS noise. It also labels the physiological peaks (cardiac and
% respiratory). 

% Rationale:
% Real-world physiological noise is not "white noise" (flat across all
% frequencies), but it's "colored noise", meaning power is heavily
% concentrated in the low-frequency range

% Compute PSD using periodogram (Signal Processing Toolbox)
[psd_hr, f_hr]   = periodogram(zscore(hr_ecg_aligned), hamming(num_time_points), [], fs_fnirs);
[psd_hbo, f_hbo] = periodogram(sc_pca1_hbo, hamming(num_time_points), [], fs_fnirs);

% Pre-compute pink noise filter amplitude spectral envelope
num_unique_pts = ceil((num_time_points+1)/2);
f_bins = (1:num_unique_pts)';
filter_amp = 1 ./ (f_bins .^ (noise_exponent / 2));

% Generate colored noise for spectral comparison (inline generation)
w_synth = randn(num_time_points, 1);
W_synth = fft(w_synth); % FFT filter
W_shaped_synth = W_synth;
W_shaped_synth(1:num_unique_pts) = W_synth(1:num_unique_pts) .* filter_amp;
if mod(num_time_points, 2) == 0
    W_shaped_synth(num_unique_pts+1:end) = conj(W_shaped_synth(num_unique_pts-1:-1:2));
else
    W_shaped_synth(num_unique_pts+1:end) = conj(W_shaped_synth(num_unique_pts:-1:2));
end
W_shaped_synth(1) = 0; % Remove DC
synthetic_noise = real(ifft(W_shaped_synth));
synthetic_noise = zscore(synthetic_noise);

[psd_pink, f_pink] = periodogram(synthetic_noise, hamming(num_time_points), [], fs_fnirs);

% Plot log-log Power Spectral Densities
figure('Name', 'Power Spectral Density (PSD) Analysis', 'Color', 'w', 'Position', [100, 100, 1100, 480]);

% Left subplot: HR (ECG) vs Shaped Noise
subplot(1, 2, 1);
h_pink1 = loglog(f_pink, psd_pink, 'Color', [0.7 0.7 0.7], 'LineWidth', 1.5);
hold on;
h_hr = loglog(f_hr, psd_hr, 'Color', [0.85 0.25 0.25], 'LineWidth', 2);
grid on;
xlim([0.01, fs_fnirs/2]);
xlabel('Frequency (Hz)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Power Density (V^2/Hz)', 'FontSize', 11, 'FontWeight', 'bold');
title(sprintf('PSD of HR (ECG) vs. 1/f^{%.1f} Noise', noise_exponent), 'FontSize', 12, 'FontWeight', 'bold');
legend([h_pink1, h_hr], {sprintf('Shaped Noise (1/f^{%.1f})', noise_exponent), 'HR (ECG)'}, 'Location', 'southwest');

% Right subplot: fNIRS SC PCA HbO vs Shaped Noise
subplot(1, 2, 2);
h_pink2 = loglog(f_pink, psd_pink, 'Color', [0.7 0.7 0.7], 'LineWidth', 1.5);
hold on;
h_hbo = loglog(f_hbo, psd_hbo, 'Color', [0.2 0.45 0.75], 'LineWidth', 2);
grid on;
xlim([0.01, fs_fnirs/2]);
xlabel('Frequency (Hz)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Power Density (V^2/Hz)', 'FontSize', 11, 'FontWeight', 'bold');
title(sprintf('PSD of SC PCA HbO vs. 1/f^{%.1f} Noise', noise_exponent), 'FontSize', 12, 'FontWeight', 'bold');
legend([h_pink2, h_hbo], {sprintf('Shaped Noise (1/f^{%.1f})', noise_exponent), 'SC PCA HbO'}, 'Location', 'southwest');

% 1. Rileva il VERO picco cardiaco direttamente sullo spettro fNIRS (0.8 - 1.8 Hz)
cardiac_range_idx = find(f_hbo >= 0.8 & f_hbo <= 1.8);
if ~isempty(cardiac_range_idx)
    [~, max_idx_card] = max(psd_hbo(cardiac_range_idx));
    cardiac_peak_freq = f_hbo(cardiac_range_idx(max_idx_card));
    xline(cardiac_peak_freq, 'r--', sprintf('Cardiac Peak (%.2f Hz)', cardiac_peak_freq), ...
          'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'center', ...
          'LineWidth', 1.3, 'FontSize', 9.5, 'HandleVisibility', 'off');
end

% 2. Rileva la frequenza respiratoria reale dalla fascia (Resp_Clean)
if ismember('Resp_Clean', physio_table.Properties.VariableNames)
    fs_physio = 1 / mean(diff(physio_table.Time_s));
    [psd_resp, f_resp] = periodogram(zscore(physio_table.Resp_Clean), hamming(height(physio_table)), [], fs_physio);
    resp_idx = find(f_resp >= 0.15 & f_resp <= 0.40);
    if ~isempty(resp_idx)
        [~, max_idx_resp] = max(psd_resp(resp_idx));
        resp_peak_freq = f_resp(resp_idx(max_idx_resp));
        xline(resp_peak_freq, 'g--', sprintf('Respiratory (%.2f Hz)', resp_peak_freq), ...
              'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'center', ...
              'LineWidth', 1.3, 'FontSize', 9.5, 'HandleVisibility', 'off');
    end
end

sgtitle(sprintf('Replication of 1/f^{%.1f} Power Law Spectrum - %s', noise_exponent, subject_name), ...
        'FontSize', 13, 'FontWeight', 'bold', 'Interpreter', 'none');


% =========================================================================
%% PART 4: MONTE CARLO SIMULATION (SPURIOUS CORRELATION)
% =========================================================================
fprintf('4. Running Monte Carlo simulation of independent shaped noise...\n');
% Rationale:
% because 1/f noise is dominated by slow, long-term drifts, 2 independent
% signals will often slide up or down together purely by chance over short
% windows -> high spurious correlations (often R > 0.4 in raw data)
% Monte Carlo simulation calculates the exact probability distribution of
% this false correlation (to not trust any raw correlation below R=0.35)

num_iterations = 1000;

% Raw distributions
corr_dist_baseline = zeros(num_iterations, 1);
corr_dist_ft       = zeros(num_iterations, 1);
corr_dist_bh       = zeros(num_iterations, 1);
corr_dist_comb     = zeros(num_iterations, 1);

% Detrended distributions
corr_dist_baseline_det = zeros(num_iterations, 1);
corr_dist_ft_det       = zeros(num_iterations, 1);
corr_dist_bh_det       = zeros(num_iterations, 1);
corr_dist_comb_det     = zeros(num_iterations, 1);

% Single isolated block distributions
idx_single_task = time_vector >= onsets_stim1(1) & time_vector <= (onsets_stim1(1) + task_duration);
corr_dist_single_task = zeros(num_iterations, 1);
corr_dist_single_task_det = zeros(num_iterations, 1);

for iter = 1:num_iterations
    % Generate two independent shaped noise signals in a 2-column matrix
    % (p1 and p2)
    w_sim = randn(num_time_points, 2); % random white noise 
    W_sim = fft(w_sim); % in the frequency domain
    
    W_shaped_sim = W_sim;
    W_shaped_sim(1:num_unique_pts, :) = W_sim(1:num_unique_pts, :) .* filter_amp;
    
    if mod(num_time_points, 2) == 0
        W_shaped_sim(num_unique_pts+1:end, :) = conj(W_shaped_sim(num_unique_pts-1:-1:2, :));
    else
        W_shaped_sim(num_unique_pts+1:end, :) = conj(W_shaped_sim(num_unique_pts:-1:2, :));
    end
    W_shaped_sim(1, :) = 0; % Remove DC component
    
    p_sim = real(ifft(W_shaped_sim));
    p1 = zscore(p_sim(:, 1));
    p2 = zscore(p_sim(:, 2));
    
    % --- 1. RAW CORRELATIONS ---
    corr_dist_baseline(iter)    = abs(corr(p1(idx_baseline), p2(idx_baseline), 'Rows', 'complete'));
    corr_dist_ft(iter)          = abs(corr(p1(idx_ft),       p2(idx_ft),       'Rows', 'complete'));
    corr_dist_bh(iter)          = abs(corr(p1(idx_bh),       p2(idx_bh),       'Rows', 'complete'));
    corr_dist_comb(iter)        = abs(corr(p1(idx_comb),     p2(idx_comb),     'Rows', 'complete'));
    corr_dist_single_task(iter) = abs(corr(p1(idx_single_task), p2(idx_single_task), 'Rows', 'complete'));
    % because the signals are independent, the true correlation should be 0
    
    % --- 2. DETRENDED CORRELATIONS ---
    corr_dist_baseline_det(iter)    = abs(corr(detrend(p1(idx_baseline)),    detrend(p2(idx_baseline)),    'Rows', 'complete'));
    corr_dist_single_task_det(iter) = abs(corr(detrend(p1(idx_single_task)), detrend(p2(idx_single_task)), 'Rows', 'complete'));
    % it detrends the segments block-by-block and calculates the correlation
    
    % FT blocks detrended individually
    p1_ft_det = []; p2_ft_det = [];
    for i = 1:length(onsets_stim1)
        idx = time_vector >= onsets_stim1(i) & time_vector <= (onsets_stim1(i) + task_duration);
        p1_ft_det = [p1_ft_det; detrend(p1(idx))];
        p2_ft_det = [p2_ft_det; detrend(p2(idx))];
    end
    corr_dist_ft_det(iter) = abs(corr(p1_ft_det, p2_ft_det, 'Rows', 'complete'));
    
    % BH blocks detrended individually
    p1_bh_det = []; p2_bh_det = [];
    for i = 1:length(onsets_stim2)
        idx = time_vector >= onsets_stim2(i) & time_vector <= (onsets_stim2(i) + task_duration);
        p1_bh_det = [p1_bh_det; detrend(p1(idx))];
        p2_bh_det = [p2_bh_det; detrend(p2(idx))];
    end
    corr_dist_bh_det(iter) = abs(corr(p1_bh_det, p2_bh_det, 'Rows', 'complete'));
    
    % COMB blocks detrended individually
    p1_comb_det = []; p2_comb_det = [];
    for i = 1:length(onsets_stim3)
        idx = time_vector >= onsets_stim3(i) & time_vector <= (onsets_stim3(i) + task_duration);
        p1_comb_det = [p1_comb_det; detrend(p1(idx))];
        p2_comb_det = [p2_comb_det; detrend(p2(idx))];
    end
    corr_dist_comb_det(iter) = abs(corr(p1_comb_det, p2_comb_det, 'Rows', 'complete'));
end

% Compute statistical means
mean_baseline     = mean(corr_dist_baseline);
mean_ft           = mean(corr_dist_ft);
mean_bh           = mean(corr_dist_bh);
mean_comb         = mean(corr_dist_comb);

mean_baseline_det = mean(corr_dist_baseline_det);
mean_ft_det       = mean(corr_dist_ft_det);
mean_bh_det       = mean(corr_dist_bh_det);
mean_comb_det     = mean(corr_dist_comb_det);

mean_single_task     = mean(corr_dist_single_task);
mean_single_task_det = mean(corr_dist_single_task_det);

fprintf('\n--- Monte Carlo Spurious Correlation Results (Mean |R|): ---\n');
fprintf('   Baseline (30s window):            Raw = %.3f | Detrended = %.3f\n', mean_baseline, mean_baseline_det);
fprintf('   Fingertapping (210s concatenated): Raw = %.3f | Detrended = %.3f\n', mean_ft, mean_ft_det);
fprintf('   Breath Holding (105s concatenated):Raw = %.3f | Detrended = %.3f\n', mean_bh, mean_bh_det);
fprintf('   Combined (105s concatenated):     Raw = %.3f | Detrended = %.3f\n', mean_comb, mean_comb_det);
fprintf('   Single Task Block (21s isolated):  Raw = %.3f | Detrended = %.3f\n', mean_single_task, mean_single_task_det);
fprintf('   (Using 1/f^%.1f noise model)\n', noise_exponent);
fprintf('------------------------------------------------------------\n');

% -------------------------------------------------------------------------
% FIGURE 2: Monte Carlo Distributions - RAW NOISE
% -------------------------------------------------------------------------
figure('Name', 'Spurious Correlation (Raw Noise)', 'Color', 'w', 'Position', [100, 100, 1000, 750]);
% Histograms plots of the 1000 correlation values to show the probability
% density of the spurious correlations

condition_dists_raw = {corr_dist_baseline, corr_dist_ft, corr_dist_bh, corr_dist_comb};
condition_names_raw = {'1. Baseline (30s Raw)', '2. Fingertapping (210s Concatenated Raw)', ...
                       '3. Breath Holding (105s Concatenated Raw)', '4. Combined (105s Concatenated Raw)'};
colors_raw = {[0.2, 0.45, 0.75], [0.85, 0.25, 0.25], [0.466, 0.674, 0.188], [0.929, 0.694, 0.125]};
means_raw = [mean_baseline, mean_ft, mean_bh, mean_comb];

for cond = 1:4
    subplot(2, 2, cond);
    h = histogram(condition_dists_raw{cond}, 'Normalization', 'pdf', 'FaceColor', colors_raw{cond}, 'EdgeColor', 'k', 'FaceAlpha', 0.6);
    hold on;
    xline(means_raw(cond), 'r--', sprintf('Spurious Baseline: %.3f', means_raw(cond)), 'LineWidth', 2.2, 'LabelVerticalAlignment', 'top', 'HandleVisibility', 'off');
    xlim([0 1]);
    xlabel('Absolute Correlation Magnitude |R|', 'FontSize', 10, 'FontWeight', 'bold');
    ylabel('Probability Density', 'FontSize', 10, 'FontWeight', 'bold');
    title(condition_names_raw{cond}, 'FontSize', 12, 'FontWeight', 'bold');
    legend(h, sprintf('Raw Noise, Mean: %.3f', means_raw(cond)), 'Location', 'northeast');
    grid on;
end

sgtitle(sprintf('Spurious Correlation Thresholds of Independent Raw 1/f^{%.1f} Noise - %s', noise_exponent, subject_name), ...
        'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'none');

% -------------------------------------------------------------------------
% FIGURE 3: Monte Carlo Distributions - DETRENDED NOISE
% -------------------------------------------------------------------------
figure('Name', 'Spurious Correlation (Detrended Noise)', 'Color', 'w', 'Position', [150, 150, 1000, 750]);

condition_dists_det = {corr_dist_baseline_det, corr_dist_ft_det, corr_dist_bh_det, corr_dist_comb_det};
condition_names_det = {'1. Baseline (30s Detrended)', '2. Fingertapping (210s Concatenated Detrended)', ...
                       '3. Breath Holding (105s Concatenated Detrended)', '4. Combined (105s Concatenated Detrended)'};
colors_det = {[0.2, 0.45, 0.75], [0.85, 0.25, 0.25], [0.466, 0.674, 0.188], [0.929, 0.694, 0.125]};
means_det = [mean_baseline_det, mean_ft_det, mean_bh_det, mean_comb_det];

for cond = 1:4
    subplot(2, 2, cond);
    h = histogram(condition_dists_det{cond}, 'Normalization', 'pdf', 'FaceColor', colors_det{cond}, 'EdgeColor', 'k', 'FaceAlpha', 0.6);
    hold on;
    xline(means_det(cond), 'r--', sprintf('Spurious Baseline: %.3f', means_det(cond)), 'LineWidth', 2.2, 'LabelVerticalAlignment', 'top', 'HandleVisibility', 'off');
    xlim([0 1]);
    xlabel('Absolute Correlation Magnitude |R|', 'FontSize', 10, 'FontWeight', 'bold');
    ylabel('Probability Density', 'FontSize', 10, 'FontWeight', 'bold');
    title(condition_names_det{cond}, 'FontSize', 12, 'FontWeight', 'bold');
    legend(h, sprintf('Detrended Noise, Mean: %.3f', means_det(cond)), 'Location', 'northeast');
    grid on;
end

sgtitle(sprintf('Spurious Correlation Thresholds of Independent Detrended 1/f^{%.1f} Noise - %s', noise_exponent, subject_name), ...
        'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'none');

% =========================================================================
%% PART 5: METHODOLOGICAL COMPARISON (DETRENDING VS. RAW NOISE)
% =========================================================================
fprintf('5. Plotting detrending vs. raw noise comparison...\n');
% To see the differences between correlation distributions of:
% - an isolated 21s single block
% - concatenated raw blocks (210s)
% - concatenated detrended blocks (210s, detrended block-by-block)

figure('Name', 'Methodological Comparison (Detrending vs. Raw Noise)', ...
       'Color', 'w', 'Position', [200, 200, 1100, 480]);

% Left subplot: Fingertapping
subplot(1, 2, 1);
h1 = histogram(corr_dist_single_task, 'Normalization', 'pdf', 'FaceColor', [0.85, 0.25, 0.25], 'EdgeColor', 'k', 'FaceAlpha', 0.3);
hold on;
h2 = histogram(corr_dist_ft, 'Normalization', 'pdf', 'FaceColor', [0.2, 0.45, 0.75], 'EdgeColor', 'k', 'FaceAlpha', 0.3);
h3 = histogram(corr_dist_ft_det, 'Normalization', 'pdf', 'FaceColor', [0.466, 0.674, 0.188], 'EdgeColor', 'k', 'FaceAlpha', 0.6);

xline(mean_single_task, 'Color', [0.85, 0.25, 0.25], 'LineStyle', '--', 'LineWidth', 1.8, 'HandleVisibility', 'off');
xline(mean_ft, 'Color', [0.2, 0.45, 0.75], 'LineStyle', '--', 'LineWidth', 1.8, 'HandleVisibility', 'off');
xline(mean_ft_det, 'Color', [0.466, 0.674, 0.188], 'LineStyle', '--', 'LineWidth', 2.2, 'HandleVisibility', 'off');
xlim([0 1]);
xlabel('Absolute Correlation Magnitude |R|', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Probability Density', 'FontSize', 11, 'FontWeight', 'bold');
title('Task Blocks: Single Raw vs. Concatenated Raw vs. Detrended', 'FontSize', 12, 'FontWeight', 'bold');
legend([h1, h2, h3], {sprintf('Single Block Raw (21s), Mean: %.3f', mean_single_task), ...
                      sprintf('Concatenated Raw (210s), Mean: %.3f', mean_ft), ...
                      sprintf('Concatenated Detrended (210s), Mean: %.3f', mean_ft_det)}, ...
       'Location', 'northeast', 'FontSize', 10);
grid on;

% Right subplot: Baseline 30s
subplot(1, 2, 2);
h4 = histogram(corr_dist_baseline, 'Normalization', 'pdf', 'FaceColor', [0.2, 0.45, 0.75], 'EdgeColor', 'k', 'FaceAlpha', 0.4);
hold on;
h5 = histogram(corr_dist_baseline_det, 'Normalization', 'pdf', 'FaceColor', [0.466, 0.674, 0.188], 'EdgeColor', 'k', 'FaceAlpha', 0.6);
xline(mean_baseline, 'Color', [0.2, 0.45, 0.75], 'LineStyle', '--', 'LineWidth', 1.8, 'HandleVisibility', 'off');
xline(mean_baseline_det, 'Color', [0.466, 0.674, 0.188], 'LineStyle', '--', 'LineWidth', 2.2, 'HandleVisibility', 'off');
xlim([0 1]);
xlabel('Absolute Correlation Magnitude |R|', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Probability Density', 'FontSize', 11, 'FontWeight', 'bold');
title('Baseline Block: Raw vs. Detrended', 'FontSize', 12, 'FontWeight', 'bold');
legend([h4, h5], {sprintf('Baseline Raw (30s), Mean: %.3f', mean_baseline), ...
                  sprintf('Baseline Detrended (30s), Mean: %.3f', mean_baseline_det)}, ...
       'Location', 'northeast', 'FontSize', 10);
grid on;

sgtitle(sprintf('Detrending & Block Length vs Spurious Correlation (1/f^{%.1f} Noise) - %s', noise_exponent, subject_name), ...
        'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'none');

% Auto-save Noise Simulation Figures (300 DPI PNG)
exportgraphics(findobj('Name', 'Power Spectral Density (PSD) Analysis'), fullfile(subj_dir, [subject_name, '_PSD_Analysis.png']), 'Resolution', 300);
exportgraphics(findobj('Name', 'Spurious Correlation (Raw Noise)'), fullfile(subj_dir, [subject_name, '_Spurious_Correlation_Raw.png']), 'Resolution', 300);
exportgraphics(findobj('Name', 'Spurious Correlation (Detrended Noise)'), fullfile(subj_dir, [subject_name, '_Spurious_Correlation_Detrended.png']), 'Resolution', 300);
exportgraphics(findobj('Name', 'Methodological Comparison (Detrending vs. Raw Noise)'), fullfile(subj_dir, [subject_name, '_Noise_Simulation_Comparison.png']), 'Resolution', 300);
fprintf('   -> Saved Noise Simulation plots (PNG 300 DPI).\n');

% Results:
% Concatenating raw blocks keeps the spurious correlation high (R=0.3)
% because the slow drift span across blocks 
% Linear detrending successfully removes the slow drift in each task
% window, shifting the entire correlation distribution back toward zero (R<0.1)

% This statistically proves that block-by-block detrending is mandatory
% before performing task-based correlation analyses to avoid reporting FP
% physiological couplings
