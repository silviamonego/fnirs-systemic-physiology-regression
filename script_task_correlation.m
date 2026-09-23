%% TASK-BASED CORRELATION ANALYSIS
% Replicating and extending the paradigm from Alexander von Lühmann
% aim: to evaluate how the correlation between different physiological
% variables and superficial scalp blood flow (short channels) changes
% dynamically across different tasks

clearvars;
close all;
clc;

% =========================================================================
% Academic Context: Master's Thesis in Bioengineering for Neuroscience
% University of Padova (DEI) | NIRx Medical Technologies LLC
% Reference: Chapter 3 (Sec. 3.5) & Chapter 4 (Sec. 4.1.2)
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
%% PART 1: LOAD fNIRS DATA AND VALID TRIALS
% =========================================================================
fprintf('\n1. Loading fNIRS data from Satori SNIRF...\n');

% Load time series and time vector
data_time_series = h5read(file_snirf, '/nirs/data1/dataTimeSeries')'; % [Time x Channels]
time_vector = h5read(file_snirf, '/nirs/data1/time');
num_time_points = length(time_vector);

% fNIRS Sampling frequency
fs_fnirs = 1 / mean(diff(time_vector));
fprintf('   -> fNIRS Sampling Rate: %.2f Hz\n', fs_fnirs);

% Load Geometry from RAW file
src_pos = h5read(file_snirf_raw, '/nirs/probe/sourcePos3D')';
det_pos = h5read(file_snirf_raw, '/nirs/probe/detectorPos3D')';

% Load ValidTrials.mat (processed by physio_general.m) containing validated triggers and tInc
% to ensure that non-compliant trials are completely ignored
file_valid_trials = fullfile(subj_dir, [subject_name, '_ValidTrials.mat']);
if ~exist(file_valid_trials, 'file')
    error('ValidTrials file not found: %s\nPlease run physio_general.m first to perform trial validation.', file_valid_trials);
end

load(file_valid_trials);
onsets_stim1 = valid_onsets_FT;
onsets_stim2 = valid_onsets_BH;
onsets_stim3 = valid_onsets_COMB;

fprintf('   -> Loaded validated onsets (FT=%d, BH=%d, COMB=%d).\n', ...
    length(onsets_stim1), length(onsets_stim2), length(onsets_stim3));

% Determine number of channels dynamically: Satori exports HbO (1..N) and HbR (N+1..2N)
% (In the study's reference montage: 54 channels, 1..54 = HbO, 55..108 = HbR)
total_columns = size(data_time_series, 2);
num_channels = total_columns / 2;

ml_src = zeros(num_channels, 1);
ml_det = zeros(num_channels, 1);

% Extract Source-Detector couples in numerical order 1..num_channels
for ch = 1:num_channels
    path_hbo = sprintf('/nirs/data1/measurementList%d', ch);
    ml_src(ch) = double(h5read(file_snirf, [path_hbo, '/sourceIndex']));
    ml_det(ch) = double(h5read(file_snirf, [path_hbo, '/detectorIndex']));
end

scale_factor_satori = 10.0; % Corrects Satori mm/cm MBLL extinction factor

% concentration matrix y: [Time x 3 x Channels] in micromolar (uM)
% column 1..num_channels = HbO | column (num_channels+1)..2*num_channels = HbR
y = zeros(num_time_points, 3, num_channels);
for ch = 1:num_channels
    y(:, 1, ch) = data_time_series(:, ch)/scale_factor_satori; % HbO (uM)
    y(:, 2, ch) = data_time_series(:, ch + num_channels)/scale_factor_satori; % HbR (uM)
    y(:, 3, ch) = y(:, 1, ch) + y(:, 2, ch); % HbT
end


% =========================================================================
%% PART 2: CALCULATE SHORT-CHANNELS PCA (SC PCA)
% =========================================================================
fprintf('2. Identifying Short Channels (< 15 mm) and calculating SC PCA...\n');

% Dynamic identification of Short Channels (distance < 15 mm)
sd_dist = zeros(num_channels, 1);
for ch = 1:num_channels
    src_pos_val = src_pos(ml_src(ch), :);
    det_pos_val = det_pos(ml_det(ch), :);
    sd_dist(ch) = norm(src_pos_val - det_pos_val);
end

short_channels = find(sd_dist < 15)'; % 8 Short Channels (< 15 mm)

if isempty(short_channels)
    error('No short channel (< 15 mm) identified in the probe geometry!');
end

fprintf('   -> Identified %d Short Channels: [%s]\n', length(short_channels), num2str(short_channels));

% --- PCA for HbO ------------------------------------
sc_data_hbo = squeeze(y(:, 1, short_channels));
R_sc_hbo = corr(sc_data_hbo, 'Rows', 'complete');

% Discard short channels with mean correlation < 0.1
mean_corr_hbo = zeros(length(short_channels), 1);
for ch = 1:length(short_channels)
    other_idx = setdiff(1:length(short_channels), ch);
    mean_corr_hbo(ch) = mean(R_sc_hbo(ch, other_idx), 'omitnan');
end
short_channels_clean_hbo = short_channels(mean_corr_hbo >= 0.1);

% Run PCA on retained clean HbO short channels
[~, score_sc_hbo, ~, ~, explained_sc_hbo] = pca(zscore(squeeze(y(:, 1, short_channels_clean_hbo))));
sc_pca1_hbo = zscore(score_sc_hbo(:, 1));

% Ensure physiological sign alignment (PC1 must correlate positively with scalp channels)
if mean(corr(sc_pca1_hbo, squeeze(y(:, 1, short_channels_clean_hbo)))) < 0
    sc_pca1_hbo = -sc_pca1_hbo;
    fprintf('   -> Inverted HbO SC PCA sign for physiological alignment (+corr with scalp channels).\n');
end

fprintf('   -> HbO Short channels used for PCA: [%s]\n', num2str(short_channels_clean_hbo));
fprintf('   -> HbO SC PCA PC1 explained variance: %.1f%%\n', explained_sc_hbo(1));

% --- PCA for HbR ------------------------------------
sc_data_hbr = squeeze(y(:, 2, short_channels));
R_sc_hbr = corr(sc_data_hbr, 'Rows', 'complete');

% Discard short channels with mean correlation < 0.1
mean_corr_hbr = zeros(length(short_channels), 1);
for ch = 1:length(short_channels)
    other_idx = setdiff(1:length(short_channels), ch);
    mean_corr_hbr(ch) = mean(R_sc_hbr(ch, other_idx), 'omitnan');
end
short_channels_clean_hbr = short_channels(mean_corr_hbr >= 0.1);

% Run PCA on retained clean HbR short channels
[~, score_sc_hbr, ~, ~, explained_sc_hbr] = pca(zscore(squeeze(y(:, 2, short_channels_clean_hbr))));
sc_pca1_hbr = zscore(score_sc_hbr(:, 1));

% Ensure physiological sign alignment (PC1 must correlate positively with scalp channels)
if mean(corr(sc_pca1_hbr, squeeze(y(:, 2, short_channels_clean_hbr)))) < 0
    sc_pca1_hbr = -sc_pca1_hbr;
    fprintf('   -> Inverted HbR SC PCA sign for physiological alignment (+corr with scalp channels).\n');
end

fprintf('   -> HbR Short channels used for PCA: [%s]\n', num2str(short_channels_clean_hbr));
fprintf('   -> HbR SC PCA PC1 explained variance: %.1f%%\n', explained_sc_hbr(1));


% =========================================================================
%% PART 3: LOAD AND ALIGN PHYSIOLOGICAL DATA
% =========================================================================
fprintf('3. Loading and aligning physiological data...\n');
physio_table = readtable(file_csv); 

time_physio = physio_table.Time_s;
resp_raw = physio_table.Resp_Clean;
eda_tonic_raw  = physio_table.EDA_Tonic;
eda_phasic_raw = physio_table.EDA_Phasic;
hr_ppg_raw     = physio_table.HR_PPG;
hr_ecg_raw     = physio_table.HR_ECG;
spo2_raw       = physio_table.SpO2_Clean;

% Interpolate to fNIRS timeline
resp_aligned = interp1(time_physio, resp_raw, time_vector, 'linear', 'extrap');
eda_tonic_aligned = interp1(time_physio, eda_tonic_raw, time_vector, 'linear', 'extrap');
eda_phasic_aligned = interp1(time_physio, eda_phasic_raw, time_vector, 'linear', 'extrap');
hr_ppg_aligned = interp1(time_physio, hr_ppg_raw, time_vector, 'linear', 'extrap');
hr_ecg_aligned = interp1(time_physio, hr_ecg_raw, time_vector, 'linear', 'extrap');
spo2_aligned = interp1(time_physio, spo2_raw, time_vector, 'linear', 'extrap');

% =========================================================================
%% PART 4: DEFINE ANALYSIS REGIONS AND STATE SEGMENTS
% =========================================================================
fprintf('4. Segmenting time series into Baseline and Task periods...\n');
% Rationale:
% to split the entire experimental timeline into the 4 distinct block design
% segments (Baseline, Fingertapping, Breath holding, Combined), it creates
% logical masks and concatenates the time points belonging to each state 

% 1. Baseline: 30s immediately preceding the first stimulus onset of the session
all_onsets = [];
if ~isempty(onsets_stim1), all_onsets = [all_onsets; onsets_stim1]; end
if ~isempty(onsets_stim2), all_onsets = [all_onsets; onsets_stim2]; end
if ~isempty(onsets_stim3), all_onsets = [all_onsets; onsets_stim3]; end
first_onset = min(all_onsets);

idx_baseline = time_vector >= (first_onset - 30) & time_vector < first_onset;

% 2. Task Periods: extract exact task duration window [onset, onset + 21s]
task_duration = 21; % seconds

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

% Diagnostics print
fprintf('   -> Baseline points: %d samples (~%.1f seconds)\n', sum(idx_baseline), sum(idx_baseline)/fs_fnirs);
fprintf('   -> Fingertapping (FT) points: %d samples (~%.1f seconds)\n', sum(idx_ft), sum(idx_ft)/fs_fnirs);
fprintf('   -> Breath Holding (BH) points: %d samples (~%.1f seconds)\n', sum(idx_bh), sum(idx_bh)/fs_fnirs);
fprintf('   -> Combined (COMB) points: %d samples (~%.1f seconds)\n', sum(idx_comb), sum(idx_comb)/fs_fnirs);

% =========================================================================
%% PART 5: COMPUTE TASK-BASED CORRELATION MATRICES
% =========================================================================
fprintf('5. Standardising signals and computing correlations...\n');

% Assemble the raw, aligned regressor matrix (ensure column vectors)
raw_regressors = [ ...
    hr_ecg_aligned(:), ...
    hr_ppg_aligned(:), ...
    resp_aligned(:), ...
    eda_tonic_aligned(:), ...
    eda_phasic_aligned(:), ...
    spo2_aligned(:), ...
    sc_pca1_hbo(:), ...
    sc_pca1_hbr(:) ...
];

% Apply Global Standardisation
regressor_matrix = zscore(raw_regressors); % each regressors with mean=0 and std=1

regressor_names = { ...
    'HR (ECG)', ...
    'HR (PPG)', ...
    'Resp', ...
    'EDA Tonic', ...
    'EDA Phasic', ...
    'SpO2', ...
    'SC PCA HbO', ...
    'SC PCA HbR' ...
};

% --- METHOD 1: RAW CONCATENATED BLOCKS ---
R_baseline_raw = corr(regressor_matrix(idx_baseline, :), 'Rows', 'complete');
R_ft_raw       = corr(regressor_matrix(idx_ft, :),       'Rows', 'complete');
R_bh_raw       = corr(regressor_matrix(idx_bh, :),       'Rows', 'complete');
R_comb_raw     = corr(regressor_matrix(idx_comb, :),     'Rows', 'complete');

% --- METHOD 2: DETRENDED BLOCKS ---
% 1. Baseline: detrending the single continuous 30s block
baseline_reg_det = detrend(regressor_matrix(idx_baseline, :));
R_baseline_det   = corr(baseline_reg_det, 'Rows', 'complete');

% 2. Fingertapping: 21s blocks detrended individually
ft_reg_det = [];
for i = 1:length(onsets_stim1)
    onset_t = onsets_stim1(i);
    idx = time_vector >= onset_t & time_vector <= (onset_t + task_duration);
    ft_reg_det = [ft_reg_det; detrend(regressor_matrix(idx, :))];
end
R_ft_det = corr(ft_reg_det, 'Rows', 'complete');

% 3. Breath Holding: 21s blocks detrended individually
bh_reg_det = [];
for i = 1:length(onsets_stim2)
    onset_t = onsets_stim2(i);
    idx = time_vector >= onset_t & time_vector <= (onset_t + task_duration);
    bh_reg_det = [bh_reg_det; detrend(regressor_matrix(idx, :))];
end
R_bh_det = corr(bh_reg_det, 'Rows', 'complete');

% 4. Combined: 21s blocks detrended individually
comb_reg_det = [];
for i = 1:length(onsets_stim3)
    onset_t = onsets_stim3(i);
    idx = time_vector >= onset_t & time_vector <= (onset_t + task_duration);
    comb_reg_det = [comb_reg_det; detrend(regressor_matrix(idx, :))];
end
R_comb_det = corr(comb_reg_det, 'Rows', 'complete');

% =========================================================================
%% PART 6: PLOT RESULTING CORRELATION MATRICES (2 SEPARATE 2x2 FIGURES)
% =========================================================================
fprintf('6. Plotting the raw and detrended correlation matrices in 2x2 layout...\n');

n_colors = 128;
c_blue  = [0.2, 0.45, 0.75]; 
c_white = [1.0, 1.0, 1.0];  
c_red   = [0.85, 0.25, 0.25]; 

r_map = [linspace(c_blue(1), c_white(1), n_colors), linspace(c_white(1), c_red(1), n_colors)]';
g_map = [linspace(c_blue(2), c_white(2), n_colors), linspace(c_white(2), c_red(2), n_colors)]';
b_map = [linspace(c_blue(3), c_white(3), n_colors), linspace(c_white(3), c_red(3), n_colors)]';
rwb_colormap = [r_map, g_map, b_map];

condition_titles_raw = { ...
    '1. Baseline (Raw)', '2. Fingertapping (Raw)', ...
    '3. Breath Holding (Raw)', '4. Combined (Raw)' ...
};

condition_titles_det = { ...
    '1. Baseline (Detrended)', '2. Fingertapping (Detrended)', ...
    '3. Breath Holding (Detrended)', '4. Combined (Detrended)' ...
};

R_raw_all = {R_baseline_raw, R_ft_raw, R_bh_raw, R_comb_raw};
R_det_all = {R_baseline_det, R_ft_det, R_bh_det, R_comb_det};

% --- FIGURA 1: RAW MATRICES (2x2) ---
fig_raw = figure('Name', 'Task-Based Correlation: Raw Signals', ...
                 'Color', 'w', 'Position', [50, 50, 950, 850]);
for cond = 1:4
    subplot(2, 2, cond);
    plot_correlation_heatmap(R_raw_all{cond}, regressor_names, ...
                             condition_titles_raw{cond}, rwb_colormap);
end
sgtitle(sprintf('Task-Based Correlation Matrices: Raw Signals - %s', subject_name), ...
        'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'none');
exportgraphics(fig_raw, fullfile(subj_dir, [subject_name, '_TaskCorrelation_Raw_2x2.png']), 'Resolution', 300);

% --- FIGURA 2: DETRENDED MATRICES (2x2) ---
fig_det = figure('Name', 'Task-Based Correlation: Detrended Signals', ...
                 'Color', 'w', 'Position', [100, 100, 950, 850]);
for cond = 1:4
    subplot(2, 2, cond);
    plot_correlation_heatmap(R_det_all{cond}, regressor_names, ...
                             condition_titles_det{cond}, rwb_colormap);
end
sgtitle(sprintf('Task-Based Correlation Matrices: Detrended Signals - %s', subject_name), ...
        'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'none');
exportgraphics(fig_det, fullfile(subj_dir, [subject_name, '_TaskCorrelation_Detrended_2x2.png']), 'Resolution', 300);

fprintf('   -> Saved 2x2 Task Correlation plots (Raw and Detrended).\n');

% =========================================================================
%% PART 7: GLOBAL SESSION-WIDE CORRELATION (ENTIRE TIMELINE)
% =========================================================================
fprintf('7. Computing and plotting Global Session-Wide Correlation...\n');

R_global_raw = corr(regressor_matrix, 'Rows', 'complete');
R_global_det = corr(detrend(regressor_matrix), 'Rows', 'complete');

fig_global = figure('Name', 'Global Session-Wide Correlation Matrix', 'Color', 'w', 'Position', [150, 150, 1100, 480]);
subplot(1, 2, 1);
plot_correlation_heatmap(R_global_raw, regressor_names, 'Full Session: Raw (No Detrending)', rwb_colormap);

subplot(1, 2, 2);
plot_correlation_heatmap(R_global_det, regressor_names, 'Full Session: Detrended', rwb_colormap);

sgtitle(sprintf('Global Session-Wide Correlation Matrix (Full Timeline) - %s', subject_name), ...
        'FontSize', 15, 'FontWeight', 'bold', 'Interpreter', 'none');

file_fig_global = fullfile(subj_dir, [subject_name, '_Global_Correlation_Matrix.png']);
exportgraphics(fig_global, file_fig_global, 'Resolution', 300);
fprintf('   -> Saved Global Correlation plot: %s\n', [subject_name, '_Global_Correlation_Matrix.png']);

% =========================================================================
%% PART 8: MULTIMODAL TIME-SERIES OVERLAY PLOT WITH SHADED TASK BLOCKS
% =========================================================================
fprintf('8. Plotting Multimodal Time-Series Overlay with Task Blocks...\n');

fig_overlay = figure('Name', 'Multimodal Time-Series', 'Color', 'w', 'Position', [50, 50, 1400, 900]);

% Data to plot in 5 stacked subplots
signals_plot = { ...
    [sc_pca1_hbo, sc_pca1_hbr], ...                      % 1. fNIRS SC PCA
    [hr_ecg_aligned, hr_ppg_aligned], ...                % 2. Heart Rate
    resp_aligned, ...                                    % 3. Respiration
    [eda_tonic_aligned, eda_phasic_aligned], ...         % 4. EDA
    spo2_aligned ...                                     % 5. SpO2
};

titles_plot = { ...
    'fNIRS Scalp SC PCA (HbO & HbR in Z-score)', ...     
    'Heart Rate (ECG & PPG in BPM)', ...                 
    'Respiration Belt Signal (Normalized A.U.)', ...       
    'Electrodermal Activity (Tonic & Phasic in \muS)', ... 
    'Arterial Oxygen Saturation (SpO2 in %)' ...          
};

legends_plot = { ...
    {'SC PCA HbO', 'SC PCA HbR'}, ...
    {'HR (ECG)', 'HR (PPG)'}, ...
    {'Respiration'}, ...
    {'EDA Tonic', 'EDA Phasic'}, ...
    {'SpO2'} ...
};

colors_plot = { ...
    {[0.2 0.45 0.75], [0.85 0.25 0.25]}, ...
    {[0.85 0.25 0.25], [0.2 0.45 0.75]}, ...
    {[0 0.45 0.74]}, ...
    {[0.466 0.674 0.188], [0.929 0.694 0.125]}, ...
    {[0.494 0.184 0.556]} ...
};

for p = 1:5
    subplot(5, 1, p);
    data_curr = signals_plot{p};
    
    % Get y-axis bounds
    y_min = min(data_curr(:)) - 0.1 * abs(min(data_curr(:)));
    y_max = max(data_curr(:)) + 0.1 * abs(max(data_curr(:)));
    if y_min == y_max, y_min = y_min - 1; y_max = y_max + 1; end
    
    hold on;
    % Draw task patches
    h_patch_ft   = [];
    h_patch_bh   = [];
    h_patch_comb = [];
    
    for i = 1:length(onsets_stim1)
        h_p = patch([onsets_stim1(i) onsets_stim1(i)+task_duration onsets_stim1(i)+task_duration onsets_stim1(i)], ...
                    [y_min y_min y_max y_max], [0.929, 0.694, 0.125], 'EdgeColor', 'none', 'FaceAlpha', 0.20, 'HandleVisibility', 'off'); % FT = Yellow
        if i == 1, h_patch_ft = h_p; end
    end
    for i = 1:length(onsets_stim2)
        h_p = patch([onsets_stim2(i) onsets_stim2(i)+task_duration onsets_stim2(i)+task_duration onsets_stim2(i)], ...
                    [y_min y_min y_max y_max], [0.2, 0.75, 0.75], 'EdgeColor', 'none', 'FaceAlpha', 0.20, 'HandleVisibility', 'off'); % BH = Aqua Green
        if i == 1, h_patch_bh = h_p; end
    end
    for i = 1:length(onsets_stim3)
        h_p = patch([onsets_stim3(i) onsets_stim3(i)+task_duration onsets_stim3(i)+task_duration onsets_stim3(i)], ...
                    [y_min y_min y_max y_max], [0.85, 0.25, 0.25], 'EdgeColor', 'none', 'FaceAlpha', 0.20, 'HandleVisibility', 'off'); % COMB = Red
        if i == 1, h_patch_comb = h_p; end
    end
    
    % Plot signal curves
    h_lines = [];
    cols = colors_plot{p};
    for col = 1:size(data_curr, 2)
        h_l = plot(time_vector, data_curr(:, col), 'LineWidth', 1.5, 'Color', cols{col});
        h_lines = [h_lines, h_l];
    end
    
    ylim([y_min, y_max]);
    title(titles_plot{p}, 'FontSize', 11, 'FontWeight', 'bold');
    
    % Add Task Block legend to the 1st subplot
    if p == 1
        % Create dummy handles for task legends
        h_leg_ft   = patch(nan, nan, [0.929, 0.694, 0.125], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
        h_leg_bh   = patch(nan, nan, [0.2, 0.75, 0.75], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
        h_leg_comb = patch(nan, nan, [0.85, 0.25, 0.25], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
        
        legend([h_lines, h_leg_ft, h_leg_bh, h_leg_comb], ...
               [legends_plot{p}, {'FT', 'BH', 'COMB'}], ...
               'Location', 'northeast', 'FontSize', 8, 'NumColumns', 5);
    else
        legend(h_lines, legends_plot{p}, 'Location', 'northeast', 'FontSize', 9);
    end
    
    grid on;
    axis tight;
    if p == 5
        xlabel('Timeline (seconds)', 'FontSize', 11, 'FontWeight', 'bold');
    end
end

sgtitle(sprintf('Multimodal Time-Series Overlay with Task Block Alignment - %s', subject_name), ...
        'FontSize', 15, 'FontWeight', 'bold', 'Interpreter', 'none');

file_fig_overlay = fullfile(subj_dir, [subject_name, '_Multimodal_Time_Series_Overlay.png']);
exportgraphics(fig_overlay, file_fig_overlay, 'Resolution', 300);
fprintf('   -> Saved Multimodal Time-Series Overlay plot: %s\n', [subject_name, '_Multimodal_Time_Series_Overlay.png']);

% =========================================================================
%% PART 9: SAVE TASK CORRELATION MATRICES (.MAT)
% =========================================================================
fprintf('\n9. Saving Task Correlation data for GLM script...\n');

% Full session continuous detrended correlation matrix
regressor_matrix_detrended = detrend(regressor_matrix);
R_detrended_full = corr(regressor_matrix_detrended, 'Rows', 'complete');
R_raw_full       = corr(regressor_matrix, 'Rows', 'complete');

file_task_corr_mat = fullfile(subj_dir, [subject_name, '_TaskCorrelation.mat']);
save(file_task_corr_mat, 'R_detrended_full', 'R_raw_full', 'R_baseline_det', ...
    'R_ft_det', 'R_bh_det', 'R_comb_det', 'regressor_names', 'regressor_matrix');
fprintf('   -> Saved Task Correlation file: %s\n', [subject_name, '_TaskCorrelation.mat']);

fprintf('\n=== TASK CORRELATION ANALYSIS COMPLETE FOR %s ===\n', subject_name);

% =========================================================================
% AUXILIARY FUNCTION: PLOT CORRELATION HEATMAP
% =========================================================================
function plot_correlation_heatmap(R, names, title_str, cmap)
    imagesc(R);
    colorbar;
    clim([-1 1]); 
    colormap(cmap);
    axis square;
    
    num_vars = length(names);
    set(gca, 'XTick', 1:num_vars, 'XTickLabel', names, ...
             'YTick', 1:num_vars, 'YTickLabel', names, ...
             'FontSize', 10, 'FontWeight', 'bold');
    xtickangle(45); % labels are rotated by 45°
    title(title_str, 'FontSize', 13, 'FontWeight', 'bold');
    
    for i = 1:num_vars
        for j = 1:num_vars
            val = R(i, j);
            if isnan(val)
                text_str = 'NaN';
            else
                if abs(val) >= 0.12 && abs(val) < 1.0  % if the corr value is higher than the noise corr value (excluding the main diagonal R=1)
                    text_str = sprintf('%+.2f*', val); % overlays the numerical correlation value (with the sign) inside each cell
                else
                    text_str = sprintf('%+.2f', val); 
                end
            end
            
            if abs(val) > 0.45
                text_color = 'w'; % color white on dark cells (R>0.45)
            else
                text_color = 'k'; % color black on light cells (R<0.45)
            end
            
            text(j, i, text_str, 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', 'Color', text_color, ...
            'FontSize', 8, 'FontWeight', 'bold'); 
        end
    end
end
