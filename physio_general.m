%% PHYSIOLOGICAL ANALYSIS SCRIPT (NIRx WINGS & SNIRF)
clearvars;
clc;
close all;

% =========================================================================
% Academic Context: Master's Thesis in Bioengineering for Neuroscience
% University of Padova (DEI) | NIRx Medical Technologies LLC
% Reference: Chapter 3 (Signal Processing and Data Analysis Pipeline, Sec. 3.3)
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
% CONFIGURATION SECTION
% =========================================================================
% Channel mapping for physiological signals in the .wings data matrix
CH_RESP = 9;   % Respiration
CH_GSR  = 10;  % Galvanic Skin Response (GSR)
CH_ECG  = 11;  % ECG (Lead II)
CH_PPG  = 12;  % PPG (Raw Pulse Wave)
CH_SPO2 = 13;  % SpO2
CH_HR   = 14;  % Heart Rate (HR)

experimental_task_duration = 21; % seconds
fs_target = 4;                   % Target sampling frequency for EDA (Hz)

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

file_wings = fullfile(subj_dir, [subject_name, '.wings']);
file_snirf = fullfile(subj_dir, [subject_name, '.snirf']);

if ~exist(file_wings, 'file')
    error('Wings file not found in subject folder: %s', file_wings);
end
if ~exist(file_snirf, 'file')
    error('Snirf file not found in subject folder: %s', file_snirf);
end

fprintf('\nSubject Directory: %s\n', subj_dir);
fprintf('Found Wings File: %s\n', [subject_name, '.wings']);
fprintf('Found SNIRF File: %s\n', [subject_name, '.snirf']);

% =========================================================================
% LOADING PHYSIOLOGICAL DATA FROM .WINGS FILE
% =========================================================================
fprintf('\nReading physiological data...\n');
raw_wings = importdata(file_wings);
data_matrix = raw_wings.data;

real_time   = data_matrix(:, 1);  
respiration = data_matrix(:, CH_RESP);  
gsr         = data_matrix(:, CH_GSR); 
ecg         = data_matrix(:, CH_ECG); 
ppg_raw     = data_matrix(:, CH_PPG); 
spo2        = data_matrix(:, CH_SPO2); 
heart_rate  = data_matrix(:, CH_HR);

% Compute dynamic original sampling rate from time vector
fs_original = 1 / mean(diff(real_time));
fprintf('   -> Original Sampling Rate: %.2f Hz\n', fs_original);

% =========================================================================
% EXTRACTING AND FILTERING EXPERIMENTAL TRIGGERS FROM .SNIRF
% =========================================================================
fprintf('Extracting and filtering triggers from .snirf file...\n');
info_snirf = h5info(file_snirf, '/nirs');
group_list = {info_snirf.Groups.Name};
stim_indices = find(contains(group_list, '/nirs/stim'));
total_conditions = length(stim_indices); 

% Structure designed to store only the valid experimental tasks
trigger_struct = struct('name', {}, 'onsets', {}, 'durations', {});
valid_count = 0;

for c = 1:total_conditions
    path_base = sprintf('/nirs/stim%d', c);
    raw_name = h5read(file_snirf, [path_base, '/name']);
    
    if iscell(raw_name)
        cond_name = strtrim(raw_name{1});
    else
        cond_name = strtrim(char(raw_name'));
    end
    
    % Inclusion filter: select only conditions 1, 2, and 3
    if strcmp(cond_name, '1') || strcmp(cond_name, '2') || strcmp(cond_name, '3')
        valid_count = valid_count + 1;
        
        stim_data = h5read(file_snirf, [path_base, '/data']);
        if size(stim_data, 2) == 3
            onsets = stim_data(:, 1);
        else
            onsets = stim_data(1, :)';
        end
        
        % Overwrite durations array with the actual protocol block length
        durations = experimental_task_duration * ones(length(onsets), 1);
        
        % Dynamic mapping to specified descriptive acronyms
        switch cond_name
            case '1'
                mapped_name = 'FT';
            case '2'
                mapped_name = 'BH';
            case '3'
                mapped_name = 'COMB';
        end
        
        trigger_struct(valid_count).name = mapped_name;
        trigger_struct(valid_count).onsets = onsets;
        trigger_struct(valid_count).durations = durations;
        
        fprintf('Included -> %s | Events: %d | Duration forced to: %d seconds\n', ...
                trigger_struct(valid_count).name, length(onsets), experimental_task_duration);
    else
        fprintf('Excluded -> Condition %s (Training or Rest block)\n', cond_name);
    end
end

% =========================================================================
% MULTI-PANEL GRAPH GENERATION (6 SUBPLOTS WITH HIGHLIGHTED BLOCKS)
% =========================================================================
figure('Name', 'NIRx WINGS2 Physiology - Experimental Conditions', ...
       'Color', [1 1 1], 'Units', 'normalized', 'Position', [0.05 0.05 0.9 0.9]);
   
signal_names = {'Respiratory Belt', 'Galvanic Skin Response (GSR)', ...
                'ECG (Lead II)', 'PPG (Raw Pulse Wave)', ...
                'Heart Rate (BPM)', 'Oxygen Saturation (SpO2)'};
            
% Clean hardware SpO2 outliers (> 100% or < 70% non-physiological spikes from motion)
spo2_clean = spo2;
spo2_clean(spo2_clean > 100) = 100;
spo2_clean(spo2_clean < 70) = NaN;
spo2_clean = fillmissing(spo2_clean, 'nearest');

signal_data = {respiration, gsr, ecg, ppg_raw, heart_rate, spo2_clean};
measurement_units = {'Amplitude (A.U.)', 'Conductance (uS)', 'Voltage (mV)', 'Absorbance (A.U.)', 'Heart Rate (BPM)', 'Oxygen Saturation (%)'};
signal_colors = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.4940 0.1840 0.5560; 0.4660 0.6740 0.1880; 0.6350 0.0780 0.1840; 0.3010 0.7450 0.9330];

% Dynamically calculate experiment timeline bounds (1 min margin before 1st task onset, 1 min margin after last task end)
if exist('trigger_struct', 'var') && ~isempty(trigger_struct)
    all_onsets_exp = [];
    all_end_times_exp = [];
    for v_tmp = 1:length(trigger_struct)
        ons_tmp = trigger_struct(v_tmp).onsets(:);
        dur_tmp = trigger_struct(v_tmp).durations(:);
        all_onsets_exp = [all_onsets_exp; ons_tmp]; %#ok<AGROW>
        all_end_times_exp = [all_end_times_exp; ons_tmp + dur_tmp]; %#ok<AGROW>
    end
    exp_start_time = max(0, min(all_onsets_exp) - 60);
    exp_end_time = min(real_time(end), max(all_end_times_exp) + 60);
else
    exp_start_time = 60;
    exp_end_time = real_time(end) - 5;
end
fprintf('Focusing graphs exclusively on Experiment Timeline: %.1f s to %.1f s\n', exp_start_time, exp_end_time);

for p = 1:6
    ax(p) = subplot(6, 1, p);
    
    % Initialization of arrays to build the legend on this axis
    legend_handles = [];
    legend_names = {};
    
    idx_valid_fig1 = (real_time >= exp_start_time) & (real_time <= exp_end_time);
    valid_data_p = signal_data{p}(idx_valid_fig1);
    y_min = min(valid_data_p);
    y_max = max(valid_data_p);
    if y_min == y_max, y_max = y_min + 1; end
    
    hold on;
    
    % Drawing shaded rectangles background blocks for valid conditions
    if ~isempty(trigger_struct)
        for v = 1:length(trigger_struct)
            onset_times = trigger_struct(v).onsets;
            task_durations = trigger_struct(v).durations;
            
            % Specific color matching based on the new mapped names
            if strcmp(trigger_struct(v).name, 'FT')
                patch_color = [1 0.92 0.5];      % Light Yellow
            elseif strcmp(trigger_struct(v).name, 'BH')
                patch_color = [0.6 0.9 0.9];     % Light Aqua Green
            else
                patch_color = [1 0.75 0.75];    % Light Red (COMB)
            end
            
            for t = 1:length(onset_times)
                t_start = onset_times(t);
                t_end = t_start + task_durations(t);
                
                h_patch = patch([t_start t_end t_end t_start], [y_min y_min y_max y_max], ...
                                patch_color, 'EdgeColor', 'none', 'FaceAlpha', 0.6);
                
                % Store only the reference handle of the very first block per condition
                if t == 1
                    legend_handles(end+1) = h_patch; %#ok<AGROW>
                    legend_names{end+1} = trigger_struct(v).name; %#ok<AGROW>
                end
            end
        end
    end
    
    % Overlay physiological signal line plot above the shaded patches
    plot(real_time, signal_data{p}, 'LineWidth', 1.2, 'Color', signal_colors(p, :));
    
    % Render the customized conditions legend exclusively inside the first panel
    if p == 1 && ~isempty(legend_handles)
        legend(legend_handles, legend_names, 'Location', 'northeast','Box', 'on','AutoUpdate', 'off');
    end
    
    title(signal_names{p}, 'FontWeight', 'bold');
    ylabel(measurement_units{p});
    grid on;
    ylim([y_min, y_max]);
    hold off;
end
linkaxes(ax, 'x');
xlim(ax(1), [exp_start_time, exp_end_time]);
xlabel('Experiment Timeline (seconds)');
fprintf('\n=== PROCESS COMPLETED ===\n');

% =========================================================================
%% RESPIRATION BELT DATA PRE-PROCESSING
% =========================================================================
% 1. Define the frame length in samples (approx 1.5 seconds, must be odd)
frame_length = round(1.5 * fs_original);
if mod(frame_length, 2) == 0
    frame_length = frame_length + 1;
end

% 2. Apply the SAVITZKY-GOLAY FILTER to the raw respiration signal
poly_order = 3; 
respiration_smoothed = sgolayfilt(respiration, poly_order, frame_length);

% 3. Figure to compare Raw vs Preprocessed signal
figure('Name', 'Respiration Preprocessing Comparison', 'Color', [1 1 1]);
plot(real_time, respiration, 'Color', [0.7 0.7 0.7], 'LineWidth', 1); % Gray raw signal
hold on;
plot(real_time, respiration_smoothed, 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5); % Blue smooth signal
title('Respiration Signal: Raw vs Savitzky-Golay Filter');
xlabel('Time (seconds)');
ylabel('Amplitude (A.U.)');
legend('Raw Respiration', 'Filtered Respiration (1.5s window)');
grid on;
axis tight;
xlim([exp_start_time, exp_end_time]);

% =========================================================================
% GRAPHICAL OUTPUT GENERATION: RAW VS FILTERED RESPIRATION WITH OVERLAID TASK BLOCKS
% =========================================================================
respiration_raw = respiration;
figure('Name', 'Respiration Comparison and Task Alignment', ...
       'Color', [1 1 1], 'Units', 'normalized', 'Position', [0.1 0.2 0.8 0.6]);

y_min = min(respiration_raw);
y_max = max(respiration_raw);
if y_min == y_max, y_max = y_min + 1; end
hold on;

handles_tasks = [];
labels_tasks = {};

if exist('trigger_struct', 'var') && ~isempty(trigger_struct)
    for v = 1:length(trigger_struct)
        onset_times = trigger_struct(v).onsets;
        task_durations = trigger_struct(v).durations;
        
        if strcmp(trigger_struct(v).name, 'FT')
            patch_color = [1 0.92 0.5];      % Light Yellow
        elseif strcmp(trigger_struct(v).name, 'BH')
            patch_color = [0.6 0.9 0.9];     % Light Aqua Green
        else
            patch_color = [1 0.75 0.75];    % Light Red (COMB)
        end
        
        for t = 1:length(onset_times)
            t_start = onset_times(t);
            t_end = t_start + task_durations(t);
            
            h_patch = patch([t_start t_end t_end t_start], [y_min y_min y_max y_max], ...
                            patch_color, 'EdgeColor', 'none', 'FaceAlpha', 0.6);
            
            if t == 1
                handles_tasks(end+1) = h_patch; %#ok<AGROW>
                labels_tasks{end+1} = trigger_struct(v).name; %#ok<AGROW>
            end
        end
    end
end

h_raw = plot(real_time, respiration_raw, 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
h_filt = plot(real_time, respiration_smoothed, 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5);

all_handles = [h_raw, h_filt, handles_tasks];
all_labels  = ['Raw Respiration', 'Filtered Respiration (1.5s window)', labels_tasks];
legend(all_handles, all_labels,'Location','northeast','Box','on','AutoUpdate','off');
title('Respiratory Belt Signal: Preprocessing & Task Alignment', 'FontWeight', 'bold');
xlabel('Experiment Timeline (seconds)');
ylabel('Amplitude (A.U.)');
grid on;
axis tight;
xlim([exp_start_time, exp_end_time]);
hold off;

% =========================================================================
%% USING RESP DATA AS A TRIAL VALIDATION TOOL FOR BH/COMB
% =========================================================================
window_sec = [-2, 2]; % Window to check around the task onset
threshold_multiplier = 1.8; 

normal_breath_amp = prctile(respiration_smoothed, 95); 

validation_tasks = {'BH', 'COMB'};
valid_onsets_struct = struct('BH', [], 'COMB', []);
rejected_onsets_struct = struct('BH', [], 'COMB', []);

tInc_physio = ones(length(real_time), 1);

for task_idx = 1:length(validation_tasks)
    task_name = validation_tasks{task_idx};
    
    v_idx = find(strcmp({trigger_struct.name}, task_name));
    
    if isempty(v_idx)
        fprintf('\nSkipping validation for %s (not found in Snirf triggers)\n', task_name);
        continue;
    end
    
    onsets_task = trigger_struct(v_idx).onsets;
    num_trials = length(onsets_task);
    valid_trials = true(num_trials, 1);
    
    fprintf('\n--- Respiration Validation: %s ---\n', task_name);
    for i = 1:num_trials
        t_start = onsets_task(i);
        idx_window = find(real_time >= (t_start + window_sec(1)) & real_time <= (t_start + window_sec(2)));
        signal_window = respiration_smoothed(idx_window);
        local_max = max(signal_window);
        
        if local_max > (normal_breath_amp * threshold_multiplier)
            valid_trials(i) = false;
            idx_trial_period = find(real_time >= t_start & real_time <= (t_start + experimental_task_duration));
            tInc_physio(idx_trial_period) = 0;
            
            fprintf('  Trial %02d (Onset: %6.1fs) -> REJECTED (Max Amp: %.0f, Normal: %.0f)\n', ...
                i, t_start, local_max, normal_breath_amp);
        else
            fprintf('  Trial %02d (Onset: %6.1fs) -> VALID\n', i, t_start);
        end
    end
    
    valid_onsets_struct.(task_name) = onsets_task(valid_trials);
    rejected_onsets_struct.(task_name) = onsets_task(~valid_trials);
    
    fprintf('Validation Complete: Kept %d/%d %s trials.\n', sum(valid_trials), num_trials, task_name);
end

valid_onsets_BH = valid_onsets_struct.BH;
rejected_onsets_BH = rejected_onsets_struct.BH;
valid_onsets_Comb = valid_onsets_struct.COMB;
rejected_onsets_Comb = rejected_onsets_struct.COMB;

v_idx_ft = find(strcmp({trigger_struct.name}, 'FT'));
if ~isempty(v_idx_ft)
    valid_onsets_FT = trigger_struct(v_idx_ft).onsets;
else
    valid_onsets_FT = [];
end
valid_onsets_COMB = valid_onsets_Comb;

valid_trials_file = fullfile(subj_dir, [subject_name, '_ValidTrials.mat']);
save(valid_trials_file, 'valid_onsets_FT', 'valid_onsets_BH', 'valid_onsets_COMB', 'tInc_physio');
fprintf('\nSaved valid onsets and exclusion vector to:\n%s\n', valid_trials_file);

% =========================================================================
%% EDA/GSR PREPROCESSING: DOWNSAMPLING, NORMALIZATION & cvxEDA
% =========================================================================
fprintf('Downsampling EDA signal from %.1f Hz to %d Hz...\n', fs_original, fs_target);
[gsr_downsampled, time_downsampled] = resample(gsr, real_time, fs_target);

fprintf('Applying Z-score normalization...\n');
gsr_normalized = (gsr_downsampled - mean(gsr_downsampled)) / std(gsr_downsampled);

delta = 1 / fs_target; % Sampling interval (0.25 s)
tau0 = 2.0;            % Slow time constant (s)
tau1 = 0.7;            % Fast time constant (s)
delta_knot = 10;       % Distance between tonic spline knots (s)
alpha = 0.0008;        % Sparsity penalty
gamma = 0.01;          % Tonic component penalty

fprintf('Running cvxEDA algorithm... (This may take a few seconds/minutes)\n');
[gsr_phasic, p_smna, gsr_tonic, l, d, e, obj] = cvxEDA(gsr_normalized, delta, tau0, tau1, delta_knot, alpha, gamma);

% =========================================================================
% GRAPHICAL OUTPUT GENERATION: EDA DECOMPOSITION INTO 2 COMPONENTS
% =========================================================================
figure('Name', 'EDA Decomposition (cvxEDA) and Task Alignment', ...
       'Color', [1 1 1], 'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);

sec_to_ignore = exp_start_time;
idx_valid_eda = (time_downsampled >= exp_start_time) & (time_downsampled <= exp_end_time);

% --- 1: Normalized Signal ---
ax1 = subplot(3, 1, 1);
plot(time_downsampled, gsr_normalized, 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2);
title('Original Downsampled EDA Signal (Z-scored, 4 Hz)', 'FontWeight', 'bold');
ylabel('Z-score');
valid_gsr_norm = gsr_normalized(idx_valid_eda);
ylim([min(valid_gsr_norm) - 0.5, max(valid_gsr_norm) + 0.5]);
grid on; 

% --- 2: Tonic Component (SCL) ---
ax2 = subplot(3, 1, 2);
plot(time_downsampled, gsr_tonic, 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.5);
title('Tonic Component (SCL)', 'FontWeight', 'bold');
ylabel('Amplitude');
valid_gsr_tonic = gsr_tonic(idx_valid_eda);
ylim([min(valid_gsr_tonic) - 0.5, max(valid_gsr_tonic) + 0.5]);
grid on; 

% --- 3: Phasic Component (SCR) with Overlaid Task Blocks ---
ax3 = subplot(3, 1, 3);
valid_gsr_phasic = gsr_phasic(idx_valid_eda);
y_min = min(valid_gsr_phasic);
y_max = max(valid_gsr_phasic);
if y_min == y_max, y_max = y_min + 1; end
y_lim_max = y_max * 1.15;
hold on;

handles_tasks = [];
labels_tasks = {};

if exist('trigger_struct', 'var') && ~isempty(trigger_struct)
    for v = 1:length(trigger_struct)
        onset_times = trigger_struct(v).onsets;
        task_durations = trigger_struct(v).durations;
        
        if strcmp(trigger_struct(v).name, 'FT')
            patch_color = [1 0.92 0.5];      % Light Yellow
        elseif strcmp(trigger_struct(v).name, 'BH')
            patch_color = [0.6 0.9 0.9];     % Light Aqua Green
        else
            patch_color = [1 0.75 0.75];     % Light Red (COMB)
        end
        
        for t = 1:length(onset_times)
            t_start = onset_times(t);
            t_end = t_start + task_durations(t);
            
            h_patch = patch([t_start t_end t_end t_start], [-1 -1 y_lim_max y_lim_max], ...
                            patch_color, 'EdgeColor', 'none', 'FaceAlpha', 0.6);
            
            if t == 1
                handles_tasks(end+1) = h_patch; %#ok<AGROW>
                labels_tasks{end+1} = trigger_struct(v).name; %#ok<AGROW>
            end
        end
    end
end

h_phasic = plot(time_downsampled, gsr_phasic, 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 1.5);

all_handles = [h_phasic, handles_tasks];
all_labels  = ['Phasic SCR', labels_tasks];
legend(all_handles, all_labels, 'Location', 'northeast', 'Box', 'on', 'AutoUpdate', 'off');

title('Phasic Component (SCR)', 'FontWeight', 'bold');
xlabel('Experiment Timeline (seconds)');
ylabel('Amplitude');
ylim([-0.2, y_lim_max]); 
grid on; 
hold off;

linkaxes([ax1, ax2, ax3], 'x');
xlim(ax1, [exp_start_time, exp_end_time]);
fprintf('\n=== EDA PREPROCESSING COMPLETED ===\n');

% =========================================================================
%% AUTOMATIC EXTRACTION OF PHASIC PEAKS (SCR Amplitude) from EDA DATA
% =========================================================================
fprintf('\n=== STATISTICAL ANALYSIS OF SKIN CONDUCTANCE RESPONSE (SCR) ===\n');

post_task_window = 10; % extended window of 10s after task
scr_results = struct();

if exist('trigger_struct', 'var') && ~isempty(trigger_struct)
    for v = 1:length(trigger_struct)
        condition_name = trigger_struct(v).name;
        onset_times = trigger_struct(v).onsets;
        task_durations = trigger_struct(v).durations;
        
        max_peaks = zeros(length(onset_times), 1);
        
        fprintf('\nCondition: %s\n', condition_name);
        fprintf('----------------------------------------\n');
        
        for t = 1:length(onset_times)
            t_start = onset_times(t);
            t_end_extended = t_start + task_durations(t) + post_task_window;
            
            idx_start = find(time_downsampled >= t_start, 1, 'first');
            idx_end = find(time_downsampled <= t_end_extended, 1, 'last');
            
            if ~isempty(idx_start) && ~isempty(idx_end)
                phasic_segment = gsr_phasic(idx_start:idx_end);
                max_peaks(t) = max(phasic_segment);
                fprintf('  Trial %d (Onset: %.1fs): Max Peak = %.3f Z-score\n', t, t_start, max_peaks(t));
            end
        end
        
        if strcmp(condition_name, 'BH')
            valid_mask = ismember(onset_times, valid_onsets_BH);
        elseif strcmp(condition_name, 'COMB')
            valid_mask = ismember(onset_times, valid_onsets_COMB);
        else
            valid_mask = true(size(onset_times));
        end
        
        scr_results(v).name = condition_name;
        scr_results(v).single_peaks = max_peaks;
        scr_results(v).mean_all = mean(max_peaks);
        scr_results(v).std_all = std(max_peaks);
        scr_results(v).mean_valid = mean(max_peaks(valid_mask));
        scr_results(v).std_valid = std(max_peaks(valid_mask));
        
        fprintf('  -> ALL TRIALS MEAN = %.3f (Std = %.3f)\n', scr_results(v).mean_all, scr_results(v).std_all);
        fprintf('  -> VALID TRIALS MEAN = %.3f (Std = %.3f)\n', scr_results(v).mean_valid, scr_results(v).std_valid);
    end
end

% Grouped Bar Chart (All Trials vs Valid Trials Only)
figure('Name', 'Sympathetic Activation Statistical Analysis', 'Color', [1 1 1], 'Position', [100, 100, 700, 480]);
names = {scr_results.name};
data_matrix_scr = [[scr_results.mean_all]', [scr_results.mean_valid]'];
error_matrix_scr = [[scr_results.std_all]', [scr_results.std_valid]'];

h_bar = bar(data_matrix_scr, 'grouped');
h_bar(1).FaceColor = [0.8500 0.3250 0.0980]; % Orange for All Trials
h_bar(2).FaceColor = [0.4660 0.6740 0.1880]; % Green for Valid Trials Only
hold on;

nbars = size(data_matrix_scr, 2);
ngroups = size(data_matrix_scr, 1);
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, data_matrix_scr(:,i), error_matrix_scr(:,i), 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 10);
end
set(gca, 'XTick', 1:ngroups, 'XTickLabel', names, 'FontSize', 12, 'FontWeight', 'bold');
title('Average Peak SCR Amplitude (All vs Valid Trials)', 'FontSize', 14);
ylabel('Max Amplitude (Z-score)', 'FontSize', 12, 'FontWeight', 'bold');
legend({'All Trials', 'Valid Trials'}, 'Location', 'northeast', 'Box', 'on');
grid on; hold off;

% =========================================================================
%% PPG PRE-PROCESSING
% =========================================================================
fprintf('\n=== PPG ANALYSIS AND HEART RATE CALCULATION ===\n');

% 1. Band-pass filtering to isolate the heartbeat
fprintf('Band-pass filtering of the PPG signal (0.5 - 5 Hz)...\n');
[b_ppg, a_ppg] = butter(3, [0.5 5]/(fs_original/2), 'bandpass');
ppg_filtered = filtfilt(b_ppg, a_ppg, ppg_raw);

% 2. Peak detection (Adaptive Moving Window Peak Detection)
fprintf('Searching for systolic peaks (Adaptive Moving Window Peak Detection)...\n');
sec_to_ignore_ppg = exp_start_time;

win_samples = round(15 * fs_original);
local_std = movstd(ppg_filtered, win_samples);
local_max_amp = movmax(ppg_filtered, win_samples);
local_thresh = 0.5 * local_std;

[all_pks, all_locs, ~, all_prom] = findpeaks(ppg_filtered, ...
    'MinPeakDistance', 0.45 * fs_original, ...
    'MinPeakHeight', 0);

valid_idx = (all_prom >= local_thresh(all_locs)) & (all_pks >= 0.25 * local_max_amp(all_locs));
pks = all_pks(valid_idx);
locs = all_locs(valid_idx);

beat_times = real_time(locs);
ibi = diff(beat_times); 
instantaneous_hr = 60 ./ ibi; 
hr_times = beat_times(2:end);

valid_hr_mask = (instantaneous_hr >= 40) & (instantaneous_hr <= 140);
hr_times_clean = hr_times(valid_hr_mask);
instantaneous_hr_clean = medfilt1(instantaneous_hr(valid_hr_mask), 5, 'truncate');

continuous_hr = interp1(hr_times_clean, instantaneous_hr_clean, real_time, 'pchip', 'extrap');

% =========================================================================
% GRAPHICAL OUTPUT GENERATION: PPG AND HR OVER TASKS
% =========================================================================
figure('Name', 'Heart Rate Analysis from Photoplethysmography (PPG)', ...
       'Color', [1 1 1], 'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.7]);

sec_to_ignore_ppg = 60;

% --- Panel 1: Filtered PPG with Detected Peaks ---
ax1 = subplot(2, 1, 1);
plot(real_time, ppg_filtered, 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
hold on;
plot(beat_times, pks, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 4);
title('Filtered PPG Signal and Systolic Peak Detection', 'FontWeight', 'bold');
ylabel('Amplitude (A.U.)');
grid on; axis tight;
xlim([exp_start_time, exp_end_time]);
hold off;

% --- Panel 2: Continuous HR Trace with Task Overlay ---
ax2 = subplot(2, 1, 2);
valid_hr = continuous_hr(real_time >= exp_start_time & real_time <= exp_end_time);
y_min_hr = max(40, min(valid_hr) - 5); 
y_max_hr = min(180, max(valid_hr) + 5); 

hold on;
handles_tasks_hr = [];
labels_tasks_hr = {};

if exist('trigger_struct', 'var') && ~isempty(trigger_struct)
    for v = 1:length(trigger_struct)
        onset_times = trigger_struct(v).onsets;
        task_durations = trigger_struct(v).durations;
        
        if strcmp(trigger_struct(v).name, 'FT')
            patch_color = [1 0.92 0.5];      % Light Yellow
        elseif strcmp(trigger_struct(v).name, 'BH')
            patch_color = [0.6 0.9 0.9];     % Light Aqua Green
        else
            patch_color = [1 0.75 0.75];    % Light Red (COMB)
        end
        
        for t = 1:length(onset_times)
            t_start = onset_times(t);
            t_end = t_start + task_durations(t);
            
            h_patch = patch([t_start t_end t_end t_start], ...
                            [y_min_hr y_min_hr y_max_hr y_max_hr], ...
                            patch_color, 'EdgeColor', 'none', 'FaceAlpha', 0.6);
            
            if t == 1
                handles_tasks_hr(end+1) = h_patch; %#ok<AGROW>
                labels_tasks_hr{end+1} = trigger_struct(v).name; %#ok<AGROW>
            end
        end
    end
end

h_hr = plot(real_time, continuous_hr, 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.5);
all_handles_hr = [h_hr, handles_tasks_hr];
all_labels_hr  = ['HR', labels_tasks_hr];
legend(all_handles_hr, all_labels_hr, 'Location', 'northeast', 'Box', 'on', 'AutoUpdate', 'off');

title('HR (PPG Derived)', 'FontWeight', 'bold');
xlabel('Experiment Timeline (seconds)');
ylabel('Heart Rate (BPM)');
ylim([y_min_hr, y_max_hr]);
grid on; axis tight;
xlim([exp_start_time, exp_end_time]);
hold off;

linkaxes([ax1, ax2], 'x');
fprintf('=== HR EXTRACTION COMPLETED ===\n\n');

% =========================================================================
% GRAPH: RAW vs FILTERED PPG COMPARISON
% =========================================================================
figure('Name', 'PPG Preprocessing Effect: Raw vs Filtered', ...
       'Color', [1 1 1], 'Units', 'normalized', 'Position', [0.15 0.15 0.7 0.7]);
ppg_raw_centered = ppg_raw - mean(ppg_raw); 

subplot(2, 1, 1);
plot(real_time, ppg_raw_centered, 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
hold on;
plot(real_time, ppg_filtered, 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 1.5);
title('PPG: Hemodynamic Baseline Removal and Centering', 'FontWeight', 'bold');
xlabel('Timeline (seconds)');
ylabel('Absorbance (Centered A.U.)');
legend('Raw PPG (Centered)', 'Filtered PPG (0.5 - 5 Hz)', 'Location', 'northeast', 'Box', 'on');
grid on; axis tight;
xlim([exp_start_time, exp_end_time]);
hold off;

zoom_start_ppg = 400;
zoom_end_ppg = 410;
idx_window = real_time >= zoom_start_ppg & real_time <= zoom_end_ppg;
offset_locale = mean(ppg_raw(idx_window));

subplot(2, 1, 2);
plot(real_time, ppg_raw - offset_locale, 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
hold on;
plot(real_time, ppg_filtered, 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 1.5);

idx_zoom_ppg = beat_times >= zoom_start_ppg & beat_times <= zoom_end_ppg;
plot(beat_times(idx_zoom_ppg), pks(idx_zoom_ppg), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6);

title(sprintf('Zoom (%d - %d s)', zoom_start_ppg, zoom_end_ppg), 'FontWeight', 'bold');
xlabel('Timeline (seconds)');
ylabel('Absorbance (Centered A.U.)');
xlim([zoom_start_ppg, zoom_end_ppg]);

ppg_zoom_window_filt = ppg_filtered(idx_window);
ylim([min(ppg_zoom_window_filt) - 0.05, max(ppg_zoom_window_filt) + 0.05]);
legend('Raw PPG (Centered)', 'Filtered PPG', 'Detected Systolic Peaks', 'Location', 'northeast', 'Box', 'on');
grid on; 
hold off;

fprintf('=== RAW/FILTERED PPG COMPARISON GRAPH GENERATED ===\n');

% =========================================================================
%% AUTOMATIC EXTRACTION OF CARDIAC DATA (HR from PPG)
% =========================================================================
fprintf('\n=== STATISTICAL ANALYSIS OF HEART RATE (HR) ===\n');

baseline_window = 5; % seconds
hr_results = struct();

if exist('trigger_struct', 'var') && ~isempty(trigger_struct)
    for v = 1:length(trigger_struct)
        condition_name = trigger_struct(v).name;
        onset_times = trigger_struct(v).onsets;
        task_durations = trigger_struct(v).durations;
        
        trial_hr_delta = zeros(length(onset_times), 1);
        fprintf('\nCondition: %s\n', condition_name);
        fprintf('----------------------------------------\n');
        
        for t = 1:length(onset_times)
            t_start = onset_times(t);
            t_end = t_start + task_durations(t);
            
            idx_base_start = find(real_time >= (t_start - baseline_window), 1, 'first');
            idx_base_end = find(real_time < t_start, 1, 'last');
            idx_task_start = find(real_time >= t_start, 1, 'first');
            idx_task_end = find(real_time <= t_end, 1, 'last');
            
            if ~isempty(idx_base_start) && ~isempty(idx_task_start)
                hr_baseline = mean(continuous_hr(idx_base_start:idx_base_end));
                hr_task = mean(continuous_hr(idx_task_start:idx_task_end));
                delta_hr = hr_task - hr_baseline;
                trial_hr_delta(t) = delta_hr;
                fprintf('  Trial %d (Onset: %.1fs): Baseline = %.1f BPM | Task = %.1f BPM -> Variation = %+.1f BPM\n', ...
                        t, t_start, hr_baseline, hr_task, delta_hr);
            end
        end
        
        if strcmp(condition_name, 'BH')
            valid_mask = ismember(onset_times, valid_onsets_BH);
        elseif strcmp(condition_name, 'COMB')
            valid_mask = ismember(onset_times, valid_onsets_COMB);
        else
            valid_mask = true(size(onset_times));
        end
        
        hr_results(v).name = condition_name;
        hr_results(v).single_deltas = trial_hr_delta;
        hr_results(v).mean_all = mean(trial_hr_delta);
        hr_results(v).std_all = std(trial_hr_delta);
        hr_results(v).mean_valid = mean(trial_hr_delta(valid_mask));
        hr_results(v).std_valid = std(trial_hr_delta(valid_mask));
        
        fprintf('  -> ALL TRIALS MEAN DELTA = %+.1f BPM (Std = %.1f)\n', hr_results(v).mean_all, hr_results(v).std_all);
        fprintf('  -> VALID TRIALS MEAN DELTA = %+.1f BPM (Std = %.1f)\n', hr_results(v).mean_valid, hr_results(v).std_valid);
    end
end

% Grouped Bar Chart (Delta HR PPG: All vs Valid)
figure('Name', 'Heart Rate Variation Statistical Analysis (PPG)', ...
       'Color', [1 1 1], 'Position', [150, 150, 700, 480]);
hr_names = {hr_results.name};
data_matrix_ppg = [[hr_results.mean_all]', [hr_results.mean_valid]'];
error_matrix_ppg = [[hr_results.std_all]', [hr_results.std_valid]'];

h_bar = bar(data_matrix_ppg, 'grouped');
h_bar(1).FaceColor = [0.8500 0.3250 0.0980]; % Orange for All Trials
h_bar(2).FaceColor = [0.4660 0.6740 0.1880]; % Green for Valid Trials Only
hold on;

nbars = size(data_matrix_ppg, 2);
ngroups = size(data_matrix_ppg, 1);
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, data_matrix_ppg(:,i), error_matrix_ppg(:,i), 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 10);
end
set(gca, 'XTick', 1:ngroups, 'XTickLabel', hr_names, 'FontSize', 12, 'FontWeight', 'bold');
yline(0, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
legend({'All Trials', 'Valid Trials'}, 'Location', 'northeast', 'Box', 'on');
grid on; hold off;

% =========================================================================
%% ECG PREPROCESSING & R-PEAK DETECTION
% =========================================================================
fprintf('\n=== ECG ANALYSIS ===\n');

% 1. Band-pass filtering [0.5 - 40 Hz]
fprintf('Band-pass filtering of the ECG signal (0.5 - 40 Hz)...\n');
[b_ecg, a_ecg] = butter(3, [0.5 40]/(fs_original/2), 'bandpass');
ecg_filtered = filtfilt(b_ecg, a_ecg, ecg);

% 2. R-Peak Detection
fprintf('Searching for R-peaks (QRS Complex)...\n');
min_peak_dist_ecg = 0.35 * fs_original;
height_threshold_ecg = max(0.3, 0.25 * max(ecg_filtered));          
prominence_threshold_ecg = max(0.3, 0.25 * max(ecg_filtered)); 

[pks_ecg, locs_ecg] = findpeaks(ecg_filtered, ...
    'MinPeakDistance', min_peak_dist_ecg, ...
    'MinPeakHeight', height_threshold_ecg, ...
    'MinPeakProminence', prominence_threshold_ecg);

beat_times_ecg = real_time(locs_ecg);
ibi_ecg = diff(beat_times_ecg); 
instantaneous_hr_ecg = 60 ./ ibi_ecg; 
hr_times_ecg = beat_times_ecg(2:end);

valid_hr_ecg_mask = (instantaneous_hr_ecg >= 40) & (instantaneous_hr_ecg <= 140);
hr_times_ecg_clean = hr_times_ecg(valid_hr_ecg_mask);
instantaneous_hr_ecg_clean = medfilt1(instantaneous_hr_ecg(valid_hr_ecg_mask), 5, 'truncate');

continuous_hr_ecg = interp1(hr_times_ecg_clean, instantaneous_hr_ecg_clean, real_time, 'pchip', 'extrap');

% =========================================================================
% GRAPHICAL OUTPUT GENERATION: ECG AND R-PEAK DETECTION
% =========================================================================
figure('Name', 'Electrocardiogram Analysis (ECG - Lead II)', ...
       'Color', [1 1 1], 'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.7]);

sec_to_ignore_ecg = exp_start_time;

% --- Panel 1: Filtered ECG with Detected R-Peaks ---
ax1_ecg = subplot(2, 1, 1);
plot(real_time, ecg_filtered, 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 1);
hold on;
plot(beat_times_ecg, pks_ecg, 'kv', 'MarkerFaceColor', 'k', 'MarkerSize', 5);
title('Filtered ECG Signal and R-Peak Detection', 'FontWeight', 'bold');
ylabel('Voltage (mV)');
grid on; axis tight;
xlim([exp_start_time, exp_end_time]);
hold off;

% --- Panel 2: Continuous HR Trace (from ECG) with Tasks ---
ax2_ecg = subplot(2, 1, 2);
valid_hr_ecg = continuous_hr_ecg(real_time >= exp_start_time & real_time <= exp_end_time);
y_min_hr_ecg = max(40, min(valid_hr_ecg) - 5); 
y_max_hr_ecg = min(180, max(valid_hr_ecg) + 5); 

hold on;
handles_tasks_ecg = [];
labels_tasks_ecg = {};

if exist('trigger_struct', 'var') && ~isempty(trigger_struct)
    for v = 1:length(trigger_struct)
        onset_times = trigger_struct(v).onsets;
        task_durations = trigger_struct(v).durations;
        
        if strcmp(trigger_struct(v).name, 'FT')
            patch_color = [1 0.92 0.5];      % Light Yellow
        elseif strcmp(trigger_struct(v).name, 'BH')
            patch_color = [0.6 0.9 0.9];     % Light Aqua Green
        else
            patch_color = [1 0.75 0.75];    % Light Red (COMB)
        end
        
        for t = 1:length(onset_times)
            t_start = onset_times(t);
            t_end = t_start + task_durations(t);
            
            h_patch = patch([t_start t_end t_end t_start], ...
                            [y_min_hr_ecg y_min_hr_ecg y_max_hr_ecg y_max_hr_ecg], ...
                            patch_color, 'EdgeColor', 'none', 'FaceAlpha', 0.6);
            if t == 1
                handles_tasks_ecg(end+1) = h_patch; %#ok<AGROW>
                labels_tasks_ecg{end+1} = trigger_struct(v).name; %#ok<AGROW>
            end
        end
    end
end

h_hr_ecg = plot(real_time, continuous_hr_ecg, 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.5);
legend([h_hr_ecg, handles_tasks_ecg], ['HR from ECG (BPM)', labels_tasks_ecg], 'Location', 'northeast', 'Box', 'on', 'AutoUpdate', 'off');
title('HR (ECG Derived)', 'FontWeight', 'bold');
xlabel('Experiment Timeline (seconds)');
ylabel('Heart Rate (BPM)');
ylim([y_min_hr_ecg, y_max_hr_ecg]);
grid on; axis tight;
xlim([exp_start_time, exp_end_time]);
hold off;

linkaxes([ax1_ecg, ax2_ecg], 'x');
fprintf('=== ECG EXTRACTION COMPLETED ===\n\n');

% =========================================================================
% GRAPH: RAW vs FILTERED ECG COMPARISON
% =========================================================================
figure('Name', 'ECG Preprocessing Effect: Raw vs Filtered', ...
       'Color', [1 1 1], 'Units', 'normalized', 'Position', [0.15 0.15 0.7 0.7]);

subplot(2, 1, 1);
plot(real_time, ecg, 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
hold on;
plot(real_time, ecg_filtered, 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 1.5);
title('ECG', 'FontWeight', 'bold');
xlabel('Timeline (seconds)');
ylabel('Voltage (mV)');
legend('Raw ECG', 'Filtered ECG (0.5 - 40 Hz)', 'Location', 'northeast', 'Box', 'on');
grid on; axis tight;
xlim([exp_start_time, exp_end_time]);
hold off;

zoom_start = 400;
zoom_end = 405;
subplot(2, 1, 2);
plot(real_time, ecg, 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
hold on;
plot(real_time, ecg_filtered, 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 1.5);

idx_zoom = beat_times_ecg >= zoom_start & beat_times_ecg <= zoom_end;
plot(beat_times_ecg(idx_zoom), pks_ecg(idx_zoom), 'kv', 'MarkerFaceColor', 'k', 'MarkerSize', 7);

title(sprintf('Zoom (%d - %d s)', zoom_start, zoom_end), 'FontWeight', 'bold');
xlabel('Timeline (seconds)');
ylabel('Voltage (mV)');
xlim([zoom_start, zoom_end]);

ecg_zoom_window = ecg(real_time >= zoom_start & real_time <= zoom_end);
ylim([min(ecg_zoom_window) - 0.1, max(ecg_zoom_window) + 0.1]);
legend('Raw ECG', 'Filtered ECG', 'Detected R-Peaks', 'Location', 'northeast', 'Box', 'on');
grid on; 
hold off;

% =========================================================================
%% AUTOMATIC EXTRACTION OF CARDIAC DATA (HR from ECG)
% =========================================================================
fprintf('\n=== STATISTICAL ANALYSIS OF HEART RATE (ECG DERIVED) ===\n');

hr_results_ecg = struct();

if exist('trigger_struct', 'var') && ~isempty(trigger_struct)
    for v = 1:length(trigger_struct)
        condition_name = trigger_struct(v).name;
        onset_times = trigger_struct(v).onsets;
        task_durations = trigger_struct(v).durations;
        
        trial_hr_delta = zeros(length(onset_times), 1);
        fprintf('\nCondition: %s\n', condition_name);
        fprintf('----------------------------------------\n');
        
        for t = 1:length(onset_times)
            t_start = onset_times(t);
            t_end = t_start + task_durations(t);
            
            idx_base_start = find(real_time >= (t_start - baseline_window), 1, 'first');
            idx_base_end = find(real_time < t_start, 1, 'last');
            idx_task_start = find(real_time >= t_start, 1, 'first');
            idx_task_end = find(real_time <= t_end, 1, 'last');
            
            if ~isempty(idx_base_start) && ~isempty(idx_task_start)
                hr_baseline = mean(continuous_hr_ecg(idx_base_start:idx_base_end));
                hr_task = mean(continuous_hr_ecg(idx_task_start:idx_task_end));
                delta_hr = hr_task - hr_baseline;
                trial_hr_delta(t) = delta_hr;
                fprintf('  Trial %d (Onset: %.1fs): Baseline = %.1f BPM | Task = %.1f BPM -> Variation = %+.1f BPM\n', ...
                        t, t_start, hr_baseline, hr_task, delta_hr);
            end
        end
        
        if strcmp(condition_name, 'BH')
            valid_mask = ismember(onset_times, valid_onsets_BH);
        elseif strcmp(condition_name, 'COMB')
            valid_mask = ismember(onset_times, valid_onsets_COMB);
        else
            valid_mask = true(size(onset_times));
        end
        
        hr_results_ecg(v).name = condition_name;
        hr_results_ecg(v).single_deltas = trial_hr_delta;
        hr_results_ecg(v).mean_all = mean(trial_hr_delta);
        hr_results_ecg(v).std_all = std(trial_hr_delta);
        hr_results_ecg(v).mean_valid = mean(trial_hr_delta(valid_mask));
        hr_results_ecg(v).std_valid = std(trial_hr_delta(valid_mask));
        
        fprintf('  -> ALL TRIALS MEAN DELTA = %+.1f BPM (Std = %.1f)\n', hr_results_ecg(v).mean_all, hr_results_ecg(v).std_all);
        fprintf('  -> VALID TRIALS MEAN DELTA = %+.1f BPM (Std = %.1f)\n', hr_results_ecg(v).mean_valid, hr_results_ecg(v).std_valid);
    end
end

% Grouped Bar Chart (Delta HR ECG: All vs Valid)
figure('Name', 'Heart Rate Variation Statistical Analysis (ECG)', ...
       'Color', [1 1 1], 'Position', [150, 150, 700, 480]);
hr_names_ecg = {hr_results_ecg.name};
data_matrix_ecg = [[hr_results_ecg.mean_all]', [hr_results_ecg.mean_valid]'];
error_matrix_ecg = [[hr_results_ecg.std_all]', [hr_results_ecg.std_valid]'];

h_bar = bar(data_matrix_ecg, 'grouped');
h_bar(1).FaceColor = [0.4940 0.1840 0.5560]; % Purple for All Trials
h_bar(2).FaceColor = [0.4660 0.6740 0.1880]; % Green for Valid Trials Only
hold on;

nbars = size(data_matrix_ecg, 2);
ngroups = size(data_matrix_ecg, 1);
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, data_matrix_ecg(:,i), error_matrix_ecg(:,i), 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 10);
end
set(gca, 'XTick', 1:ngroups, 'XTickLabel', hr_names_ecg, 'FontSize', 12, 'FontWeight', 'bold');
yline(0, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
legend({'All Trials', 'Valid Trials'}, 'Location', 'northeast', 'Box', 'on');
grid on; hold off;

% =========================================================================
% GRAPH: HR METHODOLOGICAL COMPARISON (WINGS2 vs PPG vs ECG)
% =========================================================================
figure('Name', 'Methodological Validation: HR Comparison', ...
       'Color', [1 1 1], 'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.7]);

% --- Panel 1: Full Experiment Timeline ---
ax_cmp1 = subplot(2, 1, 1);
plot(real_time, heart_rate, 'Color', [0.7 0.7 0.7], 'LineWidth', 2, 'DisplayName', 'HR WINGS2 Hardware');
hold on;
plot(real_time, continuous_hr, 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 1.5, 'DisplayName', 'HR PPG (Derived)');
plot(real_time, continuous_hr_ecg, 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 1.5, 'DisplayName', 'HR ECG (Gold Standard)');
title('Full Experiment Timeline: PPG vs ECG Heart Rate Comparison', 'FontWeight', 'bold');
xlabel('Experiment Timeline (seconds)');
ylabel('Heart Rate (BPM)');
xlim([exp_start_time, exp_end_time]);
grid on; legend('Location', 'northeast', 'Box', 'on');
hold off;

% --- Panel 2: Zoomed Window Detail ---
comparison_zoom_start = 500;
comparison_zoom_end = 650;
ax_cmp2 = subplot(2, 1, 2);
plot(real_time, heart_rate, 'Color', [0.7 0.7 0.7], 'LineWidth', 2.5, 'DisplayName', 'HR WINGS2 Hardware');
hold on;
plot(real_time, continuous_hr, 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 1.5, 'DisplayName', 'HR PPG (Derived)');
plot(real_time, continuous_hr_ecg, 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 1.5, 'DisplayName', 'HR ECG (Gold Standard)');

title(sprintf('Zoom Detail (%d - %d s)', comparison_zoom_start, comparison_zoom_end), 'FontWeight', 'bold');
xlabel('Experiment Timeline (seconds)');
ylabel('Heart Rate (BPM)');
xlim([comparison_zoom_start, comparison_zoom_end]);

hr_zoom = heart_rate(real_time >= comparison_zoom_start & real_time <= comparison_zoom_end);
ylim([min(hr_zoom)-5, max(hr_zoom)+5]);
legend('Location', 'northeast', 'Box', 'on');
grid on; 
hold off;

fprintf('=== HR METHODOLOGICAL COMPARISON GRAPH GENERATED ===\n');

% =========================================================================
%% SpO2 PRE-PROCESSING AND EXTRACTION
% =========================================================================
fprintf('\n=== STARTING SpO2 ANALYSIS ===\n');

median_window = round(5 * fs_original); 
spo2_filtered = medfilt1(spo2_clean, median_window, 'truncate');

spo2_results = struct('name', {}, 'single_deltas', {}, 'avg_delta', {}, 'std_delta', {});
baseline_sec = 5;

for v = 1:length(trigger_struct)
    condition_name = trigger_struct(v).name;
    onset_times = trigger_struct(v).onsets;
    task_durations = trigger_struct(v).durations;
    
    current_deltas = zeros(1, length(onset_times));
    
    for t = 1:length(onset_times)
        t_start_task = onset_times(t);
        t_end_task = t_start_task + task_durations(t);
        t_start_baseline = t_start_task - baseline_sec;
        t_end_baseline = t_start_task;
        
        idx_baseline = real_time >= t_start_baseline & real_time < t_end_baseline;
        idx_task = real_time >= t_start_task & real_time < t_end_task;
        
        avg_baseline = mean(spo2_filtered(idx_baseline), 'omitnan');
        avg_task = mean(spo2_filtered(idx_task), 'omitnan');
        current_deltas(t) = avg_task - avg_baseline;
    end
    
    spo2_results(v).name = condition_name;
    spo2_results(v).single_deltas = current_deltas;
    spo2_results(v).avg_delta = mean(current_deltas);
    spo2_results(v).std_delta = std(current_deltas);
end

% Plot SpO2 Cleaning
figure('Name', 'SpO2 Cleaning and Validation', 'Color', [1 1 1]);
hold on;

idx_valid_spo2 = (real_time >= exp_start_time) & (real_time <= exp_end_time);
valid_spo2_data = spo2_filtered(idx_valid_spo2);
y_min_spo2 = max(80, floor(min(valid_spo2_data) - 1));
y_max_spo2 = 100.5;

for v = 1:length(trigger_struct)
    if strcmp(trigger_struct(v).name, 'FT'), patch_color = [1 0.92 0.5];      % Light Yellow
    elseif strcmp(trigger_struct(v).name, 'BH'), patch_color = [0.6 0.9 0.9];     % Light Aqua Green
    else, patch_color = [1 0.75 0.75];    % Light Red (COMB)
    end
    
    for t = 1:length(trigger_struct(v).onsets)
        t_start = trigger_struct(v).onsets(t);
        t_end = t_start + trigger_struct(v).durations(t);
        patch([t_start t_end t_end t_start], [y_min_spo2 y_min_spo2 y_max_spo2 y_max_spo2], ...
            patch_color, 'EdgeColor', 'none', 'FaceAlpha', 0.4, 'HandleVisibility', 'off');
    end
end

plot(real_time, spo2, 'Color', [0.8 0.8 0.8], 'LineWidth', 1, 'DisplayName', 'Raw SpO2');
plot(real_time, spo2_filtered, 'Color', [0 0.4470 0.7410], 'LineWidth', 2, 'DisplayName', 'Median Filtered SpO2');

title('SpO2 Signal: Raw vs Cleaned', 'FontWeight', 'bold');
xlabel('Time (s)'); ylabel('Oxygen Saturation (%)');
xlim([exp_start_time, exp_end_time]);
ylim([y_min_spo2, y_max_spo2]);
legend('Location', 'best'); grid on; 
hold off;

fprintf('=== SpO2 ANALYSIS COMPLETED ===\n');

% =========================================================================
%% EXPORT DATA FOR SATORI (SHORT CHANNEL REGRESSION / GLM)
% =========================================================================
fprintf('\n=== EXPORT PROCESSED PHYSIO DATA FOR SATORI ===\n');

% Realignment (interp from 4Hz)
gsr_tonic_interp = interp1(time_downsampled, gsr_tonic, real_time, 'pchip', 'extrap');
gsr_phasic_interp = interp1(time_downsampled, gsr_phasic, real_time, 'pchip', 'extrap');

% Ensure column vectors (Nx1)
real_time = real_time(:);
respiration_smoothed = respiration_smoothed(:);
gsr_tonic_interp = gsr_tonic_interp(:);
gsr_phasic_interp = gsr_phasic_interp(:);
continuous_hr = continuous_hr(:);
continuous_hr_ecg = continuous_hr_ecg(:);
spo2_filtered = spo2_filtered(:);

physio_table = table(real_time, ...
                     respiration_smoothed, ...
                     gsr_tonic_interp, ...
                     gsr_phasic_interp, ...
                     continuous_hr, ...
                     continuous_hr_ecg, ...
                     spo2_filtered, ...
                     'VariableNames', {'Time_s', 'Resp_Clean', 'EDA_Tonic', 'EDA_Phasic', 'HR_PPG', 'HR_ECG', 'SpO2_Clean'});

[~, base_name, ~] = fileparts(file_wings);
output_file = fullfile(subj_dir, [base_name, '_ProcessedPhysio.csv']);
writetable(physio_table, output_file);
fprintf('File saved successfully for Satori in:\n%s\n', output_file);

% EXPORT DATA TO .SDM FORMAT FOR SATORI 
predictor_names = physio_table.Properties.VariableNames(2:end);
num_predictors = length(predictor_names);
num_datapoints = height(physio_table);
sdm_data_matrix = table2array(physio_table(:, 2:end));

% Standardization (Z-Score)
sdm_data_matrix_z = zscore(sdm_data_matrix);

sdm_output_file = fullfile(subj_dir, [base_name, '_PhysioRegressors.sdm']);
fid = fopen(sdm_output_file, 'w');

fprintf(fid, 'FileVersion:        1\n');
fprintf(fid, 'NrOfPredictors:     %d\n', num_predictors);
fprintf(fid, 'NrOfDataPoints:     %d\n', num_datapoints);
fprintf(fid, 'IncludesConstant:   0\n');
fprintf(fid, 'FirstConfoundPredictor: 1\n\n');

for i = 1:num_predictors
    fprintf(fid, '150 150 150\n'); 
end
fprintf(fid, '\n');

for i = 1:num_predictors
    fprintf(fid, '"%s" ', predictor_names{i});
end
fprintf(fid, '\n');

for row = 1:num_datapoints
    fprintf(fid, '%f ', sdm_data_matrix_z(row, :));
    fprintf(fid, '\n');
end

fclose(fid);
fprintf('SDM File created successfully and ready for Satori in:\n%s\n', sdm_output_file);

% Auto-save key Preprocessing Figures (300 DPI PNG)
fig_resp = findobj('Name', 'Respiration Comparison and Task Alignment');
if ~isempty(fig_resp), exportgraphics(fig_resp, fullfile(subj_dir, [base_name, '_Respiration_Preprocessing.png']), 'Resolution', 300); end

fig_eda = findobj('Name', 'EDA Decomposition (cvxEDA) and Task Alignment');
if ~isempty(fig_eda), exportgraphics(fig_eda, fullfile(subj_dir, [base_name, '_cvxEDA_Decomposition.png']), 'Resolution', 300); end

fig_hr = findobj('Name', 'Methodological Validation: HR Comparison');
if ~isempty(fig_hr), exportgraphics(fig_hr, fullfile(subj_dir, [base_name, '_HR_Validation.png']), 'Resolution', 300); end

fig_spo2 = findobj('Name', 'SpO2 Cleaning and Validation');
if ~isempty(fig_spo2), exportgraphics(fig_spo2, fullfile(subj_dir, [base_name, '_SpO2_Cleaning.png']), 'Resolution', 300); end

fig_scr_stat = findobj('Name', 'Sympathetic Activation Statistical Analysis');
if ~isempty(fig_scr_stat), exportgraphics(fig_scr_stat, fullfile(subj_dir, [base_name, '_EDA_Average_SCR_Peak.png']), 'Resolution', 300); end

fig_hr_ecg_stat = findobj('Name', 'Heart Rate Variation Statistical Analysis (ECG)');
if ~isempty(fig_hr_ecg_stat), exportgraphics(fig_hr_ecg_stat, fullfile(subj_dir, [base_name, '_HR_ECG_Average_Change.png']), 'Resolution', 300); end

fig_hr_ppg_stat = findobj('Name', 'Heart Rate Variation Statistical Analysis (PPG)');
if ~isempty(fig_hr_ppg_stat), exportgraphics(fig_hr_ppg_stat, fullfile(subj_dir, [base_name, '_HR_PPG_Average_Change.png']), 'Resolution', 300); end

fprintf('   -> Saved Physio Preprocessing plots (PNG 300 DPI).\n');

fprintf('\n=== ANALYSIS COMPLETE FOR SUBJECT %s ===\n', subject_name);
