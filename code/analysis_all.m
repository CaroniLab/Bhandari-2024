
%% Load Data
clear
read_Intan_RHD2000_file

fileinfo = dir(fullfile(path, 'time.dat'));
num_samples = fileinfo.bytes/4; % int32 = 4 bytes
fid = fopen(fullfile(path,'time.dat'), 'r');
t = fread(fid, num_samples, 'int32');
fclose(fid);
t = t / frequency_parameters.amplifier_sample_rate;

num_channels = length(board_adc_channels); % ADC input info from header file
fileinfo = dir(fullfile(path,'analogin.dat'));
num_samples = fileinfo.bytes/(num_channels * 2); % uint16 = 2 bytes
fid = fopen(fullfile(path,'analogin.dat'), 'r');
v = fread(fid, [num_channels, num_samples], 'uint16');
fclose(fid);
v = v * 0.000050354;

ch =0;
fileinfo = dir(fullfile(path,'digitalin.dat'));
num_samples = fileinfo.bytes/2; % uint16 = 2 bytes
fid = fopen(fullfile(path,'digitalin.dat'), 'r');
digital_word = fread(fid, num_samples, 'uint16');
fclose(fid);
digital_input_ch = (bitand(digital_word, 2^ch) > 0);



%% Camera timestamp extraction
camera_diff = diff(digital_input_ch);
camera_ts = t(camera_diff ==1);
disp(sprintf ('%d frames detected', length(camera_ts)));
disp('extracted camera timestamps')


%% Loading behavior data
mkdir(fullfile(path, 'plots')) 
mkdir(fullfile(path, 'plots_DLC'))
disp('Made plots directory')
annotfile = dir(fullfile(path, '*Annot_*.csv'));
annotfile  = annotfile.name;
events = csvread(fullfile(path, annotfile), 1, 0);

home1_start = camera_ts(events(1,1)+1);
home1_stop = camera_ts(events(1,2)+1);
FOR_start = camera_ts(events(2,1)+1);
FOR_stop = camera_ts(events(2,2)+1);
home2_start = camera_ts(events(3,1)+1);
home2_stop = camera_ts(events(3,2));
disp('Loaded Behavior data')

%% Import DLC annotations to run only first time, avoid if curation done
overwrite = 1;
thresh=70; %Set number of pixels for interaction
manual_obj = 0; %set to 1 if override DLC object location with manual label and to 0 to use DLC detected objects.
if isfile((fullfile(path, 'Object_1_approach.csv')))
    prompt = "Are you sure you wish to overwrite the annotations? 1 for yes or 0 for no: ";
    overwrite= input(prompt);
end
if overwrite == 1
    
    DLC_loading;
    disp('DLC data processed')
elseif overwrite == 0
    disp("Did not load DLC data again, kept annotations from before")
else 
    disp("Incorrect input, did not load DLC data again, kept annotations from before")
end

%% Load post curation labels
obj1_app_DLC = readtable((fullfile(path, 'Object_1_approach.csv')));
obj1_app_DLC  = camera_ts(obj1_app_DLC {:,:});
obj2_app_DLC =  readtable((fullfile(path, 'Object_2_approach.csv')));
obj2_app_DLC  = camera_ts(obj2_app_DLC {:,:});
obj_app_DLC = [obj1_app_DLC;obj2_app_DLC];
context_app_DLC = readtable((fullfile(path, 'Context_1_interact.csv')));
context_app_DLC  = camera_ts(context_app_DLC {:,:});
disp('Loaded timestamps')

%% Single unit analysis


%% load spikes:
group_to_plot = 'All';   % WIN (Wide Interneuron); NIN (Narrow Interneuron); PYR (pyramidal); All = all; INH = WIN + NIN; All 

disp('Loading spike data')
kilosort_dir = dir(fullfile(path, 'Kilosort*'));
kilosort_dir = kilosort_dir.name;
spike_times = readNPY(fullfile(path, kilosort_dir, 'spike_times.npy'));
spike_clusters = readNPY(fullfile(path, kilosort_dir, 'spike_clusters.npy'));
[~, ~, spike_info] = tsvread(fullfile(path, kilosort_dir, 'cluster_info.tsv'));
properties = spike_info(1,:);
spike_info = cell2table(spike_info(2:end,:));
spike_info.Properties.VariableNames = properties;
CellExp_cell_metrics_fn = dir(fullfile(path, '*cell_metrics.cellinfo.mat'));
CellExp_cell_metrics = load (fullfile(path, CellExp_cell_metrics_fn.name));
spike_info.ch = str2double(spike_info.ch); 


good_unit_ids_ctx = spike_info.cluster_id(strcmp(spike_info.group, 'good') & (spike_info.ch < 32));
good_unit_ids_vhc = spike_info.cluster_id(strcmp(spike_info.group, 'good') & (spike_info.ch >= 32));
put_cell_ID = CellExp_cell_metrics.cell_metrics.putativeCellType;

neurons_ctx = {};
neurons_vhc = {};
for neuron = 1:length(good_unit_ids_ctx)
    include = 0;
    if strcmp(group_to_plot, 'WIN')
        include = strcmp(put_cell_ID(neuron), 'Wide Interneuron');
    elseif strcmp(group_to_plot, 'NIN')
        include = strcmp(put_cell_ID(neuron), 'Narrow Interneuron');
    
    elseif strcmp(group_to_plot, 'INH')
        include = strcmp(put_cell_ID(neuron), 'Narrow Interneuron')  | strcmp(put_cell_ID(neuron), 'Narrow Interneuron');
    
    elseif strcmp(group_to_plot, 'PYR')
        include =  strcmp(put_cell_ID(neuron), 'Pyramidal Cell');
    
    elseif strcmp(group_to_plot, 'All')
        include = 1;
    end
    if include == 1
        cluster_id_temp = str2double(good_unit_ids_ctx(neuron)); 
        spike_times_temp = spike_times(spike_clusters == cluster_id_temp);
        neurons_ctx_av_firing_rate(neuron) = length(spike_times_temp)/t(end); 
        neurons_ctx{end+1} = t(spike_times_temp);
    end
end
for neuron = 1:length(good_unit_ids_vhc)
    include = 0;
    if strcmp(group_to_plot, 'WIN')
        include = strcmp(put_cell_ID(neuron + length(good_unit_ids_ctx)), 'Wide Interneuron');
    elseif strcmp(group_to_plot, 'NIN')
        include = strcmp(put_cell_ID(neuron + length(good_unit_ids_ctx)), 'Narrow Interneuron');
    
    elseif strcmp(group_to_plot, 'INH')
        include = strcmp(put_cell_ID(neuron + length(good_unit_ids_ctx)), 'Narrow Interneuron')  | strcmp(put_cell_ID(neuron + length(good_unit_ids_ctx)), 'Narrow Interneuron');
    
    elseif strcmp(group_to_plot, 'PYR')
        include =  strcmp(put_cell_ID(neuron + length(good_unit_ids_ctx)), 'Pyramidal Cell');
    
    elseif strcmp(group_to_plot, 'All')
        include = 1;
    end
    if include ==1
        cluster_id_temp = str2double(good_unit_ids_vhc(neuron)); 
        spike_times_temp = spike_times(spike_clusters == cluster_id_temp);
        neurons_vhc_av_firing_rate(neuron) = length(spike_times_temp)/t(end); 
        neurons_vhc{end+1} = t(spike_times_temp);
    end
end

bins_av_firingrate = 0:2.5:50;
figure
histogram(neurons_ctx_av_firing_rate, bins_av_firingrate,'Normalization','probability')
xlabel('firing rate')
ylabel('number of cortical neurons')
title('Distribution of firing rates of cortical neurons')
saveas(gcf, fullfile(path, 'plots', sprintf('Firing_rate_distribution_CTX.png')))
figure
histogram(neurons_vhc_av_firing_rate, bins_av_firingrate,'Normalization','probability')
xlabel('firing rate')
ylabel('number of VHC neurons')
title('Distribution of firing rates of VHC neurons')
saveas(gcf, fullfile(path, 'plots', sprintf('Firing_rate_distribution_VHC.png')))

disp('Done')
%% Plot Activity of neurons aligned to Context approach 
event = context_app_DLC;
save_tag = 'Context_interact';
before_context = 2;
after_context = 2;   %in seconds
sample_rate = 30000; %in Hz
n_trials = length(event);


for neuron = 1:length(good_unit_ids_ctx)
    spikes_neruon = neurons_ctx{neuron};
    figure
    for trial = 1:n_trials
       window_start_context = event(trial) - abs(before_context);
       window_end_context = event(trial) + abs(after_context);
       spikes_trial = spikes_neruon((spikes_neruon>window_start_context) & (spikes_neruon<window_end_context)) - event(trial);
       scatter(spikes_trial, repelem(trial, length(spikes_trial)), 2)
       hold on
       xline(0 , '--')
    end
    
    
    saveas(gcf, fullfile(path, 'plots', sprintf('%s_ctx_neuron_%d_spike_psth_Context_approach.png',save_tag, neuron)))
    hold off
end

for neuron = 1:length(good_unit_ids_vhc)
    spikes_neruon = neurons_vhc{neuron};
    figure
    for trial = 1:n_trials
       window_start_context = event(trial) - abs(before_context);
       window_end_context = event(trial) + abs(after_context);
       spikes_trial = spikes_neruon((spikes_neruon>window_start_context) & (spikes_neruon<window_end_context)) - event(trial);
       scatter(spikes_trial, repelem(trial, length(spikes_trial)), 2, 'filled')
       hold on
       xline(0, '--')
    end
    saveas(gcf, fullfile(path, 'plots', sprintf('%s_vhc_neuron_%d_spike_psth_Context_approach.png',save_tag, neuron)))
    hold off
end

%% Plot Activity of neurons aligned to Object approach 
%%exclusion_threshold = 0.3; %Fraction of trials with less than min_num_spikes spikes above which neurons will be excluded; if 0 no neurons are excluded
%%min_num_spikes = 1; %neurons having less than min_num_spikes spike in more than exclusion_threshold fraction of trials will be excluded

event = obj_app_DLC;
save_tag = 'Object1_approach';
before = 4;
after = 4;   %in seconds
sample_rate = 30000; %in Hz
n_trials = length(event);
for neuron = 1:length(good_unit_ids_ctx)
    spikes_neruon = neurons_ctx{neuron};
    figure
    for trial = 1:n_trials
       window_start = event(trial) - abs(before);
       window_end = event(trial) + abs(after);
       spikes_trial = spikes_neruon((spikes_neruon>window_start) & (spikes_neruon<window_end)) - event(trial);
       scatter(spikes_trial, repelem(trial, length(spikes_trial)), 2)
       %count = count + (length(spikes_trial) < min_num_spikes);
       hold on

       xline(0 , '--')

    end
    
    %bad_trials_ratio = count/n_trials;
    %if bad_trials_ratio > exclusion_threshold
     %   disp(sprintf('CTX neuron %d excluded', neuron));
    %else

        saveas(gcf, fullfile(path, 'plots', sprintf('%s_ctx_neuron_%d_spike_psth_Object_approach.png',save_tag, neuron)))
    end
    hold off
%end

for neuron = 1:length(good_unit_ids_vhc)
    spikes_neruon = neurons_vhc{neuron};
    figure
    for trial = 1:n_trials
       window_start = event(trial) - abs(before);
       window_end = event(trial) + abs(after);
       spikes_trial = spikes_neruon((spikes_neruon>window_start) & (spikes_neruon<window_end)) - event(trial);
       scatter(spikes_trial, repelem(trial, length(spikes_trial)), 2, 'filled')
       %count = count + (length(spikes_trial) < min_num_spikes);

       hold on
       xline(0, '--')
    end
     %  bad_trials_ratio = count/n_trials;
    %if bad_trials_ratio > exclusion_threshold
     %  disp(sprintf('VHC neuron %d excluded', neuron));
    %else
        saveas(gcf, fullfile(path, 'plots', sprintf('%s_vhc_neuron_%d_spike_psth_object_approach.png',save_tag, neuron)))
    end
    hold off
%end


%% Average firing rates for object approach
% New feature 07/02/2023: average of average activity of trials over time for each neuron
% in different behaviors to allow for a comparison;
% window can be selected below in time_av_window_before and time_av_window_after 
time_av_window_before = 2; % in seconds
time_av_window_after = 2; % in seconds


event = obj_app_DLC;
save_tag = 'Object_approach';
before = 4;
after = 4;   %in seconds
sample_rate = 30000; %in Hz
n_trials = length(event);
bin_width = 50;  %ms
bins = -abs(before):bin_width/1000:after;
time = bins(1: end-1) + bin_width/1000;
av_firing_rate_neuron_obj_ctx = {};
sem_firing_rate_neuron_obj_ctx = {};

for neuron = 1:length(good_unit_ids_ctx)
    spikes_neruon = neurons_ctx{neuron};
    figure
    firing_rate_trials = zeros(n_trials, length(time));
    for trial = 1:n_trials
       window_start = event(trial) - abs(before);
       window_end = event(trial) + abs(after);
       spikes_trial = spikes_neruon((spikes_neruon>window_start) & (spikes_neruon<window_end)) - event(trial);

       firing_rate_trials(trial, :) = histcounts(spikes_trial, bins)*1000/bin_width;
       plot(time, firing_rate_trials(trial, :))
       scatter(spikes_trial, repelem(trial, length(spikes_trial)), 2)
       hold on
       xline(0 , '--')
    end
    av_firing_rate_neuron_obj_ctx{neuron} = mean(firing_rate_trials, 1);
    sem_firing_rate_neuron_obj_ctx{neuron} = std(firing_rate_trials, 1)/sqrt(n_trials);


    time_av_before = (before -  time_av_window_before)/(bin_width/1000);
    time_av_after = (before +  time_av_window_after)/(bin_width/1000);
    time_av_firing_rate_neuron_obj_ctx(neuron) = mean(av_firing_rate_neuron_obj_ctx{neuron}(time_av_before:time_av_after), 2);
    
    
    saveas(gcf, fullfile(path, 'plots', sprintf('%s_ctx_neuron_%d_spike_psth_obj_approach_firing rate.png',save_tag, neuron)))
    hold off
end


av_firing_rate_neuron_obj_vhc = {};
sem_firing_rate_neuron_obj_vhc = {};
for neuron = 1:length(good_unit_ids_vhc)
    spikes_neruon = neurons_vhc{neuron};
    figure
    firing_rate_trials = zeros(n_trials, length(time));
    for trial = 1:n_trials
       window_start = event(trial) - abs(before);
       window_end = event(trial) + abs(after);
       spikes_trial = spikes_neruon((spikes_neruon>window_start) & (spikes_neruon<window_end)) - event(trial);
       firing_rate_trials(trial, :) = histcounts(spikes_trial, bins)*1000/bin_width;
       plot(time, firing_rate_trials(trial, :))
       scatter(spikes_trial, repelem(trial, length(spikes_trial)), 2)
       xline(0 , '--')
    end
    av_firing_rate_neuron_obj_vhc{neuron} = mean(firing_rate_trials, 1);
    sem_firing_rate_neuron_obj_vhc{neuron} = std(firing_rate_trials, 1)/sqrt(n_trials);

    time_av_before = (before -  time_av_window_before)/(bin_width/1000);
    time_av_after = (before +  time_av_window_after)/(bin_width/1000);
    time_av_firing_rate_neuron_obj_vhc(neuron) = mean(av_firing_rate_neuron_obj_vhc{neuron}(time_av_before:time_av_after), 2);

    saveas(gcf, fullfile(path, 'plots', sprintf('%s_ctx_neuron_%d_spike_psth_obj_approach_firing rate.png',save_tag, neuron)))
    hold off
end
% for neuron = 1:length(good_unit_ids_vhc)
%     spikes_neruon = neurons_vhc{neuron};
%     figure
%     for trial = 1:n_trials
%        window_start = event(trial) - abs(before);
%        window_end = event(trial) + abs(after);
%        spikes_trial = spikes_neruon((spikes_neruon>window_start) & (spikes_neruon<window_end)) - event(trial);
%        scatter(spikes_trial, repelem(trial, length(spikes_trial)), 2, 'filled')
%        hold on
%        xline(0, '--')
%     end
%     saveas(gcf, fullfile(path, 'plots', sprintf('%s_vhc_neuron_%d_spike_psth_obj_approach.png',save_tag, neuron)))
%     hold off
% end
%% Average firing rates for object1 approach
% New feature 07/02/2023: average of average activity of trials over time for each neuron
% in different behaviors to allow for a comparison;
% window can be selected below in time_av_window_before and time_av_window_after 
time_av_window_before = 1; % in seconds
time_av_window_after = 1; % in seconds


event = obj1_app_DLC;
save_tag = 'Object1_approach';
before = 4;
after = 4;   %in seconds
sample_rate = 30000; %in Hz
n_trials = length(event);
bin_width = 50;  %ms
bins = -abs(before):bin_width/1000:after;
time = bins(1: end-1) + bin_width/1000;
av_firing_rate_neuron_obj1_ctx = {};
for neuron = 1:length(good_unit_ids_ctx)
    spikes_neruon = neurons_ctx{neuron};
   figure
    firing_rate_trials = zeros(n_trials, length(time));
    for trial = 1:n_trials
       window_start = event(trial) - abs(before);
       window_end = event(trial) + abs(after);
       spikes_trial = spikes_neruon((spikes_neruon>window_start) & (spikes_neruon<window_end)) - event(trial);
       firing_rate_trials(trial, :) = histcounts(spikes_trial, bins)*1000/bin_width;
       plot(time, firing_rate_trials(trial, :))
       scatter(spikes_trial, repelem(trial, length(spikes_trial)), 2)
       hold on
       xline(0 , '--')
    end
    av_firing_rate_neuron_obj1_ctx{neuron} = mean(firing_rate_trials, 1);
    sem_firing_rate_neuron_obj1_ctx{neuron} = std(firing_rate_trials, 1)/sqrt(n_trials);
    
    time_av_before = (before -  time_av_window_before)/(bin_width/1000);
    time_av_after = (before +  time_av_window_after)/(bin_width/1000);
    time_av_firing_rate_neuron_obj1_ctx(neuron) = mean(av_firing_rate_neuron_obj1_ctx{neuron}(time_av_before:time_av_after), 2);




    saveas(gcf, fullfile(path, 'plots', sprintf('%s_ctx_neuron_%d_spike_psth_obj1_approach_firing rate.png',save_tag, neuron)))
    hold off
end


av_firing_rate_neuron_obj1_vhc = {};
for neuron = 1:length(good_unit_ids_vhc)
    spikes_neruon = neurons_vhc{neuron};
    figure
    firing_rate_trials = zeros(n_trials, length(time));
    for trial = 1:n_trials
       window_start = event(trial) - abs(before);
       window_end = event(trial) + abs(after);
       spikes_trial = spikes_neruon((spikes_neruon>window_start) & (spikes_neruon<window_end)) - event(trial);
       firing_rate_trials(trial, :) = histcounts(spikes_trial, bins)*1000/bin_width;
       plot(time, firing_rate_trials(trial, :))
       scatter(spikes_trial, repelem(trial, length(spikes_trial)), 2)
       hold on
       xline(0 , '--')
    end
    av_firing_rate_neuron_obj1_vhc{neuron} = mean(firing_rate_trials, 1);
    sem_firing_rate_neuron_obj1_vhc{neuron} = std(firing_rate_trials, 1)/sqrt(n_trials);


    time_av_before = (before -  time_av_window_before)/(bin_width/1000);
    time_av_after = (before +  time_av_window_after)/(bin_width/1000);
    time_av_firing_rate_neuron_obj1_vhc(neuron) = mean(av_firing_rate_neuron_obj1_vhc{neuron}(time_av_before:time_av_after), 2);


    saveas(gcf, fullfile(path, 'plots', sprintf('%s_ctx_neuron_%d_spike_psth_obj1_approach_firing rate.png',save_tag, neuron)))
    hold off
end
% for neuron = 1:length(good_unit_ids_vhc)
%     spikes_neruon = neurons_vhc{neuron};
%     figure
%     for trial = 1:n_trials
%        window_start = event(trial) - abs(before);
%        window_end = event(trial) + abs(after);
%        spikes_trial = spikes_neruon((spikes_neruon>window_start) & (spikes_neruon<window_end)) - event(trial);
%        scatter(spikes_trial, repelem(trial, length(spikes_trial)), 2, 'filled')
%        hold on
%        xline(0, '--')
%     end
%     saveas(gcf, fullfile(path, 'plots', sprintf('%s_vhc_neuron_%d_spike_psth_obj_approach.png',save_tag, neuron)))
%     hold off
% end

%% Average firing rates for object2 approach
% New feature 07/02/2023: average of average activity of trials over time for each neuron
% in different behaviors to allow for a comparison;
% window can be selected below in time_av_window_before and time_av_window_after 
time_av_window_before = 1; % in seconds
time_av_window_after = 1; % in seconds

event = obj2_app_DLC;
save_tag = 'Object2_approach';
before = 4;
after = 4;   %in seconds
sample_rate = 30000; %in Hz
n_trials = length(event);
bin_width = 50;  %ms
bins = -abs(before):bin_width/1000:after;
time = bins(1: end-1) + bin_width/1000;
av_firing_rate_neuron_obj2_ctx = {};
for neuron = 1:length(good_unit_ids_ctx)
    spikes_neruon = neurons_ctx{neuron};
     figure
    firing_rate_trials = zeros(n_trials, length(time));
    for trial = 1:n_trials
       window_start = event(trial) - abs(before);
       window_end = event(trial) + abs(after);
       spikes_trial = spikes_neruon((spikes_neruon>window_start) & (spikes_neruon<window_end)) - event(trial);
       firing_rate_trials(trial, :) = histcounts(spikes_trial, bins)*1000/bin_width;
       plot(time, firing_rate_trials(trial, :))
       scatter(spikes_trial, repelem(trial, length(spikes_trial)), 2)
        hold on
       xline(0 , '--')
    end
    av_firing_rate_neuron_obj2_ctx{neuron} = mean(firing_rate_trials, 1);
    sem_firing_rate_neuron_obj2_ctx{neuron} = std(firing_rate_trials, 1)/sqrt(n_trials);


    time_av_before = (before -  time_av_window_before)/(bin_width/1000);
    time_av_after = (before +  time_av_window_after)/(bin_width/1000);
    time_av_firing_rate_neuron_obj2_ctx(neuron) = mean(av_firing_rate_neuron_obj2_ctx{neuron}(time_av_before:time_av_after), 2);


    saveas(gcf, fullfile(path, 'plots', sprintf('%s_ctx_neuron_%d_spike_psth_obj2_approach_firing rate.png',save_tag, neuron)))
    hold off
end


av_firing_rate_neuron_obj2_vhc = {};
for neuron = 1:length(good_unit_ids_vhc)
    spikes_neruon = neurons_vhc{neuron};
    figure
    firing_rate_trials = zeros(n_trials, length(time));
    for trial = 1:n_trials
       window_start = event(trial) - abs(before);
       window_end = event(trial) + abs(after);
       spikes_trial = spikes_neruon((spikes_neruon>window_start) & (spikes_neruon<window_end)) - event(trial);
       firing_rate_trials(trial, :) = histcounts(spikes_trial, bins)*1000/bin_width;
       plot(time, firing_rate_trials(trial, :))
       scatter(spikes_trial, repelem(trial, length(spikes_trial)), 2)
       hold on
       xline(0 , '--')
    end
    av_firing_rate_neuron_obj2_vhc{neuron} = mean(firing_rate_trials, 1);
    sem_firing_rate_neuron_obj2_vhc{neuron} = std(firing_rate_trials, 1)/sqrt(n_trials);

    time_av_before = (before -  time_av_window_before)/(bin_width/1000);
    time_av_after = (before +  time_av_window_after)/(bin_width/1000);
    time_av_firing_rate_neuron_obj2_vhc(neuron) = mean(av_firing_rate_neuron_obj2_vhc{neuron}(time_av_before:time_av_after), 2);



    saveas(gcf, fullfile(path, 'plots', sprintf('%s_ctx_neuron_%d_spike_psth_obj2_approach_firing rate.png',save_tag, neuron)))
    hold off
end
% for neuron = 1:length(good_unit_ids_vhc)
%     spikes_neruon = neurons_vhc{neuron};
%     figure
%     for trial = 1:n_trials
%        window_start = event(trial) - abs(before);
%        window_end = event(trial) + abs(after);
%        spikes_trial = spikes_neruon((spikes_neruon>window_start) & (spikes_neruon<window_end)) - event(trial);
%        scatter(spikes_trial, repelem(trial, length(spikes_trial)), 2, 'filled')
%        hold on
%        xline(0, '--')
%     end
%     saveas(gcf, fullfile(path, 'plots', sprintf('%s_vhc_neuron_%d_spike_psth_obj_approach.png',save_tag, neuron)))
%     hold off
% end

%% Average firing rates for Context Interact
% New feature 07/02/2023: average of average activity of trials over time for each neuron
% in different behaviors to allow for a comparison;
% window can be selected below in time_av_window_before and time_av_window_after 
time_av_window_before = 1; % in seconds
time_av_window_after = 1; % in seconds


event = context_app_DLC;
save_tag = 'Context_interact';
before_context = 4;
after_context = 4;   %in seconds
sample_rate = 30000; %in Hz
n_trials = length(event);
bin_width = 50;  %ms
bins = -abs(before_context):bin_width/1000:after;
time = bins(1: end-1) + bin_width/1000;
av_firing_rate_neuron_ctx_context = {};
for neuron = 1:length(good_unit_ids_ctx)
    spikes_neruon = neurons_ctx{neuron};
    figure
    firing_rate_trials = zeros(n_trials, length(time));
    for trial = 1:n_trials
       window_start_context = event(trial) - abs(before_context);
       window_end_context = event(trial) + abs(after_context);
       spikes_trial = spikes_neruon((spikes_neruon>window_start_context) & (spikes_neruon<window_end_context)) - event(trial);
       firing_rate_trials(trial, :) = histcounts(spikes_trial, bins)*1000/bin_width;
       plot(time, firing_rate_trials(trial, :))
       scatter(spikes_trial, repelem(trial, length(spikes_trial)), 2)
       hold on
       xline(0 , '--')
    end
    av_firing_rate_neuron_ctx_context{neuron} = mean(firing_rate_trials, 1);
    sem_firing_rate_neuron_ctx_context{neuron} = std(firing_rate_trials, 1)/sqrt(n_trials);



    time_av_before = (before_context -  time_av_window_before)/(bin_width/1000);
    time_av_after = (before_context +  time_av_window_after)/(bin_width/1000);
    time_av_firing_rate_neuron_ctx_context(neuron) = mean(av_firing_rate_neuron_ctx_context{neuron}(time_av_before:time_av_after), 2);



    saveas(gcf, fullfile(path, 'plots', sprintf('%s_ctx_neuron_%d_spike_psth_Context_interact_firing rate.png',save_tag, neuron)))
    hold off
end


av_firing_rate_neuron_vhc_context = {};
for neuron = 1:length(good_unit_ids_vhc)
    spikes_neruon = neurons_vhc{neuron};
     figure
    firing_rate_trials = zeros(n_trials, length(time));
    for trial = 1:n_trials
       window_start_context = event(trial) - abs(before_context);
       window_end_context = event(trial) + abs(after_context);
       spikes_trial = spikes_neruon((spikes_neruon>window_start_context) & (spikes_neruon<window_end_context)) - event(trial);
       firing_rate_trials(trial, :) = histcounts(spikes_trial, bins)*1000/bin_width;
       plot(time, firing_rate_trials(trial, :))
       scatter(spikes_trial, repelem(trial, length(spikes_trial)), 2)
       hold on
       xline(0 , '--')
    end
    av_firing_rate_neuron_vhc_context{neuron} = mean(firing_rate_trials, 1);
    sem_firing_rate_neuron_vhc_context{neuron} = std(firing_rate_trials, 1)/sqrt(n_trials);


    time_av_before = (before_context -  time_av_window_before)/(bin_width/1000);
    time_av_after = (before_context +  time_av_window_after)/(bin_width/1000);
    time_av_firing_rate_neuron_vhc_context(neuron) = mean(av_firing_rate_neuron_vhc_context{neuron}(time_av_before:time_av_after), 2);
    
    saveas(gcf, fullfile(path, 'plots', sprintf('%s_ctx_neuron_%d_spike_psth_Context_interact_firing rate.png',save_tag, neuron)))
    hold off
end

%% Summary plotting of obj approach firing rates :: Obj1 vs Obj2
before = 4;
after = 4;   %in seconds
sample_rate = 30000; %in Hz
n_trials = length(event);
bin_width = 50;  %ms
bins = -abs(before):bin_width/1000:after;
time = bins(1: end-1) + bin_width/1000;
for neuron = 1: length(good_unit_ids_ctx)
    figure
    plot(time, av_firing_rate_neuron_obj_ctx{neuron}, 'g','DisplayName','Object approach')
    hold on 
    patch([time fliplr(time)], [(av_firing_rate_neuron_obj_ctx{neuron}-sem_firing_rate_neuron_obj_ctx{neuron}) fliplr((av_firing_rate_neuron_obj_ctx{neuron}+sem_firing_rate_neuron_obj_ctx{neuron}))],'g', 'FaceAlpha',.3, 'EdgeAlpha',0, 'DisplayName','' );
    
    %plot(time, av_firing_rate_neuron_obj2_ctx{neuron}, 'r','DisplayName','Novel approach')
    %patch([time fliplr(time)], [(av_firing_rate_neuron_obj2_ctx{neuron}-sem_firing_rate_neuron_obj2_ctx{neuron}) fliplr((av_firing_rate_neuron_obj2_ctx{neuron}+sem_firing_rate_neuron_obj2_ctx{neuron}))],'r', 'FaceAlpha',.3, 'EdgeAlpha',0 , 'DisplayName','');
    xline(0,'DisplayName','Object contact')
    xlabel('Time(s)')
    ylabel('Firing rate (Z-score)')
    hold off
    legend
    saveas(gcf, fullfile(path, 'plots', sprintf('ctx_neuron_%d_spike_psth_object_approach_firing_rate.png', neuron)))
end
for neuron = 1: length(good_unit_ids_vhc)
    figure
    plot(time, av_firing_rate_neuron_obj_vhc{neuron}, 'g','DisplayName','Object approach')
    hold on 
    patch([time fliplr(time)], [(av_firing_rate_neuron_obj_vhc{neuron}-sem_firing_rate_neuron_obj_vhc{neuron}) fliplr((av_firing_rate_neuron_obj_vhc{neuron}+sem_firing_rate_neuron_obj_vhc{neuron}))],'g', 'FaceAlpha',.3, 'EdgeAlpha',0, 'DisplayName','' );
    
    %plot(time, av_firing_rate_neuron_obj2_vhc{neuron}, 'r','DisplayName','Novel approach')
    %patch([time fliplr(time)], [(av_firing_rate_neuron_obj2_vhc{neuron}-sem_firing_rate_neuron_obj2_vhc{neuron}) fliplr((av_firing_rate_neuron_obj2_vhc{neuron}+sem_firing_rate_neuron_obj2_vhc{neuron}))],'r', 'FaceAlpha',.3, 'EdgeAlpha',0 , 'DisplayName','');
    xline(0,'DisplayName','Object contact')
    xlabel('Time(s)')
    ylabel('Firing rate (Z-score)')
    hold off
    legend
    saveas(gcf, fullfile(path, 'plots', sprintf('vhc_neuron_%d_spike_psth_object_approach_firing_rate.png', neuron)))
end

%% Summary plotting of obj approach firing rates :: Obj1 vs Nov Obj
before = 1;
after = 1;   %in seconds
sample_rate = 30000; %in Hz
n_trials = length(event);
bin_width = 50;  %ms
bins = -abs(before):bin_width/1000:after;
time = bins(1: end-1) + bin_width/1000;
for neuron = 1: length(good_unit_ids_ctx)
    figure
    plot(time, av_firing_rate_neuron_obj1_ctx{neuron}, 'g','DisplayName','Object1 approach')
    hold on 
    patch([time fliplr(time)], [(av_firing_rate_neuron_obj1_ctx{neuron}-sem_firing_rate_neuron_obj1_ctx{neuron}) fliplr((av_firing_rate_neuron_obj1_ctx{neuron}+sem_firing_rate_neuron_obj1_ctx{neuron}))],'g', 'FaceAlpha',.3, 'EdgeAlpha',0, 'DisplayName','' );
    
    plot(time, av_firing_rate_neuron_obj2_ctx{neuron}, 'g','DisplayName','Novel approach')
    patch([time fliplr(time)], [(av_firing_rate_neuron_obj2_ctx{neuron}-sem_firing_rate_neuron_nov_ctx{neuron}) fliplr((av_firing_rate_neuron_nov_ctx{neuron}+sem_firing_rate_neuron_nov_ctx{neuron}))],'r', 'FaceAlpha',.3, 'EdgeAlpha',0 , 'DisplayName','');
    plot(time, av_firing_rate_neuron_obj_ctx{neuron}, 'black','DisplayName','Object approach')
    patch([time fliplr(time)], [(av_firing_rate_neuron_obj_ctx{neuron}-sem_firing_rate_neuron_obj_ctx{neuron}) fliplr((av_firing_rate_neuron_obj_ctx{neuron}+sem_firing_rate_neuron_obj_ctx{neuron}))],'black', 'FaceAlpha',.3, 'EdgeAlpha',0 , 'DisplayName','');
    xline(0)
    hold off
    legend
    saveas(gcf, fullfile(path, 'plots', sprintf('ctx_neuron_%d_spike_psth_obj_approach_firing rate.png', neuron)))
end
for neuron = 1: length(good_unit_ids_vhc)
    figure
    plot(time, av_firing_rate_neuron_obj1_vhc{neuron}, 'g','DisplayName','Object1 approach')
    hold on 
    patch([time fliplr(time)], [(av_firing_rate_neuron_obj1_vhc{neuron}-sem_firing_rate_neuron_obj1_vhc{neuron}) fliplr((av_firing_rate_neuron_obj1_vhc{neuron}+sem_firing_rate_neuron_obj1_vhc{neuron}))],'g', 'FaceAlpha',.3, 'EdgeAlpha',0, 'DisplayName','' );
    
    plot(time, av_firing_rate_neuron_obj2_vhc{neuron}, 'g','DisplayName','Novel Object approach')
    patch([time fliplr(time)], [(av_firing_rate_neuron_obj2_vhc{neuron}-sem_firing_rate_neuron_nov_vhc{neuron}) fliplr((av_firing_rate_neuron_nov_vhc{neuron}+sem_firing_rate_neuron_nov_vhc{neuron}))],'r', 'FaceAlpha',.3, 'EdgeAlpha',0 , 'DisplayName','');
    plot(time, av_firing_rate_neuron_obj_vhc{neuron}, 'black','DisplayName','Object approach')
    patch([time fliplr(time)], [(av_firing_rate_neuron_obj_vhc{neuron}-sem_firing_rate_neuron_obj_vhc{neuron}) fliplr((av_firing_rate_neuron_obj_vhc{neuron}+sem_firing_rate_neuron_obj_vhc{neuron}))],'black', 'FaceAlpha',.3, 'EdgeAlpha',0 , 'DisplayName','');
    xline(0)
    hold off
    legend
    saveas(gcf, fullfile(path, 'plots', sprintf('vhc_neuron_%d_spike_psth_obj_approach_firing rate.png', neuron)))
end

%% %% Summary plotting of obj approach firing rates :: Obj1 vs Context
before = 1;
after = 1;   %in seconds
sample_rate = 30000; %in Hz
n_trials = length(event);
bin_width = 50;  %ms
bins = -abs(before):bin_width/1000:after;
time = bins(1: end-1) + bin_width/1000;
for neuron = 1: length(good_unit_ids_ctx)
    figure
    plot(time, av_firing_rate_neuron_obj_ctx{neuron}, 'g','DisplayName','Object interaction')
    hold on 
    patch([time fliplr(time)], [(av_firing_rate_neuron_obj_ctx{neuron}-sem_firing_rate_neuron_obj_ctx{neuron}) fliplr((av_firing_rate_neuron_obj_ctx{neuron}+sem_firing_rate_neuron_obj_ctx{neuron}))],'g', 'FaceAlpha',.3, 'EdgeAlpha',0, 'DisplayName','' );
    
    %plot(time, av_firing_rate_neuron_ctx_context{neuron}, 'r','DisplayName','Context interaction')
    %patch([time fliplr(time)], [(av_firing_rate_neuron_ctx_context{neuron}-sem_firing_rate_neuron_ctx_context{neuron}) fliplr((av_firing_rate_neuron_ctx_context{neuron}+sem_firing_rate_neuron_ctx_context{neuron}))],'r', 'FaceAlpha',.3, 'EdgeAlpha',0 , 'DisplayName','');
    xline(0,'DisplayName','Object contact')
    xlabel('Time(s)')
    ylabel('Firing rate (Z-score)')
    hold off
    legend
    saveas(gcf, fullfile(path, 'plots', sprintf('ctx_neuron_%d_spike_psth_object_firing_rate.png', neuron)))
end
for neuron = 1: length(good_unit_ids_vhc)
    figure
    %plot(time, av_firing_rate_neuron_obj_vhc{neuron}, 'g','DisplayName','Object interaction')
    hold on 
    %patch([time fliplr(time)], [(av_firing_rate_neuron_obj_vhc{neuron}-sem_firing_rate_neuron_obj_vhc{neuron}) fliplr((av_firing_rate_neuron_obj_vhc{neuron}+sem_firing_rate_neuron_obj_vhc{neuron}))],'g', 'FaceAlpha',.3, 'EdgeAlpha',0, 'DisplayName','' );
    
    plot(time, av_firing_rate_neuron_vhc_context{neuron}, 'r','DisplayName','Context interaction')
    patch([time fliplr(time)], [(av_firing_rate_neuron_vhc_context{neuron}-sem_firing_rate_neuron_vhc_context{neuron}) fliplr((av_firing_rate_neuron_vhc_context{neuron}+sem_firing_rate_neuron_vhc_context{neuron}))],'r', 'FaceAlpha',.3, 'EdgeAlpha',0 , 'DisplayName','');
    xline(0,'DisplayName','Object contact')
    xlabel('Time(s)')
    ylabel('Firing rate (Z-score)')
    hold off
    legend
    saveas(gcf, fullfile(path, 'plots', sprintf('vhc_neuron_%d_spike_psth_object_firing_rate.png', neuron)))
end


%% Compare average firing rates of neurons
len_home1 = home1_stop - home1_start;
len_FOR = FOR_stop - FOR_start;
len_home2 = home2_stop - home2_start;
for neuron = 1:length(neurons_ctx)
    spikes_temp = neurons_ctx{neuron};
    av_fr_ctx_home1(neuron) =  length(spikes_temp((spikes_temp>home1_start) & (spikes_temp<home1_stop)))/len_home1;
end
for neuron = 1:length(neurons_ctx)
    spikes_temp = neurons_ctx{neuron};
    av_fr_ctx_home2(neuron) =  length(spikes_temp((spikes_temp>home2_start) & (spikes_temp<home2_stop)))/len_home2;
end
for neuron = 1:length(neurons_ctx)
    spikes_temp = neurons_ctx{neuron};
    av_fr_ctx_FOR(neuron) =  length(spikes_temp((spikes_temp>FOR_start) & (spikes_temp<FOR_stop)))/len_FOR;
end

for neuron = 1:length(neurons_vhc)
    spikes_temp = neurons_vhc{neuron};
    av_fr_vhc_home1(neuron) =  length(spikes_temp((spikes_temp>home1_start) & (spikes_temp<home1_stop)))/len_home1;
end
for neuron = 1:length(neurons_vhc)
    spikes_temp = neurons_vhc{neuron};
    av_fr_vhc_home2(neuron) =  length(spikes_temp((spikes_temp>home2_start) & (spikes_temp<home2_stop)))/len_home2;
end
for neuron = 1:length(neurons_vhc)
    spikes_temp = neurons_vhc{neuron};
    av_fr_vhc_FOR(neuron) =  length(spikes_temp((spikes_temp>FOR_start) & (spikes_temp<FOR_stop)))/len_FOR;
end

map = brewermap(3,'Set1'); 

figure
histogram(av_fr_ctx_home1, bins_av_firingrate,'Normalization','probability', 'facealpha',.5,'edgecolor','none', 'facecolor',map(1,:))
hold on
histogram(av_fr_ctx_FOR, bins_av_firingrate,'Normalization','probability', 'facealpha',.5,'edgecolor','none', 'facecolor',map(2,:))
histogram(av_fr_ctx_home2, bins_av_firingrate,'Normalization','probability', 'facealpha',.5,'edgecolor','none','facecolor',map(3,:))
xlabel('firing rate')
ylabel('number of cortical neurons')
title('Distribution of firing rates of cortical neurons')
legend('Home1','FOR','Home2')

hold off

saveas(gcf, fullfile(path, 'plots', sprintf('Firing_rate_distribution_context_CTX.png')))


figure
histogram(av_fr_vhc_home1, bins_av_firingrate,'Normalization','probability', 'facealpha',.5,'edgecolor','none', 'facecolor',map(1,:))
hold on
histogram(av_fr_vhc_FOR, bins_av_firingrate,'Normalization','probability', 'facealpha',.5,'edgecolor','none', 'facecolor',map(2,:))
histogram(av_fr_vhc_home2, bins_av_firingrate,'Normalization','probability', 'facealpha',.5,'edgecolor','none', 'facecolor',map(3,:))
xlabel('firing rate')
ylabel('number of VHC neurons')
title('Distribution of firing rates of VHC neurons')
legend('Home1','FOR','Home2')
hold off
saveas(gcf, fullfile(path, 'plots', sprintf('Firing_rate_distribution_context_VHC.png')))


%% 
%%scatter plot  
home1_x = repmat( "home1" ,1, length(av_fr_ctx_home1));
home2_x = repmat( "home2" ,1, length(av_fr_ctx_home2));
FOR_x = repmat( "FOR" ,1, length(av_fr_ctx_FOR));
scatter(categorical(home1_x), av_fr_ctx_home1, [], map(1,:))
hold on 
scatter(categorical(home2_x), av_fr_ctx_home2, [], map(2,:))
scatter(categorical(FOR_x), av_fr_ctx_FOR, [], map(3,:))
barx = categorical(["home1", "home2", "FOR"]);
bary = [mean(av_fr_ctx_home1), mean(av_fr_ctx_home2), mean(av_fr_ctx_FOR)];
bar(barx, bary,  'facealpha',.5,'edgecolor','none')
xlabel('Context')
ylabel('Firing rate in Hz')
title('firing rates of CTX neurons in contexts')

saveas(gcf, fullfile(path, 'plots', sprintf('Firing_rate_context_ctx.png')))
hold off

figure
home1_x = repmat( "home1" , length(av_fr_vhc_home1));
home2_x = repmat( "home2" , length(av_fr_vhc_home2));
FOR_x = repmat( "FOR" , length(av_fr_vhc_FOR));
scatter(categorical(home1_x), av_fr_vhc_home1, [], map(1,:))
hold on 
scatter(categorical(home2_x), av_fr_vhc_home2, [], map(2,:))
scatter(categorical(FOR_x), av_fr_vhc_FOR, [], map(3,:))
barx = categorical(["home1", "home2", "FOR"]);
bary = [mean(av_fr_vhc_home1), mean(av_fr_vhc_home2), mean(av_fr_vhc_FOR)];
bar(barx, bary, 'facealpha',.5,'edgecolor','none')
xlabel('Context')
ylabel('Firing rate in Hz')
title('firing rates of VHC neurons in contexts')


saveas(gcf, fullfile(path, 'plots', sprintf('Firing_rate_context_VHC.png')))
hold off

%% Make summary table with the time average firing rates of neurons

time_ave_fr_table_ctx= table(av_fr_ctx_home1, av_fr_ctx_home2,av_fr_ctx_FOR, time_av_firing_rate_neuron_obj_ctx);
time_ave_fr_table_vhc= table(av_fr_vhc_home1, av_fr_vhc_home2,av_fr_vhc_FOR, time_av_firing_rate_neuron_obj_vhc);

writetable(time_ave_fr_table_ctx, fullfile(path, 'Average_neuron_firing_rate_ctx.csv'));
writetable(time_ave_fr_table_vhc, fullfile(path, 'Average_neuron_firing_rate_vhc.csv'));
%% LFP .lfp generated by neuroscope at 1250Hz
lfp_file = dir(fullfile(path, "*.lfp")).name;
lfp = bz_LoadBinary((fullfile(path, lfp_file)));
lfp_ch = reshape(lfp, 64, []);
lfp_ctx = lfp_ch(1:32,:);
lfp_vhc = lfp_ch(33:64,:);
t_lfp = downsample(t,24);   %downsampled by factor of 24 because 30000/1250 =24



disp('first 32 channels assigned to CTX, next 32 assigned to VHC')


%% curation of ctx raw volatage channels 
sample_range = 1:1250;
for ii = 1:size(lfp_ctx, 1)
    disp('Processing channel:')
    disp(ii)
    figure
    toPlot = lfp_ctx(ii, sample_range);
    plot(sample_range/1250, toPlot)

    title(ii)
    disp('Visualising the CTX channels, proceed to exclusion if needed')
end
%% Remove noisy channels in ctx data

exc_ctx = [9, 16];
lfp_ctx(exc_ctx, :) = [];
disp('Following channels excluded from CTX probe:')
disp(exc_ctx)

%%  curation of vhc raw volatage channels 

for ii = 1:size(lfp_vhc, 1)
    disp('Processing channel:')
    disp(ii)
    figure
    toPlot = lfp_vhc(ii, sample_range);
    plot(sample_range/30000, toPlot)

    title(ii)
    disp('Visualising the VHC channels, proceed to exclusion if needed')

end
%% Remove noisey channels in vhc data

exc_vhc = [7, 11];
lfp_vhc(exc_vhc, :) = [];
disp('Following channels excluded from vhc probe:')
disp(exc_vhc)

%% Common median referencing of data after z-scoring (based on Lemke et. al 2019)
sample_rate_LFP = 1250; %in Hz
lfp_ctx = double(lfp_ctx');
lfp_vhc = double(lfp_vhc');
ctx_v_z = zscore(lfp_ctx);
vhc_v_z = zscore(lfp_vhc);
disp('Z-scored all data, proceeding to common mode referencing')
ctx_v_CMR = ctx_v_z- median(ctx_v_z, 1);
vhc_v_CMR =  vhc_v_z - median(vhc_v_z, 1); 
ctx_v_CMR = ctx_v_CMR';
vhc_v_CMR = vhc_v_CMR';
lfp_ctx = lfp_ctx';
lfp_vhc = lfp_vhc';

disp('common mode referencing complete')

%% Split LFP into bands 

%theta 3 - 10 Hz
CTX_lfp.theta = eegfilt(ctx_v_CMR,sample_rate_LFP, 3, 10, 0, 0, 0, 'fir1', 0);
%beta 10-30Hz
CTX_lfp.beta = eegfilt(ctx_v_CMR,sample_rate_LFP, 10, 30,  0, 0, 0, 'fir1', 0);
%slow gamma 30-55Hz
CTX_lfp.s_gamma = eegfilt(ctx_v_CMR,sample_rate_LFP, 30, 55,  0, 0, 0, 'fir1', 0);
%fast gamma 55-100Hz
CTX_lfp.f_gamma = eegfilt(ctx_v_CMR,sample_rate_LFP, 55, 100,  0, 0, 0, 'fir1', 0);
CTX_lfp.lfp = ctx_v_CMR;

%theta 3 - 10 Hz
VHC_lfp.theta = eegfilt(vhc_v_CMR,sample_rate_LFP, 3, 10,  0, 0, 0, 'fir1', 0);
%beta 10-30Hz
VHC_lfp.beta = eegfilt(vhc_v_CMR,sample_rate_LFP, 10, 30,  0, 0, 0, 'fir1', 0);
%slow gamma 30-55Hz
VHC_lfp.s_gamma = eegfilt(vhc_v_CMR,sample_rate_LFP, 30, 55,  0, 0, 0, 'fir1', 0);
%fast gamma 55-100Hz
VHC_lfp.f_gamma = eegfilt(vhc_v_CMR,sample_rate_LFP, 55, 100,  0, 0, 0, 'fir1', 0);
VHC_lfp.lfp = vhc_v_CMR;


%% Visualisation of the data after CMR (CTX)
for ii = 1:size(ctx_v_CMR, 1)
    disp('Processing channel:')
    disp(ii)
    figure
    toPlot = ctx_v_CMR(ii, sample_range);
    plot(sample_range/1250, toPlot, 'Color', [0.6350 0.0780 0.1840], 'DisplayName','CMR')
    hold on
    toPlot = CTX_lfp.theta(ii, sample_range);
    plot(sample_range/1250, toPlot, 'Color', [0.4940 0.1840 0.5560],  'DisplayName','theta')
    toPlot = CTX_lfp.beta(ii, sample_range);
    plot(sample_range/1250, toPlot, 'Color', [0.9290 0.6940 0.1250],  'DisplayName','beta')
    toPlot = CTX_lfp.s_gamma(ii, sample_range);
    plot(sample_range/1250, toPlot, 'Color', [0.4660 0.6740 0.1880],  'DisplayName','slow gamma')
    toPlot = CTX_lfp.f_gamma(ii, sample_range);
    plot(sample_range/1250, toPlot, 'Color', [0.3010 0.7450 0.9330],  'DisplayName','fast gamma')
    title("CTX CMR" + ii)
    legend()
    %figure
    %toPlot = ctx_v_CMR_ds(time_ds, ii);
    %plot(time_ds/1000, toPlot, 'g')
    %title("CTX CMR ds" + ii)
    %disp('Visualising the CTX channels after Common mode referencing, sanity check')
    saveas(gcf, fullfile(path, 'plots', sprintf('%s_Ctx_channel_%d_spectrum.png', save_tag, ii )))
end


%% Visualise theta in CTX

for ii = 1:size(CTX_lfp.theta, 1)
    disp('Processing channel:')
    disp(ii)
    figure
    toPlot = CTX_lfp.theta(ii, sample_range);
    plot(sample_range/1250, toPlot, 'r')
    title("CTX theta" + ii)
    figure
    toPlot = ctx_v_CMR_ds(time_ds, ii);
    plot(time_ds/1000, toPlot, 'g')
    title("CTX CMR ds" + ii)
    disp('Visualising the CTX channels after Common mode referencing, sanity check')
    
end


%% Visualisation of the data after CMR (VHC)

for ii = 1:size(vhc_v_CMR, 1)
    disp('Processing channel:')
    disp(ii)
    figure
    toPlot = vhc_v_CMR(ii, sample_range);
    plot(sample_range/1250, toPlot, 'Color', [0.6350 0.0780 0.1840], 'DisplayName','CMR')
    hold on
    toPlot = VHC_lfp.theta(ii, sample_range);
    plot(sample_range/1250, toPlot, 'Color', [0.4940 0.1840 0.5560],  'DisplayName','theta')
    toPlot = VHC_lfp.beta(ii, sample_range);
    plot(sample_range/1250, toPlot, 'Color', [0.9290 0.6940 0.1250],  'DisplayName','beta')
    toPlot = VHC_lfp.s_gamma(ii, sample_range);
    plot(sample_range/1250, toPlot, 'Color', [0.4660 0.6740 0.1880],  'DisplayName','s gamma')
    toPlot = VHC_lfp.f_gamma(ii, sample_range);
    plot(sample_range/1250, toPlot, 'Color', [0.3010 0.7450 0.9330],  'DisplayName','f gamma')
    title("VHC CMR" + ii)
    legend()
    %figure
    %toPlot = vhc_v_CMR_ds(time_ds, ii);
    %plot(time_ds/1000, toPlot, 'g')
    %title("VHC CMR ds" + ii)
    disp('Visualising the VHC channels after Common mode referencing, sanity check')
    
end
%     hold on
%     toPlot = vhc_v_CMR_ds(time_ds, ii);
%     plot(time_ds/1000, toPlot, 'g')
%     title("VHC CMR ds" + ii)
%     disp('Visualising the VHC channels after Common mode referencing, sanity check')
    %%%downsampling is creating a shift, will try to use non downsapled and
    %%%see how it goes and if needed will adapt accordingly


%% Visualise theta in vHC

for ii = 1:size(VHC_lfp.theta, 1)
    disp('Processing channel:')
    disp(ii)
    figure
    toPlot = VHC_lfp.theta(ii, sample_range);
    plot(sample_range/1250, toPlot, 'r')
    title("VHC theta" + ii)
    %figure
    %toPlot = vhc_v_CMR_ds(time_ds, ii);
    %plot(time_ds/1000, toPlot, 'g')
    %title("VHC CMR ds" + ii)
    disp('Visualising the VHC channels after Common mode referencing, sanity check')
    
end



%% LFP in home cage vs context vs homecage
% spectopo(ctx_v_CMR_ds', 0, 1000); 
%make windows of 5 seconds with 5 seconds before and 5 seconds after
% Let's start with obj approach, irresptive of identity
before = 4;  %in seconds
after = 4;   %in seconds

%% Run only for acq day:Extracting LFP and theta band activity at object1+2 approach
% If you want for day 1 and day 2 separately, say event = obj1_app_DLC or obj2_app_DLC 
% (this is also true for coherence)
event = obj_app_DLC;
save_tag = 'Object_approach'; %IMPORTANT; Change as necessary
n_trials = length(event);

empty_array = zeros((abs(before)+after)*sample_rate_LFP, n_trials); %% Will be populated later as (frames, trials) array
fprintf('Collecting trials of voltage for object approaches on following channels from - %d to + %d s wrt approach\n', abs(before), after)
for ch = 1:size(ctx_v_CMR,1)
   disp("CTX"+ ch)
   v_ctx_obj{ch} = empty_array;   %% intitialise the cell array element to be an empty array
   theta_ctx_lfp{ch} = empty_array;  
   gamma_ctx_lfp{ch} = empty_array;
   for ii = 1:n_trials
       window_start = event(ii) - abs(before);
       window_end = event(ii) + abs(after);
       v_ctx_obj{ch}(:, ii) = ctx_v_CMR(ch, (t_lfp>window_start) & (t_lfp<=window_end));
       theta_ctx_lfp{ch}(:, ii) = CTX_lfp.theta(ch, (t_lfp>window_start) & (t_lfp<=window_end));
       gamma_ctx_lfp{ch}(:, ii) = CTX_lfp.s_gamma(ch, (t_lfp>window_start) & (t_lfp<=window_end));
   end
end

for ch = 1:size(vhc_v_CMR,1)
   disp("VHC"+ch)
   theta_vhc_lfp{ch} = empty_array;
   gamma_vhc_lfp{ch} = empty_array;
   v_vhc_obj{ch} = empty_array;  %% intitialise the cell array element to be an empty array
   for ii = 1:n_trials
       window_start = event(ii) - abs(before);
       window_end = event(ii) + abs(after);
      
       v_vhc_obj{ch}(:, ii) = vhc_v_CMR(ch, (t_lfp>window_start) & (t_lfp<=window_end));
       theta_vhc_lfp{ch}(:, ii) = VHC_lfp.theta(ch, (t_lfp>window_start) & (t_lfp<=window_end));
       gamma_vhc_lfp{ch}(:, ii) = VHC_lfp.s_gamma(ch, (t_lfp>window_start) & (t_lfp<=window_end));
   end
end

disp('Done, proceed to plotting of the power spectra etc')

%% Now lets use the newtimef function from EEGlab (Lemke et al., 2019)
% 
% 
% before = 4;
% after = 4;

ersp_vhc = zeros(length(v_ctx_obj),1);
for ch = 1:length(v_ctx_obj)
   figure
   [ersp,itc,powbase,times,freqs,erspboot,itcboot] = newtimef(v_ctx_obj{ch}, (before+after)*sample_rate_LFP, [-(abs(before)*1000) after*1000], sample_rate_LFP, 'cycles', 0, 'maxfreq', 110, 'baseline', NaN);
   saveas(gcf, fullfile(path, 'plots', sprintf('ctx_channel_%d_ersp_cut_spectrum_%d_%d.png', ch, before, after)))
   ersp_cut = ersp(:, (times > 0) & (times < 800));
   ersp_baseline = ersp(:, (times > -4000) & (times < -2500));
   mean_esrp = mean(ersp_cut,2);
   mean_baseline = mean(ersp_baseline,2);
   baseline_subtracted = mean_esrp - mean_baseline;
   esrp_ctx = table(freqs' , mean_esrp, mean_baseline, baseline_subtracted);
   writetable(esrp_ctx,fullfile(path, sprintf('LFP_mean_ctx_%d.csv', ch)))
end
disp('Finished CTX channels')
%%ersp_ctx is (nfreqs,timesout)

for ch = 1:length(v_vhc_obj)
   figure
   [ersp,itc,powbase,times,freqs,erspboot,itcboot] = newtimef(v_vhc_obj{ch}, (before+after)*sample_rate_LFP, [-(abs(before)*1000) after*1000], sample_rate_LFP, 'cycles', 0, 'maxfreq', 110, 'baseline', NaN);
   saveas(gcf, fullfile(path, 'plots', sprintf('vhc_channel_%d_ersp_cut_spectrum_%d_%d.png', ch, before, after)))
   ersp_cut = ersp(:, (times > 0) & (times < 800));
   ersp_baseline = ersp(:, (times > -4000) & (times < -2500));
   mean_esrp = mean(ersp_cut,2);
   mean_baseline = mean(ersp_baseline,2);
   baseline_subtracted = mean_esrp - mean_baseline;
   esrp_vhc = table(freqs' , mean_esrp, mean_baseline, baseline_subtracted);
   writetable(esrp_vhc,fullfile(path, sprintf('LFP_mean_vhc_%d.csv', ch)))
end

disp('Finished VHC channels')
disp('Done')

%% %% Run only for context LFP, LFP averaged across trials
event = context_app_DLC;
save_tag = 'Context_approach'; 
n_trials = length(event);

empty_array = zeros((abs(before)+after)*sample_rate_LFP, n_trials); %% Will be populated later as (frames, trials) array
fprintf('Collecting trials of voltage for Context approaches on following channels from - %d to + %d s wrt approach\n', abs(before), after)
for ch = 1:size(ctx_v_CMR,1)
   disp("CTX"+ ch)
   v_ctx_cntxt{ch} = empty_array;   %% intitialise the cell array element to be an empty array
   theta_ctx_lfp{ch} = empty_array;  
   gamma_ctx_lfp{ch} = empty_array;
   for ii = 1:n_trials
       window_start = event(ii) - abs(before);
       window_end = event(ii) + abs(after);
       v_ctx_cntxt{ch}(:, ii) = ctx_v_CMR(ch, (t_lfp>window_start) & (t_lfp<=window_end));
       theta_ctx_lfp{ch}(:, ii) = CTX_lfp.theta(ch, (t_lfp>window_start) & (t_lfp<=window_end));
       gamma_ctx_lfp{ch}(:, ii) = CTX_lfp.s_gamma(ch, (t_lfp>window_start) & (t_lfp<=window_end));
   end
end

for ch = 1:size(vhc_v_CMR,1)
   disp("VHC"+ch)
   theta_vhc_lfp{ch} = empty_array;
   gamma_vhc_lfp{ch} = empty_array;
   v_vhc_cntxt{ch} = empty_array;  %% intitialise the cell array element to be an empty array
   for ii = 1:n_trials
       window_start = event(ii) - abs(before);
       window_end = event(ii) + abs(after);
       v_vhc_cntxt{ch}(:, ii) = vhc_v_CMR(ch, (t_lfp>window_start) & (t_lfp<=window_end));
       theta_vhc_lfp{ch}(:, ii) = VHC_lfp.theta(ch, (t_lfp>window_start) & (t_lfp<=window_end));
       gamma_vhc_lfp{ch}(:, ii) = VHC_lfp.s_gamma(ch, (t_lfp>window_start) & (t_lfp<=window_end));
   end
end

disp('Done, proceed to plotting of the power spectra etc')

%% Plot Context LFP

%%ersp_vhc = zeros(length(v_ctx_cntxt),1);
for ch = 1:length(v_ctx_cntxt)
   figure
   [ersp,itc,powbase,times,freqs,erspboot,itcboot] = newtimef(v_ctx_cntxt{ch}, (before+after)*sample_rate_LFP, [-(abs(before)*1000) after*1000], sample_rate_LFP, 'cycles', 0, 'maxfreq', 110, 'baseline', NaN);
   saveas(gcf, fullfile(path, 'plots', sprintf('ctx_channel_%d_ersp_cut_Context_spectrum_%d_%d.png', ch, before, after)))
   ersp_cut = ersp(:, (times > 0) & (times < 800));
   %ersp_baseline = ersp(:, (times > -4000) & (times < -3000));
   mean_esrp = mean(ersp_cut,2);
   %mean_baseline = mean(ersp_baseline,2);
   %baseline_subtracted = mean_esrp - mean_baseline;
   esrp_ctx = table(freqs' , mean_esrp, mean_baseline, baseline_subtracted);
   writetable(esrp_ctx,fullfile(path, sprintf('ContextLFP_mean_ctx_%d.csv', ch)))
end
disp('Finished CTX channels')
%%ersp_ctx is (nfreqs,timesout)

for ch = 1:length(v_vhc_cntxt)
   figure
   [ersp,itc,powbase,times,freqs,erspboot,itcboot] = newtimef(v_vhc_cntxt{ch}, (before+after)*sample_rate_LFP, [-(abs(before)*1000) after*1000], sample_rate_LFP, 'cycles', 0, 'maxfreq', 110, 'baseline', NaN);
   saveas(gcf, fullfile(path, 'plots', sprintf('vhc_channel_%d_ersp_cut_Context_spectrum_%d_%d.png', ch, before, after)))
   ersp_cut = ersp(:, (times > 0) & (times < 800));
   %ersp_baseline = ersp(:, (times > -4000) & (times < -3000));
   mean_esrp = mean(ersp_cut,2);
   %mean_baseline = mean(ersp_baseline,2);
   %baseline_subtracted = mean_esrp - mean_baseline;
   esrp_vhc = table(freqs' , mean_esrp, mean_baseline, baseline_subtracted);
   writetable(esrp_vhc,fullfile(path, sprintf('ContextLFP_mean_vhc_%d.csv', ch)))
end

disp('Finished VHC channels')
disp('Done')


%% Run only for recall: Calculating power for Object1 on Recall day
 
event = obj1_app_DLC;
n_trials = length(obj1_app_DLC);
save_tag = 'Object_approach';

empty_array = zeros((abs(before)+after)*sample_rate_LFP, n_trials); %% Will be populated later as (frames, trials) array
fprintf('Collecting trials of voltage for object approaches on following channels from - %d to + %d s wrt approach\n', abs(before), after)
for ch = 1:size(ctx_v_CMR,1)
   disp("CTX"+ ch)
   v_ctx_obj{ch} = empty_array;   %% intitialise the cell array element to be an empty array
   theta_ctx_lfp{ch} = empty_array;  
   gamma_ctx_lfp{ch} = empty_array;
   for ii = 1:n_trials
       window_start = event(ii) - abs(before);
       window_end = event(ii) + abs(after);
       v_ctx_obj{ch}(:, ii) = ctx_v_CMR(ch, (t_lfp>window_start) & (t_lfp<=window_end));
       theta_ctx_lfp{ch}(:, ii) = CTX_lfp.theta(ch, (t_lfp>window_start) & (t_lfp<=window_end));
       gamma_ctx_lfp{ch}(:, ii) = CTX_lfp.s_gamma(ch, (t_lfp>window_start) & (t_lfp<=window_end));
   end
end

for ch = 1:size(vhc_v_CMR,1)
   disp("VHC"+ch)
   theta_vhc_lfp{ch} = empty_array;
   gamma_vhc_lfp{ch} = empty_array;
   v_vhc_obj{ch} = empty_array;  %% intitialise the cell array element to be an empty array
   for ii = 1:n_trials
       window_start = event(ii) - abs(before);
       window_end = event(ii) + abs(after);
      
       v_vhc_obj{ch}(:, ii) = vhc_v_CMR(ch, (t_lfp>window_start) & (t_lfp<=window_end));
       theta_vhc_lfp{ch}(:, ii) = VHC_lfp.theta(ch, (t_lfp>window_start) & (t_lfp<=window_end));
       gamma_vhc_lfp{ch}(:, ii) = VHC_lfp.s_gamma(ch, (t_lfp>window_start) & (t_lfp<=window_end));
   end
end

disp('Done, proceed to plotting of the power spectra etc')
%% Plot power spectra for Object1 on recall day 

parfor ch = 1:length(v_ctx_obj)
   figure
   newtimef(v_ctx_obj{ch}, (before+after)*sample_rate_LFP, [-(abs(before)*1000) after*1000], sample_rate_LFP, 'cycles', 0, 'maxfreq', 110, 'baseline', NaN);
   saveas(gcf, fullfile(path, 'plots', sprintf('ctx_channel_%d_spectrum_%d_%d_for_Object1.png', ch, before, after)))
end
disp('Finished CTX channels')
parfor ch = 1:length(v_vhc_obj)
   figure
   newtimef(v_vhc_obj{ch}, (before+after)*sample_rate_LFP, [-(abs(before)*1000) after*1000], sample_rate_LFP, 'cycles', 0, 'maxfreq', 110, 'baseline', NaN);
   saveas(gcf, fullfile(path, 'plots', sprintf('vhc_channel_%d_spectrum_%d_%d_for_object1.png', ch, before, after)))
end
disp('Finished VHC channels')
disp('Done')
%% Get mean lfp for object1
% before = 4;
% after = 4;

event = obj1_app_DLC;
ersp_vhc = zeros(length(v_ctx_obj),1);
for ch = 1:length(v_ctx_obj)
   figure
   [ersp,itc,powbase,times,freqs,erspboot,itcboot] = newtimef(v_ctx_obj{ch}, (before+after)*sample_rate_LFP, [-(abs(before)*1000) after*1000], sample_rate_LFP, 'cycles', 0, 'maxfreq', 110, 'baseline', NaN);
   saveas(gcf, fullfile(path, 'plots', sprintf('ctx_channel_%d_spectrum_%d_%d_Object1.png', ch, before, after)))
   ersp_cut = ersp(:, (times > 0) & (times < 800));
   ersp_baseline = ersp(:, (times > -4000) & (times < -3000));
   mean_esrp = mean(ersp_cut,2);
   mean_baseline = mean(ersp_baseline,2);
   baseline_subtracted = mean_esrp - mean_baseline;
   esrp_ctx = table(freqs' , mean_esrp, mean_baseline, baseline_subtracted);
   writetable(esrp_ctx,fullfile(path, sprintf('LFP_mean_ctx_%d_Object1_baseline+4:+5.csv', ch)))
end
disp('Finished CTX channels')
%ersp_ctx is (nfreqs,timesout)

for ch = 1:length(v_vhc_obj)
   figure
   [ersp,itc,powbase,times,freqs,erspboot,itcboot] = newtimef(v_vhc_obj{ch}, (before+after)*sample_rate_LFP, [-(abs(before)*1000) after*1000], sample_rate_LFP, 'cycles', 0, 'maxfreq', 110, 'baseline', NaN);
   saveas(gcf, fullfile(path, 'plots', sprintf('vhc_channel_%d_spectrum_%d_%d_Object1.png', ch, before, after)))
   ersp_cut = ersp(:, (times > 0) & (times < 800));
   ersp_baseline = ersp(:, (times > -4000) & (times < -3000));
   mean_esrp = mean(ersp_cut,2);
   mean_baseline = mean(ersp_baseline,2);
   baseline_subtracted = mean_esrp - mean_baseline;
   esrp_vhc = table(freqs' , mean_esrp, mean_baseline, baseline_subtracted);
   writetable(esrp_vhc,fullfile(path, sprintf('LFP_mean_vhc_%d_Object1_baseline+4:+5.csv', ch)))
end

disp('Finished VHC channels')
disp('Done')
%% %% run only for Recall day: Calculatting power for Object2 on Recall day
 
event = obj2_app_DLC;
n_trials = length(obj2_app_DLC);
save_tag = 'Object_approach';
empty_array = zeros((abs(before)+after)*sample_rate_LFP, n_trials); %% Will be populated later as (frames, trials) array
fprintf('Collecting trials of voltage for object approaches on following channels from - %d to + %d s wrt approach\n', abs(before), after)
for ch = 1:size(ctx_v_CMR,1)
   disp("CTX"+ ch)
   v_ctx_obj{ch} = empty_array;   %% intitialise the cell array element to be an empty array
   theta_ctx_lfp{ch} = empty_array;  
   gamma_ctx_lfp{ch} = empty_array;
   for ii = 1:n_trials
       window_start = event(ii) - abs(before);
       window_end = event(ii) + abs(after);
       v_ctx_obj{ch}(:, ii) = ctx_v_CMR(ch, (t_lfp>window_start) & (t_lfp<=window_end));
       theta_ctx_lfp{ch}(:, ii) = CTX_lfp.theta(ch, (t_lfp>window_start) & (t_lfp<=window_end));
       gamma_ctx_lfp{ch}(:, ii) = CTX_lfp.s_gamma(ch, (t_lfp>window_start) & (t_lfp<=window_end));
   end
end

for ch = 1:size(vhc_v_CMR,1)
   disp("VHC"+ch)
   theta_vhc_lfp{ch} = empty_array;
   gamma_vhc_lfp{ch} = empty_array;
   v_vhc_obj{ch} = empty_array;  %% intitialise the cell array element to be an empty array
   for ii = 1:n_trials
       window_start = event(ii) - abs(before);
       window_end = event(ii) + abs(after);
      
       v_vhc_obj{ch}(:, ii) = vhc_v_CMR(ch, (t_lfp>window_start) & (t_lfp<=window_end));
       theta_vhc_lfp{ch}(:, ii) = VHC_lfp.theta(ch, (t_lfp>window_start) & (t_lfp<=window_end));
       gamma_vhc_lfp{ch}(:, ii) = VHC_lfp.s_gamma(ch, (t_lfp>window_start) & (t_lfp<=window_end));
   end
end

disp('Done, proceed to plotting of the power spectra etc')
%% Plot power spectra for Object2 on recall day

parfor ch = 1:length(v_ctx_obj)
   figure
   newtimef(v_ctx_obj{ch}, (before+after)*sample_rate_LFP, [-(abs(before)*1000) after*1000], sample_rate_LFP, 'cycles', 0, 'maxfreq', 110, 'baseline', NaN);
   saveas(gcf, fullfile(path, 'plots', sprintf('ctx_channel_%d_spectrum_%d_%d_for_Object2.png', ch, before, after)))
end
disp('Finished CTX channels')
parfor ch = 1:length(v_vhc_obj)
   figure
   newtimef(v_vhc_obj{ch}, (before+after)*sample_rate_LFP, [-(abs(before)*1000) after*1000], sample_rate_LFP, 'cycles', 0, 'maxfreq', 110, 'baseline', NaN);
   saveas(gcf, fullfile(path, 'plots', sprintf('vhc_channel_%d_spectrum_%d_%d_for_Object2.png', ch, before, after)))
end
disp('Finished VHC channels')
disp('Done')

%% Get mean lfp for object2
% before = 4;
% after = 4;

event = obj2_app_DLC;
ersp_vhc = zeros(length(v_ctx_obj),1);
for ch = 1:length(v_ctx_obj)
   figure
   [ersp,itc,powbase,times,freqs,erspboot,itcboot] = newtimef(v_ctx_obj{ch}, (before+after)*sample_rate_LFP, [-(abs(before)*1000) after*1000], sample_rate_LFP, 'cycles', 0, 'maxfreq', 110, 'baseline', NaN);
   saveas(gcf, fullfile(path, 'plots', sprintf('ctx_channel_%d_spectrum_%d_%d_Object2.png', ch, before, after)))
   ersp_cut = ersp(:, (times > 0) & (times < 800));
   ersp_baseline = ersp(:, (times > -4000) & (times < -3000));
   mean_esrp = mean(ersp_cut,2);
   mean_baseline = mean(ersp_baseline,2);
   baseline_subtracted = mean_esrp - mean_baseline;
   esrp_ctx = table(freqs' , mean_esrp, mean_baseline, baseline_subtracted);
   writetable(esrp_ctx,fullfile(path, sprintf('LFP_mean_ctx_%d_Object2.csv', ch)))
end
disp('Finished CTX channels')
%ersp_ctx is (nfreqs,timesout)

for ch = 1:length(v_vhc_obj)
   figure
   [ersp,itc,powbase,times,freqs,erspboot,itcboot] = newtimef(v_vhc_obj{ch}, (before+after)*sample_rate_LFP, [-(abs(before)*1000) after*1000], sample_rate_LFP, 'cycles', 0, 'maxfreq', 110, 'baseline', NaN);
   saveas(gcf, fullfile(path, 'plots', sprintf('vhc_channel_%d_spectrum_%d_%d_Object2.png', ch, before, after)))
   ersp_cut = ersp(:, (times > 0) & (times < 800));
   ersp_baseline = ersp(:, (times > -4000) & (times < -3000));
   mean_esrp = mean(ersp_cut,2);
   mean_baseline = mean(ersp_baseline,2);
   baseline_subtracted = mean_esrp - mean_baseline;
   esrp_vhc = table(freqs' , mean_esrp, mean_baseline, baseline_subtracted);
   writetable(esrp_vhc,fullfile(path, sprintf('LFP_mean_vhc_%d_Object2.csv', ch)))
end

disp('Finished VHC channels')
disp('Done')

%% Power spectra with Baseline substracted
parfor ch = 1:length(v_ctx_obj)
   figure
   newtimef(v_ctx_obj{ch}, (before+after)*sample_rate_LFP, [-(abs(before)*1000) after*1000], sample_rate_LFP, 'cycles', 0, 'maxfreq', 100, 'baseline', [-4000 -3000]);
   saveas(gcf, fullfile(path, 'plots', sprintf('ctx_channel_%d_Object2_spectrum_%d_%d_baseline_subtracted(4-3s).png', ch, before, after)))
end
disp('Finished CTX channels')
parfor ch = 1:length(v_vhc_obj)
   figure
   newtimef(v_vhc_obj{ch}, (before+after)*sample_rate_LFP, [-(abs(before)*1000) after*1000], sample_rate_LFP, 'cycles', 0, 'maxfreq', 100, 'baseline', [-4000 -3000]);
   saveas(gcf, fullfile(path, 'plots', sprintf('vhc_channel_%d_Object2_spectrum_%d_%d_baseline_subtracted(4-3s).png', ch, before, after)))
end
disp('Finished VHC channels')
disp('Done')

%% %% Power spectra with Baseline substracted for Context LFP

parfor ch = 1:length(v_ctx_cntxt)
   figure
   newtimef(v_ctx_cntxt{ch}, (before+after)*sample_rate_LFP, [-(abs(before)*1000) after*1000], sample_rate_LFP, 'cycles', 0, 'maxfreq', 100, 'baseline', [-10000 -8500]);
   saveas(gcf, fullfile(path, 'plots', sprintf('ctx_channel_%d_Context_spectrum_%d_%d_baseline_subtracted.png', ch, before, after)))
end
disp('Finished CTX channels')
parfor ch = 1:length(v_vhc_cntxt)
   figure
   newtimef(v_vhc_cntxt{ch}, (before+after)*sample_rate_LFP, [-(abs(before)*1000) after*1000], sample_rate_LFP, 'cycles', 0, 'maxfreq', 100, 'baseline', [-10000 -8500]);
   saveas(gcf, fullfile(path, 'plots', sprintf('vhc_channel_%d_Context_spectrum_%d_%d_baseline_subtracted.png', ch, before, after)))
end
disp('Finished VHC channels')
disp('Done')

%% Plotting LFP band activity at trials
mkdir(fullfile(path, 'plots_raw_voltage_behavior_aligned'))

x = linspace(-abs(before),after,(before+after)*sample_rate_LFP);
for ch = 1:length(v_ctx_obj)
    figure
    lh = plot(x, v_ctx_obj{ch}, 'Color', [0.7 0.7 0.7 0.3], 'Linewidth', 1);
    %lh.Color(4) = 0.3;
    hold on
    plot(x, mean(v_ctx_obj{ch},2), 'black', 'Linewidth', 1.5)
    xline(0)
    saveas(gcf, fullfile(path, 'plots', sprintf('%s_ctx_channel_%d_voltage.png', save_tag, ch)))
    hold off
end

for ch = 1:length(v_vhc_obj)
    figure
    lh = plot(x, v_vhc_obj{ch},  'Color', [0.7 0.7 0.7 0.3], 'Linewidth', 1);
    %lh.Color(4) = 0.3;
    hold on
    plot(x, mean(v_vhc_obj{ch},2), 'black', 'Linewidth', 1.5)
    xline(0)
    saveas(gcf, fullfile(path, 'plots_raw_voltage_behavior_aligned', sprintf('%s_vhc_channel_%d_voltage.png',save_tag, ch)))
end


%% Plotting theta band activity at trials
mkdir(fullfile(path, 'plots_theta_behavior_aligned'))

x = linspace(-abs(before),after,(before+after)*sample_rate_LFP);
for ch = 1:length(theta_ctx_lfp)
    n_trials = size(v_ctx_obj{ch});
    n_trials = n_trials(2);
    figure
    lh = plot(x, theta_ctx_lfp{ch}, 'Color', [0.7 0.7 0.7 0.3], 'Linewidth', 1);
    %lh.Color(4) = 0.3;
    hold on
    plot(x, mean(theta_ctx_lfp{ch},2), 'black', 'Linewidth', 1.5)
    sem_theta = std(theta_ctx_lfp{ch}, 0, 2)/sqrt(n_trials);
    patch([x fliplr(x)], [(mean(theta_ctx_lfp{ch},2)-sem_theta)' fliplr((mean(theta_ctx_lfp{ch},2)+sem_theta)')],'r', 'FaceAlpha',.3, 'EdgeAlpha',0 , 'DisplayName','');
    xline(0)
    saveas(gcf, fullfile(path, 'plots_theta_behavior_aligned', sprintf('%s_ctx_channel_%d_theta.png', save_tag, ch)))
    hold off
end

for ch = 1:length(theta_vhc_lfp)
    figure
    lh = plot(x, theta_vhc_lfp{ch},  'Color', [0.7 0.7 0.7 0.3], 'Linewidth', 1);
    %lh.Color(4) = 0.3;lfp
    hold on
    plot(x, mean(theta_vhc_lfp{ch},2), 'black', 'Linewidth', 1.5)
    sem_theta = std(theta_vhc_lfp{ch},0,2)/sqrt(n_trials);
    patch([x fliplr(x)], [(mean(theta_vhc_lfp{ch},2)+sem_theta)' fliplr((mean(theta_vhc_lfp{ch},2)-sem_theta)')],'r', 'FaceAlpha',.3, 'EdgeAlpha',0 , 'DisplayName','');
    xline(0)
    saveas(gcf, fullfile(path, 'plots_theta_behavior_aligned', sprintf('%s_vhc_channel_%d_theta.png',save_tag, ch)))
end


%% Ploting Gamma

mkdir(fullfile(path, 'plots_gamma_behavior_aligned'))

x = linspace(-abs(before),after,(before+after)*sample_rate_LFP);
for ch = 1:length(gamma_ctx_lfp)
    n_trials = size(v_ctx_obj{ch});
    n_trials = n_trials(2);
    figure
    lh = plot(x, gamma_ctx_lfp{ch}, 'Color', [0.7 0.7 0.7 0.3], 'Linewidth', 1);
    %lh.Color(4) = 0.3;
    hold on
    plot(x, mean(gamma_ctx_lfp{ch},2), 'black', 'Linewidth', 1.5)
    sem_gamma = std(gamma_ctx_lfp{ch}, 0, 2)/sqrt(n_trials);
    patch([x fliplr(x)], [(mean(gamma_ctx_lfp{ch},2)-sem_gamma)' fliplr((mean(gamma_ctx_lfp{ch},2)+sem_gamma)')],'r', 'FaceAlpha',.3, 'EdgeAlpha',0 , 'DisplayName','');
    xline(0)
    saveas(gcf, fullfile(path, 'plots_gamma_behavior_aligned', sprintf('%s_ctx_channel_%d_gamma.png', save_tag, ch)))
    hold off
end

for ch = 1:length(gamma_vhc_lfp)
    figure
    lh = plot(x, gamma_vhc_lfp{ch},  'Color', [0.7 0.7 0.7 0.3], 'Linewidth', 1);
    %lh.Color(4) = 0.3;
    hold on
    plot(x, mean(gamma_vhc_lfp{ch},2), 'black', 'Linewidth', 1.5)
    sem_gamma = std(gamma_vhc_lfp{ch},0,2)/sqrt(n_trials);
    patch([x fliplr(x)], [(mean(gamma_vhc_lfp{ch},2)+sem_gamma)' fliplr((mean(gamma_vhc_lfp{ch},2)-sem_gamma)')],'r', 'FaceAlpha',.3, 'EdgeAlpha',0 , 'DisplayName','');
    xline(0)
    saveas(gcf, fullfile(path, 'plots_gamma_behavior_aligned', sprintf('%s_vhc_channel_%d_gamma.png',save_tag, ch)))
end


%% LFP Coherence between CTX and VHC using cohgramc from Chronux as used by Lemke et al 2019


movingwin = [0.5, 0.5];  %% moving window params 
params.tapers = [3 5]; %% scale (Nw) and number of tapers (K) k=2Nw-1 
params.pad = 0; %%amount of padding used by FFT
params.Fs = 1250;
params.fpass = [0 100];
params.err = 0;
params.trialave = 1;
coherence = {};
coherence_phi = {};
coherence_f = {};
% set window to observe average in as [before, after] for eg [0.5,0.5] in seconds:
average_window = [0.5, 0.0]; % in seconds      
%do one channel at a time to visualise heatmap
for ch_vhc = 7
     for ch_ctx = 23 
          data1 = v_ctx_obj{ch_ctx};
          data2 = v_vhc_obj{ch_vhc};
          [C,phi,S12,S1,S2,time_coherence,f_coherence] = cohgramc(data1,data2,movingwin,params);
          coherence = C; % C is of the form time x frequency  
          coherence_phi = phi;
          coherence_f = f_coherence;  
     end
     h = heatmap(time_coherence,f_coherence, coherence', 'GridVisible', 'off', 'Colormap', parula);
    
     time_labels = string(round((time_coherence-before),2));
     for ii =1:length(time_labels)
         if mod(ii,3) ~= 0
             time_labels(ii) = "";
         end
     end
    %xline(before, '--k')
    freq_labels = string(round(f_coherence,2))
    h.XDisplayLabels = time_labels;
    h.YDisplayLabels = freq_labels;
    h.XLabel = sprintf('time relative to closest point of interaction in %s (s)', save_tag);
    h.YLabel = 'Frequency(Hz)';
    saveas(gcf, fullfile(path, sprintf('%d_ctx_%d_vhc_%s_coherence_heatmap.png',ch_ctx, ch_vhc, save_tag)))
    
    
    figure
    coherence_in_window = coherence(abs(before) - abs(average_window(1)):abs(after) + abs(average_window(2)), :); % average_window x frequency
    mean_coherence = mean(coherence_in_window, 1);
    plot(f_coherence, mean_coherence);
    hold on
    xlabel('Frequency (Hz)')
    ylabel(sprintf('Average coherence in %d s to %d s around point of closest approach', average_window(1),average_window(2)))
    xline(3, '--k')
    xline(10, '--k')
    xline(30, '--k')
    xline(55, '--k')
    xline(100, '--k')
    saveas(gcf, fullfile(path, sprintf('%d_ctx_%d_vhc_%s_coherence_lineplot.png',ch_ctx, ch_vhc, save_tag)))
    hold off

    %if you only want to plot the coherence of one frequency as a line
    %plot, check the element number of that freq in f_coherence and
    %uncomment the lines below with that element number
    %element_number = 1;
    %plot(time_coherence, coherence(; , element_number))
end

%% Load raw voltage
% fileinfo = dir(fullfile(path,'amplifier.dat'));
% num_channels = length(amplifier_channels); % amplifier channel info from header file
% num_samples = fileinfo.bytes/(num_channels * 2);
% fid = fopen(fullfile(path,'amplifier.dat'), 'r');
% v_amp = fread(fid, [num_channels, num_samples], 'int16');
% fclose(fid);
% v_amp = v_amp * 0.195;
% 
% disp('Data has been loaded, proceed to extract camera timestamps')




%% function to retrive timestamps from video frames
% function frame2ts(vid_frames, camera_ts)
% frames
% for frame_id = 1:length(frames)
%     
% 
% end
