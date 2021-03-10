%% Clean start
clc;clear all;close all;

%% EEG SSVEP Processing Variables
FS = 512;                                       % sampling frequency
STI_F = ['s1+L','s2+L','s3+L','s4+L','s0+L','s1+R','s2+R','s3+R','s4+R','s0+R'];          
n_sti = 10;


%% For offline experiments exp 1, 15*5*5s

% number of trials for one stimulus frequency
Stimulus_N = 15;
% rest time in between
Rest_t = 1;
% window length
window = 8;
% initial_rest time in seconds
Stimulus_start_Time = 10;

% the sequence to be used 
Stimulus_frequencies_raw = [];
Stimulus_frequencies_easy = [];
for i = 1:Stimulus_start_Time
    Stimulus_frequencies_raw = [Stimulus_frequencies_raw,16];
end

Stimulus_space = 1:n_sti;

%step 1 conversion to array of stimulus frequency indexes
for i=1:Stimulus_N
    x=Stimulus_space(randperm(length(Stimulus_space)))
    for j=1:length(Stimulus_space)
        Stimulus_frequencies_easy = [Stimulus_frequencies_easy,x(j)];
        for k = 1:window
            Stimulus_frequencies_raw = [Stimulus_frequencies_raw,x(j)];
        end
        %adding a 1s rest time by assigning a 16
        for l = 1:Rest_t
            Stimulus_frequencies_raw = [Stimulus_frequencies_raw,16];
        end
    end
end
for i = 1:Stimulus_start_Time
    Stimulus_frequencies_raw = [Stimulus_frequencies_raw,16];
end
prompts_name = sprintf('prompts1_start_%dwindow_%drest_%d.mat',Stimulus_start_Time,window,Rest_t);
save(prompts_name,'Stimulus_frequencies_easy');
Stimulus_frequencies_raw = timeseries(Stimulus_frequencies_raw)