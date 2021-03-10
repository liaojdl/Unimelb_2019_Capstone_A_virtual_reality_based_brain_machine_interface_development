%% Clean start
clc
clear all
close all

%% CCA Processing Variables
Fs = 256;                       % sampling frequency
window_time = 4;                % length of processing window in seconds
window = window_time*Fs;        % length of processing window in samples
overlap = round(3.5*Fs);          % length of overlapping step in sample
sti_f = [10.4,11.4,12.4,13.4];             % stimulus frequencies
n_sti = length(sti_f);          % number of stimulus frequencies

%% FBCAA Processing Variables
H = 2;                          % number of harmonics in refsig
N = 9;                          % number of frequency bands 
Wn = [6,60];                    % entire usable frequency range in Hz
step_size = (Wn(2)-Wn(1))/N;    % each band's band width in Hz

% this is for storting the filter coefficients for the filter used
numx = zeros(9,13);             
denx = zeros(9,13);

% the weighting coefficients for combining the banks' CCA weighting
% check simulink block diagrams for explicit details
a = 1.25;
b = 0.5;

% this part is actually hard coded into simulink, just as reference here
for fb=1:N
    Wp = [Wn(1)+(fb-1)*step_size, Wn(2)]/(Fs/2);
    [num,den] = cheby1(6,1,Wp,'bandpass');
    num1 = num;
    den1 = den;
end

%% Generating Reference Signals
% References signal with sine-cosine square waveforms
sc1=refsig_square(sti_f(1),Fs,window_time*Fs,H);
sc2=refsig_square(sti_f(2),Fs,window_time*Fs,H);
sc3=refsig_square(sti_f(3),Fs,window_time*Fs,H);
sc4=refsig_square(sti_f(4),Fs,window_time*Fs,H);

%% For offline experiments

% number of trials for all stimulus frequencies
Stimulus_N = 15;

% initial_rest time in seconds
Stimulus_start_Time = 10;

% the sequence to be used 
Stimulus_frequencies_raw = [];
Stimulus_frequencies_easy = [];
for i = 1:Stimulus_start_Time
    Stimulus_frequencies_raw = [Stimulus_frequencies_raw,0];
end

Stimulus_frequencies_timeseries = zeros(1,Stimulus_start_Time*Fs+window*n_sti*Stimulus_N);
Stimulus_space = 1:n_sti;

%step 1 conversion to array of stimulus frequency indexes
for i=1:Stimulus_N
    x=Stimulus_space(randperm(length(Stimulus_space)));
    for j=1:n_sti
        Stimulus_frequencies_easy = [Stimulus_frequencies_easy,x(j)];
        for k = 1:window_time
            Stimulus_frequencies_raw = [Stimulus_frequencies_raw,x(j)];
        end
        Stimulus_frequencies_raw = [Stimulus_frequencies_raw,0];
    end
end
for i = 1:Stimulus_start_Time
    Stimulus_frequencies_raw = [Stimulus_frequencies_raw,0];
end
Stimulus_frequencies_raw = timeseries(Stimulus_frequencies_raw);
