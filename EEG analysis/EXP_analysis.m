%% This file is used to analyze the experimental data collected
%  on the 3rd of July for the VR-EEG test, The test subject is
%  Harfi D. Santoso. 
% Experiment type : SSVEP evoking experiments with flashing in-game lights
% Analysis: using Canonical Correlation Analysis and simple threshold 
% detecting method for comparison
%%
%% clean start
clc;
clear all;
close all;

%% FBCCA Processing Variables
Fs = 512;                       % sampling frequency
window_ssvep = 5;               % length of ssvep window in seconds
window_motor = 8;               % length of motor window in seconds
rest_time = 1;

%% exp variables
% number of trials for all stimulus frequencies for SSVEP
Stimulus_N = 75;
% number of trials for all prompts for motor imagery/action
Motor_N = 45;
% initial_rest time in seconds
Stimulus_start_Time = 10.2;

%% loading raw data
%% SSVEP data for each of the three stimulus types
sequence_flick = load('D:\VR_BCI\2019_10_03_3experiments\harfi\flick\prompts1_start_10window_5rest_1.mat');
sequence_flick = sequence_flick.Stimulus_frequencies_easy;
eeg_flick = load('D:\VR_BCI\2019_10_03_3experiments\harfi\flick\EEG_raw.mat');
eeg_flick = eeg_flick.eeg;
eeg_flick_t = eeg_flick(1,:);

sequence_shrink = load('D:\VR_BCI\2019_10_03_3experiments\harfi\shrink\prompts1_start_10window_5rest_1.mat');
sequence_shrink = sequence_shrink.Stimulus_frequencies_easy;
eeg_shrink = load('D:\VR_BCI\2019_10_03_3experiments\harfi\shrink\EEG_raw.mat');
eeg_shrink = eeg_shrink.eeg;
eeg_shrink_t = eeg_shrink(1,:);

sequence_chess = load('D:\VR_BCI\2019_10_03_3experiments\harfi\chess\prompts1_start_10window_5rest_1.mat');
sequence_chess = sequence_chess.Stimulus_frequencies_easy;
eeg_chess = load('D:\VR_BCI\2019_10_03_3experiments\harfi\chess\EEG_raw.mat');
eeg_chess = eeg_chess.eeg;
eeg_chess_t = eeg_chess(1,:);

%% Motor action and motor imagery ERS/ERD raw data
sequence_action = load('D:\VR_BCI\2019_10_03_3experiments\harfi\mi\mi_start_10window_8rest_1.mat');
sequence_action = sequence_action.Stimulus_frequencies_easy;
eeg_action = load('D:\VR_BCI\2019_10_03_3experiments\harfi\mi\EEG_raw.mat');
eeg_action = eeg_action.eeg;
eeg_action_t = eeg_action(1,:);

%% start_time and finish time arrays for each trial,
t1m_start = Stimulus_start_Time;
t1m_flick = timing_array(Stimulus_start_Time,Stimulus_N*2,window_ssvep,rest_time);
t1m_shrink = timing_array(Stimulus_start_Time,Stimulus_N*2,window_ssvep,rest_time);
t1m_chess = timing_array(Stimulus_start_Time,Stimulus_N*2,window_ssvep,rest_time);
t1m_action = timing_array(Stimulus_start_Time,Motor_N*2,window_motor,rest_time);
t1m_mi = timing_array(Stimulus_start_Time,Motor_N*2,window_motor,rest_time);
[t1mi_flick, t1mt_flick] = time2index(t1m_flick,eeg_flick_t,Fs);
[t1mi_shrink, t1mt_shrink] = time2index(t1m_shrink,eeg_shrink_t,Fs);
[t1mi_chess, t1mt_chess] = time2index(t1m_chess,eeg_chess_t,Fs);
[t1mi_action, t1mt_action] = time2index(t1m_action,eeg_action_t,Fs);

%% now need to separate the raw data for each individual channel
eegdata_flick = raw_cut2seg(eeg_flick,t1mi_flick);
eegdata_shrink = raw_cut2seg(eeg_shrink,t1mi_shrink);
eegdata_chess = raw_cut2seg(eeg_chess,t1mi_chess);
eegdata_action = raw_cut2seg(eeg_action,t1mi_action);

%% filters necessary for the task 
% 5th order butterworth band pass filter, needed just for CCA 8-20Hz
% analysis
order_n = 5;
band_wn = [12,30]/(Fs/2);
[num1,denom1] = butter(order_n,band_wn,'bandpass');
% % 2nd order IIR notch filter for 50HZ power
% W_0 = 50/(Fs/2);  
% Q_factor = 35;
% B_W = W_0/Q_factor;
% [num2,denom2] = iirnotch(W_0,B_W);
% creates the cascate filter
filter_12_30_bandpass = dfilt.df2t(num1,denom1);
% filter2 = dfilt.df2t(num2,denom2);
% filter3 = dfilt.cascade(filter1,filter2);
% fvtool(filter3);

%% applying filters to the 9 occipital lobe electrodes
% for i = 1:length(eegdata_flick)
%     eegdata_flick(i).pz = filter(filter_8_20_bandpass,eegdata_flick(i).pz);
%     eegdata_flick(i).po3 = filter(filter_8_20_bandpass,eegdata_flick(i).po3);
%     eegdata_flick(i).po4 = filter(filter_8_20_bandpass,eegdata_flick(i).po4);
%     eegdata_flick(i).po7 = filter(filter_8_20_bandpass,eegdata_flick(i).po7);
%     eegdata_flick(i).po8 = filter(filter_8_20_bandpass,eegdata_flick(i).po8);
%     eegdata_flick(i).poz = filter(filter_8_20_bandpass,eegdata_flick(i).poz);
%     eegdata_flick(i).o1 = filter(filter_8_20_bandpass,eegdata_flick(i).o1);
%     eegdata_flick(i).o2 = filter(filter_8_20_bandpass,eegdata_flick(i).o2);
%     eegdata_flick(i).oz = filter(filter_8_20_bandpass,eegdata_flick(i).oz);
%     eegdata_flick(i).c3 = filter(filter_8_20_bandpass,eegdata_flick(i).c3);
%     eegdata_flick(i).c1 = filter(filter_8_20_bandpass,eegdata_flick(i).c1);
%     eegdata_flick(i).cz = filter(filter_8_20_bandpass,eegdata_flick(i).cz);
%     eegdata_flick(i).c2 = filter(filter_8_20_bandpass,eegdata_flick(i).c2);
%     eegdata_flick(i).c4 = filter(filter_8_20_bandpass,eegdata_flick(i).c4);
%     eegdata_flick(i).fcz = filter(filter_8_20_bandpass,eegdata_flick(i).fcz);
% end

% for i = 1:length(eegdata_flick)
%     eegdata_flick(i).pz = filter(filter_8_20_bandpass,eegdata_flick(i).pz);
%     eegdata_flick(i).po3 = filter(filter_8_20_bandpass,eegdata_flick(i).po3);
%     eegdata_flick(i).po4 = filter(filter_8_20_bandpass,eegdata_flick(i).po4);
%     eegdata_flick(i).po7 = filter(filter_8_20_bandpass,eegdata_flick(i).po7);
%     eegdata_flick(i).po8 = filter(filter_8_20_bandpass,eegdata_flick(i).po8);
%     eegdata_flick(i).poz = filter(filter_8_20_bandpass,eegdata_flick(i).poz);
%     eegdata_flick(i).o1 = filter(filter_8_20_bandpass,eegdata_flick(i).o1);
%     eegdata_flick(i).o2 = filter(filter_8_20_bandpass,eegdata_flick(i).o2);
%     eegdata_flick(i).oz = filter(filter_8_20_bandpass,eegdata_flick(i).oz);
% end
% 
% for i = 1:length(eegdata_shrink)
%     eegdata_shrink(i).pz = filter(filter_8_20_bandpass,eegdata_shrink(i).pz);
%     eegdata_shrink(i).po3 = filter(filter_8_20_bandpass,eegdata_shrink(i).po3);
%     eegdata_shrink(i).po4 = filter(filter_8_20_bandpass,eegdata_shrink(i).po4);
%     eegdata_shrink(i).po7 = filter(filter_8_20_bandpass,eegdata_shrink(i).po7);
%     eegdata_shrink(i).po8 = filter(filter_8_20_bandpass,eegdata_shrink(i).po8);
%     eegdata_shrink(i).poz = filter(filter_8_20_bandpass,eegdata_shrink(i).poz);
%     eegdata_shrink(i).o1 = filter(filter_8_20_bandpass,eegdata_shrink(i).o1);
%     eegdata_shrink(i).o2 = filter(filter_8_20_bandpass,eegdata_shrink(i).o2);
%     eegdata_shrink(i).oz = filter(filter_8_20_bandpass,eegdata_shrink(i).oz);
% end

%% SSVEP analysis
%% References signal with sine-cosine waveforms
sti_f = [9,11,13,15];
n_stif = length(sti_f);
H = 2;
sc1=refsig_sin(sti_f(1),Fs,window_ssvep*Fs+1,H);
sc2=refsig_sin(sti_f(2),Fs,window_ssvep*Fs+1,H);
sc3=refsig_sin(sti_f(3),Fs,window_ssvep*Fs+1,H);
sc4=refsig_sin(sti_f(4),Fs,window_ssvep*Fs+1,H);


%% CCA
CCA_threshold = 0.2;  %threshold for null target

cca_results_flick = zeros(Stimulus_N,n_stif);
cca_results_shrink = zeros(Stimulus_N,n_stif);
cca_results_chess = zeros(Stimulus_N,n_stif);
ssvepdata_flick = zeros(Stimulus_N,9,window_ssvep*Fs+1);
ssvepdata_shrink = zeros(Stimulus_N,9,window_ssvep*Fs+1);
ssvepdata_chess = zeros(Stimulus_N,9,window_ssvep*Fs+1);
cca_windows = 0.5:0.5:5;
cca_rate_flick = zeros(1,length(cca_windows));
cca_rate_shrink = zeros(1,length(cca_windows));
cca_rate_chess = zeros(1,length(cca_windows));

for i = 1:Stimulus_N
    ssvepdata_flick(i,1,:) = eegdata_flick(i).pz;
    ssvepdata_flick(i,2,:) = eegdata_flick(i).po7;
    ssvepdata_flick(i,3,:) = eegdata_flick(i).po3;
    ssvepdata_flick(i,4,:) = eegdata_flick(i).poz;
    ssvepdata_flick(i,5,:) = eegdata_flick(i).po4;
    ssvepdata_flick(i,6,:) = eegdata_flick(i).po8;
    ssvepdata_flick(i,7,:) = eegdata_flick(i).o1;
    ssvepdata_flick(i,8,:) = eegdata_flick(i).oz;
    ssvepdata_flick(i,9,:) = eegdata_flick(i).o2;
end

for i = 1:Stimulus_N
    ssvepdata_shrink(i,1,:) = eegdata_shrink(i).pz;
    ssvepdata_shrink(i,2,:) = eegdata_shrink(i).po7;
    ssvepdata_shrink(i,3,:) = eegdata_shrink(i).po3;
    ssvepdata_shrink(i,4,:) = eegdata_shrink(i).poz;
    ssvepdata_shrink(i,5,:) = eegdata_shrink(i).po4;
    ssvepdata_shrink(i,6,:) = eegdata_shrink(i).po8;
    ssvepdata_shrink(i,7,:) = eegdata_shrink(i).o1;
    ssvepdata_shrink(i,8,:) = eegdata_shrink(i).oz;
    ssvepdata_shrink(i,9,:) = eegdata_shrink(i).o2;
end

for i = 1:Stimulus_N
    ssvepdata_chess(i,1,:) = eegdata_chess(i).pz;
    ssvepdata_chess(i,2,:) = eegdata_chess(i).po7;
    ssvepdata_chess(i,3,:) = eegdata_chess(i).po3;
    ssvepdata_chess(i,4,:) = eegdata_chess(i).poz;
    ssvepdata_chess(i,5,:) = eegdata_chess(i).po4;
    ssvepdata_chess(i,6,:) = eegdata_chess(i).po8;
    ssvepdata_chess(i,7,:) = eegdata_chess(i).o1;
    ssvepdata_chess(i,8,:) = eegdata_chess(i).oz;
    ssvepdata_chess(i,9,:) = eegdata_chess(i).o2;
end

for i = 1:Stimulus_N
    ssvep_9row = zeros (9,window_ssvep*Fs+1);
    for j = 1:9
        for k = 1: (window_ssvep*Fs+1)
            ssvep_9row(j,k) = ssvepdata_flick(i,j,k);
        end
    end
    for m = 1:length(cca_windows)
        [~,~,r1] = cca(ssvep_9row(:,1:cca_windows(m)*Fs),sc1(:,1:cca_windows(m)*Fs));
        [~,~,r2] = cca(ssvep_9row(:,1:cca_windows(m)*Fs),sc2(:,1:cca_windows(m)*Fs));
        [~,~,r3] = cca(ssvep_9row(:,1:cca_windows(m)*Fs),sc3(:,1:cca_windows(m)*Fs));
        [~,~,r4] = cca(ssvep_9row(:,1:cca_windows(m)*Fs),sc4(:,1:cca_windows(m)*Fs));
        weights = [r1,r2,r3,r4];
        [v,idx] = max([max(r1),max(r2),max(r3),max(r4)]);
        cca_results_flick(i,:) = [max(r1),max(r2),max(r3),max(r4)];
%         if v < CCA_threshold
%             idx = 0;
%         end
        if idx == sequence_flick(i)
            cca_rate_flick(m) = cca_rate_flick(m)+1;
        end
    end
end
cca_rate_flick = cca_rate_flick/(Stimulus_N-15)



for i = 1:Stimulus_N
    ssvep_9row = zeros (9,window_ssvep*Fs+1);
    for j = 1:9
        for k = 1: (window_ssvep*Fs+1)
            ssvep_9row(j,k) = ssvepdata_shrink(i,j,k);
        end
    end
    for m = 1:length(cca_windows)
        [~,~,r1] = cca(ssvep_9row(:,1:cca_windows(m)*Fs),sc1(:,1:cca_windows(m)*Fs));
        [~,~,r2] = cca(ssvep_9row(:,1:cca_windows(m)*Fs),sc2(:,1:cca_windows(m)*Fs));
        [~,~,r3] = cca(ssvep_9row(:,1:cca_windows(m)*Fs),sc3(:,1:cca_windows(m)*Fs));
        [~,~,r4] = cca(ssvep_9row(:,1:cca_windows(m)*Fs),sc4(:,1:cca_windows(m)*Fs));
        weights = [r1,r2,r3,r4];
        [v,idx] = max([max(r1),max(r2),max(r3),max(r4)]);
        cca_results_shrink(i,:) = [max(r1),max(r2),max(r3),max(r4)];
%         if v < CCA_threshold
%             idx = 0;
%         end
        if idx == sequence_shrink(i)
            cca_rate_shrink(m) = cca_rate_shrink(m)+1;
        end
    end
end
cca_rate_shirnk = cca_rate_shrink/(Stimulus_N-15)

for i = 1:Stimulus_N
    ssvep_9row = zeros (9,window_ssvep*Fs+1);
    for j = 1:9
        for k = 1: (window_ssvep*Fs+1)
            ssvep_9row(j,k) = ssvepdata_chess(i,j,k);
        end
    end
    for m = 1:length(cca_windows)
        [~,~,r1] = cca(ssvep_9row(:,1:cca_windows(m)*Fs),sc1(:,1:cca_windows(m)*Fs));
        [~,~,r2] = cca(ssvep_9row(:,1:cca_windows(m)*Fs),sc2(:,1:cca_windows(m)*Fs));
        [~,~,r3] = cca(ssvep_9row(:,1:cca_windows(m)*Fs),sc3(:,1:cca_windows(m)*Fs));
        [~,~,r4] = cca(ssvep_9row(:,1:cca_windows(m)*Fs),sc4(:,1:cca_windows(m)*Fs));
        weights = [r1,r2,r3,r4];
        [v,idx] = max([max(r1),max(r2),max(r3),max(r4)]);
        cca_results_chess(i,:) = [max(r1),max(r2),max(r3),max(r4)];
%         if v < CCA_threshold
%             idx = 0;
%         end
        if idx == sequence_chess(i)
            cca_rate_chess(m) = cca_rate_chess(m)+1;
        end
    end
end
cca_rate_shirnk = cca_rate_chess/(Stimulus_N-15)


%% FBCCA trials
FBCCA_threshold = 0.2;
fbcca_windows = 0.5:0.5:5;

fbcca_results_flick = zeros(Stimulus_N,4);
fbcca_results_shrink = zeros(Stimulus_N,4);
fbcca_results_chess = zeros(Stimulus_N,4);
fbcca_rate_flick = zeros(1,length(fbcca_windows));
fbcca_rate_shrink = zeros(1,length(fbcca_windows));
fbcca_rate_chess = zeros(1,length(fbcca_windows));

for i = 1:Stimulus_N
    ssvep_9row = zeros (9,window_ssvep*Fs+1);
    for j = 1:9
        for k = 1: (window_ssvep*Fs+1)
            ssvep_9row(j,k) = ssvepdata_flick(i,j,k);
        end
    end
    for m = 1:length(fbcca_windows)
        rb1 = fbcca(ssvep_9row(:,1:fbcca_windows(m)*Fs),sc1(:,1:fbcca_windows(m)*Fs),[6,60],9,512);
        rb2 = fbcca(ssvep_9row(:,1:fbcca_windows(m)*Fs),sc2(:,1:fbcca_windows(m)*Fs),[6,60],9,512);
        rb3 = fbcca(ssvep_9row(:,1:fbcca_windows(m)*Fs),sc3(:,1:fbcca_windows(m)*Fs),[6,60],9,512);
        rb4 = fbcca(ssvep_9row(:,1:fbcca_windows(m)*Fs),sc4(:,1:fbcca_windows(m)*Fs),[6,60],9,512);
        weights = [rb1,rb2,rb3,rb4];
        [v,idx] = max([max(rb1),max(rb2),max(rb3),max(rb4)]);
        fbcca_results_flick(i,:) = [max(rb1),max(rb2),max(rb3),max(rb4)];
%         if v < FBCCA_threshold
%             idx = 0;
%         end
        if idx == sequence_flick(i)
            fbcca_rate_flick(m) = fbcca_rate_flick(m)+1;
        end
    end
end
fbcca_rate_flick = fbcca_rate_flick/(Stimulus_N-15)

for i = 1:Stimulus_N
    ssvep_9row = zeros (9,window_ssvep*Fs+1);
    for j = 1:9
        for k = 1: (window_ssvep*Fs+1)
            ssvep_9row(j,k) = ssvepdata_shrink(i,j,k);
        end
    end
    for m = 1:length(fbcca_windows)
        rb1 = fbcca(ssvep_9row(:,1:fbcca_windows(m)*Fs),sc1(:,1:fbcca_windows(m)*Fs),[6,60],9,512);
        rb2 = fbcca(ssvep_9row(:,1:fbcca_windows(m)*Fs),sc2(:,1:fbcca_windows(m)*Fs),[6,60],9,512);
        rb3 = fbcca(ssvep_9row(:,1:fbcca_windows(m)*Fs),sc3(:,1:fbcca_windows(m)*Fs),[6,60],9,512);
        rb4 = fbcca(ssvep_9row(:,1:fbcca_windows(m)*Fs),sc4(:,1:fbcca_windows(m)*Fs),[6,60],9,512);
        weights = [rb1,rb2,rb3,rb4];
        [v,idx] = max([max(rb1),max(rb2),max(rb3),max(rb4)]);
        fbcca_results_shrink(i,:) = [max(rb1),max(rb2),max(rb3),max(rb4)];
%         if v < FBCCA_threshold
%             idx = 0;
%         end
        if idx == sequence_shrink(i)
            fbcca_rate_shrink(m) = fbcca_rate_shrink(m)+1;
        end
    end
end
fbcca_rate_shrink = fbcca_rate_shrink/(Stimulus_N-15)

for i = 1:Stimulus_N
    ssvep_9row = zeros (9,window_ssvep*Fs+1);
    for j = 1:9
        for k = 1: (window_ssvep*Fs+1)
            ssvep_9row(j,k) = ssvepdata_chess(i,j,k);
        end
    end
    for m = 1:length(fbcca_windows)
        rb1 = fbcca(ssvep_9row(:,1:fbcca_windows(m)*Fs),sc1(:,1:fbcca_windows(m)*Fs),[6,60],9,512);
        rb2 = fbcca(ssvep_9row(:,1:fbcca_windows(m)*Fs),sc2(:,1:fbcca_windows(m)*Fs),[6,60],9,512);
        rb3 = fbcca(ssvep_9row(:,1:fbcca_windows(m)*Fs),sc3(:,1:fbcca_windows(m)*Fs),[6,60],9,512);
        rb4 = fbcca(ssvep_9row(:,1:fbcca_windows(m)*Fs),sc4(:,1:fbcca_windows(m)*Fs),[6,60],9,512);
        weights = [rb1,rb2,rb3,rb4];
        [v,idx] = max([max(rb1),max(rb2),max(rb3),max(rb4)]);
        fbcca_results_shrink(i,:) = [max(rb1),max(rb2),max(rb3),max(rb4)];
%         if v < FBCCA_threshold
%             idx = 0;
%         end
        if idx == sequence_chess(i)
            fbcca_rate_chess(m) = fbcca_rate_chess(m)+1;
        end
    end
end
fbcca_rate_chess = fbcca_rate_chess/(Stimulus_N-15)


%% Motor processing
%% action
for i = 1:length(eegdata_action)
    eegdata_action(i).c3 = filter(filter_12_30_bandpass,eegdata_action(i).c3);
    eegdata_action(i).c1 = filter(filter_12_30_bandpass,eegdata_action(i).c1);
    eegdata_action(i).cz = filter(filter_12_30_bandpass,eegdata_action(i).cz);
    eegdata_action(i).c2 = filter(filter_12_30_bandpass,eegdata_action(i).c2);
    eegdata_action(i).c4 = filter(filter_12_30_bandpass,eegdata_action(i).c4);
    eegdata_action(i).fcz = filter(filter_12_30_bandpass,eegdata_action(i).fcz);
end

pl_action = [];
pr_action = [];
po_action = [];
for i = 1:Motor_N
    c3_norm = eegdata_action(i).c3-eegdata_action(i).fcz;
    [pc3,wc3] = pwelch(c3_norm,256,[],[],Fs);
    c4_norm = eegdata_action(i).c4-eegdata_action(i).fcz;
    [pc4,wc4] = pwelch(c4_norm,256,[],[],Fs);
    pc4 = trapz(wc4(12:30),pc4(12:30));
    pc3 = trapz(wc3(12:30),pc3(12:30));
    if sequence_action(i)==0
        po_action = [po_action;pc3-pc4];
    elseif sequence_action(i)==1
        pl_action = [pl_action;pc3-pc4];
    elseif sequence_action(i)==2
        pr_action = [pr_action;pc3-pc4];
    end
end
figure
plot(pl_action(:,1));
hold on
% plot(pl_action(:,2));
% legend("pc3","pc4")
% figure
% plot(po_action(:,1));
% hold on
% plot(po_action(:,2));
% legend("pc3","pc4")
% figure
% plot(pr_action(:,1));
% hold on
% plot(pr_action(:,2));
% legend("pc3","pc4")
plot(po_action);
plot(pr_action);
legend("left","nope","right")




%% functions 

%% this function takes tm the array of start-finish time arrays
%  and the the actual data recording timing array
%  convert the timings in tm to the closest index in the ta
%  make sure that time in tm is definitely within the recording time frame
% param @Mt2i is the index of the timing in ta
% param @Mt2tc is the closest timing in ta
function [Mt2i, Mt2tc]= time2index(tm,ta,frequency)
    %initiasing output
    Mt2i = tm;
    Mt2tc = tm;
    %initial seeking index at the start
    start_index = 1;
    %closest seeking distance
    eps = (1/frequency)/2;
    
    for i=1:length(tm)
        time_value = tm(i);
        for j = start_index : length(ta)
            %found the clostest one from LHS
            if ((ta(j)) >= (time_value-eps))
                Mt2i(i) = j;
                Mt2tc(i) = ta(j);
                %update start index
                start_index = j;
                break;
            end
        end
    end
end


%% this function takes in start time, number of actions,interval of action,
% interval between each action, and gives back the timing

function tm = timing_array(t_start,N,t_A,t_rest)
    tm = [1:N];
    for i=1:N
        if (i==1)
            tm(i)=t_start;
        else
            if (rem(i,2)==1)
                tm(i) = tm(i-1)+t_rest;
            else 
                tm(i) = tm(i-1)+t_A;
            end
        end
    end
end

%% cut raw exp data into segments based on tmi, 
%  and stores in a struct of segments of struct of raw_data, hopefully
%  this is not confusing lol
% format is start,stop pairs in tmi

% &param raw_segment is an array of structs of segments of the rawdata
% for instance, for 30 actions, it has 30 structs
% each struct has the following fields defined
%  t : time
%  FC4,FC2... etc etc,
% you can access by using raw_seg = raw_cut2seg(rawdata,tmi).
% then use raw_seg[1].FC4. this is the first yes action

function raw_segment = raw_cut2seg(rawdata,tmi)

    %% struct to store raw data of a segment of raw_data
    %  separated into channels first
    field0 = 'tx';
    field1 = 'pz';
    field2 = 'po7';
    field3 = 'po3';
    field4 = 'poz';
    field5 = 'po4';
    field6 = 'po8';
    field7 = 'o1';
    field8 = 'oz';
    field9 = 'o2';
    field10 = 'c3';
    field11 = 'c1';
    field12 = 'cz';
    field13 = 'c2';
    field14 = 'c4';
    field15 = 'fcz';
    
    for i = 1:(length(tmi)/2)        
        start_index = tmi(i*2-1);
        finish_index = tmi(i*2);
        
        
        %temporary storage
        tx = rawdata(1,start_index:finish_index);
        PZx = rawdata(2,start_index:finish_index);
        PO7x = rawdata(3,start_index:finish_index);
        PO3x = rawdata(4,start_index:finish_index);
        POZx = rawdata(5,start_index:finish_index);
        PO4x = rawdata(6,start_index:finish_index);
        PO8x = rawdata(7,start_index:finish_index);
        O1x = rawdata(8,start_index:finish_index);
        OZx = rawdata(9,start_index:finish_index);
        O2x = rawdata(10,start_index:finish_index);
        C3x = rawdata(11,start_index:finish_index);
        C1x = rawdata(12,start_index:finish_index);
        CZx = rawdata(13,start_index:finish_index);
        C2x = rawdata(14,start_index:finish_index);
        C4x = rawdata(15,start_index:finish_index);
        FCZx = rawdata(16,start_index:finish_index);
        
        BCIx = struct(field0,tx,field1,PZx,field2,PO7x,field3,PO3x,field4,POZx, ...
        field5,PO4x,field6,PO8x,field7,O1x,field8,OZx,field9,O2x,field10,C3x,...
        field11,C1x,field12,CZx,field13,C2x,field14,C4x,field15,FCZx);
    
        raw_segment(i) = BCIx;
    end
end




%% create fft data for each channel data segment
function [p_fft,f_fft] = channel_fft(channel,frequency)
    Y = fft(channel);
    L = length(channel);
    P2 = abs(Y/L);
    P1 = P2(1:floor(L/2)+1);
    P1(2:end-1) = 2*P1(2:end-1);
    p_fft = P1;
    f_fft = frequency*(0:(L/2))/L;
end
