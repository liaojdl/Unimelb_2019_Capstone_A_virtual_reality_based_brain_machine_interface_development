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