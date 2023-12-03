function [avgHR]=avgHR_analysis(ecg,fs)
% You may delete any of the following lines of code. You may also add as
% many lines of code as you'd like. 

    fs = double(fs);
    lower_cutoff = 0.5/(fs / 2);   
    upper_cutoff = 50/(fs / 2);  
    if upper_cutoff >=1
        upper_cutoff = 0.99;
    end
    if lower_cutoff >= 1
        lower_cutoff = 0.001;
    end

    % Normalize the cutoff frequencies to the Nyquist frequency
    Wn = [lower_cutoff, upper_cutoff]; 


    [b, a] = butter(3, Wn);% 3rd order Butterworth bandpass filter normalized to Nyquist frequency (applicable to multiple sampling rates
    average_ecg = mean(ecg,2);
    filteredECG = filter(b, a, average_ecg);

    [peaks, locs] = findpeaks(filteredECG, 'MinPeakHeight', 0.6, 'MinPeakDistance', 1/100000);

    intervals = diff(locs) / fs; 
    heartRates = 60 ./ intervals; % Convert to bpm


    avgHR = mean(heartRates);  % Compute average heart rate


end

