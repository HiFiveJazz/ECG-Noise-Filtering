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
    %filteredECG = filter(b, a, ecg);
    filteredbnECG = filter(b, a, average_ecg);
    window_size = 5; 

    % Apply the moving average filter
    smoothed_signal = movmean(filteredbnECG, window_size);
    filteredECG = (smoothed_signal - min(smoothed_signal)) / (max(smoothed_signal) - min(smoothed_signal));
    %plot(filteredECG)


    N = length(filteredECG);
    
    
    X = fft(filteredbnECG); 
    freq_resolution = fs / N;
    mag_spectrum = abs(X);
    [~, max_index] = max(mag_spectrum);
    cardiacFreq = (max_index - 1) * freq_resolution;
    
    mad_threshold = 1.4 * mad(filteredECG, 1);

    % outliers based on the MAD threshold
    outliers = abs(filteredECG - median(filteredECG)) > 3 * mad_threshold;
    heightECG = filteredECG;
    heightECG(outliers) = NaN;

    % find the max value considering outliers
    max_ecg = max(heightECG, [], 'omitnan');
    height = max_ecg;
    plot(filteredECG)

    distance = round(fs / (2));
    %height = 0.5*max(filteredECG);
    [peaks, locs] = findpeaks(filteredECG, 'MinPeakHeight', height, 'MinPeakDistance', distance);

    intervals = diff(locs) / fs; 
    heartRates = 60 ./ intervals; % Convert to bpm


    avgHR = round(mean(heartRates));  % Compute average heart rate


end

