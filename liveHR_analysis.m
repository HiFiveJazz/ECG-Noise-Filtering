function [new_peak_loc]=liveHR_analysis(ecg_avail,fs,liveBeatStart)

    % Initialize
    new_peak_loc = NaN;
    fs = double(fs);
    lower_cutoff = 0.5/(fs / 2);   
    upper_cutoff = 50/(fs / 2);  
    if upper_cutoff >=1
        upper_cutoff = 0.99999999;
    end
    if lower_cutoff >= 1
        lower_cutoff = 0.001;
    end
    last_processed_sample = 0;
    if ~isempty(liveBeatStart)
        last_processed_sample = liveBeatStart(end);
    end
    if last_processed_sample > 101
        new_data_start = last_processed_sample - 100;
    else
        new_data_start = last_processed_sample + 1;
    end 
    if new_data_start > length(ecg_avail)
        return;
    end

    Wn = [lower_cutoff, upper_cutoff]; 


    [b, a] = butter(3, Wn);% 3rd order Butterworth bandpass filter normalized to Nyquist frequency (applicable to multiple sampling rates
    average_ecg = mean(ecg_avail,2);
    %filteredECG = filter(b, a, ecg);
    filteredbnECG = filter(b, a, average_ecg);
    diffECG = diff(filteredbnECG);

    % Squaring
    squaredECG = diffECG .^ 2;
    window_size = round(0.150 * fs);

    % Apply the moving average filter
    smoothed_signal = movmean(squaredECG, window_size);
    filteredECG = (smoothed_signal - min(smoothed_signal)) / (max(smoothed_signal) - min(smoothed_signal));
    plot(filteredECG)


    N = length(filteredECG);
    X = fft(filteredbnECG); 
    freq_resolution = fs / N;
    mag_spectrum = abs(X);
    normalized_spectrum = mag_spectrum /max(mag_spectrum);
    [p,l]=findpeaks(normalized_spectrum);
    [~, max_index] = max(p);
    dominantFrequencyIndex = l(max_index);
    cardiacFreq = (dominantFrequencyIndex - 1) * (fs / length(filteredbnECG));
    distance = round(fs / (2*cardiacFreq));
    if distance < fs/2
        distance = fs/2;
    end
    
    mad_threshold = 1.4 * mad(filteredECG, 1);

    % outliers based on the MAD threshold
    outliers = abs(filteredECG - median(filteredECG)) > 3 * mad_threshold;
    heightECG = filteredECG;
    heightECG(outliers) = NaN;

    % find the max value considering outliers
    max_ecg = max(heightECG, [], 'omitnan');
    height = max_ecg;
   % plot(filteredECG);

    [QRS_Peaks, peakindices] = findpeaks(filteredECG, 'MinPeakHeight', height, 'MinPeakDistance', distance);

    if ~isempty(peakindices)
        stored_locations = peakindices + new_data_start - 1;
        %new_peak_loc = stored_locations;
    else
        stored_locations = [];
    end
    new_locations = setdiff(stored_locations, liveBeatStart);
    new_peak_loc = new_locations;
   

    % Loop through stored locations
    % for i = 1:length(stored_locations)
    % % Check if the current location has been processed or exists in liveBeatStart
    %     if last_processed_sample == stored_locations(i) || any(stored_locations(i) == liveBeatStart)
    %         continue;
    %     elseif any(abs(stored_locations(i) - liveBeatStart) <= 200)
    %     % Check if the stored location is within 50 time points of any last processed samples
    %         continue;
    %     else
    %     % Set new_peak_loc, update last_processed_sample, and break the loop
    %         new_peak_loc = stored_locations(i);
    %         break;
    %     end
    % end
end