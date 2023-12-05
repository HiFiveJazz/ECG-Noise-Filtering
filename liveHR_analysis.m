function [new_peak_loc]=liveHR_analysis(ecg_avail,fs,liveBeatStart)

    % Initialize
    new_peak_loc = NaN;
   
    last_processed_sample = 0;
    if ~isempty(liveBeatStart)
        last_processed_sample = liveBeatStart(end);
    end
    if last_processed_sample > (fs/2)+1
        new_data_start = last_processed_sample - (fs/2);
    else
        new_data_start = last_processed_sample + 1;
    end 
    if new_data_start > length(ecg_avail)
        return;
    end

    % Extract new data for analysis
    new_data = ecg_avail(new_data_start:end,:);
    average_ecg = mean(new_data,2);

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
    
    mad_threshold = 1.4 * mad(filteredECG, 1);

    % outliers based on the MAD threshold
    outliers = abs(filteredECG - median(filteredECG)) > 3 * mad_threshold;
    heightECG = filteredECG;
    heightECG(outliers) = NaN;

    % find the max value considering outliers
    max_ecg = max(heightECG, [], 'omitnan');
    height = max_ecg;
    plot(filteredECG)
    distance = round(fs/2);

    [QRS_Peaks, peakindices] = findpeaks(filteredECG, 'MinPeakHeight', height, 'MinPeakDistance', distance);

    if ~isempty(peakindices)
        stored_locations = peakindices + new_data_start - 1;
    else
        stored_locations = [];
    end
    % Loop through stored locations
    for i = 1:length(stored_locations)
    % Check if the current location has been processed or exists in liveBeatStart
        if last_processed_sample == stored_locations(i) || any(stored_locations(i) == liveBeatStart)
            continue;
        elseif any(abs(stored_locations(i) - liveBeatStart) <= 200)
        % Check if the stored location is within 50 time points of any last processed samples
            continue;
        else
        % Set new_peak_loc, update last_processed_sample, and break the loop
            new_peak_loc = stored_locations(i);
            break;
        end
    end
end
