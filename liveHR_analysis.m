function [new_peak_loc]=liveHR_analysis(ecg_avail,fs,liveBeatStart);
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
    

    % Initialize
    new_peak_loc = NaN;
    % Adaptive threshold for peak detection
    threshold = 0.6; %determined by looking at the signal characteristics
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

    % Extract new data for analysis
    new_data = ecg_avail(new_data_start:end,:);
    average_ecg = mean(new_data, 2);

    [b, a] = butter(3, Wn);% 3rd order Butterworth bandpass filter normalized to Nyquist frequency (applicable to multiple sampling rates
    filtered_ecg = filter(b, a, average_ecg);


    % Find peaks
    if length(filtered_ecg) < 3
        peaks = [0];
        locs = [0];
    else
        [peaks, locs] = findpeaks(filtered_ecg, 'MinPeakDistance', 1 / 1000);
    end

    if ~isempty(peaks)
    % Find peaks above the threshold
    peaks_above_threshold = find(peaks > threshold);
        if ~isempty(peaks_above_threshold)
            % Adjust and store locations
            stored_locations = locs(peaks_above_threshold);
            stored_locations = stored_locations + new_data_start - 1;
        else
            stored_locations = [];
        end
        % Loop through stored locations
        for i = 1:length(stored_locations)
            % Check if the current location has been processed
            if last_processed_sample == stored_locations(i)
                continue;
            else
                % Set new_peak_loc and break the loop
                new_peak_loc = stored_locations(i);
                break;
            end
        end
    end
end