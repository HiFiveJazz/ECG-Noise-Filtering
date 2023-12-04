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
    [b, a] = butter(3, Wn);
    

    % Initialize
    new_peak_loc = NaN;
   
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
    average_ecg = mean(new_data,2);
    %filteredECG = filter(b, a, ecg);
    filteredbnECG = filter(b, a, average_ecg);
    window_size = 5; % Adjust this value based on your requirements

    % Apply the moving average filter
    smoothed_signal = movmean(filteredbnECG, window_size);
    filteredECG = (smoothed_signal - min(smoothed_signal)) / (max(smoothed_signal) - min(smoothed_signal));

    mad_threshold = 1.4 * mad(filteredECG, 1);

    % outliers based on MAD threshold
    outliers = abs(filteredECG - median(filteredECG)) > 3 * mad_threshold;
    heightECG = filteredECG;
    heightECG(outliers) = NaN;

    % find the max value considering outliers
    max_ecg = max(heightECG, [], 'omitnan');
    height = max_ecg;

    %distance = round(fs / (2));

    % Find peaks
    if length(filteredECG) < 3
        peaks = [0];
        locs = [0];
    else
        [peaks, locs] = findpeaks(filteredECG, 'MinPeakHeight', height);
    end

    if ~isempty(peaks)
    % Find peaks above the threshold
    peaks_above_threshold = find(peaks >= height);
        if ~isempty(peaks_above_threshold)
            % Adjust and store locations
            stored_locations = locs(peaks_above_threshold);
            stored_locations = stored_locations + new_data_start - 1;
        else
            stored_locations = [];
        end
        % Loop through stored locations
        index_last_processed = find(stored_locations == last_processed_sample);
        if ~isempty(index_last_processed)
            for i = index_last_processed:length(stored_locations)
                if last_processed_sample == stored_locations(i)
                    continue;
                else
            % Set new_peak_loc and break the loop
                    new_peak_loc = stored_locations(i);
                    break;
                end
            end
        else
            new_peak_loc = stored_locations(1);
        end

    end