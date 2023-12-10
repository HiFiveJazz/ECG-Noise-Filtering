liveBeatStart = [];
liveBeatLag = [];

for j=5*fs:100:size(ecg,1)
        
        ecg_avail=ecg(1:j,:);  
    
        % analyze the ecg and add to the log if you find a peak
    
        [new_peak_loc]=liveHR_analysis(ecg_avail,fs,liveBeatStart);
    
        if ~isnan(new_peak_loc)
            liveBeatStart=[liveBeatStart; new_peak_loc];
            liveBeatLag = [liveBeatLag; (j-(new_peak_loc))*1000/(fs)]; %in milliseconds
        end
end
t = (0:length(filteredECG)-1) / fs;
    figure;
    plot(t, filteredECG);
    hold on;

    % Mark the QRS peaks on top of the filtered ECG signal
    plot(liveBeatStart/fs, filteredECG(liveBeatStart), 'go', 'MarkerSize', 10);