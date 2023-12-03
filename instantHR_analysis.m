function [instantHR,beatStart]=instantHR_analysis(ecg,fs)
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
    

    %filteredECG = ecg;
    %period = 1/fs;
    %time = [0:1:length(filteredECG)-1].*period;
    N = length(filteredECG);
    
    %Normalizing signal from 0 to 1
    ecgmin = min(filteredECG); ecgmax = max(filteredECG);
    range = ecgmax-ecgmin;
    ecgnorm = (filteredECG-ecgmin)/range;
    
    N_min = min(ecgnorm); N_max = max(ecgnorm); N_range = N_max-N_min;
    
    
    X = fft(ecgnorm);
    freq = (0:N-1)*fs/N;%converting to frequency axis/domain from the time
    %Only looking at frequencies between 0.5-40 Hz
    freq = freq(find(freq>=lower_cutoff*(fs/2)));
    X = X(length(X)-length(freq)+1:end);
    freq = freq(find(freq<=(upper_cutoff*(fs/2))));
    X = X(1:length(freq));
    
    [peak,index] = max(abs(X)); %finding the biggest peak in frequency spectrum and its index
    cardiacFreq = freq(index); %the frequency where the biggest peak is located
    
    height = 0.6*N_range+N_min;
    low = -N_max + N_range*0.4;
    
    distance = 1/(2*cardiacFreq);   %adding 40% error to the frequency, so we set MinPeakDistance 40% less to capture all important peaks
    
    [QRS_Peaks, peakindices] = findpeaks(filteredECG, 'MinPeakHeight', height, 'MinPeakDistance', distance);
    beatStart = peakindices.';
    numPeaks = length(QRS_Peaks);
    for i=1:length(peakindices)-1
        peak_to_peak_distance(i) = peakindices(i+1)-peakindices(i);
    end
    instantHR = floor(60./peak_to_peak_distance);

end