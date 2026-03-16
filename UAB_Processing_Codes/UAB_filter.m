function [filtered_data] = UAB_filter(out,Arts,base,dirname,out_path)
%% Code to filter at 60 Hz
path = string(append(out_path,'/',dirname));
cd(path);
filtered_data = base;
fs = out.fs;
stimlist = fieldnames(out.removedart_trials); 
stim = string(stimlist(1));
resplist = fieldnames(out.epoch_data.(stim)); 
for i = 1:length(stimlist)%each stim
    S = string(stimlist(i));
    for j = 1:length(resplist)%each response
        R = string(resplist(j));
        if Arts.(S)(j) == 1
            F = ceil(fs*1.4);
            filtered_data.filteredmeans.(S).(R) = NaN(1,F);
            continue;
        else
            %concatenate epochs
            A = [];
            for k = 2:(size(out.removedart_trials.(S).(R),1)-1)
                if k == size(out.removedart_trials.(S).(R),1)-1
                    trial = out.removedart_trials.(S).(R)(k,1:end);
                else
                    trial = out.removedart_trials.(S).(R)(k,1:(fs-1));
                end
                    A = [A trial];
            end
            B = A; 
            n = 500;
            Nyq = fs/2;
            lowerbound = (59/Nyq);
            upperbound = (61/Nyq);
            f1 = 0;
            f2 = lowerbound - (0.1*lowerbound);
            f3 = lowerbound;
            f4 = upperbound;
            f5 = upperbound + (0.1*upperbound);
            f6 = 1;
            f = [f1 f2 f3 f4 f5 f6];
            a = [1 1 0 0 1 1];
            a2 = 1;
            b = firls(n,f,a);
            filteredsignal = filtfilt(b,a2,A);
%                 %low pass @150Hz 
%                 f1 = 0;
%                 f2 = 150/Nyq;
%                 f3 = 165/Nyq;
%                 f4 = 1;
%                 f = [f1 f2 f3 f4];
%                 a = [1 1 0 0];
%                 a2 = 1;
%                 b = firls(n,f,a);
%                 filteredsignal = filtfilt(b,a2,filteredsignal);
%                 %high pass @2 Hz
%                 f1 = 0;
%                 f2 = 1/Nyq;
%                 f3 = 2/Nyq;
%                 f4 = 1;
%                 f = [f1 f2 f3 f4];
%                 a = [0 0 1 1];
%                 a2 = 1;
%                 b = firls(n,f,a);
%                 filteredsignal = filtfilt(b,a2,filteredsignal);
%                 filteredsignalB = filteredsignal;
            %re-epoching
            epochs = zeros(k,round(1.4*fs)+1);
            for k = 1:((size(epochs,1))-1)%losing first trial for all stims 
                epoch = filteredsignal(1:round(1.4*fs)+1);
                baseline1 = round(.05*fs);
                baseline2 = round(.450*fs);
                baseline = epoch(baseline1:baseline2);
                MeanBaseline = mean(baseline); 
                epoch = epoch - MeanBaseline; 
                epochs(k,:) = epoch;
                filteredsignal(1:fs-1) = [];
            end
            checkforzeros = mean(epochs,2);
            cz = ismember(checkforzeros,0);
            zeroidx = find(cz == 1);
            if ~isempty(zeroidx)
                    epochs(zeroidx,:) = [];
                meanfiltsig = mean(epochs);
            else
                meanfiltsig = mean(epochs);
            end
            filtered_data.filteredmeans.(S).(R) = meanfiltsig;
            filtered_data.filteredepochs.(S).(R) = epochs;
        end
        filtered_data.filteredmeans.(S).(R) = meanfiltsig; %only at 60 Hz
        filtered_data.filteredepochs.(S).(R) = epochs;%only at 60 Hz
        clearvars zeroidx
    end%each response
end%each stim
save("filtered_data","filtered_data",'-v7.3');
