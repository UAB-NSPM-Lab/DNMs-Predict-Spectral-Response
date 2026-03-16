function [out, base, Arts] = UAB_removeartifacts(out)
%code to remove artifacts from NYU data
fs = out.fs;
stimlist = fieldnames(out.epoch_data); 
stim = string(stimlist(1));
resplist = string(fieldnames(out.epoch_data.(stim)));
stimmedPairChannels = zeros(size(resplist));
artChannels = zeros(size(resplist));
ignoreChannels = zeros(size(resplist));
for o = 1:length(resplist)%for each stim channel
    rlabel = resplist(o);
    if contains(rlabel,'_C')
        resplist(o) = "nan";
        ignoreChannels(o) = 1;
    elseif contains(rlabel,'_DC')
        resplist(o) = "nan";
        ignoreChannels(o) = 1;
    elseif contains(rlabel,'_SKIP')
        resplist(o) = "nan";
        ignoreChannels(o) = 1;
    elseif contains(rlabel,'_ECG')
        resplist(o) = "nan";
        ignoreChannels(o) = 1;
    elseif contains(rlabel,'_EKG')
        resplist(o) = "nan";
        ignoreChannels(o) = 1;
    end
end
for i = 1:length(stimlist)%each stim
    S = string(stimlist(i));
    S1_2 = strsplit(S,'_');
    S1 = S1_2(1);
    S2 = S1_2(2);
    SMChannels = zeros(size(resplist,1),1);
    meanBaseArt = zeros(size(resplist,1),1);
    residStimArt = zeros(size(resplist,1),1);
    stdBaseArt = zeros(size(resplist,1),1);
    stdAll = zeros(size(resplist,1),1);
    for j = 1:length(resplist)%each response
        R = string(resplist(j));
        if contains(R,"nan")
            continue;
        end
        %remove stim artifacts on trial by trial basis
        [mean_data_trial,datatrials] = removestimart_getmeanresponse(out.epoch_data.(S).(R),fs); 
        out.means.(S).(R) = mean_data_trial;
        out.removedart_trials.(S).(R) = datatrials;
        % Check for wonky baselines
        mean0 = mean(mean_data_trial(round(0.05*fs):round(0.45*fs))); %from 50ms to 450 ms
        mean1 = mean(mean_data_trial(round(end-0.45*fs):round(end-0.05*fs))); %from 400 ms before end to end
        if mean1 > (mean0 + 200) || mean1 < (mean0 - 200)
            meanBaseArt(j) = 1;
        else
            meanBaseArt(j) = 0;
        end
        %checking for residual stim artifacts
        maxVal = max(abs(mean_data_trial)); % - mean(mean_data_trial));
        if maxVal > 1500%3*mean_data_trial_std %trying to make this less arbitrary; was '8000'
            residStimArt(j) = 1;
        end
        %check std
        std_trial = std(out.epoch_data.(S).(R)(2:end,:));
        std_trial(round(0.494*fs):round(0.515*fs)) = 1; %accounts for stim artifact
        stdMed = median(std_trial);
        if stdMed>800
            SMChannels(j) = 1;
        end
        std0 = mean(std_trial(round(0.05*fs):round(0.45*fs))); %from 50ms to 450 ms
        std1 = mean(std_trial(round(end-0.45*fs):round(end-0.05*fs))); %from 400 ms before end to end
        if std1 > (std0 + 200) || std1 < (std0 - 200)
            stdBaseArt(j) = 1;
        else
            stdBaseArt(j) = 0;
        end
        largeStd0 = find(std_trial(1:round(0.494*fs))>1000);
        largeStd1 = find(std_trial(round(0.6*fs):end)>1000);
        largeStd = [largeStd0 largeStd1];
        if size(largeStd, 2) > 10
            stdAll(j) = 1;
        end
    end %each response
    for l = 1:length(resplist)
        R = string(resplist(l));
        R1_2 = strsplit(R,'_');
        if sum(ismember(R1_2,S1)) >= 1
            stimmedPairChannels(l) = 1;
        elseif sum(ismember(R1_2,S2)) >=1
            stimmedPairChannels(l) = 1;
        else
            continue;
        end
    end
    artChannels = artChannels + stimmedPairChannels;
    artChannels = artChannels + ignoreChannels;
    artChannels = artChannels + SMChannels;
    artChannels = artChannels + meanBaseArt;
    artChannels = artChannels + residStimArt;
    artChannels = artChannels + stdBaseArt;
    artChannels =   artChannels + stdAll;
    artChannels(artChannels>=1) = 1;
    artChannels = logical(artChannels);
    Arts.(S) = artChannels;
    artChannels = zeros(size(resplist));
    stimmedPairChannels = zeros(size(resplist));
end %each stim
% %%Plot artifactual channels
% name = string(out.name);
% S = fieldnames(Arts);
% M = zeros(length(S),length(Arts.(string(S(1)))));
% for i = 1:length(S)
%     M(i,:) = Arts.(string(S(i)));
% end
% figure; imagesc(M); title((name));
%% Save
save('ArtifactualChannels',"Arts",'-v7.3');
save('out',"out",'-v7.3');
base.name = out.name;
base.bipolarresplabels = out.bipolarresplabels;
base.fs = out.fs;
base.stimlabels = out.stimlabels;
save("base","base",'-v7.3');