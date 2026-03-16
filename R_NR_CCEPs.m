%% Responsive vs non-responsive CCEP identifier
clearvars 
cd /data/project/NSPMlab/hbriny99/UAB_CCEPdata/CCSR_Resonance_Comparison
load("PercentFreqMatchAllRand.mat");
pt_list = string(fieldnames(PercentFreqMatchAllRand));
thresh = 100;
for k = 1:length(pt_list)
    pt = pt_list(k);
    cd(append('/data/project/NSPMlab/hbriny99/UAB_CCEPdata/Output/',pt));
    load(append(pt,'_short.mat'));
    load("ArtifactualChannels.mat")
    fs = short.fs;
    stims = string(fieldnames(short.average_waveforms));
    resps = string(short.resplabels);
    CCEPs.(pt).index = zeros(length(stims),length(resps));
    for i = 1:length(stims)
        S = stims(i);
        for j = 1:length(resps)
            if Arts.(S)(j) == 1
                CCEPs.(pt).index(i,j) = nan;
            else
                ccep = short.average_waveforms.(S)(:,j);
                baselineccep = ccep(ceil(0.05*fs):ceil(0.45*fs));
                poststimccep = ccep(ceil(0.515*fs):fs);
                baselineamp = mean(baselineccep);
                poststimamp = max(poststimccep);
                if abs(poststimamp-baselineamp) > thresh
                    CCEPs.(pt).(S)(:,j) = ccep;
                    CCEPs.(pt).index(i,j) = 1;
                end
            end
        end
    end
    rccepstims = string(fieldnames(CCEPs.(pt)));
    sq = ceil(sqrt(length(rccepstims)));
    figure;
    for i = 1:length(rccepstims)
        S = string(rccepstims(i));
        for j = 1:size(CCEPs.(pt).(S),2)
            ccep = CCEPs.(pt).(S)(:,j);
            if mean(ccep) == 0
                continue
            else
                subplot(sq,sq,i);
                hold on
                plot(ccep)
                title(S)
            end
        end
    end
    sgtitle(string(short.name));
    clearvars -except pt_list thresh CCEPs
end
cd /data/project/NSPMlab/hbriny99/UAB_CCEPdata/CCSR_Resonance_Comparison
save('CCEPs',"CCEPs");