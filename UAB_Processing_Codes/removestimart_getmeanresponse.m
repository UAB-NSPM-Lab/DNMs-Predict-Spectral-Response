%function to remove stim arts for each trial and to get mean response

function [mean_data_trial,datatrials] = removestimart_getmeanresponse(stimepochdata,fs)

% figure;
% datatrials = zeros(size(stimepochdata));
% removeidx = zeros(size(stimepochdata,1),1);
kk = 0;
for k = 2:size(stimepochdata,1) %each trial, skipping first in case there are zeros: P0130 work around
    data_trial = stimepochdata(k,:); %x
%     plot(data_trial); hold on;
    y2 = data_trial;%15 ms
%     if sum(data_trial) == 0 %Changed for P0130 only
% %         removeidx(k) = 1;
% %         kk = k+1;
%         continue;
%     end
    kk = kk+1;
%     if kk == 1
%         continue;
%     end
    stim_onset = floor(0.5*fs)-3; %dp at 500 ms %start and 3 dp before for better measure
    %% y2
    lapse = ceil(0.015*fs)+3;%dp %duration, 15 ms 
    stim_end = stim_onset+lapse;
    pre_art_zone = y2(stim_onset-lapse:stim_onset-1); %get outer edges
    pre_art_zone_flip = fliplr(pre_art_zone); %flip them
    pre_art_zone_linspace = pre_art_zone_flip.*linspace(1,0,lapse); %weight them
    post_art_zone = y2(stim_end+1:stim_end+lapse);
    post_art_zone_flip = fliplr(post_art_zone);
    post_art_zone_flip_linspace = post_art_zone_flip.*linspace(0,1,lapse);
    y2(stim_onset:stim_end-1) = pre_art_zone_linspace + post_art_zone_flip_linspace;
    datatrials(kk,:) = y2;
end
% get mean 
% adjust_factor = 0; 
% for i = 1:size(datatrials,1) 
%     if removeidx(i) == 1
%         datatrials(i-adjust_factor,:) = []; 
%         adjust_factor = adjust_factor + 1;
%     end 
% end
% datatrials(removeidx,:) = [];
if exist("datatrials")
    mean_data_trial = mean(datatrials(1:end,:));
else
    mean_data_trial = nan(size(stimepochdata,2),1);
end