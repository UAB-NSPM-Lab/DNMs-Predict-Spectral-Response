function UAB_Bode(pt, avg_waveforms,code_path,out_path)

%SEEG only, will write additional script for GRID
fs = avg_waveforms.fs;
dirCode = append(code_path,'/TFM_codes'); %path to modeling scripts
StimNames = fieldnames(avg_waveforms.average_waveforms);%stims 
CC_array = nan(size(StimNames,1),length(avg_waveforms.resplabels));%for calculating correlation coefficients to determine model performance
for q = 1:size(StimNames,1) %each stim
    disp(q);%displays stim number to track progress 
    stim_data = avg_waveforms.average_waveforms.(StimNames{q});
    nanidx = isnan(stim_data(1,:));
    stim_data(:,nanidx) = [];%%%
    dataSEEG = stim_data';%%%
    cd(dirCode);
    stim_onset = ceil(fs/2);
    [A,B,C] = calc_SSmodel(dataSEEG,stim_onset,2);
    savePlace = append(out_path,'/',pt,'/TFMsBode_figures/ModelParams');
    if q == 1
        mkdir(savePlace);
    end
    cd(savePlace);
    save([StimNames{q} '_modelParams_SEEG'],'A','B','C','fs');
    
    cd(dirCode);
    [~,~,freqs_log,singValues_log,~,~,bode_matrix,freq_matrix] = calc_BodePlotParams(A,B,C,fs,-1,3,100);
    
    %This plots the Sigma plot
    plotBode(freqs_log,singValues_log,'k');
    title(string(append(StimNames{q},'_Sigma')),'Interpreter','none');
    savePlace = append(out_path,'/',pt,'/TFMsBode_figures/Sigma');%'C:\Users\rache\OneDrive - Johns Hopkins\Documents\UAB\Research\CCEPs Code\OneDrive_2023-09-20\Code for Erik Gommer\P0123_figures\';
    if q == 1
        mkdir(savePlace);
    end
    cd(savePlace)
    savefig([StimNames{q} '_sigma_SEEG']);
    close
    cd(dirCode)
    
    %This plots the True Bode plots (one plot for every stim-response pair)
    n_subplots = ceil(sqrt(size(A,1)));
    good_response_names = avg_waveforms.resplabels(~nanidx);
    figure;
    for i = 1:size(A,1) %every response
        subplot(n_subplots,n_subplots,i);
        plot(log(freq_matrix{i}),log(bode_matrix{i}));
        title(string(good_response_names(i)),'Interpreter','none');
    end
    sgtitle(string(append(StimNames{q},'_Bode')),'Interpreter','none');
    savePlace = append(out_path,'/',pt,'/TFMsBode_figures/Bode');
    if q == 1
        mkdir(savePlace);
    end
    cd(savePlace);
    savefig([StimNames{q} '_Bodes_SEEG']);
    save([StimNames{q} '_Bode_Matrices'],'bode_matrix','freq_matrix');
    cd(dirCode);
    close

    % Plot Model Reconstructions 
    u = zeros(1,size(dataSEEG,2));
    u(stim_onset:stim_onset+1) = 1; %this is so it is centered on the "stimulus onset point"
    %Intitialize
    x_hat = nan(size(dataSEEG));
    x_hat_iter = zeros(size(dataSEEG,1),1);
    for n = 1:size(dataSEEG,2) %loop through every time point in the window
         x_hat_temp = A*x_hat_iter+B*u(n); 
         x_hat(:,n) = x_hat_temp.*diag(C); %assign predicted value to time series
         x_hat_iter = x_hat_temp; %reset "initial condition" to current time point
    end
    %To plot output:
    figure;
    nChans = size(A,1);
    nSubPlots = ceil(sqrt(nChans));
    for j = 1:size(dataSEEG,1)
        meanTrialIter = dataSEEG(j,:);
        meanmean = mean(meanTrialIter);
        stdTrialIter = ones(size(meanTrialIter));
        subplot(nSubPlots,nSubPlots,j);            
        fill([1/fs:1/fs:1 1:-1/fs:1/fs 1/fs],[(meanTrialIter(1:1*fs)-meanmean)+(stdTrialIter(1:1*fs)-meanmean) fliplr((meanTrialIter(1:1*fs)-meanmean)-(stdTrialIter(1:1*fs)-meanmean)) (meanTrialIter(1)-meanmean)+(stdTrialIter(1)-meanmean)],'y');
% fill([1/fs:1/fs:2 2:-1/fs:1/fs 1/fs],[(meanTrialIter(1:2*fs)-meanmean)+(stdTrialIter(1:2*fs)-meanmean) fliplr((meanTrialIter(1:2*fs)-meanmean)-(stdTrialIter(1:2*fs)-meanmean)) (meanTrialIter(1)-meanmean)+(stdTrialIter(1)-meanmean)],'y');
        hold on;
        plot((1:1*fs)./fs,meanTrialIter(1:1*fs)-meanmean,'k','LineWidth',3);
        
        
        cc1 = meanTrialIter(1:1*fs)-meanmean; 
% plot((1:2*fs)./fs,meanTrialIter(1:2*fs)-meanmean,'k','LineWidth',3);
    %     title(parameters.ChannelNames.Value{j});
    %     ylim([-3000 3000]);
        hold on;
        x_hat_vec = x_hat(j,:);
        plot((1:1*fs)./fs,x_hat_vec(1:1*fs),'b','LineWidth',3); 
        cc2 = x_hat_vec(1:1*fs);
% plot((1:2*fs)./fs,x_hat_vec(1:2*fs),'b','LineWidth',3);
        set(gca,'xtick',[]);
%             title(['Channel: ' parameters.ChannelNames.Value{j}]);
%             xlabel('Time (sec)');
%             ylabel('Amp (uV)');
        CC = corrcoef(cc1(ceil(end/2):end),cc2(ceil(end/2):end));
        resp_idx = ismember(string(avg_waveforms.resplabels),string(good_response_names(j)));
        CC_array(q,resp_idx) = CC(2,1);%for each stim, vector of response channel corrcoef to model recon
        title(append(string(good_response_names(j)),' - ',string(CC(2,1))),'Interpreter','none');
    end
    sgtitle(string(append(StimNames{q},'_ModelReconstructions')),'Interpreter','none');
    savePlace = append(out_path,'/',pt,'/TFMsBode_figures/ModelReconstructions');
    if q == 1
        mkdir(savePlace);
    end
    cd(savePlace);
    savefig([StimNames{q} '_reconstruct_SEEG']);
    cd(dirCode);
    close
end
% figure; imagesc(CC_array);colorbar;
T = append(pt,' - Model Reconstruction Correlation Coefficient Matrix'); 
ylabel('Stimulated Channels'); xlabel('Response Channels');
title((T))
cd(savePlace);
save('CC_array',"CC_array");
save('StimNames','StimNames');
% figure; boxplot(CC_array');