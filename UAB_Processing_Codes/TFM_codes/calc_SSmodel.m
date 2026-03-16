function [A,B,C] = calc_SSmodel_UAB(data,sampBeforeStim,stimLength)
% INPUTS: data = an n x m matrix where n is the number of response channels and m is the number
% of time points. The matrix is the average evoked response response
% for each channel to a given stimulating channel. You need to take out all
% the NaNs beforehand!
%         sampBeforeStim = the number of samples from the start to the stimulus
% sample (e.g. 0.5*fs)
%         stimLength = the number of samples you want your model to use to
% simulate stimulation (I have always used 2).

        % create stimulation vector, using input stimLength
        u = zeros(1,size(data,2));
        u(sampBeforeStim:sampBeforeStim+stimLength-1) = 1; %this is so it is centered on the "stimulus onset point"

        bMatrix = data(:,2:end); %b is the "predicted" value that we are going to estimate A from. 
        b = sparse(double(bMatrix(:))); %make sparse (I don't know why this is necessary, but it is.)

        % build sparse H matrix directly %I build this using the method from Adam's 2017 paper
        % tic;
        N = size(data,1);
        T = size(data,2);
        chanStart = 1:N:(T-1)*N;
        H_x = sparse([],[],[],(T-1)*N,N^2);
        for j=1:N 
            rows = repmat(chanStart+(j-1),N,1); 
            cols = repmat((j-1)*N+1:j*N,(T-1),1)';
            dataTemp = data(:,1:end-1);
            H_x_temp = sparse(rows(:),cols(:),dataTemp(:),(T-1)*N,N^2);
            H_x = H_x + H_x_temp;
%             disp(j);
        end
        % toc;
        % This comes from the method in the Proctor et al. 2016 paper. See equation
        % 3.14 where [A B] = X' [X
        %                        U]; 
        inputData = sparse([],[],[],size(H_x,1),N);
        for j = 1:T-1
            for k = 1:N
                row = (j-1)*N+k;
                col = k;
                inputData_temp = sparse(row,col,u(j),size(H_x,1),N);
                inputData = inputData+inputData_temp;
            end
        end

        H = [H_x inputData];

        clear H_x inputData

        AB = H\b; %this is where you do the least squares estimation

        Atotal = reshape(AB(1:end-N),[N,N])'; %re-separate A and B
        Btotal = AB(end-N+1:end);
            
        A = full(Atotal);
        B = full(Btotal);
        
        x_hat = nan(size(data)); %initialize
        x_hat_iter = zeros(size(data,1),1);
        for j = 1:size(data,2) %loop through every time point in the window
            x_hat_temp = A*x_hat_iter+B*u(j); %in first iteration, A is multiplied by x0, but then for every successive iteration, A is multiplied by the previous time point.
            x_hat(:,j) = x_hat_temp; %assign predicted value to time series
            x_hat_iter = x_hat_temp; %reset "initial condition" to be the current time point.
        end
        
        scaleFactor = nan(size(data,1),1);
        for j = 1:size(data,1)
            scaleFactor(j) = (max(data(j,:))-min(data(j,:)))/(max(x_hat(j,:))-min(x_hat(j,:)));
        end
        C = diag(scaleFactor);
       
end