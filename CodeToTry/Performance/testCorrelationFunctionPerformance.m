try
    parpool('local')
catch
    disp('Pool already running')
end


for i = 1:1
    data = load(['sphericalHarmonics' num2str(i) '.mat']);
    sphericalHarmonic = data.sphericalHarmonic;
    for counter = 1:6
        sphericalHarmonic = [sphericalHarmonic; sphericalHarmonic];
    end
    clearvars -except i tocCase1 tocCase2 tocCase3 sphericalHarmonic tocFftSH tocFftCorrFunc tocCorrFunc tocSumCorrFunc tocConj
    disp(['Dataset: ' num2str(i)])
    lags = 500001;
    
    [~,timeSteps] = size(sphericalHarmonic);
    zeroPaddingLength = 2^(nextpow2(timeSteps)+1);
    
    memoryLinux();
    
    disp('CASE 1:')
    ticCase1 = tic;

    fftSH = fft(sphericalHarmonic,zeroPaddingLength,2);

    conjFftSH = conj(fftSH);
    
    fftCorrelationFunction = fftSH.*conjFftSH;
    clearvars fftSH conjFftSH
    
    correlationFunction = ifft(fftCorrelationFunction,[],2);
    clearvars fftCorrelationFunction
    
    sumCorrelationFunction1 = sum(correlationFunction(:,1:lags),1)/timeSteps;
    clearvars correlationFunction
    
    tocCase1(i) = toc(ticCase1);
    
    disp(['Overall time: ' num2str(tocCase1(i))])
    
    lags = 500001;
    [~,timeSteps] = size(sphericalHarmonic);
    zeroPaddingLength = 2^(nextpow2(timeSteps)+1);
    
    disp('CASE 2:')
    disp('parallel clearing')
    ticCase2 = tic;
    fftSH = fft(sphericalHarmonic,zeroPaddingLength,2);
    correlationFunction = ifft(fftSH.*conj(fftSH),[],2);
    clearvars fftSH
    sumCorrelationFunction2 = sum(correlationFunction(:,1:lags),1)/timeSteps;
    clearvars correlationFunction
    tocCase2(i) = toc(ticCase2);
    disp(['Overall time: ' num2str(tocCase2(i))])
      
    lags = 500001;
    [~,timeSteps] = size(sphericalHarmonic);
    zeroPaddingLength = 2^(nextpow2(timeSteps)+1);
    
    disp('CASE 3:')
    disp('variable clearing in background')
    ticCase3 = tic;
    % run with parallel toolbox
    memoryLinux();
    fftSH = fft(sphericalHarmonic,zeroPaddingLength,2);
    memoryLinux();
    correlationFunction = ifft(fftSH.*conj(fftSH),[],2);
    memoryLinux();
    parfeval(@clearvars,0,'fftSH');
    memoryLinux();
    sumCorrelationFunction3 = sum(correlationFunction(:,1:lags),1)/timeSteps;
    memoryLinux();
    parfeval(@clearvars,0,'correlationFunction');
    memoryLinux();
    tocCase3(i) = toc(ticCase3);
    disp(['Overall time: ' num2str(tocCase3(i))])
    
    disp('Differences between Cases: ')
    disp(['    Correlation Function Difference Case 1 and Case 2: ' ...
        num2str(sum(sum(abs(sumCorrelationFunction1) ...
        -abs(sumCorrelationFunction2))))])
    disp(['    Correlation Function Difference Case 1 and Case 3: ' ...
        num2str(sum(sum(abs(sumCorrelationFunction1) ...
        -abs(sumCorrelationFunction3))))])
    
    disp('===================================================')
    
    % no parallel fft with parfeval bacause parfeval runs in background
    
end
disp('ENDRESULT')
% disp('Average times needed: ');
% disp(['    FT of SH: ' num2str(mean(tocFftSH))]);
% % disp(['    Conj: ' num2str(mean(tocConj))]);
% disp(['    Multiplication fftSH fftSH*: ' num2str(mean(tocFftCorrFunc))]);
% disp(['    ifft of fft Correlation function: ' ...
%     num2str(mean(tocCorrFunc))]);
% disp(['    summation of correlation functions_ ' ...
%     num2str(mean(tocSumCorrFunc))]);

disp('Overall differences between cases: ');
disp(['    Case 1: ' num2str(mean(tocCase1)) ]);
disp(['    Case 2: ' num2str(mean(tocCase2)) ]);
disp(['    Case 3: ' num2str(mean(tocCase3)) ]);
