function [powers, maxPower, powerVsTime] = slmg_LFP_analysis(conf)
% Data extraction and LFP analysis around grooming events

% Experiments root path
geneData_path = 'Data';
% Number of animals
nAnimals = length(conf.animalList);
% Empty and zero initializations
sumMiceOFC = [];
nPowerOFC = 0;
powerVsTime.ofc = [];
powerVsFreq.ofc = [];
maxPower = struct();

% Loop through animal(s)
fprintf('> LFP analysis < \n')
for iAnimal = 1:nAnimals
    
    animal = conf.animalList{iAnimal};
    % Retrieves configuration for the current animal
    mouseConf = conf.mouseConf(animal);
    % Data for the current animal is loaded from the Data directory
    fprintf('Loading data for animal %s...', animal)
    %load([geneData_path '\extractedData\' animal '\' 'datasource_DE_' animal '.mat']);
    load([geneData_path '\extractedData\' animal '\' 'datasource_' animal '.mat']);
    fprintf('done!\n')
    
    % Initilization of variables
    data = [];
    totEvents = 0;
    minDuration = +Inf;
    clear detail;
    detail(1)=struct();
    
    %%  Data events extraction
    
    % Loop through each recording session for the current animal
    for i = 1:length(datasource.recording)
        % Check if the "session" is required
        if isfield(mouseConf, 'sessions')
            session_found=false;
            for is=1:length(mouseConf.sessions)
                s = mouseConf.sessions(is);
                if strfind(datasource.recording(i).directory, s)
                    session_found=true;
                end
            end
            if ~session_found
                continue;
            end
        end
        
        % Events around grooming onset are extracted
        switch conf.typeAnalysis
            case 'Gr'
                recEvents = datasource.recording(i).groomingEvents.event;
                evtName = 'Grooming';
            otherwise
                return
        end
        
        % Number of events in nEvents
        nEvents = length(recEvents);
        fprintf('Extracting %d events for this mouse. \n', nEvents)
        
        % Loop through events within each recording ression
        for iEv = 1:nEvents
            
            % Minimum duration of all events is tracked
            prevMinDuration = minDuration;
            
            % Data from events is concatenated, ensuring that all events have the same length
            minDuration = min(minDuration, size(recEvents(iEv).data, 2));
            if minDuration < prevMinDuration && prevMinDuration ~= +Inf
                data = [data(1:minDuration, :) recEvents(iEv).data(:, 1:minDuration)'];
            else
                data = [data recEvents(iEv).data(:, 1:minDuration)'];
            end
            
            % Tracking information relative to each event is saved in the detail structure
            detail(totEvents + iEv).recNum = i;
            P = strsplit(datasource.recording(i).directory, '\');
            dec = 0;
            
            if strcmp(P{end},'')
                dec = -1;
            end
            
            detail(totEvents + iEv).session = P{end + dec};
            detail(totEvents + iEv).animal = P{end-2 + dec};
            detail(totEvents + iEv).sessionLength = nEvents;
            detail(totEvents + iEv).iEvent = iEv;
            
        end
        
        totEvents = totEvents + nEvents;
    end
    
    %%  Wavelet Analysis
    
    % > Setting up the analysis environment <
    
    % Define the range of events to analyze
    firstEvent = 1;
    maxEvents = +Inf;
    % Sampling frequency
    params.Fs = datasource.recording(1).fs;
    
    % Analysis of grooming events
    switch conf.typeAnalysis
        % Retrieves grooming events data from datasource
        case 'Gr'
            recEvents = datasource.recording(1).groomingEvents.event;
            evtName = 'Grooming';
    end
    
    nChannels = size(recEvents(1).data, 1);
    % Time before and after an event
    t_before_datasource = datasource.config.times.before;
    t_after_datasource = datasource.config.times.after;
    
    % Adjust time window for wavelet representation
    t_before = min(datasource.config.times.before, conf.wavelet.t_before);
    t_after = min(datasource.config.times.after, conf.wavelet.t_after);
    
    analysisName = [animal ' (' evtName ')'];
    fprintf([analysisName ' - analysis start, per event\n']);
    
    % > Filter design: low pass <
    
    % Cutoff frequency
    fcut_hz = conf.fpass(2);
    
    %Target sampling rate (for LFP, the sf is usually 1kHz)
    newSamplingRate = fcut_hz*100;  % sampling rate after downsampling
    downSamplingCoef = ceil(params.Fs/newSamplingRate); % downsamplig coef
    newSamplingRate = round(params.Fs/downSamplingCoef);
    
    % Low pass Butterworth filter of order 10
    fnorm =fcut_hz/(newSamplingRate/2); % normalized cutoff frequency
    [b,a] = butter(10,fnorm,'low');
    
    % Initialization of result variables
    wholeMeanOFC = [];
    wholePowerOFC = [];
    
    countOFC = 0;
    % Iteration over each event
    for iEv = max(1, firstEvent):min(totEvents,maxEvents + max(1, firstEvent) -1)
        
        % Prints current event index and total number of events
        fprintf('Event %d of %d\n', iEv, totEvents);
        
        % Data extraction from all channels
        eventData = data(:, ((iEv-1)*nChannels+1):(iEv*nChannels));
        
        % > Downsampling  <
        
        % The extracted data is downsampled
        dataW = downsample(eventData, downSamplingCoef);
        
        % > Filtering  <
        
        data_filt1 = filtfilt(b, a, dataW);
        
        data_filt = data_filt1;

        % > Wavelet analysis: single event all channels <

        if conf.singleEvents
            figure; % New figure for each event
        end
        
        switch conf.typePlot
            case 'SELECTION_NORMALIZED'
                
                OFC_channels = data_filt(:, mouseConf.electrodes_OFC);
                mean_cfs = [];
                
                % Loop through channels for cwt calculation
                for iChannel=1:size(OFC_channels,2)
                    
                    % Continuous wavelet transform
                    [cfs, frq, coi] = cwt(OFC_channels(:, iChannel), conf.wavelet.type, newSamplingRate);
                    
                    %  The resulting wavelet coefficients are normalized and processed to keep only the
                    % relevant frequencies within the specified frequency range 
                    
                    % Coefficients outside the specified passband are
                    % excluded
                    cfs(frq>conf.fpass(2), :)=0;
                    cfs(frq<conf.fpass(1), :)=0;
                    
                    % Keep only time-window of interest (from -t_before to t_after)
                    fulltms = (0:size(data_filt, 1)-1)/newSamplingRate - t_before_datasource;
                    toKeep = (fulltms>=-t_before) & (fulltms<=t_after);
                    coi = coi(toKeep);
                    cfs = cfs(:, toKeep);
                    
                    % Normalizing the wavelet power across multiple
                    % channels
                    absCfs = abs(cfs).^2;
                    normCfs = absCfs/max(max(absCfs)); 
                    if numel(mean_cfs)>0
                        mean_cfs = mean_cfs+normCfs;
                    else
                        mean_cfs = normCfs;
                    end
                end
                
                % Mean across channels
                meanOFC=mean_cfs/size(OFC_channels,2); 
       
                % Normalization by event
                meanOFC = meanOFC/max(max(meanOFC));
                
                % Optionally the power within a specified window and
                % frequency band is computed 
                if (conf.powerConf.enable==1) && (size(OFC_channels,2)>0)
                    times=-t_before:(t_after_datasource+t_before)/(size(meanOFC,2)-1):t_after;
                    t_keep = ((times>=conf.powerConf.window(1)) .* (times<=conf.powerConf.window(2)))==1;
                    f_keep = ((frq>=conf.powerConf.bandPower(1))) .* (frq<=conf.powerConf.bandPower(2))==1;
                    evtPower = mean(mean(meanOFC(f_keep, t_keep)));
                end
                
                % Accumulates the mean wavelet power across grooming events
                % (wholeMeanOFC) and optionally accumulates the total power
                % across events (wholePowerOFC).
                if numel(wholeMeanOFC)>0
                    
                    % Checks if number of time points are different and if
                    % so, it ensures compatibility across events
                    oldSz = size(wholeMeanOFC,2);
                    newSz = size(meanOFC,2);
                  
                    if (oldSz~=newSz)
                        meanOFC = (resample(meanOFC', oldSz, newSz))';
                    end
                    
                    % Accumulates the mean wavelet power across events 
                    wholeMeanOFC = wholeMeanOFC+meanOFC;
                    % Increments the count of grooming events processed
                    countOFC=countOFC+1;
                    
                    if (conf.powerConf.enable) && (size(OFC_channels,2)>0)
                        wholePowerOFC = wholePowerOFC + evtPower;
                    end
                else
                    wholeMeanOFC = meanOFC;
                    countOFC=countOFC+1;
                    
                    if (conf.powerConf.enable) && (size(OFC_channels,2)>0)
                        wholePowerOFC = evtPower;
                    end
                end
                
                % Calculating and storing the maximum wavelet power within
                % a specifies frequency range for each event. 
                if conf.maxPower.enable && (size(OFC_channels,2)>0)
                    
                    % Time vector
                    times=-t_before:(t_after+t_before)/size(meanOFC,2):t_after;
                    % Focus on wavelet power values inside the specified
                    % frequency range to get maximum value
                    meanToAnalyse = meanOFC;
                    meanToAnalyse((frq<conf.maxPower.flim(1)) | (frq>conf.maxPower.flim(2)), :) = -inf;
                    [maxP, maxCoord] = maxOf(meanToAnalyse);
                    
                    maxPower.animal(iAnimal).animalName = animal;
                    maxPower.animal(iAnimal).event(iEv).ofc.value = maxP;
                    maxPower.animal(iAnimal).event(iEv).ofc.t = times(maxCoord(2));
                    maxPower.animal(iAnimal).event(iEv).ofc.f = frq(maxCoord(1));
                end
                
                tms = (0:size([meanOFC], 2)-1)/newSamplingRate - t_before;
                
                % > Wavelet Contour plot <
                nPlot = 0;
                iPlot=1;
                if numel(meanOFC)
                    nPlot=nPlot+1;
                end
                
                if numel(meanOFC) && conf.singleEvents
                    plotWaveletContour(meanOFC.^0.5, frq, tms, coi, t_before, t_after, conf.fpass, conf.powDb, 1);
                    title('OFC');
                    iPlot=iPlot+1;
                end
                
                % If other type of anaylsis, add here another case
            otherwise
        end
    end
    
    figure;
    
    nPlots = 0;
    if (numel(wholeMeanOFC))
        nPlots=nPlots+1;
    end
    
    if (numel(wholePowerOFC))
        nPowerOFC = nPowerOFC+1;
        powers.OFC(nPowerOFC).animal = animal;
        powers.OFC(nPowerOFC).value = wholePowerOFC/countOFC;
    end
    
    iPlot=1;
    
    if numel(wholeMeanOFC)
        subplot(nPlots,1,iPlot);
        if strcmp(conf.typePlot,'SELECTION_NORMALIZED')
            wholeMeanOFC = wholeMeanOFC/countOFC;
            wholeMeanOFC=wholeMeanOFC/max(max(wholeMeanOFC));
            plotWaveletContour((wholeMeanOFC).^0.5, frq, tms, coi, t_before, t_after, conf.fpass, conf.powDb, 1);
        else
            plotWaveletContour((wholeMeanOFC/countOFC).^0.5, frq, tms, coi, t_before, t_after, conf.fpass, conf.powDb);
        end
        title([animal ' OFC']);
        iPlot=iPlot+1;
        
        
        if numel(sumMiceOFC)
            
            oldSz = size(wholeMeanOFC,2);
            newSz = size(sumMiceOFC,2);
            if (oldSz~=newSz)
                sumMiceOFC = (resample(sumMiceOFC', oldSz, newSz))';
            end
            oldSz = size(wholeMeanOFC,1);
            newSz = size(sumMiceOFC,1);
            if (oldSz~=newSz)
                sumMiceOFC = (resample(sumMiceOFC, oldSz, newSz));
            end
            sumMiceOFC = sumMiceOFC + wholeMeanOFC;
            nMiceOFC = nMiceOFC + 1;
            
            matrix_wholeMeanOFC{1, nMiceOFC} = wholeMeanOFC;
        else
            sumMiceOFC = wholeMeanOFC;
            nMiceOFC = 1;
            matrix_wholeMeanOFC = {wholeMeanOFC};
        end
        
        % calculate the power vs time curve for the current mouse and
        % add it to the list of power vs time curves
        if (conf.powerVsTime.enable)
            fToKeep = (frq>=conf.powerVsTime.fband(1)) & (frq<=conf.powerVsTime.fband(2));
            currentPowerVsTime = 1/sum(fToKeep) * sum(wholeMeanOFC(fToKeep,:), 1);
            
            if (numel(powerVsTime.ofc))
                oldSz = size(wholeMeanOFC,2);
                newSz = size(powerVsTime.ofc,2);
                if (oldSz~=newSz)
                    powerVsTime.ofc = (resample(powerVsTime.ofc', oldSz, newSz))';
                end
                
                powerVsTime.ofc = [powerVsTime.ofc; currentPowerVsTime];
            else
                powerVsTime.ofc = currentPowerVsTime;
            end
        end
        
        % calculate the power vs freq curve for the current mouse and
        % add it to the list of power vs freq curves
        if (conf.powerVsFreq.enable)
            
            fToKeep = (frq>=conf.powerVsFreq.flim(1)) & (frq<=conf.powerVsFreq.flim(2));
            if isfield(conf.powerVsFreq, 'tval')
                tToKeep = [];
                for i=1:length(conf.powerVsFreq.tval)
                    tval = conf.powerVsFreq.tval(i);
                    [~, tm_tToKeep] = min(abs(tms - tval));
                    tToKeep = [tToKeep tm_tToKeep];
                end
                currentPowerVsFreq = wholeMeanOFC(fToKeep, tToKeep);
            else
                tToKeep = (tms>=conf.powerVsFreq.tlim(1)) & (tms<=conf.powerVsFreq.tlim(2));
                currentPowerVsFreq = 1/length(tToKeep) * sum(wholeMeanOFC(fToKeep, tToKeep), 2);
            end
            
            if (numel(powerVsFreq.ofc))
                powerVsFreq.ofc = [powerVsFreq.ofc; currentPowerVsFreq'];
            else
                powerVsFreq.ofc = currentPowerVsFreq';
            end
        end
        
        if (conf.powerConf.enable==1)
            % keep only the power during conf.powerConf.window
            times=-t_before:(t_after+t_before)/(size(wholeMeanOFC,2)-1):t_after;
            t_keep = ((times>=conf.powerConf.window(1)) .* (times<=conf.powerConf.window(2)))==1;
            f_keep = ((frq>=conf.powerConf.bandPower(1))) .* (frq<=conf.powerConf.bandPower(2))==1;
            powers.OFC(nPowerOFC).valueFromAverage = mean(mean(wholeMeanOFC(f_keep, t_keep)));
        end
        
        if conf.maxPower.enable
            % get the max value and max coordinates
            times=-t_before:(t_after+t_before)/size(wholeMeanOFC,2):t_after;
            % set power out of frequency band expected to -inf
            meanToAnalyse = wholeMeanOFC;
            meanToAnalyse((frq<conf.maxPower.flim(1)) | (frq>conf.maxPower.flim(2)), :) = -inf;
            [maxP, maxCoord] = maxOf(meanToAnalyse);
            
            maxPower.animal(iAnimal).animalName = animal;
            maxPower.animal(iAnimal).ofc.value = maxP;
            maxPower.animal(iAnimal).ofc.t = times(maxCoord(2));
            maxPower.animal(iAnimal).ofc.f = frq(maxCoord(1));
            
            if conf.maxPower.show
                plot(maxPower.animal(iAnimal).ofc.t,maxPower.animal(iAnimal).ofc.f,'k+');
            end
        end
    end
    
end

%% Display figures
% Generates multiple plots representing different analyses of data related to
% the orbitofrontal cortex (OFC), including wavelet contour plots, power vs time plots,
% and power vs frequency plots.

figure;

% Count the number of plots to be created based on nb. subjects
nPlots = 0;
if (numel(sumMiceOFC))
    nPlots=nPlots+1;
end

% Plots wavelet contour and finds maximum power point within a frequency
% band
if numel(sumMiceOFC)
    meanMiceOFC = sumMiceOFC/nMiceOFC;
    
    if strcmp(conf.typePlot,'SELECTION_NORMALIZED')
        meanMiceOFC=meanMiceOFC/max(max(meanMiceOFC));
        plotWaveletContour((meanMiceOFC).^0.5, frq, tms, coi, t_before, t_after, conf.fpass, conf.powDb, 1);
    else
        plotWaveletContour(meanMiceOFC.^0.5, frq, tms, coi, t_before, t_after, conf.fpass, conf.powDb);
    end
    
    % Computes maximum power along with its corresponding time and
    % frequency coordinates
    if conf.maxPower.enable
        % generates time vector
        times=-t_before:(t_after+t_before)/size(meanMiceOFC,2):t_after;
        % set power out of frequency band expected to -inf
        meanToAnalyse = meanMiceOFC;
        meanToAnalyse((frq<conf.maxPower.flim(1)) | (frq>conf.maxPower.flim(2)), :) = -inf;
        [maxP, maxCoord] = maxOf(meanToAnalyse);
        
        maxPower.wholeMean.ofc.value = maxP;
        maxPower.wholeMean.ofc.t = times(maxCoord(2));
        maxPower.wholeMean.ofc.f = frq(maxCoord(1));
        
    end
    
    title('OFC - All Mice');
end

% Plot power vs time
if (numel(powerVsTime.ofc)>0)
    figure;
    
    if conf.powDb
        powerVsTime.ofc = pow2db(powerVsTime.ofc);
    end
    
    if (numel(powerVsTime.ofc))
        
        % calculate the mean SNR ration over events
        mean_pvt = mean(powerVsTime.ofc, 1);
        
        % calculate the standard deviation of SNR over events
        se_mean_pvt = std(powerVsTime.ofc, 1,  1)/sqrt(size(powerVsTime.ofc,1));
        
        % smooth the values
        mean_pvt_smooth = smooth(mean_pvt, 0.01,'rloess' );
        sd_mean_pvt_smooth = smooth(se_mean_pvt, 0.01,'rloess');
        
        boundedline(tms, mean_pvt_smooth, sd_mean_pvt_smooth);
        
        % set title
        title(['Power vs time - OFC']);
        
        xlabel('t (s)');
        
        axOfc = gca;
    end
    
end

% Plot power vs frequency
if (numel(powerVsFreq.ofc) >0)
    % Frequency range of interest
    f = frq((frq>=conf.powerVsFreq.flim(1)) & (frq<=conf.powerVsFreq.flim(2)));
    
    if conf.powDb
        powerVsFreq.ofc = pow2db(powerVsFreq.ofc);
    end
    
    if isfield(conf.powerVsFreq, 'tval')
        nPlotGroup = length(conf.powerVsFreq.tval);
        plotGroupOFC = mod(0:size(powerVsFreq.ofc,1)-1, nPlotGroup) + 1;
    else
        nPlotGroup = 1;
        plotGroupOFC = ones(1, size(powerVsFreq.ofc,1));
    end
    
    for iPlotGroup = 1:nPlotGroup
        figure;
        
        if isfield(conf.powerVsFreq, 'tval')
            titleComplement = [' @' num2str(conf.powerVsFreq.tval(iPlotGroup)) 's'];
        else
            titleComplement = '';
        end
        
        if (numel(powerVsFreq.ofc))
            powerVsFreq_ofc_data = powerVsFreq.ofc(plotGroupOFC==iPlotGroup, :);
            
            % Mean and std error of the power values
            mean_pvf = mean(powerVsFreq_ofc_data, 1);
            se_mean_pvf = std(powerVsFreq_ofc_data, 1,  1)/sqrt(size(powerVsFreq_ofc_data,1));
            
            % Smooths using a robust local regression methos
            mean_pvf_smooth = smooth(mean_pvf, 0.01,'rloess' );
            sd_mean_pvf_smooth = smooth(se_mean_pvf, 0.01,'rloess');
            
            % Plots the mean power vs frequency along with shaded error
            % bars
            boundedline(f, mean_pvf_smooth, sd_mean_pvf_smooth);
            
            % set title
            title(['Power vs frequency - OFC' titleComplement]);
            
            xlabel('Frequency (Hz)');
            
        end
        
    end
    
end

end

%% Functions

function [] = plotWaveletContour(cfs, frq, tms, coi, t_before, t_after, fpass, powDb, maxCAxis)
% Creates a contour plot of continuous wavelet transform (cwt) coefficients
% over time and frewuency.
% cfs - CWT coefficients
% frq - frequency vector
% tms - time vector
% fpass - frequency passband
% powDB - flag for decibels
% maxCAxis - max value for the color axis

if powDb
    cfs = pow2db(cfs);
end

% Opengl renderer to improve rendering quality
set(gcf,'renderer','opengl');
% Creates a surfer plot of the CWT coefficients over time and frequency
helperCWTTimeFreqPlot(cfs,tms,(frq), 'surf');

% Modify the x-axis tick labels to display time (s) ranging from t_before -
% t_after
get(gca,'XTickLabel');
xticks(-t_before:t_after);
set(gca,'XTickLabel',num2str((-t_before:t_after)')) % x labels every seconds from -t-before to t_after_datasource
xlabel('Time (s)');

% Vertical line at the beginning of the grooming bout
l=vline (0,'m');
set(l, 'LineWidth', 1.5);   % set the linewidth

colormap(jet);

hold on;

% Plots the cone of influence
plot(tms,coi,'k--','linewidth',2);

% Sets the y-axis to the frequency passband specified by fpass
ylim(fpass);
% Sets the x-axis time limits
xlim([-t_before t_after]);

% Checks if maxCAxis is provided
if nargin>8
    % Sets the color axis limits
    caxis([0 maxCAxis]);
end

end

function [val, coord] = maxOf(A)

[m, idx] = max(A);
[val, idx2] = max(m);

coord = [idx(idx2(1)), idx2(1)];
end
