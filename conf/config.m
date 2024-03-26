% Flags and variable definition used in the wavelet analysis function
% (slmg_LFP_analysis.m)

% Pre-grooming
confWV.filePrefix = '';

% Type of events to analyse (Gr = grooming)
confWV.typeAnalysis = 'Gr';
% Filtering
confWV.fpass = [1.5 10];
% Wavelet configuration
confWV.wavelet.type = 'morse';

% For pre-grooming:
confWV.wavelet.t_before = 2;
confWV.wavelet.t_after = 2;

% Power calculation configuration
confWV.powerConf.enable = 1;    % 1 to enable
confWV.powerConf.bandPower = [1.5 3];  

% For pre-grooming:
confWV.powerConf.window = [-1 0];

% configurations specific to each mouse
confWV.mouseConf = containers.Map('KeyType','char','ValueType','any');

% Power in dB
confWV.powDb = 0;

% Power vs. time
confWV.powerVsTime.enable = 1;
confWV.powerVsTime.fband = [1.5 4];  

% Power vs freq
confWV.powerVsFreq.enable = 1;
confWV.powerVsFreq.flim = [1.5 10];

% For pre-grooming:
confWV.powerVsFreq.tlim = [-1 0];

% Get max power value
confWV.maxPower.enable = 1;
confWV.maxPower.flim = [1.5 10];
confWV.maxPower.show = 1;


% List of mice
% ***********************************
confWV.animalList = { 'M1'};

% selection of type of plot that we want
confWV.typePlot = 'SELECTION_NORMALIZED';

confWV.singleEvents = 0; % If 1: single events + mean, If 0: only mean for single animal and meann for all animals


% Animals configuration

% ************** M1 *****************
confWV.mouseConf('M1')=struct( ...
    'electrodes_OFC', [5 9 21 30], ...    % electrodes selection in OFC
    'electrodes_', [], ...       % discarded electrodes
    'all_OFC', [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 ], ...    % all electrodes in OFC
    'all_', []);  % % discarded electrodes

% ************** M2 *****************
confWV.mouseConf('M2')=struct( ...
    'electrodes_OFC', [1 6 21], ...    % electrodes selection in OFC 
    'electrodes_', [], ...       % discarded electrodes
    'all_OFC', [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 ], ...    % all electrodes in OFC
    'all_', []);  % % discarded electrodes

% ************** M3  *****************
confWV.mouseConf('M3')=struct( ...
    'electrodes_OFC', [1 7 10 13 18 23], ...    % electrodes selection in OFC
    'electrodes_', [], ...       % discarded electrodes
    'all_OFC', [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 ], ...    % all electrodes in OFC
    'all_', []);  % % discarded electrodes

% ************** M4  *****************
confWV.mouseConf('M4')=struct( ...
    'electrodes_OFC', [8 14 18 21], ...    % electrodes selection in OFC
    'electrodes_', [], ...      % discarded electrodes
    'all_OFC', [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 ], ...    % all electrodes in OFC
    'all_', []);  % discarded electrodes

% ************* M5 ******************
confWV.mouseConf('M5')=struct( ...
    'electrodes_OFC', [1 5 20], ...    % all electrodes in OFC
    'electrodes_', [ ], ...        % discarded electrodes
    'all_OFC', [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20], ...    % all electrodes in OFC
    'all_', [], ... % discarded electrodes
    'sessions', ["S2" "S3" "S6" "S7" "S10" "S8"]); 

% ************** M6  *****************
confWV.mouseConf('M6')=struct( ...
    'electrodes_OFC', [1 8 15 20 32], ...    % electrodes selection in OFC
    'electrodes_', [], ...       % discarded electrodes
    'all_OFC', [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 ], ...    % all electrodes in OFC
    'all_', []);  %  discarded electrodes


% ************** M7  *****************
confWV.mouseConf('M7')=struct( ...
    'electrodes_OFC', [5 9 22], ...    % electrodes selection in OFC
    'electrodes_', [], ...       % discarded electrodes
    'all_OFC', [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 ], ...    % all electrodes in OFC
    'all_', [], ...% discarded electrodes
    'sessions', ["S1" "S2" "S3_1" "S3_2" "S4" "S5" "T1"]);

% ************** M8  *****************
confWV.mouseConf('M8')=struct( ...
    'electrodes_OFC', [2 10 14 21], ...    % electrodes selection in OFC (removed 18)
    'electrodes_', [], ...       % discarded electrodes
    'all_OFC', [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 ], ...    % all electrodes in OFC
    'all_', []); % discarded electrodes
% 'sessions', ["S3" "S4" ]);

% ************** M9  *****************
confWV.mouseConf('M9')=struct( ...
    'electrodes_OFC', [1 20 21 22 29 ], ...    % electrodes selection in OFC
    'electrodes_', [], ...       % discarded electrodes
    'all_OFC', [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 ], ...    % all electrodes in OFC
    'all_', []);  % discarded electrodes

% *************** M10 ****************
confWV.mouseConf('M10')=struct( ...
    'electrodes_OFC', [1 5 15 ], ...    % electrodes selection in OFC
    'electrodes_', [], ...       % discarded electrodes
    'all_OFC', [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 ], ...    % all electrodes in OFC
    'all_', [], ...  % discarded electrodes
    'sessions', ["S11" "S12" "S14"]);
