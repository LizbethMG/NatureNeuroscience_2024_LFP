% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>|
% LFP_Analysis
% INPUT
%   config.m configuration file of the experiment(s)
%   Data/ data files (can be downloaded in the corresponding OSF folder)
% DATA PROCESSING: 
    % Loads data for each experiment, and prepares it for further analysis. It ensures that the data is standardized across events and sessions to facilitate consistent analysis.
    % Low pass filtering
    % Wavelet analysis
    % Compute power within a specified frequency band
    % Compute maximum wavelet power within specified frequency band
    % Visualization
    %   plotWaveletContour Creates a contour plot of continuous wavelet transform (cwt) coefficients over time and frequency.
% OUTPUT
%   Generates different plots within a frequency band:
%       Wavelet contour
%       Power vs time
%       Power vs frequency
%       maxPower Maximum power and its corresponding time and frequency values
% Subfunctions: Folders including fcn_dataExtraction and fcn_thridparty
% MAT-files required: none
%
% +++++++++++++++++++++++++++++++++
% Author: SLMG
% GitHub: LizbethMG

%  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<|

close all;
clc;
clear all;

%% 1) Load the configuration for Data Extraction (confDE)

% additional duration in second to have no useful data outside the cone of influence
% of the wavelet
confDE.additionalDuration = 1.5; 
config %configuration file
%% 2) Analysis and figures

[powers, maxP, powerVsTime]= slmg_LFP_analysis(confWV);
