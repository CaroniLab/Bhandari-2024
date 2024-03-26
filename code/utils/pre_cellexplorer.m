clear
clc
%% load session / only run first time
dir = uigetdir();

%session = sessionTemplate('W:\nobackup\garber\kanohars\HA\043_02\CellExplorer\43_02');
session = sessionTemplate(dir,  'showGUI',true);
%% Run cell metrics
cell_metrics = ProcessCellMetrics('session', session,'showGUI',true);
%% Launch Cell Explorer
CellExplorer('metrics',cell_metrics)
%% Launch Neuroscope
NeuroScope2('session',session)