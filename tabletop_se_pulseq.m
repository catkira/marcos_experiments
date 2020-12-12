addpath('../pulseq/matlab')
close all
clear
gamma = 42.57E6;
sequencerRasterTime = 7E-9; % make sure all times are a multiple of sequencer raster time

TE = 12E-3;
Nx = 128;
samplingDuration = 2.5e-3;

% set system limits
rfDeadTime = 500e-6; % [us], minicircuits PA needs 500 us to turn on
adcDeadTime = 0;
sys = mr.opts('rfDeadTime', round(rfDeadTime/sequencerRasterTime)*sequencerRasterTime, ...
    'adcDeadTime', round(adcDeadTime/sequencerRasterTime)*sequencerRasterTime, ...
    'rfRasterTime', 1.003e-6, 'gradRasterTime',1.003e-6);
seq=mr.Sequence(sys);              % Create a new sequence object

% Create HF pulses
rf90duration=0.08e-3
rf90 = mr.makeBlockPulse(pi/2, 'duration', rf90duration, 'sys', sys);
rf180 = mr.makeBlockPulse(pi, 'duration', rf90duration*2, 'sys',sys);

% ADC event
adc = mr.makeAdc(Nx,'Duration',round(samplingDuration/Nx/sequencerRasterTime)*Nx*sequencerRasterTime, ...
    'Delay',0,'sys',sys);

% Calculate timing
delayTE1 = ceil((TE/2 - (mr.calcDuration(rf90) - rf90.delay)/2 ...
    - (mr.calcDuration(rf180) - rf180.delay)/2)/seq.gradRasterTime)*seq.gradRasterTime;
delayTE2 = ceil((TE - (mr.calcDuration(rf90)  - rf90.delay)/2 - delayTE1 ...
    - mr.calcDuration(rf180))/seq.gradRasterTime)*seq.gradRasterTime;

% Loop over phase encodes and define sequence blocks
seq.addBlock(rf90);
seq.addBlock(mr.makeDelay(delayTE1));
seq.addBlock(rf180);
seq.addBlock(mr.makeDelay(delayTE2));
seq.addBlock(adc);


%% prepare sequence export
seq.setDefinition('Name', 'se');
seq.setDefinition('TE [s]', TE);
seq.setDefinition('Nx', Nx);
seq.setDefinition('Bandwidth [Hz]', 1/adc.dwell);

seq.plot();

seq.write('tabletop_se_pulseq.seq')       % Write to pulseq file
parsemr('tabletop_se_pulseq.seq');
