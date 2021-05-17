addpath('../pulseq/matlab')
close all
clear
gamma = 42.57E6;
sequencerRasterTime = 7E-9; % make sure all times are a multiple of sequencer raster time
tx_t = 1E-6;
grad_t = 10E-6;
Nx = 256;
samplingDuration = 6e-3;

% set system limits
rfDeadTime = 500e-6; % [us], minicircuits PA needs 500 us to turn on
adcDeadTime = 0;
sys = mr.opts('rfDeadTime', round(rfDeadTime/sequencerRasterTime)*sequencerRasterTime, ...
    'adcDeadTime', round(adcDeadTime/sequencerRasterTime)*sequencerRasterTime, ...
    'rfRasterTime', tx_t);
seq=mr.Sequence(sys);              % Create a new sequence object

% Create HF pulses
rf90duration=0.08e-3
rf90 = mr.makeBlockPulse(pi/2, 'duration', rf90duration,...
    'PhaseOffset', 0, 'sys', sys);

% ADC event
adc = mr.makeAdc(Nx,'Duration',round(samplingDuration/Nx/sequencerRasterTime)*Nx*sequencerRasterTime, ...
    'Delay',0,'sys',sys);

% Calculate timing
delayFID = 0.01e-3;

% Loop over phase encodes and define sequence blocks
seq.addBlock(rf90);
%seq.addBlock(mr.makeDelay(delayFID));
seq.addBlock(adc);


%% prepare sequence export
seq.setDefinition('Name', 'fid');
seq.setDefinition('Nx', Nx);
seq.setDefinition('tx_t', tx_t);
seq.setDefinition('grad_t', grad_t);

seq.plot();

seq.write('tabletop_fid_pulseq.seq')       % Write to pulseq file
parsemr('tabletop_fid_pulseq.seq');
