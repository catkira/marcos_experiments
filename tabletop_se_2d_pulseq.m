addpath('../pulseq/matlab')
close all
clear
gamma = 42.57E6;
sequencerRasterTime = 7E-9; % make sure all times are a multiple of sequencer raster time

fov=10e-3; Nx=128; Ny=51;       % Define FOV and resolution
TE=12e-3;
TR=5; % not used               

gxFlatTime = 4e-3;

% set system limits
maxGrad = 125; % [mT/m], value for tabletop coils and gpa fhdo
rfDeadTime = 500e-6; % [us], minicircuits PA needs 500 us to turn on
adcDeadTime = 0;
sys = mr.opts('MaxGrad', maxGrad, 'GradUnit', 'mT/m', ...
    'MaxSlew', 300, 'SlewUnit', 'T/m/s', ...
    'rfDeadTime', rfDeadTime, 'adcDeadTime', adcDeadTime, ...
    'rfRasterTime', 1.003e-6, 'gradRasterTime',10.003e-6); % TODO: try shorter gradRasterTime
seq=mr.Sequence(sys);              % Create a new sequence object

% Create HF pulses, 500 us delay for tx gate
rf90duration=0.08e-3;
rf90 = mr.makeBlockPulse(pi/2, 'duration', rf90duration,...
    'PhaseOffset', 0, 'sys', sys);
rf180 = mr.makeBlockPulse(pi, 'duration', rf90duration*2,...
    'PhaseOffset', pi/2, 'sys',sys);


% Define other gradients and ADC events
deltak=1/fov;
kWidth=deltak*Nx;
kHeight=deltak*Ny;
gx = mr.makeTrapezoid('x','FlatArea',kWidth,'FlatTime',gxFlatTime,'sys',sys);
fprintf('Sequence bandwidth: %.3f Hz\n',gx.amplitude*1E-3*fov);
fprintf('Pixelbandwidth: %.3f Hz\n',gx.amplitude*1E-3*fov/Nx);
gx.delay = 0; % assumes rfDeadTime > gx.riseTime !!
gxPre = mr.makeTrapezoid('x','Area',gx.area/2,'Duration',gx.flatTime/2,'sys',sys);
%gy_area = 1.3e+03; % TODO: calculate this value
gy_area = kHeight/2; % TODO: this value doesnt work, why?
gy = mr.makeTrapezoid('y','Area',gy_area,'Duration',gx.flatTime/2,'sys',sys);
oversamplingFactor = 1.3;
adc = mr.makeAdc(round(oversamplingFactor*Nx),'Duration',gx.flatTime,'Delay',gx.riseTime,'sys',sys);

% Calculate timing
delayTE1 = ceil((TE/2 - (mr.calcDuration(rf90)-rf90.delay)/2 ...
    - mr.calcDuration(gxPre)...
    - rf180.delay - (mr.calcDuration(rf180)-rf180.delay)/2)/seq.gradRasterTime)*seq.gradRasterTime;
delayTE2 = ceil((TE/2 - (mr.calcDuration(rf180) - rf180.delay)/2 ...
    - mr.calcDuration(gx)/2 )/seq.gradRasterTime)*seq.gradRasterTime;
fprintf('delay1: %.3f ms \ndelay2: %.3f ms \n',delayTE1*1E3,delayTE2*1E3)

% Phase looping is done in python script
seq.addBlock(rf90);
seq.addBlock(mr.makeDelay(delayTE1));
seq.addBlock(gxPre,gy);
seq.addBlock(rf180);
seq.addBlock(mr.makeDelay(delayTE2));
seq.addBlock(gx,adc);


%% prepare sequence export
seq.setDefinition('Name', 'se_2d');
seq.setDefinition('FOV', [fov fov]);
seq.setDefinition('TE [s]', TE);
seq.setDefinition('TR', TR);
seq.setDefinition('Nx', Nx);
seq.setDefinition('Ny', Ny);
seq.setDefinition('Bandwidth [Hz]', 1/adc.dwell);

seq.plot();

seq.write('tabletop_se_2d_pulseq.seq')       % Write to pulseq file
parsemr('tabletop_se_2d_pulseq.seq');
