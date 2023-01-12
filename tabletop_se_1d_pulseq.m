addpath('../pulseq/matlab')
close all
clear
gamma = 42.57E6;

fov=8e-3; Nx=128; Ny=1;       % Define FOV and resolution
TE=12e-3;                     % echo time

gradFlatTime = 3e-3;

% set system limits
maxGrad = 400;       % [mT/m], value for tabletop coils and gpa fhdo
rfDeadTime = 500e-6; % [us], minicircuits PA needs 500 us to turn on
adcDeadTime = 0;
sys = mr.opts('MaxGrad', maxGrad, 'GradUnit', 'mT/m', ...
    'MaxSlew', 800, 'SlewUnit', 'T/m/s', ...
    'rfDeadTime', rfDeadTime, 'adcDeadTime', adcDeadTime);
seq=mr.Sequence(sys);              % Create a new sequence object

% Create HF pulses
rf90duration=0.10e-3;
rf90 = mr.makeBlockPulse(pi/2, 'duration', rf90duration,...
    'PhaseOffset', 0, 'sys', sys);
rf180 = mr.makeBlockPulse(pi, 'duration', rf90duration*2,...
    'PhaseOffset', pi/2, 'sys',sys);

% Define other gradients and ADC events
deltak=1/fov;
grad = mr.makeTrapezoid('y','FlatArea',Nx*deltak,'FlatTime',gradFlatTime,'sys',sys);
fprintf('Sequence bandwidth: %.3f Hz\n',grad.amplitude*1E-3*fov);
fprintf('Pixelbandwidth: %.3f Hz\n',grad.amplitude*1E-3*fov/Nx);
grad.delay = 0; % assumes rfDeadTime > gx.riseTime !!
gradPre = mr.makeTrapezoid('y','Area',grad.area/2,'Duration',grad.flatTime/2,'sys',sys);
oversamplingFactor = 2;
adc = mr.makeAdc(oversamplingFactor*Nx,'Duration',grad.flatTime,'Delay',grad.riseTime,'sys',sys);

% Calculate timing
delayTE1_2 = 1e-3;
delayTE1 = ceil((TE/2 - (mr.calcDuration(rf90)-rf90.delay)/2 ...
    - mr.calcDuration(gradPre) - rf180.delay ...
    - (mr.calcDuration(rf180)-rf180.delay)/2 ...
    - delayTE1_2)/seq.gradRasterTime)*seq.gradRasterTime;
delayTE2 = ceil((TE/2 - (mr.calcDuration(rf180) - rf180.delay)/2 ...
    - mr.calcDuration(grad)/2)/seq.gradRasterTime)*seq.gradRasterTime;
fprintf('delay1: %.3f ms \ndelay2: %.3f ms \n',delayTE1*1E3,delayTE2*1E3)

seq.addBlock(rf90);
seq.addBlock(mr.makeDelay(delayTE1));
seq.addBlock(gradPre);
seq.addBlock(mr.makeDelay(delayTE1_2));
seq.addBlock(rf180);
seq.addBlock(mr.makeDelay(delayTE2));
seq.addBlock(grad,adc);


%% prepare sequence export
seq.setDefinition('Name', 'se');
seq.setDefinition('FOV', [fov fov]);
seq.setDefinition('TE [s]', TE);
seq.setDefinition('Nx', Nx);
seq.setDefinition('Bandwidth [Hz]', 1/adc.dwell);
seq.setDefinition('grad_t', 10);  % gradient raster time in us

seq.plot();

seq.write('tabletop_se_1d_pulseq.seq')       % Write to pulseq file
parsemr('tabletop_se_1d_pulseq.seq');
