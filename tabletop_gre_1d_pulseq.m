addpath('../pulseq/matlab')
close all
clear
gamma = 42.57E6;

fov=10e-3; Nx=128; Ny=1;       % Define FOV and resolution
TE=2.2e-3;

gxFlatTime = 2e-3;

% set system limits
maxGrad = 400; % [mT/m], value for tabletop coils and gpa fhdo
rfDeadTime = 500e-6; % [us], minicircuits PA needs 500 us to turn on
adcDeadTime = 0;
sys = mr.opts('MaxGrad', maxGrad, 'GradUnit', 'mT/m', ...
    'MaxSlew', 1200, 'SlewUnit', 'T/m/s', ...
    'rfDeadTime', rfDeadTime, 'adcDeadTime', adcDeadTime, ...
    'rfRasterTime', 1e-6, 'gradRasterTime',10e-6);
seq=mr.Sequence(sys);              % Create a new sequence object

% Create HF pulses, 500 us delay for tx gate
rf90duration=0.10e-3;
rf90 = mr.makeBlockPulse(pi/2, 'duration', rf90duration,...
    'PhaseOffset', 0, 'sys', sys);

% Define other gradients and ADC events
deltak=1/fov;
gx = mr.makeTrapezoid('x','FlatArea',Nx*deltak,'FlatTime',gxFlatTime,'sys',sys);
fprintf('Sequence bandwidth: %.3f Hz\n',gx.amplitude*1E-3*fov);
fprintf('Pixelbandwidth: %.3f Hz\n',gx.amplitude*1E-3*fov/Nx);
gx.delay = 0; % assumes rfDeadTime > gx.riseTime !!
gxPre = mr.makeTrapezoid('x','Area',-gx.area/2,'Duration',gx.flatTime/2,'sys',sys);
oversamplingFactor = 10;
adc = mr.makeAdc(oversamplingFactor*Nx,'Duration',gx.flatTime,'Delay',gx.riseTime,'sys',sys);

% Calculate timing
delayTE = (TE - (mr.calcDuration(rf90)-rf90.delay)/2 - mr.calcDuration(gxPre) ...
    - mr.calcDuration(gx)/2);
fprintf('delay1: %.3f ms \n',delayTE*1E3)

seq.addBlock(rf90);
seq.addBlock(gxPre);
seq.addBlock(mr.makeDelay(delayTE));
seq.addBlock(gx,adc);
seq.addBlock(mr.makeDelay(delayTE));


%% prepare sequence export
seq.setDefinition('Name', 'se');
seq.setDefinition('FOV', [fov fov]);
seq.setDefinition('TE [s]', TE);
seq.setDefinition('Nx', Nx);
seq.setDefinition('Bandwidth [Hz]', 1/adc.dwell);

seq.plot();

seq.write('tabletop_gre_1d_pulseq.seq')       % Write to pulseq file
parsemr('tabletop_gre_1d_pulseq.seq');
