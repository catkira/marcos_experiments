addpath('../pulseq/matlab')
close all
clear
gamma = 42.57E6;

fov=10e-3; Nx=128; Ny=1;       % Define FOV and resolution
TE=12e-3;

gxFlatTime = 4e-3;

% set system limits
maxGrad = 125; % [mT/m], value for tabletop coils and gpa fhdo
spA=1000; % spoiler area in 1/m (=Hz/m*s)
rfDeadTime = 500e-6; % [us], minicircuits PA needs 500 us to turn on
adcDeadTime = 0;
sys = mr.opts('MaxGrad', maxGrad, 'GradUnit', 'mT/m', ...
    'MaxSlew', 500, 'SlewUnit', 'T/m/s', ...
    'rfDeadTime', rfDeadTime, 'adcDeadTime', adcDeadTime, ...
    'rfRasterTime', 1e-6, 'gradRasterTime',10e-6);
seq=mr.Sequence(sys);              % Create a new sequence object

% Create HF pulses, 500 us delay for tx gate
rf90duration=0.10e-3;
rf90 = mr.makeBlockPulse(pi/2, 'duration', rf90duration,...
    'PhaseOffset', 0, 'sys', sys);
rf180 = mr.makeBlockPulse(pi, 'duration', rf90duration*2,...
    'PhaseOffset', pi/2, 'sys',sys);

% Define other gradients and ADC events
deltak=1/fov;
gx = mr.makeTrapezoid('x','FlatArea',Nx*deltak,'FlatTime',gxFlatTime,'sys',sys);
fprintf('Sequence bandwidth: %.3f Hz\n',gx.amplitude*1E-3*fov);
fprintf('Pixelbandwidth: %.3f Hz\n',gx.amplitude*1E-3*fov/Nx);
gx.delay = 0; % assumes rfDeadTime > gx.riseTime !!
gxPre = mr.makeTrapezoid('x','Area',gx.area/2,'Duration',gx.flatTime/2,'sys',sys);
g_sp = mr.makeTrapezoid('x','Area',spA,'Duration',0.5e-3,'system',sys);
oversamplingFactor = 1;
adc = mr.makeAdc(oversamplingFactor*Nx,'Duration',gx.flatTime,'Delay',gx.riseTime,'sys',sys);

% Calculate timing
delayTE1_2 = 1e-3;
delayTE1 = ceil((TE/2 - (mr.calcDuration(rf90)-rf90.delay)/2 ...
    - mr.calcDuration(gxPre) -  mr.calcDuration(g_sp)...
    - rf180.delay - (mr.calcDuration(rf180)-rf180.delay)/2 - delayTE1_2)/seq.gradRasterTime)*seq.gradRasterTime;
delayTE2 = ceil((TE/2 - (mr.calcDuration(rf180) - rf180.delay)/2 ...
    - mr.calcDuration(gx)/2  -  mr.calcDuration(g_sp))/seq.gradRasterTime)*seq.gradRasterTime;
fprintf('delay1: %.3f ms \ndelay2: %.3f ms \n',delayTE1*1E3,delayTE2*1E3)

% Loop over phase encodes and define sequence blocks
for i=1:Ny
    for c=1:length(TE)    
        seq.addBlock(rf90);
        seq.addBlock(mr.makeDelay(delayTE1(c)));
        seq.addBlock(gxPre);
        seq.addBlock(mr.makeDelay(delayTE1_2(c)));
        seq.addBlock(g_sp);
        seq.addBlock(rf180);
        seq.addBlock(g_sp);
        seq.addBlock(mr.makeDelay(delayTE2(c)));
        seq.addBlock(gx,adc);
    end
end


%% prepare sequence export
seq.setDefinition('Name', 'se');
seq.setDefinition('FOV', [fov fov]);
seq.setDefinition('TE [s]', TE);
seq.setDefinition('Nx', Nx);
seq.setDefinition('Bandwidth [Hz]', 1/adc.dwell);

seq.plot();

seq.write('tabletop_se_1d_pulseq.seq')       % Write to pulseq file
parsemr('tabletop_se_1d_pulseq.seq');
