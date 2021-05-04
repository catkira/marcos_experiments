addpath('../pulseq/matlab')
close all
clear
gamma = 42.57E6;
sequencerRasterTime = 7E-9; % make sure all times are a multiple of sequencer raster time
grad_interval = 10E-6;
rf_interval = 1E-6;

fov=12e-3; Nx=200; Ny=63;   % Define FOV and resolution
TE=12e-3; % [s]
TR=5; % [s]     
readoutOversamplingFactor = 4;
sliceThickness = 10;
use_slice = 0;

gxFlatTime = 4e-3;
spA=1000; % spoiler area in 1/m (=Hz/m*s)

% set system limits
maxGrad = 125; % [mT/m], value for tabletop coils and gpa fhdo
rfDeadTime = 500e-6; % [us], minicircuits PA needs 500 us to turn on
adcDeadTime = 0;
sys = mr.opts('MaxGrad', maxGrad, 'GradUnit', 'mT/m', ...
    'MaxSlew', 500, 'SlewUnit', 'T/m/s', ...
    'rfDeadTime', rfDeadTime, 'adcDeadTime', adcDeadTime, ...
    'rfRasterTime', rf_interval, 'gradRasterTime', grad_interval);
seq=mr.Sequence(sys);              % Create a new sequence object

% Create HF pulses, 500 us delay for tx gate
rf90duration=0.1e-3;
if use_slice == 1
    [rf90, gs] = mr.makeBlockPulse(pi/2, 'duration', rf90duration,...
        'PhaseOffset', 0, 'sys', sys, 'SliceThickness', sliceThickness);
    gs.channel='z'; % change it to X because we want sagittal orientation
else
    rf90 = mr.makeBlockPulse(pi/2, 'duration', rf90duration,...
        'PhaseOffset', 0, 'sys', sys);    
end
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
g_sp = mr.makeTrapezoid('x','Area',spA,'Duration',0.5e-3,'system',sys);
gy_area = kHeight/2;
gy = mr.makeTrapezoid('y','Area',gy_area,'Duration',gx.flatTime/2,'sys',sys);
Ny = round(Ny);
adc = mr.makeAdc(round(readoutOversamplingFactor*Nx),'Duration',gx.flatTime,'Delay',gx.riseTime,'sys',sys);

% Calculate timing
delayTE1_2 = 1e-3;
delayTE1 = ceil((TE/2 - (mr.calcDuration(rf90)-rf90.delay)/2 ...
    - mr.calcDuration(gxPre) -  mr.calcDuration(g_sp)...
    - rf180.delay - (mr.calcDuration(rf180)-rf180.delay)/2 - delayTE1_2)/seq.gradRasterTime)*seq.gradRasterTime;
delayTE2 = ceil((TE/2 - (mr.calcDuration(rf180) - rf180.delay)/2 ...
    - mr.calcDuration(gx)/2  -  mr.calcDuration(g_sp))/seq.gradRasterTime)*seq.gradRasterTime;
delayTR = TR - TE -rf90.delay -mr.calcDuration(rf90)/2 - mr.calcDuration(gx)/2;
fprintf('delay1: %.3f ms \ndelay2: %.3f ms \n',delayTE1*1E3,delayTE2*1E3)

phase_factor = linspace(-1,1,Ny);
for n=1:Ny
    if use_slice == 1
        seq.addBlock(rf90, gs);
    else
        seq.addBlock(rf90);
    end
    seq.addBlock(mr.makeDelay(delayTE1));
    gy = mr.makeTrapezoid('y','Area',gy_area*phase_factor(n),'Duration',gx.flatTime/2,'sys',sys);    
    seq.addBlock(gxPre,gy);
    seq.addBlock(mr.makeDelay(delayTE1_2));
    seq.addBlock(g_sp);    
    seq.addBlock(rf180);
    seq.addBlock(g_sp);    
    seq.addBlock(mr.makeDelay(delayTE2));
    seq.addBlock(gx,adc);
    seq.addBlock(mr.makeDelay(delayTR));
end


%% prepare sequence export
seq.setDefinition('Name', 'se_2d');
seq.setDefinition('FOV', [fov fov]);
seq.setDefinition('TE [s]', TE);
seq.setDefinition('TR', TR);
seq.setDefinition('Nx', Nx);
seq.setDefinition('Ny', Ny);
seq.setDefinition('Bandwidth [Hz]', 1/adc.dwell);
seq.setDefinition('grad_interval]', grad_interval);
seq.setDefinition('rf_interval]', rf_interval);

seq.plot();

seq.write('tabletop_se_v2_2d_pulseq.seq')       % Write to pulseq file
parsemr('tabletop_se_v2_2d_pulseq.seq');
