close all
clear
gamma = 42.57E6;
sequencerRasterTime = 7E-9; % make sure all times are a multiple of sequencer raster time
grad_interval = 10E-6;
rf_interval = 1E-6;

fov=10e-3; Nx=200; Nspokes=80;       % Define FOV and resolution
Ndummy=2;                            % number of dummy scans
TE=12e-3; % [s]
TR=5; % [s]     
angle=360;
oversamplingFactor = 4;

gxFlatTime = 3e-3;  % = adc read time
rf90duration=0.1e-3;
rf180duration=2*rf90duration;
sliceThickness = 10;
use_slice = 0;

% set system limits
maxGrad = 500; % [mT/m], value for tabletop coils and gpa fhdo
rfDeadTime = 500e-6; % [us], minicircuits PA needs 500 us to turn on
adcDeadTime = 0;
sys = mr.opts('MaxGrad', maxGrad, 'GradUnit', 'mT/m', ...
    'MaxSlew', 500, 'SlewUnit', 'T/m/s', ...
    'rfDeadTime', rfDeadTime, 'adcDeadTime', adcDeadTime, ...
    'rfRasterTime', rf_interval, 'gradRasterTime', grad_interval);
seq=mr.Sequence(sys);              % Create a new sequence object

% Create HF pulses, 500 us delay for tx gate
if use_slice == 1
    [rf90, gs] = mr.makeBlockPulse(pi/2, 'duration', rf90duration,...
        'PhaseOffset', 0, 'sys', sys, 'SliceThickness', sliceThickness);
    gs.channel='z'; % change it to X because we want sagittal orientation
else
    rf90 = mr.makeBlockPulse(pi/2, 'duration', rf90duration,...
        'PhaseOffset', 0, 'sys', sys); 
    sliceThickness = 0;
end

rf180 = mr.makeBlockPulse(pi, 'duration', rf180duration,...
    'PhaseOffset', pi/2, 'sys', sys);


% Define other gradients and ADC events
deltak=1/fov;
kWidth=deltak*Nx;
gx = mr.makeTrapezoid('x','FlatArea',kWidth,'FlatTime',gxFlatTime,'sys',sys);
gy = mr.makeTrapezoid('y','FlatArea',kWidth,'FlatTime',gxFlatTime,'sys',sys);
fprintf('Sequence bandwidth: %.3f Hz\n',gx.amplitude*1E-3*fov);
fprintf('Pixelbandwidth: %.3f Hz\n',gx.amplitude*1E-3*fov/Nx);
gx.delay = 0; % assumes rfDeadTime > gx.riseTime !!
gxPre = mr.makeTrapezoid('x','Area',gx.area/2,'Duration',gx.flatTime/2,'sys',sys);
gyPre = mr.makeTrapezoid('y','Area',gx.area/2,'Duration',gx.flatTime/2,'sys',sys);
adc = mr.makeAdc(round(oversamplingFactor*Nx),'Duration',gx.flatTime,'Delay',gx.riseTime,'sys',sys);

% Calculate timing
delayTE1 = ceil((TE/2 - (mr.calcDuration(rf90)-rf90.delay)/2 ...
    - mr.calcDuration(gxPre)...
    - rf180.delay - (mr.calcDuration(rf180)-rf180.delay)/2)/seq.gradRasterTime)*seq.gradRasterTime;
delayTE2 = ceil((TE/2 - (mr.calcDuration(rf180) - rf180.delay)/2 ...
    - mr.calcDuration(gx)/2 )/seq.gradRasterTime)*seq.gradRasterTime;
delayTR = TR - TE -rf90.delay -mr.calcDuration(rf90)/2 - mr.calcDuration(gx)/2;
fprintf('delay1: %.3f ms \ndelay2: %.3f ms \n',delayTE1*1E3,delayTE2*1E3)

extra_delay = 1e-3
delta = angle/360*2*pi / Nspokes;            % angular increment;
for i=(1-Ndummy):Nspokes
    if use_slice == 1
        seq.addBlock(rf90, gs);
    else
        seq.addBlock(rf90);
    end
    seq.addBlock(mr.makeDelay(delayTE1 - extra_delay));
    seq.addBlock(mr.rotate('z',delta*(i-1),gyPre)); 
    seq.addBlock(mr.makeDelay(extra_delay));
    seq.addBlock(rf180);
    seq.addBlock(mr.makeDelay(delayTE2));
    if (i>0)
        seq.addBlock(mr.rotate('z',delta*(i-1),gy,adc));
    else
        seq.addBlock(mr.rotate('z',delta*(i-1),gy));   
    end
    seq.addBlock(mr.makeDelay(delayTR));
end


%% prepare sequence export
seq.setDefinition('Name', 'se_2d');
seq.setDefinition('FOV', [fov fov]);
seq.setDefinition('SliceThickness', sliceThickness);
seq.setDefinition('TE [s]', TE);
seq.setDefinition('TR', TR);
seq.setDefinition('Nx', Nx);
seq.setDefinition('Nspokes', Nspokes);
seq.setDefinition('Bandwidth [Hz]', 1/adc.dwell);
seq.setDefinition('grad_interval', grad_interval);
seq.setDefinition('rf_interval', rf_interval);
seq.setDefinition('angle', angle);

seq.plot();

seq.write('tabletop_radial_v2_2d_pulseq.seq')       % Write to pulseq file
parsemr('tabletop_radial_v2_2d_pulseq.seq');
