close all
clear
%gamma = 42.57E6;
sequencerRasterTime = 1/(122.88E6); % make sure all times are a multiple of sequencer raster time
grad_interval = ceil(10E-6/sequencerRasterTime)*sequencerRasterTime;
rf_interval = ceil(1E-6/sequencerRasterTime)*sequencerRasterTime;

fov=10e-3; Nx=200; Nspokes=80;       % Define FOV and resolution
Ndummy=0;                            % number of dummy scans
TE=ceil(12e-3/sequencerRasterTime)*sequencerRasterTime; % [s]
TR=5; % [s]     
angle=360;
oversampling_factor = 4;
sliceThickness = 10;
use_slice = 0;
sp_amplitude = 2000; % spoiler area in 1/m (=Hz/m*s)
sp_duration = ceil(0.5E-3/sequencerRasterTime)*sequencerRasterTime;


requested_gxFlatTime = 3e-3;  % = adc read time [s]
requested_rf90duration = 0.1e-3;
% put the dwell time and rf90duration on a 2*sequencerRasterTime raster, so that
% mr.calcDuration(gx)/2 and rf90duration/2 are still on the sequencer raster 
dwellTime = ceil(requested_gxFlatTime/(Nx*oversampling_factor)/(2*sequencerRasterTime))*(2*sequencerRasterTime);
gxFlatTime = dwellTime * Nx * oversampling_factor;
rf90duration=ceil(requested_rf90duration/(2*sequencerRasterTime))*(2*sequencerRasterTime);
rf180duration=2*rf90duration;

% set system limits
maxGrad = 400; % [mT/m], value for tabletop coils and gpa fhdo
rfDeadTime = ceil(500e-6/sequencerRasterTime)*sequencerRasterTime; % [us], minicircuits PA needs 500 us to turn on
adcDeadTime = 0;
sys = mr.opts('MaxGrad', maxGrad, 'GradUnit', 'mT/m', ...
    'MaxSlew', 800, 'SlewUnit', 'T/m/s', ...
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
    'PhaseOffset', pi/2, 'use','refocusing', 'sys', sys);


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
g_sp = mr.makeTrapezoid('x','Area',sp_amplitude,'Duration',sp_duration,'system',sys);
adc = mr.makeAdc(round(oversampling_factor*Nx),'Duration',gx.flatTime,'Delay',gx.riseTime,'sys',sys);

% Calculate timing
delayTE1 = round((TE/2 - (mr.calcDuration(rf90)-rf90.delay)/2 ...
    - mr.calcDuration(gxPre) -  mr.calcDuration(g_sp)...
    - rf180.delay - (mr.calcDuration(rf180)-rf180.delay)/2)/sequencerRasterTime)*sequencerRasterTime;
delayTE2 = round((TE/2 - (mr.calcDuration(rf180) - rf180.delay)/2 ...
    - mr.calcDuration(gx)/2 -  mr.calcDuration(g_sp))/sequencerRasterTime)*sequencerRasterTime;
delayTR = round((TR - TE -rf90.delay -(mr.calcDuration(rf90) - rf90.delay)/2 - mr.calcDuration(gx)/2)/sequencerRasterTime)*sequencerRasterTime;
fprintf('delay1: %.3f ms \ndelay2: %.3f ms \n',delayTE1*1E3,delayTE2*1E3)

extra_delay = ceil(1e-3/sequencerRasterTime)*sequencerRasterTime;
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
    seq.addBlock(g_sp);    
    seq.addBlock(rf180);
    seq.addBlock(g_sp);
    seq.addBlock(mr.makeDelay(delayTE2));
    if (i>0)
        seq.addBlock(mr.rotate('z',delta*(i-1),gy,adc));
    else
        seq.addBlock(mr.rotate('z',delta*(i-1),gy));   
    end
    seq.addBlock(mr.makeDelay(delayTR));
end

% some checks
calculated_t_ref = (mr.calcDuration(rf90) - rf90.delay)/2 + mr.calcDuration(g_sp) ...
    + delayTE1 + mr.calcDuration(gxPre) + rf180.delay + (mr.calcDuration(rf180) - rf180.delay)/2;
assert(abs(calculated_t_ref - TE/2) < sequencerRasterTime)
calulatedTE = (mr.calcDuration(rf90) - rf90.delay)/2 + 2*mr.calcDuration(g_sp)...
    + delayTE1 + mr.calcDuration(gxPre) + mr.calcDuration(rf180) + delayTE2 ...
    + mr.calcDuration(gx)/2;
assert(abs(calulatedTE - TE) < sequencerRasterTime)
calulatedTE2 = (mr.calcDuration(rf90) - rf90.delay)/2 + 2*mr.calcDuration(g_sp) ...
    + delayTE1 + mr.calcDuration(gxPre) + mr.calcDuration(rf180) + delayTE2 ...
    + adc.delay + adc.dwell*adc.numSamples/2;
assert(abs(calulatedTE2 - TE) < sequencerRasterTime)

%% prepare sequence export
seq.setDefinition('Name', 'se_2d');
seq.setDefinition('FOV', [fov fov]);
seq.setDefinition('SliceThickness', sliceThickness);
seq.setDefinition('TE [s]', TE);
seq.setDefinition('TR', TR);
seq.setDefinition('Nx', Nx);
seq.setDefinition('Nspokes', Nspokes);
seq.setDefinition('Bandwidth [Hz]', 1/adc.dwell);
seq.setDefinition('grad_t', grad_interval*1E6);
seq.setDefinition('tx_t', rf_interval*1E6);
seq.setDefinition('angle', angle);

seq.plot();
seq.write('tabletop_radial_v2_2d_pulseq.seq')       % Write to pulseq file


% some more checks
parsemr('tabletop_radial_v2_2d_pulseq.seq');

% calculate k-space but only use it to check timing
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();
if Ndummy==0
    assert(abs(t_refocusing(1)-t_excitation(1)-TE/2)<1e-6); % check that the refocusing happens at the 1/2 of TE
    assert(abs(t_adc(Nx*oversampling_factor/2)-t_excitation(1)-TE)<adc.dwell); % check that the echo happens as close as possible to the middle of the ADC elent
end
