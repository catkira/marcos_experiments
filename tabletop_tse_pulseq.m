close all
clear
gamma = 42.57E6;
sequencerRasterTime = 1/(122.88E6); % make sure all times are a multiple of sequencer raster time
grad_interval = 10E-6;
rf_interval = 1E-6;

fov=10e-3; Nx=200; Ny=128;   % Define FOV and resolution
TR=5; % [s]     
ETL=8;
ESP=8e-3;
Ndummy = 1; % excitations
oversampling_factor = 4;
sliceThickness = 5;
use_slice = 1;
sp_amplitude = 0; % spoiler area in 1/m (=Hz/m*s)
sp_duration = 0.5E-3;

assert(mod(Ny,ETL) == 0)

gxFlatTime = 3e-3;  % = adc read time [s]
rf90duration = 0.1e-3;
dwellTime = gxFlatTime/(Nx*oversampling_factor);
rf180duration=2*rf90duration;

% set system limits
maxGrad = 400; % [mT/m], value for tabletop coils and gpa fhdo
rfDeadTime = 200e-6; % [us], minicircuits PA needs 500 us to turn on, Motorola amp needs 200 us
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
end
rf180 = mr.makeBlockPulse(pi, 'duration', rf180duration,...
    'PhaseOffset', pi/2, 'use', 'refocusing', 'sys',sys);


%% Define other gradients and ADC events
deltak=1/fov;
kWidth=deltak*Nx;
kHeight=deltak*Ny;
gx = mr.makeTrapezoid('x','FlatArea',kWidth,'FlatTime',gxFlatTime,'sys',sys);
fprintf('Sequence bandwidth: %.3f Hz\n',gx.amplitude*1E-3*fov);
fprintf('Pixelbandwidth: %.3f Hz\n',gx.amplitude*1E-3*fov/Nx);
gx.delay = 0; % assumes rfDeadTime > gx.riseTime !!
gxPre = mr.makeTrapezoid('x','Area',gx.area/2,'Duration',mr.calcDuration(gx),'sys',sys);
g_sp = mr.makeTrapezoid('x','Area',sp_amplitude,'Duration',sp_duration,'system',sys);
if sp_amplitude == 0
    g_sp.riseTime = 0;
    g_sp.fallTime = 0;
    g_sp.flatTime = 0;
end
gy_area = kHeight/2;
gy = mr.makeTrapezoid('y','Area',gy_area,'Duration',gx.flatTime/2,'sys',sys);
adc_duration = dwellTime * Nx * oversampling_factor;
adc_delay = (mr.calcDuration(gx) - adc_duration)/2 - 1/2*grad_interval; % - 1/2*grad_interval is needed because first sample is already in the ramp
adc = mr.makeAdc(round(oversampling_factor*Nx),'Duration',adc_duration,'Delay',adc_delay,'sys',sys);

%% Calculate delays
delayTE = ESP/2 - mr.calcRfCenter(rf90) ...
    - mr.calcDuration(gxPre) - mr.calcDuration(g_sp) ...
    - rf180.delay - mr.calcRfCenter(rf180);
delayTE1 = ESP/2 -  mr.calcDuration(gx)/2 - gx.flatTime/2 ...
    - mr.calcDuration(g_sp) - rf180.delay - mr.calcRfCenter(rf180);
delayTE2 = ESP/2 - mr.calcRfCenter(rf180) ...
    - mr.calcDuration(gx)/2  -  mr.calcDuration(g_sp) -  mr.calcDuration(gy);
delayTR = TR - ETL*ESP - mr.calcDuration(rf90) - mr.calcDuration(gx)/2;
fprintf('delay1: %.3f ms \ndelay2: %.3f ms \n',delayTE1*1E3,delayTE2*1E3)

%% assemble sequence
nex = Ny / ETL;
pe_steps=(1:(ETL*nex))-0.5*ETL*nex-1;
pe_steps_interleaved = reshape(pe_steps, [nex, ETL]); 
phase_areas = pe_steps_interleaved * deltak;
for n = 1-Ndummy:nex
    if use_slice == 1
        seq.addBlock(rf90, gs);
    else
        seq.addBlock(rf90);
    end
    seq.addBlock(mr.makeDelay(delayTE));
    seq.addBlock(gxPre);
    for m=1:ETL  
        if n > 0
            gy = mr.makeTrapezoid('y','Area',phase_areas(n,m),'Duration',gx.flatTime/2,'sys',sys);    
            gy_rev = mr.makeTrapezoid('y','Area',-phase_areas(n,m),'Duration',gx.flatTime/2,'sys',sys);    
        else
            gy = mr.makeTrapezoid('y','Area',phase_areas(1,m),'Duration',gx.flatTime/2,'sys',sys);    
            gy_rev = mr.makeTrapezoid('y','Area',-phase_areas(1,m),'Duration',gx.flatTime/2,'sys',sys);    
        end
        if m ~= 1
            seq.addBlock(mr.makeDelay(delayTE1));
        end
        if sp_amplitude > 0
            seq.addBlock(g_sp);    
            seq.addBlock(rf180);
            seq.addBlock(g_sp);    
        else
            seq.addBlock(rf180);
        end
        seq.addBlock(mr.makeDelay(delayTE2));
        seq.addBlock(gy);      
        if n > 0
            seq.addBlock(gx,adc);
        else
            seq.addBlock(gx,adc);
        end
        seq.addBlock(gy_rev);
    end
    if n < nex
        seq.addBlock(mr.makeDelay(delayTR));
    end    
end

%% some checks
% time of first echo should be at ESP + rf90.delay + mr.calcRfCenter(rf90)
t_1st_echo = mr.calcDuration(rf90) + delayTE + mr.calcDuration(gxPre) + mr.calcDuration(g_sp) ...
    + mr.calcDuration(rf180) +  mr.calcDuration(g_sp) + delayTE2 + mr.calcDuration(gy) ...
    + mr.calcDuration(gx)/2;
assert(abs(t_1st_echo - (ESP + rf90.delay + mr.calcRfCenter(rf90))) < sequencerRasterTime) 
% center of 2nd 180deg pulse should be at t_first_echo + ESP/2
t_2nd_rf180 = t_1st_echo + mr.calcDuration(gx)/2 + mr.calcDuration(gy_rev) + delayTE1 ...
    + mr.calcDuration(g_sp) + rf180.delay + mr.calcRfCenter(rf180);
assert(abs(t_2nd_rf180 - (t_1st_echo + ESP/2)) < sequencerRasterTime)
% center of 2nd echo should be at t_1st_echo + ESP
t_2nd_echo = t_2nd_rf180 + mr.calcRfCenter(rf180) + mr.calcDuration(g_sp) ...
    + delayTE2 + mr.calcDuration(gy) + mr.calcDuration(gx)/2;
assert(abs(t_2nd_echo - (t_1st_echo + ESP)) < sequencerRasterTime)

%% prepare sequence export
seq.setDefinition('Name', 'tse_2d');
seq.setDefinition('FOV', [fov fov]);
seq.setDefinition('ESP [s]', ESP);
seq.setDefinition('ETL', ETL);
seq.setDefinition('TR', TR);
seq.setDefinition('Nx', Nx);
seq.setDefinition('Ny', Ny);
seq.setDefinition('Ndummy', Ndummy);
seq.setDefinition('Bandwidth [Hz]', 1/adc.dwell);
seq.setDefinition('grad_t', grad_interval*1E6);
seq.setDefinition('tx_t', rf_interval*1E6);
seq.setDefinition('SliceThickness', sliceThickness);

seq.plot();

seq.write('tabletop_tse_pulseq.seq')       % Write to pulseq file
%parsemr('tabletop_tse_pulseq.seq');

%% plot trajectory
%[ktraj_adc, ktraj, t_excitation, t_refocusing] = seq.calculateKspace();
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();

% plot k-spaces

%figure; plot(ktraj'); % plot the entire k-space trajectory
figure; plot(ktraj(1,:),ktraj(2,:),'b',...
             ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % a 2D plot
axis('equal'); % enforce aspect ratio for the correct trajectory display
