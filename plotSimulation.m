%% plotSimulation.m
% 
% An example script that runs a simulation using a structure from the 
% modelNeuron.m function.  This periodically calls the function stepTime() 
% to advance the time simulation and plots some variables. We find it
% should be run with time steps (dT) of .1 ms or smaller, in order to 
% reproduce the key findings.
% 
% Note that plotting in MATLAB is slow, so we advance time through
% plotInterval steps, storing chunks of time data in a buffer, between each
% plot.
% 
% In our hands, depending on the input, the simulation takes at least 100
% seconds of simulated time to converge on a stable post-synaptic rate and
% distribution of conductances, and sometimes substantially longer.  Tip:
% Increasing the strength of STDP and the presynaptic rates will speed
% convergence.
% 
% In this script we make an example model neuron (which we will call
% aNeuron) by calling the modelNeuron function:
%   aNeuron = modelNeuron();
% 
% aNeuron currently has the default properties set in the modelNeuron
% script. We can change these default properties in the following way:
% 
%     aNeuron.exSynapses.rate = 15;
% 
% This has changed the excitatory presynaptic rate of action potentials
% from 25 to 15.
% 
% We then advance simulation time by calling the stepTime(neuron, dT)
% function.  This function takes the example neuron we've generate and the
% desired time step (in sec.) as input arguments and returns an updated
% version of our example neuron as output:
% 
%     aNeuron = stepTime(aNeuron,stepSize);
% 
% - JSB & AEB 2/2013 - AVB & SLH 2/2016

%% Note this will *clear all* variables from your workspace
clear all
close all

%% Simulation parameters
simLength    =  999;         % Time to simulate, (sec)
stepSize     = .0001;        % Time step resolution, (sec)
plotInterval =  5000;        % Period of plot updates, (time steps)
plotWindow   =     1;        % Window to plot Vm over (sec)

%% Housekeeping variables

nTimePoints = round(simLength/stepSize);  % # of points in simulation

% Figures for display
traceFig = figure();         % Displays trace of Vm and conductances
histFig  = figure();         % Displays histogram of cond. strength
psthFig  = figure();         % Displays firing rate over time
blahFig = figure();

% Stores instantaneous firing rates, averaged over the plotWindow, over
% time.  A new point is added every plotInterval.
psth    = [];   % # of spikes per plotInterval, over time.
psthN   = 1;    % # of time points in the psth

% Buffers for storing data between plot updates
tChunk   = zeros(plotInterval+1,1); % Time (sec)
VmChunk  = zeros(plotInterval+1,1); % Vm (mV)
gExChunk = zeros(plotInterval+1,1); % Excitatory conductance (1/Rin)
gInChunk = zeros(plotInterval+1,1); % Inhibitory conductance (1/Rin)
rasChunk = zeros(plotInterval+1,1); % Post-synaptic spike raster

%% Make a modelNeuron structure
aNeuron = modelNeuron();
% Change some of the default parameters of the simulation
aNeuron.exSynapses.rate = 15;
aNeuron.inSynapses.rate = 10;
aNeuron.exSynapses.Aplus = .020;
aNeuron.exSynapses.Aminus = 1.05*aNeuron.exSynapses.Aplus;

ras = [];
preras = [];
prespikes = zeros(round(0.100 / stepSize), 1);

%% Simulation loop over points in simulated time
%  n is the step # in the simulation
for n=1:nTimePoints
    
    aNeuron = stepTime(aNeuron,stepSize);   % Advance the sim time by 1 step
    plotIdx = mod(n,plotInterval)+2;        % Index into the buffers, ranges from [2,plotInterval+1]
    
    tChunk(plotIdx)   = n*stepSize;    % Store the time in buffer
    VmChunk(plotIdx)  = aNeuron.Vm;    % Store the Vm in buffer
    gExChunk(plotIdx) = aNeuron.gEx;   % Store excitatory conductance
    gInChunk(plotIdx) = aNeuron.gIn;   % Store inhibitory conductance
    rasChunk(plotIdx) = aNeuron.spike; % Store the raster
   % ras = [ras; aNeuron.spike];
   % preras = [preras; sum(aNeuron.exSynapses.preSynapticSpike)];
    
    % Every 10 plotWindows, clear the traceFig to prevent unseen Vm
    % trace from eating up memory.
    if mod(n,round(plotWindow*10/stepSize)) == 0
        figure(traceFig); clf; hold off;
    end

   % if n >= 10 / stepSize && ras(end - round(length(prespikes) / 2))
    %    prespikes = prespikes + preras(end - length(prespikes) + 1 : end);
   % end
    
    % If the buffers are full, plot the traces to figures
    if (plotIdx == plotInterval + 1)
        
       % spikes = find(ras);
       % intervals = spikes(2:end) - spikes(1:end-1);
       % cv = std(intervals) / mean(intervals);
       % disp(cv);
        
       % figure(blahFig);
        %plot(prespikes');        
        
        % Plot traces of Vm
        figure(traceFig); subplot(3,1,1:2);
        plot(tChunk,VmChunk,'b'); hold on;
        xlim([tChunk(end)-plotWindow, tChunk(end)+plotWindow*.05]);
        ylim([aNeuron.Ein aNeuron.Vthresh*.9]); ylabel('Vm (mV)');
        
        % Plot traces of gEx (blue) and gIn (red)
        figure(traceFig); subplot(3,1,3);
        plot(tChunk,gExChunk,'b'); hold on;
        plot(tChunk,gInChunk,'r');
        xlim([tChunk(end)-plotWindow, tChunk(end)+plotWindow*.05]);
        xlabel('Time (s)'); ylabel('gEx (blue) gIn (red)');
        
        % Plot histograms
        gAs = [aNeuron.exSynapses.gA];
        gMax = aNeuron.exSynapses.gMax;
        nBins = 20;
        figure(histFig);
        hist(gAs/gMax,[1/(2*nBins):1/nBins:(1-1/(2*nBins))]);
        axis tight; xlabel('gA/gMax'); ylabel('N synapses');
        
        % Plot PSTH
        figure(psthFig);
        psth(psthN) = nnz(rasChunk);
        plot(([1:psthN]-1)*plotInterval*stepSize,psth./(plotInterval*stepSize));
        psthN = psthN + 1;
        xlabel('Time (s)'); ylabel('Post-synaptic Rate (Hz)');
        
        % Clear plot buffers, storing the last element of the
        % current buffer into the first element of the next buffer.
        tChunk(1)   = tChunk(plotInterval+1);
        VmChunk(1)  = VmChunk(plotInterval+1);
        gExChunk(1) = gExChunk(plotInterval+1);
        gInChunk(1) = gInChunk(plotInterval+1);
        rasChunk(1) = rasChunk(plotInterval+1);
        
        pause(.01); % Pause to allow plots to update to screen
    end % End if (buffers are full)
end % End for all points in time
