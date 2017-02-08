%% stepTime.m
% 
% This function advances the time by one step of length dT. It takes a
% model neuron as input and returns an updated model neuron as output.
% 
% - JSB & AEB 2/2013
% - AVB & SLH 2/2016

function neuron = stepTime(neuron, dT)

logistic = @(Vthresh, lambda, V) 1 ./ (1 + exp(- (V - Vthresh) / lambda));

% lambda = 0.00001   CV = 1.4674
% lambda = 0.1   CV = 1.5
% lambda = 1      CV =  1.5852

lambda = 2;
P_spike = logistic(neuron.Vthresh, lambda, neuron.Vm);
do_spike = P_spike > rand;

% Determine if the neuron will spike, or not
%if (neuron.Vm > neuron.Vthresh)
if (do_spike)
    neuron.spike = true;        % Note there's a spike
    % Reset Vm to the reset voltage
    neuron.Vm = neuron.Vreset;
    
    % Update learning rule M
    neuron.exSynapses.M = neuron.exSynapses.M - neuron.exSynapses.Aminus;
    % Update conductances as a result of the learning rule applied
    % to post-synaptic spikes.
    neuron.exSynapses.gA = neuron.exSynapses.gA + neuron.exSynapses.Pa.*neuron.exSynapses.gMax;
    % Don't allow conductances out of the range [0,gMax]
    neuron.exSynapses.gA = neuron.exSynapses.gA - ...
        (neuron.exSynapses.gA > neuron.exSynapses.gMax).*...
        (neuron.exSynapses.gA - neuron.exSynapses.gMax);
    neuron.exSynapses.gA = neuron.exSynapses.gA - ...
        (neuron.exSynapses.gA < 0).*...
        (neuron.exSynapses.gA);
else % If it doesn't spike...
    neuron.spike = false;       % Note there's no spike
    % Update membrane voltage based on conductances
    dV = (dT/neuron.tauM)*(neuron.Vrest - neuron.Vm ...
        + neuron.gEx*(neuron.Eex - neuron.Vm) ...
        + neuron.gIn*(neuron.Ein - neuron.Vm) );
    neuron.Vm = neuron.Vm + dV;
end

% Allow conductances to decay exponentially
dgEx = -neuron.gEx*dT/neuron.tauEx;
dgIn = -neuron.gIn*dT/neuron.tauIn;
neuron.gEx = neuron.gEx + dgEx;
neuron.gIn = neuron.gIn + dgIn;

% Generate Poisson presynaptic spikes, 1 for spike, 0 for none
neuron.exSynapses.preSynapticSpike = (rand(neuron.Nex,1) < dT*neuron.exSynapses.rate);
neuron.inSynapses.preSynapticSpike = (rand(neuron.Nin,1) < dT*neuron.inSynapses.rate);

% Presynaptic spikes generate conductances in the post-synaptic
% cell
exCond = neuron.exSynapses.preSynapticSpike.*neuron.exSynapses.gA;
inCond = neuron.inSynapses.preSynapticSpike.*neuron.inSynapses.gA;
neuron.gEx = neuron.gEx + sum(exCond);
neuron.gIn = neuron.gIn + sum(inCond);

% Update learning rule: Pa increases conductances on
% post-synaptic spiking, and is incremented on pre-synaptic
% spiking.
neuron.exSynapses.Pa = neuron.exSynapses.Pa +...
    neuron.exSynapses.preSynapticSpike.*neuron.exSynapses.Aplus;

% Update the conductances as a result of the learning rule
% applied to pre-synaptic spikes.
neuron.exSynapses.gA = neuron.exSynapses.gA + ...
    neuron.exSynapses.preSynapticSpike.*neuron.exSynapses.M*neuron.exSynapses.gMax;
% Don't allow conductances out of the range [0,gMax]
neuron.exSynapses.gA = neuron.exSynapses.gA - ...
    (neuron.exSynapses.gA > neuron.exSynapses.gMax).*...
    (neuron.exSynapses.gA - neuron.exSynapses.gMax);
neuron.exSynapses.gA = neuron.exSynapses.gA - ...
    (neuron.exSynapses.gA < 0).*...
    (neuron.exSynapses.gA);

% The learning rule functions M and Pa decay exponentially
dM = -neuron.exSynapses.M*dT/neuron.exSynapses.tauMinus;
neuron.exSynapses.M = neuron.exSynapses.M + dM;
dPa = -neuron.exSynapses.Pa.*dT/neuron.exSynapses.tauPlus;
neuron.exSynapses.Pa = neuron.exSynapses.Pa + dPa;

end % End stepTime()