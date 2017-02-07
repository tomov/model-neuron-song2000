%% modelNeuron.m
%
% This function defines a model neuron with STDP, as described in
% Song, et. al., 2000.  
% 
% The model neuron is a structure which we call 'neuron' inside this
% function. Note that when you call this function your structure will be
% called whatever you set the output of the function as.  For example, if
% you call:
%   
%         aNeuron = modelNeuron() 
% 
% your structure will be called aNeuron. 
% 
%
% This structure has many fields.  In this function we set the default
% values for these fields in this way: 
% 
%         neuron.Vm = -58;    % Membrane voltage, mV
%
% This sets the 'default' membrane voltage to -58mV. 
% 
% Note: you should not change the default values for a field in this
% document. Instead you should first create a default model neuron by
% calling aNeuron = modelNeuron(). Then in a separate script you should
% change the values of your example neuron with code that looks like this: 
% 
%         aNeuron.Vm = -65mV; 
% 
% This has changed the membrane voltage from the default of -58mV to -65mV.
% The script plotSimulation has an example of us doing this.
% 
% 
% - JSB & AEB 2/2013
% - AVB & SLH 2/2016

function neuron = modelNeuron()

% Variables - these change over the simulation
neuron.Vm      =   -58;    % Membrane voltage, mV
neuron.gEx     =     0;    % Total excitatory conductance
neuron.gIn     =     0;    % Total inhibitory conductance
neuron.spike   = false;    % 'true' if the neuron is spiking, else 'false'

% Neuron properties
neuron.tauM    =  .020;    % Membrane time constant, sec
neuron.Vrest   =   -70;    % Resting membrane voltage, mV
neuron.Eex     =     0;    % Excitatory reversal potential, mV
neuron.Ein     =   -70;    % Inhibitory reversal potential, mV
neuron.tauEx   =  .005;    % Time constant of excitatory conductances
neuron.tauIn   =  .005;    % Time constant of inhibitory conductances
neuron.Vthresh =   -54;    % Spike threshold voltage, mV
neuron.Vreset  =   -60;    % Post-spike reset voltage, mV

%% Define excitatory synapses
neuron.Nex                 = 1000; % # of excitatory synapses
neuron.exSynapses.rate     =   25; % Rate of pre-synaptic APs
% Maximum peak excitatory conductance
neuron.exSynapses.gMax     = .015;
% Peak excitatory conductance for each synapse a, these are
% initially set to gMax
neuron.exSynapses.gA       = neuron.exSynapses.gMax.*ones(neuron.Nex,1);
% A vector listing which presynaptic synapses are firing
neuron.exSynapses.preSynapticSpike = zeros(neuron.Nex,1);
% M and Pa are house keeping variables for implementing STDP
% rule.  Both decay exponentially toward zero.
%    Pa is positive, used to increase the strength of
%       synapses, and is incremented on pre-synaptic spiking.
%       On postsynaptic spiking, gA -> gA + P*gMax
%    M is negative, used to decrease the strength of synapses,
%       and is decremented on post-synaptic spiking.
%       On presynaptic spiking,  gA -> gA + M*gMax
neuron.exSynapses.Pa       = zeros(neuron.Nex,1);
neuron.exSynapses.M        = 0;
% STDP parameters for excitatory synapses
neuron.exSynapses.Aplus    =      .005; % Magnitude of synapse strengthening
neuron.exSynapses.Aminus   = 1.05*.005; % Magnitude of synapse weakening
neuron.exSynapses.tauPlus  =      .020; % Time constant of strengthening (sec)
neuron.exSynapses.tauMinus =      .020; % Time constant of weakening (sec)


%% Define inhibitory synapses (nb. These are not subject to STDP)
neuron.Nin                 = 200;  % # of inhibitory synapses
neuron.inSynapses.rate     =  10;  % Rate of presynaptic APs
% A vector listing which presynaptic synapses are firing
neuron.inSynapses.preSynapticSpike = zeros(neuron.Nin,1);
% Peak excitatory conductance for each synapse a
neuron.inSynapses.gA               = .050*ones(neuron.Nin,1);

end % End modelNeuron()