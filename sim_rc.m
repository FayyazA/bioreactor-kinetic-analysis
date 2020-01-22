function [data] = sim_rc(Vin, varargin)

%% Generate a circuit of the form in lab (resistor and capacitor in series)

%Resistor, in ohms
R = 10000;
%Possible capacitor values, in farads
possible_C = [.00047, .001, .0022, .0047, .01, .022]*10^(-6);
%Pick a random value from the 5 options, or assign the defined value
if isempty(varargin)
    C = randsample(possible_C,1);
else
    C = varargin{1};
end


%% Feed the input voltage step function through the circuit
%Measure voltage across capacitor as the output, which we will call Vout
%Assume Vout(0)=0, uncharged capacitor

%Create a 10ms time vector, in ms, sampled at 10MHz
%time steps = 1/100MHz = 1*10^(-4)ms/sample
dt = 1e-4;
t=0:dt:10;
Vout = zeros(size(Vin));

%Create a generic function handle that will generate an output at each time
%point given the previous voltage time since last sample, dt
rc_fxn = @(Vin, Vout_0, dt, R, C)Vin+(Vout_0-Vin)*exp(-(dt/1000)/(R*C));

for k = 2:length(Vin)
    Vout(k) = rc_fxn(Vin(k), Vout(k-1), dt ,R, C);
end

%% Plot response if the capacitor value is unknown
if isempty(varargin)
    figure; plot(t,Vout)
    ylabel('Voltage across C1 (V)')
    xlabel('Time (ms)')
    axis([0,10, -0.5,4.5])
end

%% Log data in output struct
data = struct('R',R, 'C',C ,'Vout',Vout);


