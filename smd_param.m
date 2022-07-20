% SMD_parameters

% This function calculates seismic parameters from an acceleration time
% series. Specifically, it calculates velocity vs time, displacement vs 
% time, peak ground acceleration (PGA), peak ground velocity (PGV), peak
% ground displacement (PGD), Arias intensity vs time, total Arias intensity
% (Ia), time between when 5% and 75% of Ia has occurred (significant
% duration D5-75), time between when 5% and 95% of Ia has occurred
% (significant duration D5-95), mean period (Tm), pseudo-acceleration
% response spectrum (Sa), pseudo-velocity response spectrum (Sv),
% displacement response spectrum (Sd), and the Fourier amplitude spectrum 
% (FAS).
%
% BY RASHID SHAMS (dt-4-April-2021)
%
% SYNTAX
% [param]=smd_param(time,acc,damp,LUF,HUF)
%
% MANDATORY INPUTS
% acc = acceleration vector in g
% time = time vector in seconds, must be the same length as acc
%
% OPTIONAL INPUTS
% damp = damping ratio for response spectra
%       default = 0.05
% LUF = lowest usable frequency of response spectra
%       default = 0.10
% HUF = highest usable frequency of response spectra
%       default = (1/dt)/2 (Nyquist frequency)
%
% OUTPUT
% param = MATLAB structure with the following fields
%
% param.vel = velocity time series in cm/s
% param.disp = displacement time series in cm
% param.PGA = peak ground acceleration in g
% param.PGV = peak ground velocity in cm/s
% param.PGD = peak ground displacement in cm
% param.aint2 = cumulative fraction of Arias intensity occurring with time
% param.arias = total arias intensity at end of time series in m/s
% param.D_5_75 = time between when 5% and 75% of Ia has occurred (significant
%       duration D5-75) in seconds
% param.t_5_75 = time when 5% and 75% of Ia has occurred in seconds
% param.D_5_95 = time between when 5% and 95% of Ia has occurred (significant
%       duration D5-95) in seconds
% param.t_5_95 = time when 5% and 95% of Ia has occurred in seconds
% param.Tm = mean period in seconds according to Rathje et al (2004)
% param.Period = periods for response spectra
% param.Sa = pseudo-acceleration response spectrum in g
% param.Sv = pseudo-velocity response spectrum in cm/s
% param.Sd = displacement response spectrum in cm
% param.FAS = Fourier amplitude spectrum in g
% param.freq = frequencies for Fourier amplitude spectrum in Hz

function [param]=smd_param(time,acc,damp,LUF,HUF)
acc = acc(:);
time = time(:);
g = 9.81;
A   = acc*g; %convert to m/s^2
dt = time(2)-time(1);
% Check what variables are specified by the user, if a variable is not
% specified, then assign the default value
if exist('damp','var') == 0;
    damp = 0.05;
end
if exist('LUF','var') == 0;
    LUF = 0.10;
end
if exist('HUF','var') == 0;
    HUF = (1/dt)/2;
    if HUF > 100;
        HUF = 100;
    end
end
% TIME SERIES
param.vel = cumsum(A)*dt*100; 
param.disp = cumsum(param.vel)*dt; 
% PEAK RESPONSES
param.PGA = max(abs(A))/g; 
param.PGV = max(abs(param.vel)); 
param.PGD = max(abs(param.disp)); 
% ARIAS INTENSITY
aint2 = cumsum(A.^2)*pi*dt/(2*g);
arias = aint2(end);
param.aint2 = aint2/arias; 
param.arias = arias; 
% DURATION
timed = time(aint2>=0.05*arias & aint2<=0.75*arias);
param.t_5_75 = [timed(1),timed(end)];
param.D_5_75 = timed(end)-timed(1);
timed = time(aint2>=0.05*arias & aint2<=0.95*arias);
param.t_5_95 = [timed(1),timed(end)];
param.D_5_95 = timed(end)-timed(1);
% RESPONSE SPECTRA
[Sa,Sv,Sd,T]=rs(acc,dt,damp,LUF,HUF);
param.Sd=Sd(:);
param.Sv=Sv(:);
param.Sa=Sa(:);
param.Period = T(:);
% FOURIER AMPLITUDE SPECTRUM
[f,U]=FAS(dt,acc);
param.FAS = U;
param.freq = f;
% MEAN PERIOD (Rathje et al, 2004)
fi = f(f>0.25 & f<20);
Ci = U(f>0.25 & f<20);
Tm = ((Ci(:)'.^2)*(1./fi(:)))/(Ci(:)'*Ci(:));
param.Tm = Tm;
function[Sa,Sv,Sd,T]=rs(acc,dt,damp,LUF,HUF)
Acccms=acc*981;%convert from g to cm/s^2
if dt > .005;
     beta = .25;
else beta = 1/6;
end
gamma= 0.5; %parameters for Newmark's method
%average acceleration method gamma = 0.5, beta = .25, linear acceleration 
%method gamma = 0.5, beta = 1/6.  Average acceleration method is
%unconditionally stable, but less accurate.  Linear acceleration method is
%stable for dt/T < 0.551 but more accurate (Chopra, 2011)
Tlong = LUF^-1; %lowest usable frequency = 1/max period
Tshort = HUF^-1; %highest usable frequency = 1/min period
T = 10.^linspace(log10(Tshort),log10(Tlong),150); %150 points
umax = zeros(1,length(T));
for j=1:length(T)
    wn = 2*pi/T(j);
    m = 1;%then c and k are in terms of damping and natural period
    k = wn^2;
    c = 2*wn*damp;
    khat = k+gamma/beta/dt*c+m/beta/dt^2;
    a = m/beta/dt+gamma*c/beta;
    b = 1/2/beta*m+dt*(gamma/2/beta-1)*c;
    u = zeros(length(Acccms),1); %oscillator starting from rest
    udot = zeros(length(Acccms),1);%pre-allocate for speed
    uddot = zeros(length(Acccms),1);
    du = zeros(length(Acccms)-1,1);
    dudot = zeros(length(Acccms)-1,1);
    duddot = zeros(length(Acccms)-1,1);
    for i = 1:length(Acccms)-1
       du(i) = (Acccms(i+1)-Acccms(i)+a*udot(i)+b*uddot(i))/khat;
       u(i+1) = u(i)+du(i);
       dudot(i) = gamma*du(i)/beta/dt-gamma*udot(i)/beta+dt*(1-gamma/2/beta)*uddot(i);
       udot(i+1) = udot(i)+dudot(i);
       duddot(i) = du(i)/beta/dt^2-udot(i)/beta/dt-uddot(i)/2/beta;
       uddot(i+1) = uddot(i)+duddot(i);
    end
    umax(j) = max(abs(u));%max displacement for every period T (cm) 
end
Sd = umax; %displacement in cm
Sv=2*pi*Sd./T;%pseudo velocity in cm/s
Sa=2*pi*Sv./T/981;%pseudo acceleration in g
function[f,U]=FAS(dt,acc)
Ny = (1/dt)/2; %Nyquist frequency (highest frequency)
L  = length(acc); %number of points in acc
NFFT = 2^nextpow2(L); % Next power of 2 from length of acc
df = 1/(NFFT*dt); %frequency spacing
U = abs(fft(acc,NFFT))*dt; %Fourier amplitudes 
U = U(2:Ny/df+1); %single sided FAS
f = linspace(df,Ny,Ny/df)'; %[small, large, number] frequencies