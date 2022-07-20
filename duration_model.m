% DURATION MODEL

dt = time(2)-time(1);g=9.81;
% ARIAS INTENSITY
%accn=accn.*100;   % in m/s2
aint2 = cumsum(accn.^2)*pi*dt/(2*g);
arias = aint2(end);  
aint2 = aint2/arias;

% DURATION (Effective Duration)
time=time(1:end-1);

[~,idx]=min(abs(aint2-0.05));
one=aint2(idx);
[~,idx]=min(abs(aint2-0.95));
two=aint2(idx);

[~,idx]=min(abs(time-one*arias));
timed(1) = time(idx);

[~,idx]=min(abs(time-two*arias));
timed(2) = time(idx);

t_5_95 = [timed(1),timed(end)];
D_5_95 = timed(end)-timed(1);