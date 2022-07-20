%% Working on Stong Motion Data 

%% Load the SMD 
  EQ_smd=load('Cape Mendocino.dat'); 

  accn_g=EQ_smd(:,2); Time_s=EQ_smd(:,1);
  accn_g=accn_g./9.81;    % convert m/s2 to g 
  
%% Plot EQ time-history
   plot(Time_s,accn_g,'b');
   
%% Plot the fourier spectrum of the given eq time history data

n = length(accn_g); 
Ts = Time_s(2)-Time_s(1);        % Sample Time
Fs = 1/Ts;                       % Sampling Frequency
NFFT = 2^nextpow2(n);            % Next power of 2 from length of data
Y = fft(accn_g,NFFT)/n;
f = Fs/2*linspace(0,1,NFFT/2+1);
Iv = 1:length(f);                 % Index Vector
plot(f,abs(Y(Iv)))
xlabel('Frequency (Hz)')
ylabel('Amplitude (g)')
title('Fourier Spectrum of the EQ Time History')
%%
for i=1:1500
    if i==1
        icr(i)=ans(i);
    else 
        icr(i)=ans(i)-ans(i-1);
    end
end
%%
for i=1:length(Time_s)
    ans(i)=cumsum((sqrt(cumtrapz((accn_g(i)).^2)))./Time_s(i));
end

        


