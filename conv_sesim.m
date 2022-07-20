% Code to convert seismogram data from one type to another

% Variables:-
%   
%   datain       =   input waveform data of type datain_type
%
%   dataout      =   output waveform data of type dataout_type
%
%                    1 - Displacement
%   datain_type  =   2 - Acceleration
%                    3 - velocity 
%                   
%                    1 - Displacement
%   dataout_type =   2 - Velocity
%                    3 - Acceleration

%vel=VDATA(:,1);time=t; clearvars -except vel time;
datain_type=input('Enter the input data type');
dataout_type=input('Enter the output data type');

% INPUT DISPLACEMENT DATA AS displ, ACCELERATION DATA AS acc and TIME DATA
% IN SECONDS AS time.

% INTEGARTION 
if datain_type==2
  
    if  dataout_type==2
        velocity=cumtrapz(time,acc);
        data_out=velocity;
    elseif dataout_type==1
        v=cumtrapz(time,acc);
        displ=cumtrapz(time,v);
        data_out=displ;
    elseif dataout_type==3
        disp('The input is not correct');
    end

% DIFFERENTIATION
elseif  datain_type==1
    
    if     dataout_type==2
           s=displ;
           t=time;
           ds=diff(s);dt=diff(t);
           velocity=ds./dt;
           data_out=velocity;
           
    elseif dataout_type==3
           s=displ;
           t=time;
           ds=diff(s);dt=diff(t);
           v=ds./dt;
           dv=diff(v);
           accn=dv./dt(2:end);
           data_out=accn;
    elseif dataout_type==1
           disp('The input is wrong');
    end
    
elseif datain_type==3
       
       if     dataout_type==3
              v=vel;
              t=time;
              dv=diff(v);dt=diff(t);
              accn=dv./dt;
              data_out=accn;
           
       elseif dataout_type==1
              v=vel;t=time;
              disp=cumtrapz(t,v);
              data_out=disp;
       end    
            
end

accn_g=accn./981;Time_s=t(2:end);

        