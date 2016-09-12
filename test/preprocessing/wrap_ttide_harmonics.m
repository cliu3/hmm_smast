% given elevation and time step, calculate harmonics for the folloing comps 
function [amp,phase] = wrap_ttide_harmonics(eta,dt_minutes,start_time)

%[NAME,FREQ,TIDECON,XOUT] = t_tide_quiet(eta,'interval',dt_minutes/60.,'start time',start_time,'output','none');
[NAME,FREQ,TIDECON,XOUT] = t_tide(eta,'interval',dt_minutes/60.,'start time',start_time,'output','none');
M2 = 0;
N2 = 0;
S2 = 0;
O1 = 0;
K1 = 0;
K2 = 0;
P1 = 0;
Q1 = 0;

dims = size(NAME);
for i=1:dims(1)
  if(sum(NAME(i,:) == 'M2  ')==4); M2 = i; end;
  if(sum(NAME(i,:) == 'N2  ')==4); N2 = i; end;
  if(sum(NAME(i,:) == 'S2  ')==4); S2 = i; end;
  if(sum(NAME(i,:) == 'O1  ')==4); O1 = i; end;
  if(sum(NAME(i,:) == 'K1  ')==4); K1 = i; end;
  if(sum(NAME(i,:) == 'K2  ')==4); K2 = i; end;
  if(sum(NAME(i,:) == 'P1  ')==4); P1 = i; end;
  if(sum(NAME(i,:) == 'Q1  ')==4); Q1 = i; end;
end;

comps = [M2,N2,S2,O1,K1,K2,P1,Q1]; 
for j=1:length(comps)
  i = comps(j);
  if(i ~= 0)
%    fprintf('%s %f %f \n',NAME(i,:),TIDECON(i,1),TIDECON(i,3));
    amp(j) = TIDECON(i,1);
    phase(j) = TIDECON(i,3);
  else
    amp(j) = 0.0;
    phase(j) = 0.0;
  end;
end;



