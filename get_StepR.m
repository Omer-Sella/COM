function result=get_StepR(ir,param,cb_step,ZT)
%ir = impulse response
%t_base=time array with equal time steps
%samp_UI = number of samples per UI for ir
%   result.SBR
% t for debug
t=(1/param.fb)/param.samples_per_ui*(0:length(ir)-1);

if cb_step
    Ag=1;
    dt=1/param.fb/param.samples_per_ui;
    edge_time=param.TR_TDR*1e-9;
    fedge=1/edge_time;
    tedge=0:dt:edge_time*2;
    %
    edge=Ag*(2*cos(2*pi*(tedge)*fedge/16-pi/4).^2-1);
    drive_pulse=[edge ones(1,param.samples_per_ui)];
    %pulse=filter(UI_ones,1,ir);
    
    pulse=filter(drive_pulse,1,ir);
else
    pulse=cumsum(ir);
end
TDR_response=(1+pulse)./(1-pulse)*ZT*2;
result.ZSR=TDR_response;
result.pulse=pulse;

