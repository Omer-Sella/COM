function [results] = TP0V_example(varargin)
% Example for TP0v determination of Vp, Vp/Vf, and ERL
% TP0V_eample(com_COM versions 3.4 and later
% Author: Richard Mellitz
% Date: 12-08-2021 
% Upadated: 10-26-2022 RIM
% Code snipits by Adam Gregory
% Syntax:
%    my_results=TP0V_example('com_m_file','config_file','fixture_file','s2p_measurement_file')
% com_m_file is without the .m extension
% example:
%    my_results=TP0V_example('com_ieee8023_93a_390','config_com_ieee8023_93a=3ck_SA _TP0V_08_17_2022.xlsx',{fixture s4p}, {s2p_measurement_file})
% 
%
%%

[com_m_file,varargin]  =varargin_extractor(varargin{:});
[config_file,varargin] =varargin_extractor(varargin{:});
[fixture_file,varargin]=varargin_extractor(varargin{:});
[s2p_measurement_file,varargin]=varargin_extractor(varargin{:});
% [fixture_file,varargin]=varargin_extractor(varargin{:});
% [s2p_RL_meas_file,varargin]=varargin_extractor(varargin{:});
% add com arguments
varargin2=varargin;
varargin=[varargin {'OP.pkg_len_select'} , {'[ 2 ]'} ]; % typically this will be the 2nd package
varargin=[varargin {'OP.AUTO_TFX'} , {'1'} ]; % Finds Tfx from the fixture parameters
%% very important assign the fixture ports correctly [ +BGA -BGA +SMA -SMA ]   <-------------------------
varargin=[varargin {'param.snpPortsOrder'},{'[2 4 1 3 ]'} ]; % Make sure Tx side is listed as the BGA port
varargin2= [varargin2 {'OP.ERL_ONLY'} , {'1'} ]; % varargin2 is for processing the s2p file
varargin2=[varargin2 {'OP.ERL'} , {'2'} ] ; % this selects s2p file as the RL file
varargin2=[varargin2 {'OP.TDR_W_TXPKG'} , {'0'} ] ;
varargin2=[varargin2 {'OP.RUNTAG'} , {'DIFF_RL_MEAS'} ] ;
%%
com_call=str2func(com_m_file);
results{1}=com_call(config_file,0,0,fixture_file,varargin{:});
%%
[~,ncases] = size(results{1});
if ncases ==1
    tfx_str=sprintf('%.10e', results{1, 1}.tfx_estimate);
    ERL_ref=results{1, 1}.ERL;
    Vf_ref=results{1, 1}.steady_state_voltage_mV;
    Pmax_by_Vf_ref= results{1, 1}.Pmax_by_Vf_est;
else
    tfx_str=sprintf('%.10e', results{1, 1}{1, 1}.tfx_estimate);
    ERL_ref=results{1, 1}{1, 1}.ERL;
    Vf_ref=results{1, 1}{1, 1}.steady_state_voltage_mV;
    Pmax_by_Vf_ref= results{1, 1}{1, 1}.Pmax_by_Vf_est;
end
display(tfx_str)
varargin2=[varargin2 'param.tfx', {tfx_str}];
if ~isempty(s2p_measurement_file)
    results{2}=com_call(config_file,0,0,s2p_measurement_file,varargin2{:});
end

dERL=results{1, 2}.ERL-ERL_ref;

if dERL > -3
    fprintf('<strong> \n  Results ... \n \t dERL = %.3g \n \t V^(ref)_f  = %.4g \n \t R^(ref)_peak = %.3g \n </strong>', dERL, Vf_ref, Pmax_by_Vf_ref )
else
    fprintf(2,'<strong> \n  Results ... \n \t dERL = %.3g \n \t V^(ref)_f = %.4g \n \t R^(ref)_peak = %.3g \n </strong>', dERL, Vf_ref, Pmax_by_Vf_ref )
end


end
%%
function [out_var,varg_out]=varargin_extractor(varargin)
% peels away first argurment.
if isempty(varargin)
    out_var=[];
    varg_out={};
else
    out_var=varargin{1};
    varg_out=varargin;
    varg_out(1)=[];
end

end