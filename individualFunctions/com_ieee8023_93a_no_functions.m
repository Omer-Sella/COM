function results=com_ieee8023_93a(varargin)
% This is NOT an official IEEE document.
%% Implementation example of annex 93A IEEE Std 802.3™
%   http://www.ieee802.org/3/ck/public/adhoc/index.html
% result=com_ieee8023_93a(config_file, num_fext, num_next [, <s4p files>])
% - config_file: xls, xls, mat file which contains configuration settings
% - num_fext: number of FEXT s4p files in the listfigure(300+package_testcase_i);
% - num_next: number of NEXT s4p files  in the list
% - <s4p_files>: (1+num_fext+num_nefxt) file names. If not supplied, program
%   will ask for each of the files interactively.opupu
%
% This program is intended for the development of standard specifications
% and reflects activity of IEEE P802.3bj, .3by, .3bm, .3bs, .3cd, .3ck
% found in Annex 93A IEEE Std 802.3™ and project =updates
% original proposal for COM may be found at
% http://www.ieee802.org/3/bj/public/jul12/mellitz_01_0712.pdf in July 2012
% from the following authors and affiliations in 2012.
%      Richard Mellitz, Intel Corporation
%      Charles Moore, Avago Technologies
%      Mike Dudek, QLogic Corporation
%      Mike Peng Li, Altera Corporation
%      Adee Ran, Intel Corporation
%
% Some of the authors and Contributors:
%   Adee Ran
%   Richard Mellitz
%   Yasuo Hidaka
%   John Ewen
%   Bill Kirkland
%   Adam Gregory
%   Howard Heck
%   Jingbo Li
%   Adam Healey
%   Matt Brown
%   Sameh Elnagar
%   Hossein Shakiba
zzz_list_of_changes()

%% Opening Preface
% acquire parsing command string and set up OP control structure. Then read in files
close(findall(0, 'tag', 'TMWWaitbar', '-or', 'tag', 'COM'));
try %  version number at end of call string
    cmdfile=mfilename;
    hindx=strfind(mfilename,'_');
    ver=cmdfile(hindx(end)+1:end);
    output_args.code_revision = [ver(1), '.',ver(2:end)];
catch
    output_args.code_revision ='';
end
teststr='';
OP.TESTING=0;
if OP.TESTING == 1 % set to 1 or pre release
    teststr='testing';
    testmsg=sprintf('Evaluation Copy: COM%s%s\n',output_args.code_revision,teststr);
    htest = msgbox(testmsg);
    set(htest,'Color','y', 'tag', 'COM');movegui(htest,'northeast');
end
display('This is NOT an official IEEE document.')
fprintf('Revision:<strong> %s%s </strong>This is a computation example for exploring COM and ERL  \n',output_args.code_revision,teststr)
disp(' for projects like IEEE P802.3bj/b/bs/cd/ck with some exploratory extensions and is not normative or official')
t0=tic;
set(0,'defaulttextinterpreter','none') % prevents subscripting character in displayed messages
% reset to tex on exit
%% file_setup
%%
% need to see what happens for version 8
if verLessThan('matlab', '7.4.1')
    error('Matlab version 7.4 or higher required')
end

results=[];

%% New Command Line parser
[config_file,num_fext,num_next,Remember_keyword,OP,varargin]=COM_CommandLine_Parse(OP,varargin{:});


%% get the first 3 arguments and allow for interactive input.
if isempty(config_file)
    config_file=input('Enter config XLS file or return will just pop a window to ask for the XLS file]:  ','s');
    if isempty(config_file)
        [config_file, config_file_path] = uigetfile([{ '*.xls;*.xlsx'} ; {'*.mat'}],'INPUT CONFIG FILE .xls');
    else
        [config_file_path,cname,cext]=fileparts(config_file);
        config_file=[cname cext];
    end
    if config_file==0
        % cancel - exit gracefully
        return;
    end
    config_file = fullfile(config_file_path, config_file);
end
output_args.config_file = config_file;
OP.SAVE_KEYWORD_FILE=0;
if OP.SAVE_KEYWORD_FILE
    if exist('keyworklog.mat','file')==2
        delete('keyworklog.mat');
    end
end
[param, OP] = read_ParamConfigFile(config_file,OP);
if OP.CONFIG2MAT_ONLY
    return;
end
if isempty(num_fext)
    if OP.RX_CALIBRATION
        num_fext=1;
        display('First prompt is for the measured test thru channel and following prompt is for Rx noise path channel')
    else
        if param.tfx(1) == -1 && OP.ERL_ONLY && OP.ERL == 2 % OP.ERL=2 is for package ERL, param.tfx = 1 means fixture time not set and need to be determinind for test fixture channel
            num_fext=1;
            display('First prompt is for the s2p measured data following prompt is for s4p of of the test fixtrure channel')
        elseif ~OP.ERL_ONLY
            num_fext=input('How many FEXT channels are to be entered? [return means no FEXT] ');
        else
            num_fext=0;
        end
    end
    if isempty(num_fext)==1, num_fext=0; end
end
if isempty(num_next)
    if param.tfx(1) == -1 && OP.ERL_ONLY && OP.ERL == 2 % OP.ERL=2 is for package ERL, param.tfx = 1 means fixture time not set and need to be determinind for test fixture channel
        num_next=0;
    else
        if OP.RX_CALIBRATION
            num_next=0;
        elseif ~OP.ERL_ONLY
            num_next=input('How many NEXT channels are to be entered? [return means no NEXT] ');
        else
            num_next=0;
        end
    end
    if isempty(num_next)==1, num_next=0; end
end
% Allow string inputs for running compiled version from OS command-line
if ischar(num_fext), num_fext=str2double(num_fext); end
if ischar(num_next), num_next=str2double(num_next); end
xtk=num_fext+num_next; % total number of crosstalk aggressors
param.num_next=num_next;
param.num_fext=num_fext;
param.num_s4p_files=num_fext+num_next+1;
% checking for data when running for rx compliance BBN calibration
if OP.RX_CALIBRATION == 1
    if num_fext ~=1
        h = msgbox('One and only noise path channel is required'); set(h,'Color',[1 .85 0]);
        movegui(h,'northwest')
        set(h,'Tag','COM') % RIM 06-13-2022 ... tak msg box for closing later
        if OP.DEBUG ~= 1
            return
        end
    end
    h = msgbox('Please make sure the measured "sigma_RJ", A_DD, and SNR_TX" fields in the config xls file have been modified from the Tx measurement. '); set(h,'Color',[0 1 1]);
    movegui(h,'southeast')
    set(h,'Tag','COM') % RIM 06-13-2022 ... tak msg box for closing later
end

if   param.tfx(1) == -1 && OP.ERL_ONLY && OP.ERL == 2
    if num_fext ~=1
        h = msgbox('One and only test channel is required'); set(h,'Color',[1 .85 0]);
         set(h,'Tag','COM') % RIM 06-13-2022 ... tak msg box for closing later
        movegui(h,'northwest')
        if OP.DEBUG ~= 1
            return
        end
    end
    h = msgbox('The test fixture file is use to gate measurements '); set(h,'Color',[0 1 1]);
    movegui(h,'southeast')
     set(h,'Tag','COM') % RIM 06-13-2022 ... tak msg box for closing later
end


% create result directory if needed
if ~exist(OP.RESULT_DIR,'dir'); mkdir(OP.RESULT_DIR); end
% allow finite impulse response input rather that s-parameters with
% OP.EXTERNAL = true. However the use_external_IR function is not provided
if ~isempty(varargin) % process case where file names are passed in function call
    if strfind(upper(char(varargin(1))),'EXTERNAL_IR') ~= 0
        error('External IR mode is no longer supported');
        %OP.EXTERNAL = true;
        %OP.GET_FD = 0;
        %ir1a= varargin(2);
        %ex_var = varargin(3);
        %[chdata OP param ]  = use_external_IR(param, OP ,num_fext,num_next,0,ir1a,ex_var);
    else
        if OP.TDMODE
            OP.GET_FD=false;
        end
        if length(varargin) < xtk +1 % check that number of varargin arguments passed is at least number of crosstalk files+1 (thru)
            error('files must include next + fext + a thru');
        end
        %% eveluate any extra arguments as possible modifications of parameters
        extra_args = varargin(xtk+2:end);
        for k=1:2:floor(length(extra_args)/2)*2
            try
                orig_value_is_str = 1;
                orig_value=eval(extra_args{k});
                if ~ischar(orig_value)
                    orig_value_is_str = 0;
                    orig_value=mat2str(orig_value);
                end
            catch eval_err
                if isequal(eval_err.identifier, 'MATLAB:nonExistentField')
                    % trying to modify a nonexistent parameter - probably a
                    % typo. save the user from his error.
                    error('COM:BadExtraParameter', 'Attempted override of a non-existing parameter %s.', extra_args{k});
                else
                    % unexpected condition
                    rethrow(eval_err);
                end
            end
            try
                if orig_value_is_str
                    mod_string = sprintf('%s = ''%s'';', extra_args{k}, extra_args{k+1});
                else
                    mod_string = sprintf('%s = %s;', extra_args{k}, extra_args{k+1});
                end
                eval(mod_string);
                fname=['mod_str' num2str(k)];
                % begin yasuo patch 2/11/2018
                % output_args.(fname)=mod_string;
                % If mod_string contains a comma, enclose it by double quotes to avoid misaligned column in the CSV output.
                
                % re-patch yasuo 3/18/2019
                % v2.56 if contains(mod_string,',')
                % v2.57 if isempty(strfind(mod_string,','))
                % Here, if-condition was inverted by the change of function from 'contains()' to 'isempty()'.
                % So, it is changed back by adding an '~' operator.
                % if isempty(strfind(mod_string,','))
                if ~isempty(strfind(mod_string,','))
                    output_args.(fname)=['"' mod_string '"'];
                else
                    output_args.(fname)=mod_string;
                end
                fprintf('Applied parameter modification: %s (override %s)\n', mod_string, orig_value);
            catch eval_err
                error(eval_err.identifier, 'Error evaluating "%s".', mod_string);
            end
        end
    end
end
%% Parameters computationally defined by values from the settings files
param.ui=1/param.fb;
param.sample_dt = param.ui/param.samples_per_ui;
param.sigma_X=sqrt( (param.levels^2-1)/ (3*(param.levels-1)^2) );
factor_3db=0.473037;
param.fb_BT_cutoff=factor_3db*param.f_r;
param.fb_BW_cutoff=param.f_r;
param.Tx_rd_sel=1;
param.Rx_rd_sel=2;
if isempty(param.snpPortsOrder) || any(isnan(param.snpPortsOrder))
    param.snpPortsOrder = [1 3 2 4]; % default order normally used.
end
%% size adjust vector parameters which may be entered as one element
param=parameter_size_adjustment(param,OP);

%% get input models
param.FLAG.S2P=0;
if OP.TDMODE
    OP.FIXTURE_CALIBRATION= 0;
    [chdata, param] = get_TD_files(param, OP, num_fext, num_next, varargin);
else
    OP.FIXTURE_CALIBRATION= 0;
    [chdata, param] = get_s4p_files(param, OP, num_fext, num_next, varargin);
    if any(strcmpi({chdata.ext},'.s2p'))
        param.FLAG.S2P=1;
    end
end

OP.SAVE_CMD_STR=1;
if OP.SAVE_CMD_STR
    cmd_str = save_cmd_line([Remember_keyword ''',''' config_file], chdata, num_fext,num_next,mfilename );
    setappdata(0,'cmd_str',cmd_str);
end
%% from here on, multiple package test cases are run. results will be saved separately.
results = cell(size(OP.pkg_len_select));
COM = inf;
min_COM=inf; % reset COM prior to calibration
% min_VEO = inf;
min_VEO_mV = inf;
max_VEC_dB = -inf;
threshold_DER=inf;
% begin yasuo patch 3/18/2019
threshold_DER_max = 0;	% reset worst DER
% end yasuo patch
sigma_bn=0;
DO_ONCE=true;
low_COM_found = 0;
% at this point only the impulse responses are needed. However vestiges of FD may be intermingled
while (OP.RX_CALIBRATION==1 || DO_ONCE==true)
    if ~DO_ONCE
        if abs(min_COM - param.pass_threshold)<0.1 || (sigma_bn==0 && min_COM < param.pass_threshold)
            break;
        elseif min_COM > param.pass_threshold
            % increase noise level linearly until low COM found; then perform binary search.
            if low_COM_found
                if OP.sigma_bn_STEP>0 % previous increase too small
                    OP.sigma_bn_STEP = OP.sigma_bn_STEP/2; % gearshift
                else % previously decrease too large
                    OP.sigma_bn_STEP = -OP.sigma_bn_STEP/2; % gearshift and change direction
                end
            end
        else % binary searchparam.Pkg_len_TX
            low_COM_found=1;
            if OP.sigma_bn_STEP>0 % previous increase too large
                OP.sigma_bn_STEP = -OP.sigma_bn_STEP/2; % gearshift and change direction
            else % previously decrease too small
                OP.sigma_bn_STEP = OP.sigma_bn_STEP/2; % gearshift
            end
        end
        min_COM = inf; % ignore previous iterations
        min_VEO_mV = inf;
        max_VEC_dB = -inf;
        sigma_bn = sigma_bn + OP.sigma_bn_STEP;
    end
    msgctr=1;
    for package_testcase_i = 1:length(OP.pkg_len_select)
        CSV_FILE=sprintf('%s%s_case%d_results.csv', OP.RESULT_DIR, chdata(1).base, package_testcase_i);
        package_testcase=OP.pkg_len_select(package_testcase_i);
        param.Pkg_len_TX = param.z_p_tx_cases(package_testcase,:);
        param.Pkg_len_NEXT = param.z_p_next_cases(package_testcase,:);
        param.Pkg_len_FEXT = param.z_p_fext_cases(package_testcase,:);
        param.Pkg_len_RX = param.z_p_rx_cases(package_testcase,:);
        param.AC_CM_RMS_TX= param.AC_CM_RMS(package_testcase);
        if param.PKG_Tx_FFE_preset ~=0
            param.Pkg_TXFFE_preset= param.PKG_Tx_FFE_preset(package_testcase,:);
        else
            param.Pkg_TXFFE_preset=0;
        end
        %         ki=package_testcase;
        % %         param.Pkg_Zc=[ param.pkg_Z_c(ki,1); param.pkg_Z_c(ki,2) ];
        %         param.Pkg_Zc=[ param.pkg_Z_c(ki,:) ];SDDp2p
        param.Pkg_Zc= param.pkg_Z_c;
        [cmele,centries] = size(param.Pkg_Zc);
        [mele, ncases] = size(param.Pkg_len_TX);
        if cmele ~=1 && centries ~=2 && mele ~= 1
            param.Pkg_Zc=reshape(param.Pkg_Zc,2,4);
        end
        param.package_testcase_i = package_testcase_i;
        
        %% Fill in chdata
        if OP.TDMODE
            [chdata, param ] = read_PR_files(param, OP, chdata);
            [chdata, param, SDDch, SDDp2p ] = TD_FD_fillin(param, OP, chdata); % fill in fd data to keep rest of SW happy
        else
            %fill in chada with s-parameters
            [chdata, SDDch, SDDp2p ] = read_s4p_files(param, OP, chdata);
            [chdata, param] = process_sxp(param, OP, chdata, SDDch);
        end
        if OP.BREAD_CRUMBS
            output_args.RL.f=chdata(1).faxis;  % RIM 07/19/2019 only use the first index
            output_args.RL.rl1=chdata(1).sdd11_raw; % RIM 07/19/2019 only use the first index
            if isfield(chdata(1),'sdd22_raw')% RIM 10/15/2019 only use the first index
                output_args.RL.r22=chdata(1).sdd22_raw; % RIM 07/19/2019 only use the first index
            end
            if isfield(chdata(1),'TX_RL')% RIM 10/09/2020 report Tx RL with RD
                output_args.RL.TXRL=chdata(1).TX_RL; %R IM 10/09/2020 report Tx RL with RD
            end
        end
        if param.FLAG.S2P, OP.ERL_ONLY =1;end
        
        %% Process TDR & ERL
        [output_args,ERL,min_ERL]=TDR_ERL_Processing(output_args,OP,package_testcase_i,chdata,param);
        if OP.ERL_ONLY
            results = cell(1);
            results{1} = output_args;
            rt=toc(t0);
            output_args.rtmin=rt/60;
            if 1
                fprintf('run time = %g min\n',output_args.rtmin)
            end
            if OP.CSV_REPORT ==1
                Write_CSV(output_args,CSV_FILE);
            end
            break;
        end
        
        %% FD processing s-parameter
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % at this point sdd21 responses and faxis (frequency) array are defined
        %most operations now wrapped into FD_Processing function
        param.number_of_s4p_files=length(chdata);
        %ICN=0;
        output_args.ICN_mV=0;
        output_args.MDNEXT_ICN_92_46_mV=0;
        output_args.MDFEXT_ICN_92_47_mV=0;
        if OP.WC_PORTZ
            param.SNR_TX=param.SNDR(param.Tx_rd_sel);
        else
            param.SNR_TX=param.SNDR(package_testcase);
        end
        
        %TD Mode now also calls FD_Processing but skips the main parts
        [chdata,output_args]=FD_Processing(chdata,output_args,param,OP,SDDp2p,DO_ONCE);
        
        %% Convert from Frequency Domain to Time Domain
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if DO_ONCE
            if ~OP.TDMODE
                chdata=COM_FD_to_TD(chdata,param,OP);
                output_args.VMC_HF_mV=chdata(1).VCM_HF_struct.DCn*1000;
                output_args.SCMR_dB=chdata(1).SCMR;
            end
        end
        
        %% Determine equalization settings
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %---------------------
        do_C2M=0;
        if param.T_O~=0 && param.Min_VEO_Test~=0
            do_C2M=1;
        end
        fom_result = optimize_fom(OP,param, chdata, sigma_bn,do_C2M);
        if fom_result.eq_failed ; return; end % RIM 12-20-2023
        
        %% Apply Equalization (returns pulse response with CTLE, TXLE, RXFFE)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        A_s=abs(fom_result.A_s); % this is the "s" in SNR (PAM4 gain in handled in last sections
        param.use_bmax=fom_result.best_bmax.';
        %AJG021820
        param.use_bmin=fom_result.best_bmin.';
        % Recommended Delta_y no larger than As/1000 or 0.01 mV
        param.current_ffegain=fom_result.best_current_ffegain;
        if OP.force_pdf_bin_size
            param.delta_y = OP.BinSize;
        else
            param.delta_y = min(A_s/1000, OP.BinSize);
        end
        % the pdf for PAM4 uses the full swing SBR but assigns voltage for the PDF accordingly
        if OP.RX_CALIBRATION, param.number_of_s4p_files=1; end
        
        chdata=Apply_EQ(param,fom_result,chdata,OP);
        
        PSD_results=[];
        if (strcmp(OP.FFE_OPT_METHOD,'MMSE') && OP.RxFFE)
            OP.WO_TXFFE=1;
            %  chdata(1).eq_pulse_response includes rx FFE from Apply_EQ
            PSD_results = get_PSDs(PSD_results,chdata(1).eq_pulse_response,fom_result.t_s, fom_result.txffe,param.ctle_gdc_values(fom_result.ctle),param.g_DC_HP_values(fom_result.best_G_high_pass),param,chdata,OP);
            OP.WO_TXFFE=0;
            PSD_results = get_PSDs(PSD_results,chdata(1).eq_pulse_response,fom_result.t_s, fom_result.txffe,param.ctle_gdc_values(fom_result.ctle),param.g_DC_HP_values(fom_result.best_G_high_pass),param,chdata,OP);
        end
        %% Create ISI PDF & Individual Crosstalk PDFs
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~OP.DISPLAY_WINDOW, fprintf('processing COM PDF '); end
        for i=1:param.number_of_s4p_files
            if ~OP.DISPLAY_WINDOW, fprintf('%d ', i); end
            if strcmp(OP.FFE_OPT_METHOD,'MMSE') && OP.RxFFE

                pdf = get_pdf(chdata(i), param.delta_y, fom_result.t_s, param, OP,fom_result.PSD_results.iphase(i)) ;
            else
                pdf = get_pdf(chdata(i), param.delta_y, fom_result.t_s, param, OP,[]);
            end
            if OP.DEBUG && OP.DISPLAY_WINDOW && ~OP.RX_CALIBRATION
                figure(150+package_testcase_i);set(gcf,'Tag','COM');
                subplot(2,1,2);
                pdf0.x=pdf.x(pdf.y~=0);
                pdf0.y=pdf.y(pdf.y~=0);
                semilogy(pdf0.x, pdf0.y,'disp', chdata(i).base);
                current_ylim=ylim; ylim([param.specBER/100, current_ylim(2)]);
                hold on; title('PDF')
                recolor_plots(gca);
            end
            
            chdata(i).pdfr=pdf;
            % reporting
            a=find(cumsum(chdata(i).pdfr.y) >1e-12,1,'first');
            chdata(i).maxquickpdf=(chdata(i).pdfr.y(a));
            
        end
        if ~OP.DISPLAY_WINDOW, fprintf('\n'); end
        
        %% Return final PDF & CDF and Package all noise parameters in Noise_Struct
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [PDF,CDF,Noise_Struct]=Create_Noise_PDF(A_s,param,fom_result,chdata,OP,sigma_bn,PSD_results);
        combined_interference_and_noise_pdf=PDF;
        combined_interference_and_noise_cdf=CDF;

              
        %% Calculate COM and other associated outputs
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % The noise and interference amplitude, A_ni, is the magnitude of the value of y0
        % that satisfies the relationship P(y0) = DER_0
        A_ni_ix=find(combined_interference_and_noise_cdf>param.specBER, 1, 'first');
        A_ni = abs(combined_interference_and_noise_pdf.x(A_ni_ix));

        % begin yasuo patch 3/18/2019
        % estimate DER at threshold COM
        threshold_ix=find(combined_interference_and_noise_pdf.x>-A_s/(10^(param.pass_threshold/20)),1);
        threshold_DER=combined_interference_and_noise_cdf(threshold_ix);
        threshold_DER_max = max(threshold_DER_max, threshold_DER);
        % end yasuo patch
        
        if OP.RX_CALIBRATION ==0 && OP.EW == 1
            [Left_EW,Right_EW,eye_contour,EH_T_C2M,EH_B_C2M]=COM_eye_width(chdata,param.delta_y,fom_result,param,OP,Noise_Struct,0);
            EW_UI=floor((Left_EW+Right_EW))/param.samples_for_C2M;
            if OP.DISPLAY_WINDOW && OP.DEBUG
                figure_name =  'Eye at DER0 estimate';
                fig=findobj('Name', figure_name);
                if isempty(fig), fig=figure('Name', figure_name); end
                figure(fig);set(gcf,'Tag','COM');
                movegui(fig,'southwest')
                plot(eye_contour)
                xlabel('UI %')
                ylabel('V')
            end
            
        else
            EW_UI=0;
            eye_contour=[];
        end
        if OP.MLSE==0
            if param.T_O ~=0
                eye_opening=EH_T_C2M-EH_B_C2M;
                A_ni=2*A_s-eye_opening;
                %eq 124E-4
                vec_arg=2*A_s/eye_opening;
                if vec_arg<eps
                    vec_arg=eps;
                end
                VEC_dB = 20*log10(vec_arg);
                COM=20*log10(2*A_s/A_ni);
                VEO_mV=eye_opening*1000;
                min_COM = min(min_COM, COM);
                min_VEO_mV = min(min_VEO_mV,VEO_mV);
                max_VEC_dB = max(max_VEC_dB, VEC_dB);
            else
                VEO_mV = 1000*(A_s-A_ni)*2;
                vec_arg=(A_s-A_ni)/A_s;
                if vec_arg<eps
                    vec_arg=eps;
                end
                VEC_dB = -20*log10(vec_arg);
                COM=20*log10(A_s/A_ni);
                min_COM = min(min_COM, COM);
                min_VEO_mV = min(min_VEO_mV,VEO_mV);
                max_VEC_dB = max(max_VEC_dB, VEC_dB);
            end
            MLSE_results=struct;
        else % MLSE case
            if OP.MLSE==1 
                [MLSE_results] =  MLSE(param,fom_result.DFE_taps(1),A_s,A_ni,PDF,CDF);
            elseif OP.MLSE==2
                [MLSE_results] =  MLSE_instu(param,fom_result.DFE_taps(1),A_s,A_ni,PDF,CDF);
            else
                warning('unsuported MLSE option')
            end
            if param.T_O ~=0
                eye_opening=EH_T_C2M-EH_B_C2M;
                A_ni=2*A_s-eye_opening;
                %eq 124E-4
                vec_arg=2*A_s/eye_opening;
                if vec_arg<eps
                    vec_arg=eps;
                end
                VEC_dB_orig = 20*log10(vec_arg); % was negative in 400 beta1 ... Fixed 2-2-23
                VEC_dB=MLSE_results.delta_com_CDF-20*log10( (10^(MLSE_results.delta_com_CDF/20)-1)*10^(VEC_dB_orig/20)+1)+VEC_dB_orig;
                COM_orig=20*log10(2*A_s/A_ni);
                COM=MLSE_results.COM_CDF;
                VEO_mV=eye_opening*1000;
                min_COM = min(min_COM, COM);
                min_VEO_mV = min(min_VEO_mV,VEO_mV);
                max_VEC_dB = max(max_VEC_dB, VEC_dB);
                output_args.delta_COM=MLSE_results.delta_com_CDF;
                output_args.delta_VEC=MLSE_results.delta_com_CDF-20*log10( (10^(MLSE_results.delta_com_CDF/20)-1)*10^(VEC_dB_orig/20)+1);
            else
                VEO_mV = 1000*(A_s-A_ni)*2;
                vec_arg=(A_s-A_ni)/A_s;
                if vec_arg<eps
                    vec_arg=eps;
                end
                VEC_dB_orig = -20*log10(vec_arg);
                VEC_dB=MLSE_results.delta_com_CDF-20*log10( (10^(MLSE_results.delta_com_CDF/20)-1)*10^(VEC_dB_orig/20)+1)+VEC_dB_orig;
                COM_orig=20*log10(A_s/A_ni);
                COM=MLSE_results.COM_CDF;
                min_COM = min(min_COM, COM);
                min_VEO_mV = min(min_VEO_mV,VEO_mV);
                max_VEC_dB = max(max_VEC_dB, VEC_dB);
                output_args.delta_COM=MLSE_results.delta_com_CDF;
            end
        end
        
        %% Create COM_SNR_Struct to hold the main COM outputs
        COM_SNR_Struct.A_s=A_s;
        COM_SNR_Struct.A_ni=A_ni;
        COM_SNR_Struct.threshold_DER=threshold_DER;
        COM_SNR_Struct.EW_UI=EW_UI;
        COM_SNR_Struct.COM=COM;
        COM_SNR_Struct.VEC_dB=VEC_dB;
        if OP.MLSE == 0
            COM_SNR_Struct.COM_orig=[];
            COM_SNR_Struct.VEC_dB_orig=[];
        else
            COM_SNR_Struct.COM_orig=COM_orig;
            COM_SNR_Struct.VEC_dB_orig=VEC_dB_orig;
        end
        COM_SNR_Struct.VEO_mV=VEO_mV;
        COM_SNR_Struct.combined_interference_and_noise_pdf=combined_interference_and_noise_pdf;
        COM_SNR_Struct.combined_interference_and_noise_cdf=combined_interference_and_noise_cdf;
        COM_SNR_Struct.eye_contour=eye_contour;
        
        
        %% Save TD
        if OP.SAVE_TD
            sbr=timeseries(fom_result.sbr,fom_result.t);
            if ~OP.TDMODE
                fir=timeseries(fom_result.IR,fom_result.t);
            end
            for i=1:param.number_of_s4p_files
                Pulses(i).uneq_responce= timeseries(chdata(i).uneq_pulse_response, chdata(i).t );
                Pulses(i).eq_responce= timeseries(chdata(i).eq_pulse_response, chdata(i).t );
                if ~OP.TDMODE
                    FIR(i).uneq_imp_response=  timeseries(chdata(i).uneq_imp_response, chdata(i).t );
                    FIR(i).eq_imp_response=  timeseries(chdata(i).eq_imp_response, chdata(i).t );
                end
            end
            if OP.TDMODE
                save( [OP.RESULT_DIR 'sbr_fir_' param.base '.mat'],'sbr','Pulses');
            else
                save( [OP.RESULT_DIR 'sbr_fir_' param.base '.mat'],'sbr','fir','Pulses','FIR')
            end
        end
        
        %% Bathtub/Contribution Plot
        if OP.DISPLAY_WINDOW && ~OP.RX_CALIBRATION
            Bathtub_Contribution_Wrapper(COM_SNR_Struct,Noise_Struct,param,chdata,OP);
        end
        
        %% Msg management
        if ~exist('msg','var')
            msg=[];
        end
        if OP.DEBUG
            [ncases, mele]=size(param.z_p_tx_cases);
            switch param.flex
                case 4
                    msg = sprintf('%s: Case %g: z_p=(%g:%g:%g:%g, %g:%g:%g:%g, %g:%g:%g:%g, %g:%g:%g:%g) (TX, RX, NEXT, FEXT):\n', ...
                        msg,package_testcase_i, param.Pkg_len_TX, param.Pkg_len_RX, param.Pkg_len_NEXT, param.Pkg_len_FEXT ...
                        );
                case 2
                    msg = sprintf('%s: Case %g: z_p=(%g:%g, %g:%g, %g:%g, %g:%g) (TX, RX, NEXT, FEXT):\n', ...
                        msg,package_testcase_i, param.Pkg_len_TX(1:2), param.Pkg_len_RX(1:2), param.Pkg_len_NEXT(1:2), param.Pkg_len_FEXT(1:2) ...
                        );
                otherwise
                    msg = sprintf('%s: Case %g: z_p=(%g, %g, %g, %g) (TX, RX, NEXT, FEXT):', ...
                        msg, package_testcase_i, param.Pkg_len_TX, param.Pkg_len_RX, param.Pkg_len_NEXT, param.Pkg_len_FEXT ...
                        );
                    
            end
        else
            msg = sprintf('Case %d:', package_testcase_i );
        end
        
        if OP.TDMODE
            min_ERL=inf;
            ERL=[inf inf];
        end
        [msg] = end_display_control(msg,param,OP,output_args,COM,min_ERL,ERL, VEO_mV,VEC_dB,threshold_DER,OP.DISPLAY_WINDOW); % {} forces no ERL print
        
        
        %% Output Args
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        output_args=Output_Arg_Fill(output_args,sigma_bn,Noise_Struct,COM_SNR_Struct,param,chdata,fom_result,OP);
        rt=toc(t0);
        output_args.rtmin=rt/60;
        
        if OP.BREAD_CRUMBS
            output_args.OP=OP;
            output_args.param=param;
            output_args.chdata=chdata;
            output_args.fom_result = fom_result;
            output_args.PDF=PDF; % for exploration
            output_args.CDF=CDF; % for exploration
            output_args.MLSE_results=MLSE_results;
        end
        % results{package_testcase_i} = output_args;% moved RIM 04-14-2023
        
        %% making csv file
        if OP.CSV_REPORT ==1
            Write_CSV(output_args,CSV_FILE);
        end
        %% making mat file
        if(OP.DEBUG)
            save (sprintf('%s%s_case%d_results.mat', OP.RESULT_DIR, chdata(1).base, package_testcase_i), ...
                'output_args','param','OP');
        end
        if 1
            fprintf(' Die to die loss = %g dB \n',output_args.IL_db_die_to_die_at_Fnq)
            fprintf('run time = %g min \n',output_args.rtmin)
        end
        
        if nargout==0
            fprintf('<strong>--- Testcase %d results ---</strong>\n', package_testcase_i);
            disp(output_args)
        end
        
        if OP.BREAD_CRUMBS
            [my_path,rootname]=fileparts(chdata(1).filename);
            if ~isempty(OP.BREAD_CRUMBS_FIELDS)
                %Attempt to reduce the size of output_args.chdata by removing certain fields
                try
                    output_args.chdata=Bread_Crumb_Chdata_Reduction(output_args.chdata,OP.BREAD_CRUMBS_FIELDS);
                catch
                    fprintf('Failed to reduce output_args.chdata\n');
                end
            end
            save (sprintf('%s%s_case%d_results.mat', OP.RESULT_DIR, rootname, package_testcase), ...
                'output_args','param','OP');
        end
        
        results{package_testcase_i} = output_args; % moved to after chdata field reduction RIM 04-14-2023
    end
    [tmp] = end_display_control('WC All cases',param,OP,output_args,min_COM,min_ERL,ERL,min_VEO_mV,max_VEC_dB,threshold_DER,0);
    %%

    if OP.RX_CALIBRATION ==1
        sigma_hp= Noise_Struct.sigma_hp; % added for clause 162 else sigma_bn = sigma_hp (RIM 09-30-2022)
        display ([' LOOP with [sigma_bn sigma_hp] = [' num2str(sigma_bn) ' ' num2str(sigma_hp) '] performed with COM = ' num2str(min_COM) ])
    end
    DO_ONCE=false;
end

%% Final cleanup
if OP.DISPLAY_WINDOW
    savefigs(param, OP);
    set(0,'defaulttextinterpreter','tex'); % reset defaut text interpreter to tex
end

if OP.RX_CALIBRATION==1 % updated for clause 162 else sigma_bn = sigma_hp (RIM 09-30-2022)
    if ~param.f_hp==0
        fprintf ('Set Tx calibration noise(sigma_hp) rms voltage to %g mV\n', sigma_hp*1000);
        if OP.DISPLAY_WINDOW
            message=sprintf('Set Tx calibration noise (sigma_hp) rms voltage to %g mV.',sigma_hp*1000);
            hlast = msgbox(message,'sigma_hp','help');
            set(hlast,'Color','y', 'tag', 'COM');
        end
    else
        fprintf ('Set calibration noise (sigma_bn)rms voltage to %g mV\n', sigma_bn*1000);
        if OP.DISPLAY_WINDOW
            message=sprintf('Set calibration noise rms (sigma_bn) voltage to %g mV.',sigma_bn*1000);
            hlast = msgbox(message,'sigma_bn','help');
            set(hlast,'Color','y', 'tag', 'COM');
        end
    end
end

if length(results)==1, results = results{1}; end
redo_cmd_str=' redo string is: eval([''My_var_0 = '' getappdata(0,''cmd_str'')])';
disp(redo_cmd_str);
if isdeployed
    if OP.exit_if_deployed
        quit
    end
end
