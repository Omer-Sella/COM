function results=com_ieee8023_93a(config_file, num_fext, num_next, varargin)
sprintf('%d',length(varargin))
disp(varargin)
disp(varargin{1})
% Implementation of Annex 93A of IEEE 802.3
%
% Usage:
% result=com_ieee8023_93a(config_file, num_fext, num_next [, <s4p files>])
% - config_file: XLS file which contains configuration settings (samples in
%   http://www.ieee802.org/3/bj/public/tools/ran_3bj_com_d2p2_02_0813.zip)
% - num_fext: number of FEXT s4p files in the list
% - num_next: number of NEXT s4p files in the list
% - <s4p_files>: (1+num_fext+num_next) file names. If not supplied, program
%   will ask for each of the files interactively.
%
% This program is intended for the development of channel specifications
% and reflects the work of IEEE P802.3bj, Annex 93A.
% It is not an official IEEE document.


% --- Internal comments - do not remove the empty line above this line ---
% $Id: com_ieee8023_93a.m.rca 1.50 Sat Apr 12 18:33:50 2014 aran Experimental $
%
% Original authors:
%   Adee Ran (adee.ran@intel.com)
%   Richard Mellitz (richard.mellitz@intel.com)

% structures:
% chdata(i)
%           i= 1 --> THRU index
%           i= 2, num_fext+1 --> FEXT channel index
%           i= num_fext+2, num_next+num_fext+1
%       base: name of THRU file
%          A: amplitude
%       type: 'THRU', 'NEXT', or 'FEXT'
%        ftr: Rise time frequency
%      fmaxi: max number of frequency points
%      faxis: frequency array [Hz]
%      sdd21: (Htot) rewritten as product of ,vtf based on pkg RL ,tx filter, Rx filter
%      sdd22: differential RL
%      sdd11: differential RL
%     sdd21p: vtf based on pkg RL , this an interim parameter and set to sdd21
%     sdd21f: raw differential IL not filtered use for FD plots

%% file_setup
% acquire parsing command string and set up OP control structure. Then read in files
close all; close(findall(0, 'tag', 'TMWWaitbar', '-or', 'tag', 'COM'));
display('COM tool for P802.3bj Draft 3.2/P802.3bm Draft 2.2 ($Revision: 1.50 $)')
display('This is not an official IEEE document.')
display('Annex 93A is normative for some implementations. This code is a sample implemention of Annex 93A and is not normative.')

% need to see what happens for version 8

if verLessThan('matlab', '7.4.1')
    error('Matlab version 7.4 or higher required')
end
% get the first 3 arguments and allow for interactive input.
if ~exist('config_file','var')
    config_file=input('Enter config XLS file or return will just pop a window to ask for the XLS file]:  ','s');
    if isempty(config_file)
        [config_file, config_file_path] = uigetfile({'*.xls'},'INPUT CONFIG FILE .xls');
    end
    if config_file==0
        % cancel - exit gracefully
        return;
    end
    config_file = fullfile(config_file_path, config_file);
end

[param OP] = read_ParamConfigFile(config_file);

if ~exist('num_fext','var')
    if OP.RX_CALIBRATION
        num_fext=1;
        display('First prompt is for the measured test thru channel and following prompt is for Rx noise path channel')
    else
        num_fext=input('How many FEXT channels are to be entered? [return is means no FEXT] ');
    end
    if isempty(num_fext)==1, num_fext=0; end
end
if ~exist('num_next','var')
    if OP.RX_CALIBRATION
        num_next=0;
    else
        num_next=input('How many NEXT channels are to be entered? [return is means no NEXT] ');
    end
    if isempty(num_next)==1, num_next=0; end
end
xtk=num_fext+num_next; % total number of crosstalk aggressors
param.num_next=num_next;
param.num_fext=num_fext;
param.num_s4p_files=num_fext+num_next+1;

% checking for data when running for rx compliance BBN calibration
if OP.RX_CALIBRATION == 1
    if num_fext ~=1
        h = msgbox('One and only noise path channel is required'); set(h,'Color',[1 .85 0]);
        movegui(h,'northwest')
        if OP.DEBUG ~= 1
            return
        end
    end
    h = msgbox('Please make sure the measured "sigma_RJ", A_DD, and SNR_TX" fields in the config xls file have been modified from the Tx measurement. '); set(h,'Color',[0 1 1]);
    movegui(h,'southeast')
end

% create result directory if needed
if ~exist(OP.RESULT_DIR,'dir'); mkdir(OP.RESULT_DIR); end

% allow finite impulse response input rather that s-parameters with
% OP.EXTERNAL = true. However the use_external_IR function is not provided
if ~isempty(varargin) % process case where file names are passed in function call
    if strfind(upper(char(varargin(1))),'EXTERNAL_IR') ~= 0
        OP.EXTERNAL = true;
        OP.GET_FD = 0;
        ir1a= varargin(2);
        ex_var = varargin(3);
        [chdata OP param ]  = use_external_IR(param, OP ,num_fext,num_next,0,ir1a,ex_var);
    else
        OP.EXTERNAL = false;
        if length(varargin) < xtk +1 % check that number of varargin arguments passed is at least number of crosstalk files+1 (thru)
            error('files must include next + fext + a thru');
        end

        %% eveluate any extra arguments as possible modifications of parameters
        extra_args = varargin(xtk+2:end);
        for k=1:2:floor(length(extra_args)/2)*2
            try
                eval(extra_args{k});
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
                mod_string = sprintf('%s = %s;', extra_args{k}, extra_args{k+1});
                eval(mod_string);
                fprintf('Applied parameter modification: %s\n', mod_string);
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
%% size adjust vector parameters which may be entered as one element
if numel(param.C_pkg_board)==1, param.C_pkg_board = param.C_pkg_board*[1 1]; end
if numel(param.C_diepad)==1, param.C_diepad = param.C_diepad*[1 1]; end
if numel(param.R_diepad)==1, param.R_diepad = param.R_diepad*[1 1]; end
if numel(param.CTLE_fp1)==1, param.CTLE_fp1 = param.CTLE_fp1*ones(size(param.ctle_gdc_values)); end
if numel(param.CTLE_fp2)==1, param.CTLE_fp2 = param.CTLE_fp2*ones(size(param.ctle_gdc_values)); end
if numel(param.CTLE_fz)==1, param.CTLE_fz = param.CTLE_fz*ones(size(param.ctle_gdc_values)); end


if ~OP.EXTERNAL
    [chdata, param] = get_s4p_files(param, OP, num_fext, num_next, varargin);
end

%% from here on, multiple package test cases are run. results will be saved separately.
results = cell(size(OP.pkg_len_select));
COM = inf;
min_COM=inf; % reset COM prior to calibration
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
        else % binary search
            low_COM_found=1;
            if OP.sigma_bn_STEP>0 % previous increase too large
                OP.sigma_bn_STEP = -OP.sigma_bn_STEP/2; % gearshift and change direction
            else % previously decrease too small
                OP.sigma_bn_STEP = OP.sigma_bn_STEP/2; % gearshift
            end
        end
        min_COM = inf; % ignore previous iterations
        sigma_bn = sigma_bn + OP.sigma_bn_STEP;

    end
    for package_testcase_i = 1:length(OP.pkg_len_select)
        package_testcase=OP.pkg_len_select(package_testcase_i);
        param.Pkg_len_TX = param.z_p_tx_cases(package_testcase);
        param.Pkg_len_NEXT = param.z_p_next_cases(package_testcase);
        param.Pkg_len_FEXT = param.z_p_fext_cases(package_testcase);
        param.Pkg_len_RX = param.z_p_rx_cases(package_testcase);
        param.package_testcase_i = package_testcase_i;

        if ~OP.EXTERNAL
            %fill in chada with s-parameters
            chdata = read_s4p_files(param, OP, chdata);
        end
        number_of_s4p_files=length(chdata);

        %% FD processing s-parameter
        % at this point sdd21 responses and faxis (frequency) array are defined
        ICN=0;
        if OP.GET_FD && DO_ONCE
            f2=param.fb;
            for i=1:number_of_s4p_files
                if i == 2, PSXT(1:length(chdata(i).sdd21f))=0;end
                a=find(chdata(i).faxis(:)>f2);
                if isempty(a)
                    f2=chdata(i).faxis(end);
                    index_f2=length(chdata(i).faxis);
                else
                    index_f2=a(1);
                end
                % R is the frequency dependent parameter for the sinc function use in the
                % PWF for ICN
                temp_angle=(param.samples_per_ui*param.sample_dt)*pi.*chdata(i).faxis;
                if(chdata(i).faxis(1)==0)
                    temp_angle(1)=1e-20;% we don't want to divide by zero
                end
                SINC = sin(temp_angle)./temp_angle;
                PWF_data=SINC.^2;
                PWF_trf=(1+(chdata(i).faxis/chdata(i).ftr).^4).^-1;
                %// bw1=2.613126; bw2=3.4142136; bw3=2.613126;
                fr=param.f_r*param.fb;
                PWF_rx=(1+(chdata(i).faxis/fr).^8).^-1;
                %// PWF_rx=1./(1+(chdata(i).faxis./(param.fb*.75)).^8); % receiver filter
                PWF_highpass=1;
                PWF=PWF_data.*PWF_trf.*PWF_rx.*PWF_highpass; % power weight function
                % freq delta for integration
                chdata(i).delta_f=chdata(i).faxis(11)-chdata(i).faxis(10);
                % from ba spec, this is basically ICN
                chdata(i).sigma_FD=sqrt(2*chdata(i).delta_f/param.fb*sum( PWF(1:index_f2).*abs(chdata(i).sdd21f(1:index_f2)).^2))/2;
                faxis_GHz = chdata(i).faxis/1e9;

                if isequal(chdata(i).type, 'THRU')
                    [chdata(i).iln chdata(i).fit] = get_ILN(chdata(i).sdd21f(1:index_f2), chdata(i).faxis(1:index_f2));
                    nq_indx=find(chdata(i).faxis >= 1/param.ui/2,1, 'first');
                    chdata(i).fit_ILatNq=-20*log10(abs(chdata(i).fit(nq_indx)));
                    chdata(i).ILatNq=-20*log10(abs(chdata(i).sdd21f(nq_indx)));
                    fit_loss=chdata(i).fit_ILatNq;
                    Nq_loss=chdata(i).ILatNq;
                    ILD=20*log10(abs(chdata(i).fit))-20*log10(abs(chdata(i).sdd21f(1:index_f2)));
                    ILD_RMS=sqrt(2*chdata(i).delta_f/param.fb*sum( PWF(1:index_f2).*ILD(1:index_f2).^2));
                    if OP.DEBUG
                        if OP.DISPLAY_WINDOW
                            figure(300+package_testcase_i);
                            screen_size=get(0,'ScreenSize');
                            pos = get(gcf, 'OuterPosition');
                            set(gcf, 'Name', 'Raw frequency-domain data', 'OuterPosition', ...
                                screen_size([3 4 3 4]).*[0 1 0 0] + pos([3 4 3 4]).*[0 -2 1 2] ...
                                - (package_testcase_i-1)*[0 20 0 0]);

                            subplot(3,1,1)
                            title('Losses')
                            plot(faxis_GHz, 20*log10(abs(chdata(i).sdd21f)), 'b', 'LineWidth', 3, 'Disp','IL w/o packages')
                            hold on
                            plot(faxis_GHz, 20*log10(abs(chdata(i).sdd21)), 'b--', 'Disp','IL with packages')
                            ylim(get(gca, 'ylim'));
                            plot(faxis_GHz(1:index_f2), 20*log10(abs(chdata(i).fit(1:index_f2))),'g','Disp','IL fit')
                            plot(faxis_GHz, 20*log10(abs(chdata(i).sdd11)),'c','Disp','RL11')
                            plot(faxis_GHz, 20*log10(abs(chdata(i).sdd22)),'m','Disp','RL22')
                            subplot(3,1,3)
                            plot(faxis_GHz(1:index_f2), ILD,'Disp','ILD')
                        else
                            display(['Insertion Loss at Nyquist = ', num2str(chdata(i).ILatNq)])
                        end
                    end
                else % NEXT or FEXT
                    PSXT=sqrt((abs(chdata(i).sdd21f)).^2+PSXT.^2); % power sum xtk
                    ICN=sqrt(2*chdata(i).delta_f/param.fb*sum( PWF(1:index_f2).*abs(PSXT(1:index_f2)).^2))/2;
                end
            end  % for loop
            if OP.DEBUG && OP.DISPLAY_WINDOW
                figure(300+package_testcase_i);
                if number_of_s4p_files > 1
                    scale=param.a_fext/param.a_thru;% added for crosstalk scaling
                    subplot(3,1,1)
                    hold on
                    plot(faxis_GHz, 20*log10(abs(PSXT*scale)),'r','Disp','PSXTK')
                    icrxi=find(chdata(i).faxis >=param.fb/2,1,'first');
                    subplot(3,1,2)
                    grid on
                    ILtemp=20*log10(abs(chdata(1).sdd21f));
                    IL4ICR=interp1(chdata(1).faxis,ILtemp,chdata(i).faxis);
                    scale=param.a_fext/param.a_thru;% added for crosstalk scaling
                    ICR=-20*log10(abs(PSXT*scale))+IL4ICR;
                    semilogx(faxis_GHz, ICR,'Disp', 'ICR')
                    hold on
                    stem(faxis_GHz(icrxi), ICR(icrxi),'g', 'disp', 'f_{Baud}/2')
                end
                subplot(3,1,1)
                title([param.base ' Losses']); ylabel('dB'); xlabel('GHz')
                grid on;  legend show
                subplot(3,1,2)
                title([param.base ' ICR']); ylabel('dB'); xlabel('GHz')
                ylim([0 80])
                grid on; legend show
                subplot(3,1,3)
                title([param.base ' ILD']); ylabel('dB'); xlabel('GHz')
                ylim([-3 3])
                grid on; legend show
            end
        end % get_FD && DO_ONCE




        if DO_ONCE
            % get impulse responses which in interim step between equation for X(f) and
            % H^(k)(t) without TX FFE or CTLE. These will we added later.
            for i=1:number_of_s4p_files
                if OP.INCLUDE_FILTER % apply RX filter
                    %% Equation 93A-20 %%
                    H_r = 1./polyval([1 2.613126 3.414214 2.613126 1], 1i*chdata(i).faxis./(param.f_r*param.fb));
                    chdata(i).sdd21=chdata(i).sdd21.*H_r;
                end
                % if i==1 and else are the same. just do it
                    if OP.EXTERNAL == false
                        [chdata(i).uneq_imp_response, ...
                            chdata(i).t, ...
                            chdata(i).causality_correction_dB, ...
                            chdata(i).truncation_dB] = s21_to_impulse_DC(chdata(i).sdd21 ,chdata(i).faxis, param.sample_dt, OP) ;
                        chdata(i).uneq_imp_response=chdata(i).uneq_imp_response*chdata(i).A; % adjust IRx for amplitude
                    end
                %------------------------------------------------------------
                % next find Pulse response (SBR) for each channel h^(k)(t)
                if ~OP.DISPLAY_WINDOW && i==1, fprintf('processing COM PDF '); end

                chdata(i).uneq_pulse_response=filter(ones(1, param.samples_per_ui), 1, chdata(i).uneq_imp_response);

                if OP.DEBUG
                    if OP.DISPLAY_WINDOW && ~OP.RX_CALIBRATION
                        figure(150+package_testcase_i);
                        screen_size=get(0,'ScreenSize');
                        pos = get(gcf, 'OuterPosition');
                        set(gcf, 'OuterPosition', ...
                            screen_size([3 4 3 4]).*[1 1 0 0] + pos([3 4 3 4]).*[-1 -1 1 1] ...
                            - (package_testcase_i-1)*[0 20 0 0]);
                        %movegui(gcf,'northeast')

                        set(gcf, 'Name', sprintf('Case %d PR & PDF - %s', package_testcase_i, chdata(i).base));
                        subplot(2,1,1);  hold on; % all plots on the same axes
                        hp=plot(chdata(i).t, chdata(i).uneq_pulse_response,'Disp', chdata(i).base);
                        % hide thru PR in order to show xtalk in a reasonable
                        % scale. thru is shown in another plot.
                        if isequal(chdata(i).type, 'THRU' ) || OP.RX_CALIBRATION
                            set(hp, 'visible', 'off');
                            set(get(get(hp,'Annotation'),'LegendInformation'), 'IconDisplayStyle','off');
                        end
                        title('Pulse response')
                        ylabel('Volts')
                        xlabel('seconds')
                        recolor_plots(gca);
                    else
                        if OP.GET_FD
                            fprintf('%s\tFD rms = %.1f mV\n', chdata(i).base, 1000*chdata(i).sigma_FD);
                        end
                        if param.ndfe~=0
                            fprintf('%s\tUnequalized pulse peak = %.1f mV\n', chdata(i).base, 1000*max(abs(chdata(i).uneq_pulse_response)));
                        end
                    end
                    if OP.ENFORCE_CAUSALITY
                        fprintf('%s\tCausality correction = %.1f dB\n', chdata(i).base, chdata(i).causality_correction_dB);
                    else
                        fprintf('%s\tCausality correction = %.1f dB (not applied)\n', chdata(i).base, chdata(i).causality_correction_dB);
                    end
                    fprintf('%s\tTruncation ratio = %.1f dB\n', chdata(i).base, chdata(i).truncation_dB);
                end
            end
        end %DO_ONCE
        % Determine equalization settings
        fom_result = optimize_fom(OP,param, chdata, sigma_bn);

        chdata(1).IRx=fom_result.IR; % IRx has the TX FFE and CTLE applied
        A_s=abs(fom_result.A_s); % this is the "s" in SNR (PAM4 gain in handled in last sections)
        %% Recommended Delta_y no larger than As/1000 or 0.01 mV
        if OP.force_pdf_bin_size
            delta_y = OP.BinSize;
        else
            delta_y = min(A_s/1000, OP.BinSize);
        end
        % the pdf for PAM4 uses the full swing SBR but assigns voltage for the PDF accordingly
        if OP.RX_CALIBRATION, number_of_s4p_files=1; end
        for i=1:number_of_s4p_files
            % get quick PDF  results but only for THRU when in Rx calibration
            if isequal(chdata(i).type, 'THRU')
                if OP.GET_FD ~=0, FD_rms_sci=chdata(i).sigma_FD; end % decide later whether to use FD or TD rms. pdfr has an rms field.
            end
            if OP.INCLUDE_CTLE==1
                % apply the CTLE that was selected in FOM optimization.
                chdata(i).eq_imp_response = TD_CTLE(chdata(i).uneq_imp_response, param.fb ...
                    , param.CTLE_fz(fom_result.ctle), param.CTLE_fp1(fom_result.ctle) ...
                    , param.CTLE_fp2(fom_result.ctle), param.ctle_gdc_values(fom_result.ctle) ...
                    , param.samples_per_ui);
            else
                chdata(i).eq_imp_response=chdata(i).uneq_imp_response;
            end

            if isequal(chdata(i).type, 'FEXT') || isequal(chdata(i).type, 'THRU')
                % apply the TXFFE that was selected in FOM optimization.
                upsampled_txffe = zeros(1, param.samples_per_ui*2+1); % start with zeros everywhere
                upsampled_txffe(1+(0:2)*param.samples_per_ui) = fom_result.txffe; % plant the coefficients in the desired locations
                chdata(i).eq_imp_response = filter(upsampled_txffe, 1, chdata(i).eq_imp_response);
            end
            chdata(i).eq_pulse_response=filter(ones(1, param.samples_per_ui), 1, chdata(i).eq_imp_response);
            %chdata(i).result=get_SBR(chdata(i).eq_imp_response, param.samples_per_ui);

            if ~OP.DISPLAY_WINDOW, fprintf('%d ', i); end
            pdf = get_pdf(chdata(i), delta_y, fom_result.t_s, param, OP) ;
            if OP.DEBUG && OP.DISPLAY_WINDOW && ~OP.RX_CALIBRATION
                figure(150+package_testcase_i);
                subplot(2,1,2);
                semilogy(pdf.x, pdf.y,'disp', chdata(i).base);
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

        if OP.RX_CALIBRATION
            H_ctf2 = (10.^(param.ctle_gdc_values(fom_result.ctle)/20) + 1i*chdata(2).faxis/param.CTLE_fz(fom_result.ctle)) ./ ...
                ((1+1i*chdata(2).faxis/param.CTLE_fp1(fom_result.ctle)).*(1+1i*chdata(2).faxis/param.CTLE_fp2(fom_result.ctle)));
            sigma_ne = get_sigma_noise( H_ctf2,  param, chdata, sigma_bn );
        else
            sigma_ne=0;
        end
        sigma_N = fom_result.sigma_N;
        %% Equation 93A-30 %%
        sigma_TX = (param.levels-1)*A_s/param.R_LM*10^(-param.SNDR/20);

        %% Equation 93A-41 %%
        sigma_G = norm([param.sigma_RJ*param.sigma_X*norm(fom_result.h_J), sigma_N, sigma_TX]);

        %% Equation 93A-42 %%
        % number of sigmas needed depends on the required BER.
        ber_q = sqrt(2)*erfcinv(2*param.specBER);
        gaussian_noise_pdf = normal_dist(sigma_G, ber_q, delta_y);
        % enable overriding the Q factor of the BBN instrument.
        if OP.force_BBN_Q_factor
            ne_noise_pdf = normal_dist(sigma_ne, OP.BBN_Q_factor, delta_y);
        else
            ne_noise_pdf = normal_dist(sigma_ne, ber_q, delta_y);
        end
        gaussian_noise_pdf = conv_fct(gaussian_noise_pdf, ne_noise_pdf);

        %% p_DD is computed using the procedure defined in 93A.1.7.1 with h(n)=A_DD*h_J(n)
        p_DD = get_pdf_from_sampled_signal(param.A_DD*fom_result.h_J, param.levels, delta_y);

        %% Equation 93A-43
        noise_pdf=conv_fct(gaussian_noise_pdf, p_DD);

        %% Implementation of 93A.1.7.3 combination procedure
        %%  (effectively Equation 93A-44) %%

        % Self-Channel Interference is thru residual result
        sci_pdf = chdata(1).pdfr;
        sci_mxi=find(cumsum(sci_pdf.y)>=param.specBER, 1, 'first');
        thru_peak_interference_at_BER=abs(sci_pdf.x(sci_mxi));
        sci_msi=find(cumsum(sci_pdf.y)>=1e-12, 1, 'first');
        sci_sigma=abs(sci_pdf.x(sci_msi)/(erfcinv(2*1e-12)*sqrt(2)));
        if OP.RX_CALIBRATION ==0
            % Co-Channel Interference PDFs (for information only):
            % initialize to deltas
            MDNEXT_cci_pdf = d_cpdf(delta_y, 0, 1);
            MDFEXT_cci_pdf = d_cpdf(delta_y, 0, 1);
            % serially convolve FEXT/NEXT PDFs
            for k=2:number_of_s4p_files
                if isequal(chdata(k).type, 'NEXT')
                    MDNEXT_cci_pdf = conv_fct(MDNEXT_cci_pdf, chdata(k).pdfr);
                else % ... must be FEXT
                    MDFEXT_cci_pdf = conv_fct(MDFEXT_cci_pdf, chdata(k).pdfr);
                end
            end

            % find "peaks" of MDNEXT/MDFEXT for reporting
            mdnxi=find(cumsum(MDNEXT_cci_pdf.y)>=param.specBER, 1, 'first');
            MDNEXT_peak_interference=abs(MDNEXT_cci_pdf.x(mdnxi));
            mdfxi=find(cumsum(MDFEXT_cci_pdf.y)>=param.specBER, 1, 'first');
            MDFEXT_peak_interference=abs(MDFEXT_cci_pdf.x(mdfxi));

            % Combined crosstalk effect
            cci_pdf = conv_fct(MDFEXT_cci_pdf, MDNEXT_cci_pdf);
            cci_mxi=find(cumsum(cci_pdf.y)>=param.specBER, 1, 'first');
            cci_msi=find(cumsum(cci_pdf.y)>=1e-12, 1, 'first');
            cci_sigma=abs(cci_pdf.x(cci_msi)/(erfcinv(2*1e-12)*sqrt(2)));
            crosstalk_peak_interference_at_BER=abs(cci_pdf.x(cci_mxi));
            % combine cci and sci
            isi_and_xtalk_pdf = conv_fct(sci_pdf, cci_pdf);
        else
            % for calibration there is no cci
            isi_and_xtalk_pdf=sci_pdf;
        end

        mxi=find(cumsum(isi_and_xtalk_pdf.y)>=param.specBER, 1, 'first');
        peak_interference_at_BER=abs(isi_and_xtalk_pdf.x(mxi));


        %% Equation 93A-45
        combined_interference_and_noise_pdf = conv_fct(isi_and_xtalk_pdf, noise_pdf);

        %% Equation 93A-37
        combined_interference_and_noise_cdf=cumsum(combined_interference_and_noise_pdf.y);

        %% The noise and interference amplitude, A_ni, is the magnitude of the value of y0
        %% that satisfies the relationship P(y0) = DER_0
        A_ni_ix=find(combined_interference_and_noise_cdf>param.specBER, 1, 'first');
        A_ni = abs(combined_interference_and_noise_pdf.x(A_ni_ix));

        %% Equation 93A-1
        COM=20*log10(A_s/A_ni);
        min_COM = min(min_COM, COM);

        %% Calculation of error propagation and burst probability
        if OP.nburst>0
            % an error burst of length N will cause each of the first N taps tap to mis-correct and create a PAM (2 or
            % 4) noise term - depending on the N'th previous symbol, with double the tap voltage. From this we calculate
            % the probability of staying in error state, i.e. burst of length N+1.

            % initialize loop with uncorrelated noise and BER
            error_propagation_noise_pdf{1}=combined_interference_and_noise_pdf; %#ok<AGROW> % PDF for burst of length 1 is the uncorrelated PDF

            % Assume an error will occur if the noise excceds the available signal
            % reduced by some dB. reduction is by COM threshold minus Error
            % propagation margin (a positive EP margin reduces uncorrelated error propability
            % below target BER).
            error_threshold =  A_s./10^((param.pass_threshold-OP.COM_EP_margin)/20);
            % Find the probability of this event by integration of the PDF. Use 1e-20 as a floor probabilty if noise PDF isn't wide enough.
            x_error_propagation = find(error_propagation_noise_pdf{1}.x >= error_threshold, 1, 'first');
            if isempty(x_error_propagation)
                p_error_propagation(1) = 1e-20;
            else
                p_error_propagation(1) = sum(error_propagation_noise_pdf{1}.y(x_error_propagation:end));  % uncorrelated BER
            end

            sorted_abs_dfe_taps = sort(abs(result.DFE_taps), 'descend');
            for k=2:min(param.ndfe, OP.nburst)
                % (arrays kept to allow tracking during development, though not really needed)
                if OP.use_simple_EP_model
                    post_error_dfe_noise_pdf{k} = get_pdf_from_sampled_signal( 2*A_s*max(sorted_abs_dfe_taps), param.levels, delta_y ); %#ok<AGROW>
                    error_propagation_noise_pdf{k} = conv_fct(error_propagation_noise_pdf{1}, post_error_dfe_noise_pdf{k}); %#ok<AGROW>
                else
                    post_error_dfe_noise_pdf{k} = get_pdf_from_sampled_signal( 2*A_s*sorted_abs_dfe_taps(k-1), param.levels, delta_y ); %#ok<AGROW>
                    error_propagation_noise_pdf{k} = conv_fct(error_propagation_noise_pdf{k-1}, post_error_dfe_noise_pdf{k}); %#ok<AGROW>
                end

                % Assume an error will propagate if this noise exceeds the threshold defined above
                x_error_propagation = find(error_propagation_noise_pdf{k}.x >= error_threshold, 1, 'first');
                if isempty(x_error_propagation)
                    p_error_propagation(k) = 1e-20; %#ok<AGROW>
                else
                    p_error_propagation(k) = sum(error_propagation_noise_pdf{k}.y(x_error_propagation:end)); %#ok<AGROW>
                end
            end

            % Assume an uncorrelated error will occur if the original noise exceeds
            % the available signal reduced by pass_threhsold dB. Find the probability
            % of this event by partial sum of the PDF.
            % p_uncorrelated_error_i = find(combined_interference_and_noise_pdf.x >= A_s./10^(param.pass_threshold/20), 1, 'first');
            % p_uncorrelated_error = sum(combined_interference_and_noise_pdf.y(p_uncorrelated_error_i:end));

            % probability of bursts of different lengths
            p_burst = cumprod(p_error_propagation);
        end



        %% reporting
        output_args.code_revision=regexp('$Revision: 1.50 $', '1\.(\d+)', 'match', 'once');
        output_args.config_file = config_file;
        fileset_str=str2csv({chdata.base});
        output_args.file_names=sprintf('"%s"', fileset_str);
        for pkg_params = {'levels', 'Pkg_len_TX', 'Pkg_len_NEXT', 'Pkg_len_FEXT', 'Pkg_len_RX'}
            output_args.(pkg_params{:})= param.(pkg_params{:});
        end
        output_args.baud_rate_GHz=param.fb/1e9;
        output_args.f_Nyquist_GHz = param.fb/2e9;

        output_args.channel_operating_margin_dB=COM;
        if OP.RX_CALIBRATION== 1, output_args.sigma_bn=sigma_bn; end
        output_args.peak_interference_mV=1000*A_ni;
        output_args.peak_channel_interference_mV=1000*peak_interference_at_BER;
        output_args.peak_ISI_mV=1000*thru_peak_interference_at_BER;
        if OP.RX_CALIBRATION == 0
            output_args.peak_MDXTK_interference_mV=1000*crosstalk_peak_interference_at_BER;
            output_args.peak_MDNEXT_interference_mV=1000*MDNEXT_peak_interference;
            output_args.peak_MDFEXT_interference_mV=1000*MDFEXT_peak_interference;
        end
        output_args.available_signal_after_eq_mV=1000*A_s;
        output_args.steady_state_voltage_mV = 1000*sum(fom_result.IR);
        output_args.VEO_mV = 1000*(A_s-A_ni)*2;
        output_args.VEO_normalized = (A_s-A_ni)/A_s;
        output_args.VEC_dB = -20*log10(output_args.VEO_normalized);
        if OP.GET_FD
            output_args.fit_loss_dB_at_Fnq = fit_loss;
            output_args.IL_dB_at_Fnq=Nq_loss;
            output_args.ILD_RMS=ILD_RMS;
            output_args.ICN_mV=ICN*1000;
        end

        if OP.DEBUG
            msg = sprintf('Case %d: z_p=(%d, %d, %d, %d) (TX, RX, NEXT, FEXT):', ...
                package_testcase_i, param.Pkg_len_TX, param.Pkg_len_RX, param.Pkg_len_NEXT, param.Pkg_len_FEXT ...
                );

            if OP.DISPLAY_WINDOW && ~OP.RX_CALIBRATION
                % display bathtub curves in one axis per test case.
                figure_name =  'Voltage bathtub curves';
                fig=findobj('Name', figure_name);
                if isempty(fig), fig=figure('Name', figure_name); end
                figure(fig);
                movegui(fig,'south')

                hax = subplot(length(OP.pkg_len_select), 1, package_testcase_i);
                plot_bathtub_curves( hax ...
                    , A_s ...
                    , sci_pdf ...
                    , cci_pdf ...
                    , isi_and_xtalk_pdf ...
                    , noise_pdf ...
                    , combined_interference_and_noise_pdf ...
                    , delta_y ...
                    );
                set(hax, 'tag', 'BTC');
                title(hax, sprintf('case %d VBC: %s ', package_testcase_i, regexprep([chdata(1).base,' '],'_',' ')));
                ylim(hax, [param.specBER/10 1]);
                % show BER target line
                hp=plot(get(hax, 'xlim'), param.specBER*[1 1], 'r:');
                set(get(get(hp,'Annotation'),'LegendInformation'), 'IconDisplayStyle','off');

                h=findall(0, 'Name', 'COM results');
                if isempty(h)
                    msgtext = cell(1, length(OP.pkg_len_select));
                    msgcolor = 'g';
                else
                    msgtext=get(findobj(h, 'type', 'text'), 'string');
                    msgcolor = get(h, 'color');
                    close(h); % will be recreated
                end
                if COM >= param.pass_threshold
                    msgtext{package_testcase_i}=sprintf('%s: COM = %.3f dB (pass)\n', ...
                        msg, COM);
                else
                    msgtext{package_testcase_i}=sprintf('%s: COM = %.3f dB (FAIL)\n', ...
                        msg, COM);
                    msgcolor = 'r';
                end
                close(findall(0, 'tag', 'COM'));
                h=msgbox(msgtext, 'COM results');
                set(h, 'color', msgcolor, 'tag', 'COM');
                movegui(h, 'center');
            else % no windows
                display(['max noise at BER = ' num2str(peak_interference_at_BER)])
                display(['signal after eq = ' num2str(A_s/(param.levels-1))])
                if COM >= param.pass_threshold
                    fprintf('%s <strong> PASS ... COM = %.3f dB</strong>\n', msg, COM);
                else
                    fprintf('%s <strong> FAIL ... COM = %.3f dB</strong>\n',  msg, COM);
                end
            end

        end

        if OP.GET_FD ~=0
            output_args.equivalent_ISI_ICN=sci_sigma;
            output_args.sci_noise_FD_RMS=FD_rms_sci;
        else
            output_args.equivalent_ISI_ICN=0;
            output_args.sci_noise_FD_RMS=0;
        end
        output_args.CTLE_zero_poles=[param.CTLE_fz(fom_result.ctle) param.CTLE_fp2(fom_result.ctle) param.CTLE_fp1(fom_result.ctle)];
        output_args.CTLE_DC_gain_dB=param.ctle_gdc_values(fom_result.ctle);
        output_args.TXLE_taps=fom_result.txffe;
        output_args.DFE_taps=fom_result.DFE_taps;
        if xtk>0 && OP.RX_CALIBRATION ==0
            output_args.cci_noise_TD_BER=crosstalk_peak_interference_at_BER;
            output_args.equivalent_ISI_ICN=cci_sigma;
        else
            output_args.cci_noise_TD_BER=0;
            output_args.equivalent_ISI_ICN=0;
        end
        output_args.peak_interference_at_BER=peak_interference_at_BER;
        output_args.FOM = fom_result.FOM;
        output_args.DFE4_RSS=norm(fom_result.DFE_taps(4:end));
        if OP.nburst>0
            output_args.error_propagation_probability = p_error_propagation;
            output_args.burst_probabilities = p_burst;
        end
        %% making mat file
        if ~OP.EXTERNAL
            if(OP.DEBUG)
                save (sprintf('%s%s_case%d_results.mat', OP.RESULT_DIR, chdata(1).base, package_testcase_i), ...
                    'output_args');
            end
        else %??
            output_args.baud_rate_GHz=param.fb/1e9;
            output_args.sigma_rj=param.sigma_RJ;
            output_args.Add_noise=param.A_DD;
        end

        %% making csv file
        if OP.CSV_REPORT ==1
            items = fieldnames(output_args);
            item_value_strings = cell(size(items));
            for field_id=1:length(items)
                field_name=items{field_id};
                field_value=output_args.(field_name);
                if ischar(field_value)
                    item_value_strings{field_id}=field_value;
                elseif numel(field_value)==1
                    item_value_strings{field_id}=num2str(field_value);
                else
                    item_value_strings{field_id}=sprintf('"%s"', mat2str(field_value));
                end
            end

            header_string = str2csv(items);
            data_string = str2csv(item_value_strings);
            if OP.EXTERNAL== false
                fid = fopen(sprintf('%s%s_case%d_results.csv', OP.RESULT_DIR, chdata(1).base, package_testcase_i),'w');
                fprintf(fid,'%s\n', header_string);
                fprintf(fid,'%s\n', data_string);
                fclose(fid);
            end
        end
        results{package_testcase_i} = output_args;
        if nargout==0
            fprintf('<strong>--- Testcase %d results ---</strong>\n', package_testcase_i);
            disp(output_args)
        end
        % package test cases

        if OP.DEBUG && OP.DISPLAY_WINDOW && OP.RX_CALIBRATION==0
            btc_axes = findobj('tag', 'BTC');
            if ~isempty(btc_axes), linkaxes(btc_axes, 'x'); end
            eqe_axes = findobj('tag', 'EQE');
            if ~isempty(eqe_axes), linkaxes(eqe_axes, 'xy'); end
        end
    end
    if ~OP.DEBUG
        % display bottom line results, minimum COM across test cases
        if OP.DISPLAY_WINDOW
            close(findall(0, 'tag', 'COM'));
            if min_COM >= param.pass_threshold
                h = msgbox(sprintf('PASS ... COM = %.3f dB', min_COM));
                set(h,'Color','g', 'tag', 'COM');
            else
                h = msgbox(sprintf('FAIL ... COM = %.3f dB', min_COM));
                set(h,'Color','r', 'tag', 'COM');
            end
        else
            if min_COM >= param.pass_threshold
                fprintf('<strong> PASS ... COM = %.3f dB</strong>\n', min_COM);
            else
                fprintf('FAIL ... COM = %.3f dB\n', min_COM);
            end
        end
    end
    %%
    if OP.RX_CALIBRATION ==1
        display ([' LOOP with sigma_bn = ' num2str(sigma_bn) ' performed with COM = ' num2str(min_COM) ])
    end
    DO_ONCE=false;
end

if OP.DISPLAY_WINDOW
    savefigs(param, OP);
end

if OP.RX_CALIBRATION==1
    fprintf ('Set Rx calibration noise rms voltage to %g mV\n', sigma_bn*1000);
    if OP.DISPLAY_WINDOW
        message=sprintf('Set Rx calibration noise rms voltage to %g mV.',sigma_bn*1000);
        hlast = msgbox(message,'sigma_bn','help');
        set(hlast,'Color','y', 'tag', 'COM');
        htmp = findobj( hlast, 'Type', 'Text');
        set(htmp,  'FontWeight', 'bold');
%        htmp = findobj( hlast, 'Type', 'uicontrol'); % OK button
%        delete(htmp); % better not do this
    end
end


if length(results)==1, results = results{1}; end
%%
%--------------------------------------------------------------------------
%--------------- Helper functions -----------------------------------------
%--------------------------------------------------------------------------
function [voltage, t_base, causality_correction_dB, truncation_dB] = ...
    s21_to_impulse_DC(IL, freq_array, time_step, OP)
% Creates a time-domain impulse response from frequency-domain IL data.
% IL does not need to have DC but a corresponding frequency array
% (freq_array) is required.
%
% Causality is imposed using the Alternating Projections Method. See also:
% Quatieri and Oppenheim, "Iterative Techniques for Minimum Phase Signal
% Reconstruction from Phase or Magnitude", IEEE Trans. ASSP-29, December
% 1981 (http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=1163714)

ILin=IL;
fmax=1/time_step/2;
freq_step=(freq_array(3)-freq_array(2))/1;
fout=0:1/round(fmax/freq_step)*fmax:fmax;
IL=interp_Sparam(ILin,freq_array,fout, ...
    OP.interp_sparam_mag, OP.interp_sparam_phase);
IL_nan = find(isnan(IL));
for in=IL_nan
    IL(in)=IL(in-1);
end

IL = IL(:);
% add padding for time steps
IL_symmetric = [IL(1:end-1);0; flipud(conj(IL(2:end-1)))];
impulse_response = real(ifft(IL_symmetric));
L = length(impulse_response);
t_base = (0:L-1)/(freq_step*L);

original_impulse_response=impulse_response;
% Correct non-causal effects frequently caused by extrapolation of IL
% Assumption: peak of impulse_response is in the first half, i.e. not anti-causal
abs_ir=abs(impulse_response);
a = find(abs_ir(1:L/2) > max(abs_ir(1:L/2))*OP.EC_PULSE_TOL);
start_ind = a(1);

err=inf;
while ~all(impulse_response==0)
    impulse_response(1:start_ind)=0;
    impulse_response(floor(L/2):end)=0;
    IL_modified=abs(IL_symmetric).*exp(1j*angle(fft(impulse_response)));
    ir_modified = real(ifft(IL_modified));
    delta = abs(impulse_response-ir_modified);

    err_prev = err;
    err=max(delta)/max(impulse_response);
    if err<OP.EC_REL_TOL || abs(err_prev-err)<OP.EC_DIFF_TOL
        break;
    end

    impulse_response=ir_modified;
end

causality_correction_dB=20*log10(norm(impulse_response-original_impulse_response)/norm(impulse_response));

if ~OP.ENFORCE_CAUSALITY
    impulse_response = original_impulse_response;
end
% truncate final samples smaller than 1e-3 of the peak
ir_peak = max(abs(impulse_response));
ir_last  = find(abs(impulse_response)>ir_peak*OP.impulse_response_truncation_threshold, 1, 'last');

voltage = impulse_response(1:ir_last);
t_base = t_base(1:ir_last);

truncation_dB=20*log10(norm(impulse_response(ir_last+1:end))/norm(voltage));


function [Sout] = interp_Sparam(Sin,fin,fout, ...
    opt_interp_Sparam_mag, opt_interp_Sparam_phase)
% Sout = interp_Sparam(Sin,fin,fout)
%
% Interpolate S-parameters Sin from frequency grid fin to frequency grid
% fout.

if ( fin(end)<fout(end) )
    %    warning('Channel high frequencies extrapolation might be inaccurate!');
end

H_mag = abs(Sin);
H_mag(H_mag<eps)=eps; % handle ill cases...
H_ph = unwrap(angle(Sin));
% For long delay channels, the result can turn anti-causal if frequency step is too coarse. Don't let the
% user ignore that.
if mean(diff(H_ph))>0
    error('Anti-causal response found. Finer frequency step is required for this channel');
end

%opt_interp_Sparam_mag='linear_trend_to_DC';
switch opt_interp_Sparam_mag
    case 'linear_trend_to_DC'
        fin_x=fin;
        H_mag_x=H_mag(:);
        if fin(1)>0
            p=polyfit(fin(1:10), H_mag(1:10), 1);
            dc_trend_val=polyval(p, 0);
            fin_x=[0, fin_x];
            H_mag_x = [dc_trend_val; H_mag_x];
        end
        if fin(end)<fout(end)
            warn_state=warning('off', 'MATLAB:polyfit:RepeatedPointsOrRescale');
            mid_freq_ind=round(length(fin)/2);
            p=polyfit(fin(mid_freq_ind:end), H_mag(mid_freq_ind:end), 1);
            warning(warn_state);
            hf_trend_val=polyval(p, fout(end));
            if hf_trend_val>H_mag(end)
                hf_trend_val=H_mag(end);
            elseif hf_trend_val<eps
                hf_trend_val=eps;
            end
            fin_x=[fin_x, fout(end)];
            H_mag_x = [H_mag_x; hf_trend_val];
        end
        H_mag_i = interp1(fin_x, H_mag_x, fout, 'linear', 0);
    case 'trend_to_DC'
        % extrapolate to trend value at DC.
        fin_x=fin;
        H_mag_x=H_mag;
        if fin(1)>0
            p=polyfit(fin(1:10)', log10(H_mag(1:10)), 1);
            dc_trend_val=10^polyval(p, 0);
            fin_x=[0, fin_x];
            H_mag_x = [dc_trend_val; H_mag_x];
        end
        if fin(end)<fout(end)
            warn_state=warning('off', 'MATLAB:polyfit:RepeatedPointsOrRescale');
            mid_freq_ind=round(length(fin)/2);
            p=polyfit(fin(mid_freq_ind:end)', log10(H_mag(mid_freq_ind:end)), 1);
            warning(warn_state);
            hf_trend_val=10^polyval(p, fout(end));
            if hf_trend_val>H_mag(end)
                hf_trend_val=H_mag(end);
            end
            fin_x=[fin_x, fout(end)];
            H_mag_x = [H_mag_x; hf_trend_val];
        end
        H_mag_i = 10.^interp1(fin_x,log10(H_mag_x),fout,'linear', 'extrap');
    case 'extrap_to_DC_or_zero'
        % same as extrap_to_DC but detect AC-coupled channels and
        % extrapolate them to 0.
        if fin(1)>0 && 20*log10(H_mag(1))<-20
            % assume AC coupling, with 0 at DC
            H_mag_i = 10.^interp1([0, fin],[-100; log10(H_mag)],fout(fout<=fin(end)),'linear', 'extrap');
        else
            H_mag_i = 10.^interp1(fin, log10(H_mag), fout(fout<=fin(end)),'linear', 'extrap');
        end
        H_mag_i(fout>fin(end)) = H_mag(end);
    case 'extrap_to_DC'
        % first extrapolate down to DC, then use highest available frequency
        % for higher frequencies
        H_mag_i = 10.^interp1(fin,log10(H_mag),fout(fout<=fin(end)),'linear', 'extrap');
        H_mag_i(fout>fin(end)) = H_mag(end);
    case 'old'
        H_mag_i = interp1(fin,H_mag,fout,'linear','extrap');
    otherwise
        error('COM:Extrap:InvalidOption', 'opt_interp_Sparam_mag valid values are "old", "extrap_to_DC"');
end

H_ph_i = interp1(fin,H_ph,fout,'linear', 0);

%opt_interp_Sparam_phase='trend_and_shift_to_DC';
%opt_interp_Sparam_phase='interp_cubic_to_dc_linear_to_inf';
switch opt_interp_Sparam_phase
    case 'old'
        H_ph_i = H_ph_i-H_ph_i(1);
    case 'zero_DC'
        H_ph_i(1) = 0;
    case 'interp_to_DC'
        if fin(1) ~= 0
            H_ph_i = interp1([0; fin(:)], [0; H_ph(:)], fout, 'linear', 'extrap');
        end
    case 'extrap_cubic_to_dc_linear_to_inf'
        if fin(1) ~= 0
            % estimate low frequency group delay
            group_delay = -diff(H_ph(:))./diff(fin(:));
            low_freq_gd = group_delay(1:50);
            %  calculate trend, throwing away outliers
            m = median(low_freq_gd); sigma = std(low_freq_gd);
            lf_trend = mean(low_freq_gd(abs(low_freq_gd-m)<sigma));
            % correct outliers in first 10 phase samples
            for k=10:-1:1
                H_ph(k) = H_ph(k+1) + lf_trend*(fin(k+1)-fin(k));
            end
            H_ph_cubic = interp1(fin, H_ph, fout, 'cubic', 'extrap');
            H_ph_linear = interp1(fin, H_ph, fout, 'linear', 'extrap');
            % modification - trend to inf
            if (1)
                high_freq_gd = group_delay(end-50:end);
                %  calculate trend, throwing away outliers
                m = median(high_freq_gd); sigma = std(high_freq_gd);
                hf_trend = -mean(high_freq_gd(abs(high_freq_gd-m)<sigma));
                hf_extrap_range = find(fout>fin(end));
                last_data_sample = hf_extrap_range(1)-1;
                H_ph_linear(hf_extrap_range) = H_ph_linear(last_data_sample) + (fout(hf_extrap_range)-fout(last_data_sample))*hf_trend;
%                for k=hf_range
%                    H_ph_linear(k) = H_ph_linear(k-1) + hf_trend*(fout(k)-fout(k-1));
%                end
            end

            [UNUSED_OUTPUT, indx] = min(abs(H_ph_cubic-H_ph_linear)); %#ok<ASGLU>
            H_ph_i=H_ph_cubic;
            H_ph_i(indx:end) = H_ph_linear(indx:end);
        end
    case 'interp_and_shift_to_DC'
        if fin(1) ~= 0
            dc_phase_trend = H_ph(1)-(H_ph(2)-H_ph(1))/(fin(2)-fin(1))*fin(1);
            H_ph_i = interp1([0; fin(:)], [0; H_ph(:)-dc_phase_trend], fout, 'linear', 'extrap');
        end
    case 'trend_and_shift_to_DC'
        % estimate low frequency group delay
        group_delay = -diff(H_ph(:))./diff(fin(:));
        low_freq_gd = group_delay(1:50);
        %  calculate trend, throwing away outliers
        m = median(low_freq_gd); sigma = std(low_freq_gd);
        lf_trend = mean(low_freq_gd(abs(low_freq_gd-m)<sigma));
        fin_x=fin;
        H_ph_x=H_ph(:);
        if fin(1) ~= 0
            % correct outliers in first 10 phase samples
            for k=10:-1:1
                H_ph(k) = H_ph(k+1) + lf_trend*(fin(k+1)-fin(k));
            end

            % shift all phase data so that DC extrapolation to 0 follows trend
            dc_phase_trend = H_ph(1)+lf_trend*(fin(1)-0);
            fin_x=[0, fin_x];
            H_ph_x=[0; H_ph(:)-dc_phase_trend];
        end
        % Modification: extrapolate using trend. (interp1 with "extrap" extrapolates using just
        % the last two samples, so noise can create an inverted slope and
        % non-causal response).
        if fout(end)>fin(end)
            group_delay = -diff(H_ph_x(:))./diff(fin_x(:));
            % p=polyfit(fin_x', H_ph_x, 1);
            hf_phase_trend = H_ph_x(end)-median(group_delay)*(max(fout)-max(fin_x));
            % hf_phase_trend=polyval(p,max(fout));
            fin_x=[fin_x, fout(end)];
            H_ph_x=[H_ph_x; hf_phase_trend];
        end
        H_ph_i = interp1(fin_x, H_ph_x, fout, 'linear', 'extrap');

    otherwise
        error('COM:Extrap:InvalidOption', ...
            'debug_interp_Sparam valid values are "old", "zero_DC", "interp_to_DC", "interp_and_shift_to_DC", "trend_and_shift_to_DC", "interp_cubic_to_dc_linear_to_inf"');
end
H_i = H_mag_i.*exp(1j*H_ph_i);
Sout=H_i;

function [data, SDD, SDC] = read_p4_s4params(infile, plot_ini_s_params, plot_dif_s_params, ports,OP)
%% FUNCTION :: read_sp4_sparams
%
% Description
%   Read the fid of single-ended 4-port complex S-parameters
%   in Touchstone format 'file' and convert to the internal
%   format using the port transform 'ports'
%
%   Created by Mike Y. He
%   April 22, 2005
%
%   Reused some code from
%   Anthony Sanders, Alex Deas, Bob Davidov (24 January 2005)
%   for touchstone 4-port S-matrix import.
%
%   Modified (2012-July-27) by Ken Young to match current indexing scheme and
%   optimized for quicker parameter matching and parsing. also, separated out
%   the plotting algorithms into their own sub-function routines
%
% Input Variables (required)
%   file                -- The s4p file to be read and converted
%   plot_ini_sparams    -- Plot the initial s-parameter information. For debugging purposes
%   plot_diff_sparams   -- Plot the differential s-parameter information. For debugging purposes
%   swap                -- Re-order the port layout
%
% Output/Return Variables
%   sch         -- the data matrix contain the network parameter data points
%   schFreq     -- the frequency vector of the network parameter data points
%   sdc         -- the differential in/common-mode out s-parameter data matrix
%   sdd         -- the differential in/differential out s-parameter data matrix
%
% Local Variables (scoped locally to current function only)
%   fid                 -- the opened file being processed
%   fileDataLine        -- the current line in the fid
%   optionLine          -- the options line in the s4p file
%   ports               -- the order of the ports being processed
%   idxNetParams        -- the position (line) of the network parameter data points in the s4p file
%   netParamDataPoints  -- the network parameter data points of the current line of the fid
%   schFreq             -- the frequency points in the s4p file stored as a separate vector
%   sch                 -- the network parameter data points from the s4p file converted in matrix form
%   S                   --
%   T                   --
%   W                   --
%   D                   --
%   sdd                 -- the differential in/differential out form of network parameter data points
%   sdc                 -- the differential in/common-mode out form of network parameter data points
%   plot_ini_sparams    -- Plot the initial s-parameter information. For debugging purposes
%   plot_diff_sparams   -- Plot the differential s-parameter information. For debugging purposes
%

% open the fid
fid = fopen(infile, 'r');

data_file = textscan(fid, '%s', 'delimiter', sprintf('\n'));
data_file = data_file{1,1}; % removes a dimension layer of cell array depth

fclose(fid);

% backwards compatibility settings. can be removed in updated code.
if ~exist('OP', 'var'); OP.DISPLAY_WINDOW = true; end
if isempty(ports); ports = [1 3 2 4]; end % default order normally used.

% adjust ports to maintain the meaning [in1, in2 , out1, out2] when one
% pair is reversed.
ports_adj=ports; for k=1:4, ports_adj(k)=find(ports==k); end; ports=ports_adj;

if OP.DISPLAY_WINDOW
    set(0,'defaulttextinterpreter','none') % prevents subscripting character in displayed messages
    hMsgBox = msgbox(infile, 'Reading S-Parameter File'); % display a progress bar for reading the s-parameter file(s)
end

% parse the file array until the s-parameters options line is found
idxLine = 1;
while true
    if regexp(data_file{idxLine}, '^#', 'once')
        % tests the option fileDataLine cell array to see if it contains...
        if ~isempty(regexpi(data_file{idxLine}, '\sGHz')); freq_scale = 1e9; end
        if ~isempty(regexpi(data_file{idxLine}, '\sMHz')); freq_scale = 1e6; end
        if ~isempty(regexpi(data_file{idxLine}, '\sKHz')); freq_scale = 1e3; end
        if ~isempty(regexpi(data_file{idxLine}, '\sHz'));  freq_scale = 1e0; end
        if ~isempty(regexpi(data_file{idxLine}, 'RI')) && ~isempty(regexpi(data_file{idxLine}, 'R'))
            format = 'complex';
        end
        if ~isempty(regexpi(data_file{idxLine}, 'MA')) && ~isempty(regexpi(data_file{idxLine}, 'R'))
            format = 'linpolar';
        end
        if ~isempty(regexpi(data_file{idxLine}, 'DB')) && ~isempty(regexpi(data_file{idxLine}, 'R'))
            format = 'dBpolar';
        end

        % if the format or freq_scale are not set then exit with error
        if isempty(freq_scale)
            error( 'read_s4p:MissingOptionsParameter', ['\n\t THE FREQUENCY SCALE PARAMETER INFORMATION IS MISSING' ...
                'FROM THE OPTIONS LINE IN THE S4P FILE']);
        elseif isempty(format)
            error( 'read_s4p:MissingOptionsParameter', ['\n\t THE FORMAT PARAMETER INFORMATION IS MISSING' ...
                'FROM THE OPTIONS LINE IN THE S4P FILE']);
        end

        idxLine = idxLine+1; %step one more line before exiting loop
        break
    end
    idxLine = idxLine+1;
end

col_1=1; col_2=2; col_3=3; col_4=4; % This is to maintain current indexing compatibility as well as readability
schFreqAxis = zeros(1,int8((length(data_file)-idxLine)/4)); % preallocated array length for memory and time purposes
idxNetParams = 1; % each 4-port network parameter data point set will receive an index

freqCounter = 0; % variable that will be used to capture number of freq points.

while true
    if isempty(data_file{idxLine}) || ~isempty(regexp(data_file{idxLine}, '^!', 'once'))
        idxLine = idxLine+1;
    else
        for row=1:4
            dataPointsLine = eval (['[' data_file{idxLine} ']']);
            if mod(length(dataPointsLine),8)==1
                schFreqAxis(idxNetParams) = dataPointsLine(1) * freq_scale;
                dataPointsLine(1) = []; % removes the frequency value from the current data fileDataLine cell array
                % so that only the network parameter data is in the array
            end
            if mod(length(dataPointsLine),8)~=0  % sanity check to make sure the network parameter data point line structure is correct
                msg = sprintf(['\n'...
                    '\n\t ******************************************************' ...
                    '\n\t ** There is an error in the s4p network parameter   **' ...
                    '\n\t ** data point information structure at line %G.'        ...
                    '\n\t **              ENDING PROGRAM                      **' ...
                    '\n\t ******************************************************' ...
                    '\n\n\n\t Press any key to exit (an error will be displayed). \n\n']);
                disp(msg)
                pause on;
                pause;
                pause off;
                error( 'COM:reads_4p:SparamFileDataError', 'INCORRECT S4P NETWORK PARAMETER DATA POINT STRUCTURE');
            else
                switch format
                    case 'complex'
                        sch(idxNetParams, ports(row), ports(col_1)) = dataPointsLine(1) + 1i*dataPointsLine(2); %#ok<AGROW>
                        sch(idxNetParams, ports(row), ports(col_2)) = dataPointsLine(3) + 1i*dataPointsLine(4); %#ok<AGROW>
                        sch(idxNetParams, ports(row), ports(col_3)) = dataPointsLine(5) + 1i*dataPointsLine(6); %#ok<AGROW>
                        sch(idxNetParams, ports(row), ports(col_4)) = dataPointsLine(7) + 1i*dataPointsLine(8); %#ok<AGROW>
                    case 'linpolar'
                        sch(idxNetParams, ports(row), ports(col_1)) = dataPointsLine(1) * exp(1i*dataPointsLine(2)*(pi/180)); %#ok<AGROW>
                        sch(idxNetParams, ports(row), ports(col_2)) = dataPointsLine(3) * exp(1i*dataPointsLine(4)*(pi/180)); %#ok<AGROW>
                        sch(idxNetParams, ports(row), ports(col_3)) = dataPointsLine(5) * exp(1i*dataPointsLine(6)*(pi/180)); %#ok<AGROW>
                        sch(idxNetParams, ports(row), ports(col_4)) = dataPointsLine(7) * exp(1i*dataPointsLine(8)*(pi/180)); %#ok<AGROW>
                    case 'dBpolar'
                        sch(idxNetParams, ports(row), ports(col_1)) = 10^(dataPointsLine(1)/20) * exp(1i*dataPointsLine(2)*(pi/180)); %#ok<AGROW>
                        sch(idxNetParams, ports(row), ports(col_2)) = 10^(dataPointsLine(3)/20) * exp(1i*dataPointsLine(4)*(pi/180)); %#ok<AGROW>
                        sch(idxNetParams, ports(row), ports(col_3)) = 10^(dataPointsLine(5)/20) * exp(1i*dataPointsLine(6)*(pi/180)); %#ok<AGROW>
                        sch(idxNetParams, ports(row), ports(col_4)) = 10^(dataPointsLine(7)/20) * exp(1i*dataPointsLine(8)*(pi/180)); %#ok<AGROW>
                    otherwise
                end
            end
            idxLine = idxLine+1;
        end

        freqCounter = freqCounter+1; % capturing number of freq points for preallocation.
        idxNetParams = idxNetParams + 1;
    end
    if idxLine > length(data_file)
        break
    end
end

D=NaN(size(sch));
% calculate differential s parameter matrix from single ended
for i=1:size(sch,1)
    S(:,:) = sch(i,:,:);
    T = [1 1 0 0 ; 1 -1 0 0 ; 0 0 1 1 ; 0 0 1 -1];
    W = T * (S / T);
    D(i,:,:) = W(:,:);
end

% D matrix should be
% Scc11 Scd11 Scc12 Scd21
% Sdc11 Sdd11 Sdc12 Sdd12
% Scc21 Scd21 Scc22 Scd22
% Sdc21 Sdd21 Sdc22 Sdd22

% proper values
SDD(:,1,1) = D(:,2,2);
SDD(:,2,2) = D(:,4,4);
SDD(:,1,2) = D(:,2,4);
SDD(:,2,1) = D(:,4,2);

SDC(:,1,1) = D(:,2,1);
SDC(:,2,2) = D(:,4,3);
SDC(:,1,2) = D(:,2,3);
SDC(:,2,1) = D(:,4,1);

% backwards compatibility output variables
data.m = sch;
schFreqAxis=schFreqAxis(1:freqCounter); % truncating preallocated array to number of freq points.
data.freq = schFreqAxis;
colors = 'rgbk';

if (plot_ini_s_params == 1)
    figure('name', 'Single-ended s-parameters');
    for mj=1:4
        subplot(2,2,mj);
        for mi=1:4
            plot(data.freq, 20*log10(abs(data.m(:,mj,mi)+1.0e-15)), ...
                colors(mi), 'linewidth', 2, 'disp', sprintf('S%d%d', mj, mi));
            hold on
        end
        xlabel('Frequency (Hz)');
        ylabel('Magnitude (dB)');
        legend show
        grid on
        title(sprintf('Output port %d', mj));
    end
end

if (plot_dif_s_params == 1)
    figure('name', 'Mixed-mode s-parameters');
    subplot(2,1,1);
    for mj=1:2
        for mi=1:2
            plot(data.freq, 20*log10(abs(SDD(:,mj,mi))), ...
                colors((mj-1)*2+mi), 'linewidth',2, 'disp', sprintf('SDD%d%d', mj, mi));
            hold on
        end
    end
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB)');
    legend show
    grid on
    title(infile);

    subplot(2,1,2);
    for mj=1:2
        for mi=1:2
            plot(data.freq, 20*log10(abs(SDC(:,mj,mi))+1.0e-15), ...
                colors((mj-1)*2+mi), 'linewidth',2, 'disp', sprintf('SDC%d%d', mj, mi));
            hold on
        end
    end
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB)');
    legend show
    grid on
end

if OP.DISPLAY_WINDOW, close(hMsgBox); end
% end read_sp4_sparam

function result = readdataSnPx(filename, nport)
%function [freq, cs] = readdataSnPx(filename, nport)
% [freq, cs] = readdataSnP(filename, nport, format, nheader)
%
% Read Touchstone file with frequencies in units of Hertz
%
% Input:
% ======
% filename: Name of the Touchstone/SnP file
% nport: Number of ports
% format: 'RI' for real/imag, 'MA' for mag/angle (check option line in the
%         Touchstone file)
% nheader: Number of header lines (comment lines plus option line in the
%          Touchstone file)
%
% Output:
% =======
% freq: Vector of frequencies [Hz]
% cs: 3-D array of complex-valued S parameters where cs(i,j,k) is S(i,j)
%     at frequency freq(k)
%
% Note: If frequency unit is not Hertz (but GHz, MHz etc.) simply scale
% frequencies appropriately after reading the data.
%
% Ref.: Touchstone(R) File Format Specification, Rev.1.1,
% EIA/IBIS Open Forum, 2002.
%
% Written by Henning Braunisch, September 2004.
% Updated by Steven Krooswyk, April 2006.


fid = fopen(filename, 'r');


% Skip header lines
str = ' ';
n = 0;
while ~strcmp(str(1),'#')
    str = fgetl(fid);
    if isempty(str)
        str=' ' ;
        if n > 1000
            display('error: could not find config line (#)')
            break
        end
    end
    n = n + 1;
end

% parse configuration line
A=sscanf(str,'%1s %2s %1s %2s %1s %2s',[1,inf]);
p = find(A=='S');           %position of 'S'
units = lower(A(2:p-1));    %units before 'S'
format = A(p+1:p+2);        %format after 'S'

% skip any more header lines
%while ~str

nk = 0; % frequency counter
while 1

    [temp, count] = fscanf(fid, '%f', 1);
    if count == 0
        temp2 = fscanf(fid, '%s', 1);
        if ~isempty(temp2), fgetl(fid); continue, end;
        break
    end
    nk = nk+1; freq(1,nk) = temp;  %#ok<AGROW>
    for ni = 1:nport
        for nj = 1:nport
            switch lower(format)
                case 'ma'
                    mag = fscanf(fid, '%f', 1); ang = fscanf(fid, '%f', 1);
                    cs(ni,nj,nk) = mag * exp(1i*ang*pi/180); %#ok<AGROW>
                case 'ri'
                    re = fscanf(fid, '%f', 1); im = fscanf(fid, '%f', 1);
                    cs(ni,nj,nk) = complex(re, im); %#ok<AGROW>
                case 'db'
                    db = fscanf(fid, '%f', 1); ang = fscanf(fid, '%f', 1);
                    M = 10^(db/20);
                    %re = M*cos(ang);
                    %im = M*sin(ang);
                    re =  M*cos(ang * pi / 180);
                    im =  M*sin(ang * pi / 180);
                    cs(ni,nj,nk) = complex(re, im); %#ok<AGROW>
                otherwise
                    error('readdataSnP: Unknown data format');
            end
        end
    end
end

fclose(fid);

% If 2-port then swap S_12 and S_21 per Touchstone spec
if nport == 2
    temp = cs(2,1,:);
    cs(2,1,:) = cs(1,2,:);
    cs(1,2,:) = temp;
end

% Update freq units to Hz
switch lower(units)
    case 'hz'

    case 'khz'
        freq=freq.*1e3;
    case 'mhz'
        freq=freq.*1e6;
    case 'ghz'
        freq=freq.*1e9;
end

% passivity check
result.freq = freq;
result.cs    = cs;

function result = reduce(var1)
% --- Reduce 1x1xn array to 1xn (aka squeeze)
out = zeros(1,length(var1));
out(1,:) = var1(1,1,:);
result=out;

function pdf=get_pdf(chdata, delta_y, t_s, param, OP)
SBR=chdata.eq_pulse_response(:)'; % row vector
type=chdata.type;
samp_UI=param.samples_per_ui;
residual_response = SBR;

if isequal(type, 'THRU')
    % for thru pulse response:
    % remove the cursor and the DFE postcursors (up to their limit), since
    % we only care about the residuals.
    ideal_cancelled_cursors = SBR(t_s+param.samples_per_ui*(0:param.ndfe));
    effective_cancelled_cursors = sign(ideal_cancelled_cursors) .* ...
        min(abs(ideal_cancelled_cursors), residual_response(t_s)*[1, param.bmax]);
    effective_cancellation_samples = kron(effective_cancelled_cursors, ones(1, param.samples_per_ui));

    % Apply a constant DFE coefficient 1/2 UI before and after each postcursor. Not
    % really needed for COM, but helps debugging. May be factored out in future revisions.
    start_cancel = t_s-param.samples_per_ui/2;
    end_cancel = t_s+(1/2+param.ndfe)*param.samples_per_ui - 1;
    residual_response(start_cancel:end_cancel) = ...
        residual_response(start_cancel:end_cancel) - effective_cancellation_samples;
    %else
    % for crosstalk pulse responses, nothing is cancelled, and all phases
    % are equally important.
end

nui=round(length(residual_response)/param.samples_per_ui);

vs=zeros(nui-2, param.samples_per_ui);
for i=1:param.samples_per_ui
    vs(:,i)=residual_response(param.samples_per_ui*(1:nui-2)+i);
end

if OP.DISPLAY_WINDOW,
    hwaitbar=waitbar(0);
end

% determine which pdf to use
if isequal(type, 'THRU')
    % one phase is interesting for thru
    phases = mod(t_s,param.samples_per_ui);
    if phases==0, phases = param.samples_per_ui; end
else
    phases=1:samp_UI;
end

mxV = zeros(size(phases));
for k=phases
    pdf_samples(k)=get_pdf_from_sampled_signal(vs(:,k), param.levels, delta_y); %#ok<AGROW>
    mxV(k)=sqrt(sum( pdf_samples(k).x.^2.*pdf_samples(k).y)); % standard deviation of PDF
    progress = k/length(phases);
    if OP.DISPLAY_WINDOW, waitbar(progress, hwaitbar, ['processing COM pdf ' chdata.base ] ); figure(hwaitbar); drawnow; end
end
[UNUSED_OUTPUT pxi]=max(mxV); %#ok<ASGLU>
pdf=pdf_samples(pxi);

if OP.DISPLAY_WINDOW,
    close(hwaitbar);
end

function [ILN, FIT]= get_ILN(sdd21,faxis_f2)
% used for FD IL fitting
% sdd21 us a complex insertion loss
fmbg=[ones(length(faxis_f2),1).*transpose(sdd21)  transpose(sqrt(faxis_f2)).*transpose(sdd21)  transpose(faxis_f2).*transpose(sdd21) transpose(faxis_f2.^2).*transpose(sdd21) ];
warning('off','MATLAB:nearlySingularMatrix');
unwraplog=log(abs(sdd21))+1i*unwrap(angle(sdd21));
LGw=transpose(sdd21.*unwraplog);
alpha = ((fmbg'*fmbg)^-1)*fmbg'*LGw;
efit=(alpha(1)+alpha(2).*sqrt(faxis_f2)+alpha(3).*faxis_f2 +faxis_f2.^2.*alpha(4)   );
FIT=transpose(exp(transpose(efit)));
ILN = sdd21-FIT;

function result=optimize_fom(OP, param, chdata, sigma_bn)
%% input
% chdata(1).uneq_imp_response is the impulse response input expected to be normalized to At, peak drive voltage
% baud_rate - baud rate in seconds
% param.samples_per_ui = samples per UI of IR
% param.max_ctle - maximum ac to dc gain in dB
% param.tx_ffe(1) - maximum pre cursor (positive value)
% param.tx_ffe(2) - maximum post cursor (positive value)
% param.tx_ffe_step - sweep step size for tx pre and post taps
% param.ndfe - number of reference dfe taps
% output
% result.eq.txle - [ precusor curosr postcursor]: pre and post are negative
% result.eq.ctle - index of CTLE parameters in table
% result.IR - impulse response
% result.avail_signal - maximum signal after equalization
% result.avail_sig_index - index in result.IR of max signal
% result.best_FOM - best raw ISI

min_number_of_UI_in_response=40;
baud_rate=1/param.ui;

H_r = 1./polyval([1 2.613126 3.414214 2.613126 1], 1i*chdata(1).faxis/(param.f_r*param.fb));

cm1_values = param.tx_ffe_cm1_values;
cp1_values = param.tx_ffe_cp1_values;
gdc_values = param.ctle_gdc_values;
best_ctle = [];
best_FOM = -inf;
best_txffe = [];
delta_sbr = [];
pxi=0;
if OP.DISPLAY_WINDOW
    hwaitbar=waitbar(0);
else
    fprintf('FOM search');
end
for ctle_index=1:length(gdc_values)
    g_dc = gdc_values(ctle_index);
    kacdc = 10^(g_dc/20);
    CTLE_fp1 = param.CTLE_fp1(ctle_index);
    CTLE_fp2 = param.CTLE_fp2(ctle_index);
    CTLE_fz = param.CTLE_fz(ctle_index);

    for k=1:param.num_s4p_files
        chdata(k).ctle_imp_response = TD_CTLE(chdata(k).uneq_imp_response, baud_rate ...
            , CTLE_fz, CTLE_fp1, CTLE_fp2, g_dc, param.samples_per_ui);
    end

    %% Equation 93A-22 %%
    H_ctf = (kacdc + 1i*chdata(1).faxis/CTLE_fz) ./ ...
        ((1+1i*chdata(1).faxis/CTLE_fp1).*(1+1i*chdata(1).faxis/CTLE_fp2));
    if OP.RX_CALIBRATION
        H_ctf2 = (kacdc + 1i*chdata(2).faxis/CTLE_fz) ./ ...
            ((1+1i*chdata(2).faxis/CTLE_fp1).*(1+1i*chdata(2).faxis/CTLE_fp2));
    end
    %% Equation 93A-35 - independent of FFE setting %%
    sigma_N = sqrt(param.eta_0*sum(abs(H_r(2:end) .* H_ctf(2:end)).^2 .* diff(chdata(1).faxis)/1e9));
    if OP.RX_CALIBRATION
        sigma_ne = get_sigma_noise( H_ctf2, param, chdata, sigma_bn); %% Equation 93A-48 %%
        sigma_NEXT=0;
    else
        %% Equations 93A-33 and 93A-34  for NEXT - independent of TXFFE setting %%
        sigma_NEXT =  get_xtlk_noise( [0 1 0], 'NEXT', param, chdata );
        sigma_ne=0;
    end

    for k_cm1=1:length(cm1_values)
        cm1=cm1_values(k_cm1);
        for k_cp1=1:length(cp1_values)
            cp1=cp1_values(k_cp1);
            pxi=pxi+1;
            progress = pxi/( length(cm1_values)*length(cp1_values)*length(gdc_values) );
            txffe = [cm1, 1-abs(cm1)-abs(cp1), cp1];
            % Skip combinations with small values of c(0), not guaranteed to be supported by all transmitters.
            if txffe(2)<param.tx_ffe_c0_min
                continue;
            end
            upsampled_txffe = zeros(1, param.samples_per_ui*2+1); % start with zeros everywhere
            upsampled_txffe(1+(0:2)*param.samples_per_ui) = txffe; % "plant" the coefficients in the desired locations
            effective_channel = filter(upsampled_txffe, 1, chdata(1).ctle_imp_response);
            sbr=filter(ones(param.samples_per_ui, 1), 1, effective_channel);
            [UNUSED_OUTPUT, sbr_peak_i]=max(abs(sbr)); %#ok<ASGLU>
            triple_transit_time = round(sbr_peak_i*2/param.samples_per_ui)+20;
            if min_number_of_UI_in_response < triple_transit_time
                min_number_of_UI_in_response = triple_transit_time;
            end

            %% initial guess at cursor location (t_s)  - based on approximate zero crossing
            zxi = find(diff(sign(sbr-.01*max(sbr)))>=1);
            zxi = zxi(zxi<sbr_peak_i);
            zxi = zxi(sbr_peak_i - zxi < 4*param.samples_per_ui);
            if isempty(zxi)
                continue;
            elseif length(zxi)>1
                zxi=zxi(end);
            end

            if param.ndfe==0
                max_dfe1=0;
            else
                max_dfe1=param.bmax(1);
            end
            %% adjust cursor_i to Solve equation 93A-25 %%
            % Muller-Mueller criterion with DFE
            mm_range = zxi+(0:2*param.samples_per_ui);
            mm_metric = ...
                abs(sbr(mm_range-param.samples_per_ui) - max(sbr(mm_range+param.samples_per_ui)-max_dfe1*sbr(mm_range), 0));
            [UNUSED_OUTPUT, mm_cursor_offset] = min(mm_metric); %#ok<ASGLU>
            cursor_i = zxi+mm_cursor_offset-1;

            cursor = sbr(cursor_i);
            %% 93A.1.6 step c defines A_s %%
            A_s = param.R_LM*cursor/(param.levels-1);
            if isempty(delta_sbr)
                delta_sbr = sbr;
            end
            sbr=sbr(:);

            %% Equation 93A-27 "otherwise" case %%
            far_cursors = sbr(cursor_i+param.samples_per_ui*(param.ndfe+1):param.samples_per_ui:end);
            precursors = sbr(cursor_i-param.samples_per_ui:-param.samples_per_ui:1);
            precursors = precursors(end:-1:1);

            % Error message if the sbr is not long enough for the specified range of Nb
            if length(sbr) < cursor_i+param.samples_per_ui*(param.ndfe+1)
                close(hwaitbar);
                error('Pulse Response contains %d samples after the cursor. Specified Nb requires %d samples after the cursor.' ...
                    , length(sbr)-cursor_i, param.samples_per_ui*(param.ndfe+1));
            end

            %% Equation 93A-27, when 1<=n<=N_b
            dfecursors=sbr(cursor_i+param.samples_per_ui*(1):param.samples_per_ui:cursor_i+param.samples_per_ui*(param.ndfe));
            excess_dfe_cursors = dfecursors - ...
                sign(dfecursors) .*min(sbr(cursor_i)*param.bmax(:), abs(dfecursors) );

            %% Eq. 93A-28 %%
            sampling_offset = mod(cursor_i, param.samples_per_ui);
            %ensure we can take early sample
            if sampling_offset<=1
                sampling_offset=sampling_offset+param.samples_per_ui;
            end
            if (OP.LIMIT_JITTER_CONTRIB_TO_DFE_SPAN)
                cursors_early_sample = sbr(cursor_i-1+param.samples_per_ui*(-1:param.ndfe));
                cursors_late_sample = sbr(cursor_i+1+param.samples_per_ui*(-1:param.ndfe));
            else
                cursors_early_sample = sbr(sampling_offset-1:param.samples_per_ui:end);
                cursors_late_sample = sbr(sampling_offset+1:param.samples_per_ui:end);
            end
            % ensure lengths are equal
            cursors_early_sample = cursors_early_sample(1:length(cursors_late_sample));
            h_J = (cursors_late_sample-cursors_early_sample)/2*param.samples_per_ui;

            %% Equation 93A-30 %%
            % since A_s = param.R_LM*cursor/(param.levels-1), cursor=(param.levels-1)*A_s/param.R_LM
            sigma_TX = (param.levels-1)*A_s/param.R_LM*10^(-param.SNDR/20);

            %% Equation 93A-31 %%
            sigma_ISI = param.sigma_X*norm([precursors; excess_dfe_cursors; far_cursors]);

            %% Equation 93A-32 %%
            sigma_J = norm([param.A_DD param.sigma_RJ])*param.sigma_X*norm(h_J);

            %% Equations 93A-33 and 93A-34 for FEXT (depends on TXFFE setting) %%
            if OP.RX_CALIBRATION
                sigma_XT=0;
            else
                sigma_FEXT =  get_xtlk_noise( upsampled_txffe, 'FEXT', param, chdata );
                sigma_XT = norm([sigma_NEXT sigma_FEXT]);
            end
            %% Equation 93A-36 denominator (actually its sqrt)
            total_noise_rms = norm([sigma_ISI sigma_J sigma_XT sigma_N sigma_TX sigma_ne]);

            %% Equation 93A-36 (note log argument is voltage rather than power ratio)
            FOM = 20*log10(A_s/total_noise_rms);

            if (FOM > best_FOM)
                best_txffe = txffe;
                best_sbr = sbr;
                best_ctle = ctle_index;
                best_FOM = FOM;
                best_cursor_i = cursor_i;
                best_IR=effective_channel;
                best_sigma_N = sigma_N;
                best_h_J = h_J;
                best_A_s=A_s;
            end
        end
        if OP.DISPLAY_WINDOW
            waitbar(progress, hwaitbar, 'Linear equalization tuning'); figure(hwaitbar); drawnow;
        else
            if ~mod(progress*100,20), fprintf('%i%% ', progress*100 );end
        end

    end
end

if ~exist('best_cursor_i', 'var')% take last setting
    best_txffe = txffe;
    best_sbr = sbr;
    best_ctle = ctle_index;
    best_FOM = FOM;
    best_cursor_i = cursor_i;
    best_IR=effective_channel;
    best_sigma_N = sigma_N;
    best_h_J = h_J;
end

best_cursor = best_sbr(best_cursor_i);
% report during debug
PRin=filter(ones(param.samples_per_ui, 1),1, chdata(1).uneq_imp_response);
f=1e8:1e8:100e9;

ctle_gain = (10^(gdc_values(best_ctle)/20) + 1i*f/param.CTLE_fz(best_ctle)) ./ ...
    ((1+1i*f/param.CTLE_fp1(best_ctle)).*(1+1i*f/param.CTLE_fp2(best_ctle)));

lsbr=length(sbr);
t=0:param.ui/param.samples_per_ui:(lsbr-1)*param.ui/param.samples_per_ui;

sampled_best_sbr_precursors_t = (best_cursor_i/param.samples_per_ui:-1:1/param.samples_per_ui)*param.ui;
sampled_best_sbr_precursors_t = sampled_best_sbr_precursors_t(end:-1:2); % exclude cursor
sampled_best_sbr_precursors = best_sbr(round(sampled_best_sbr_precursors_t/param.ui*param.samples_per_ui));
sampled_best_sbr_postcursors_t = (best_cursor_i:param.samples_per_ui:lsbr)/param.samples_per_ui*param.ui;
sampled_best_sbr_postcursors_t = sampled_best_sbr_postcursors_t(2:end); % exclude cursor
sampled_best_sbr_postcursors = best_sbr(round(sampled_best_sbr_postcursors_t/param.ui*param.samples_per_ui));
sampled_best_sbr_dfecursors_t = (best_cursor_i/param.samples_per_ui+(1:param.ndfe))*param.ui;
% apply max tap value constraint
dfe_cursors = sampled_best_sbr_postcursors(1:param.ndfe);
DFE_taps_mV = sign(dfe_cursors).*min(best_cursor*param.bmax(:), abs(dfe_cursors) );
sampled_best_sbr_postcursors(1:param.ndfe) = dfe_cursors-DFE_taps_mV;

if OP.DEBUG ~=0
    display(['FOM:                ' ,num2str(best_FOM, 2),' dB']);
    display(['TXFFE coefficients: ' ,mat2str(best_txffe) ] );
    display(['CTLE DC gain:       ' ,num2str(gdc_values(best_ctle)), ' dB']);
    display(['CTLE peaking gain:  ' ,num2str(20*log10(max(abs(ctle_gain))), 2), ' dB']);
    display(['Available signal:   ' ,num2str(best_cursor)]);

    if OP.DISPLAY_WINDOW && ~OP.RX_CALIBRATION
        figure(110);
        set(gcf, 'Name', 'CTLE selection');
        movegui(gcf, 'southeast');

        semilogx(f,20*log10(abs(ctle_gain)), 'disp', sprintf('Case %d', param.package_testcase_i));
        hold on
        fbaud_tick=find(f >= baud_rate, 1);
        fnq_tick=find(f >= baud_rate/2, 1);
        stem(f(fnq_tick),20*log10(abs(ctle_gain(fnq_tick))),'g', 'handlevisibility', 'off');
        stem(f(fbaud_tick),20*log10(abs(ctle_gain(fbaud_tick))),'g', 'handlevisibility', 'off');
        recolor_plots(gca);
        title('CTLE response')
        ylabel('dB')
        xlabel('Hz')

        % display pulse responses in one axis per test case.
        figure_name = sprintf('Equalization effect: %s', param.base);
        fig=findobj('Name', figure_name);
        if isempty(fig), fig=figure('Name', figure_name); end
        figure(fig);
        movegui(fig,'north')

        hax = subplot(length(OP.pkg_len_select), 1, param.package_testcase_i);

        plot(t,best_sbr,'disp', 'Equalized end-to-end PR');
        hold on
        PRplt(1:param.samples_per_ui+5)=PRin(1); % line up with the ffe introduced delay
        PRplt(param.samples_per_ui+6:length(t))=PRin(1:length(t)-param.samples_per_ui-5);
        plot(t,PRplt,'r','disp', 'Unequalized end-to-end PR');
        stem(t(best_cursor_i),best_sbr(best_cursor_i),'g','disp','Cursor (sample point)');
        title(sprintf('Case %d', param.package_testcase_i));
        ylabel('volts')
        xlabel('seconds')

        plot(sampled_best_sbr_precursors_t-param.ui/param.samples_per_ui, sampled_best_sbr_precursors', 'kx', 'disp','Pre cursors');
        plot(sampled_best_sbr_postcursors_t-param.ui/param.samples_per_ui, sampled_best_sbr_postcursors', 'ko', 'disp','Post cursors');
        stem(sampled_best_sbr_dfecursors_t-param.ui/param.samples_per_ui,dfe_cursors','m', 'LineWidth',2,'disp','DFE-canceled cursors');
        grid on
        legend show
        zoom xon
        set(hax, 'tag', 'EQE');
    end
end
if OP.DISPLAY_WINDOW
    close(hwaitbar);
else
    fprintf('\n');
end

% % eq_data
result.txffe = best_txffe;
result.ctle = best_ctle;
result.DFE_taps = DFE_taps_mV/best_cursor; %relative
result.DFE_taps_i = best_cursor_i+(1:param.ndfe)*param.samples_per_ui;

result.IR = best_IR;
result.A_s = best_A_s;

result.t_s = best_cursor_i;
result.sigma_N = best_sigma_N;
result.h_J = best_h_J;
result.FOM = best_FOM;

function [impulse_response, p1_ctle, p2_ctle, z_ctle] = TD_CTLE(ir_in, fb, f_z, f_p1, f_p2, kacdc_dB, oversampling)
%% Equation 93A-22 implemented in z-space and applied to the impulse response.
p1_ctle = -2*pi*f_p1;
p2_ctle = -2*pi*f_p2;
z_ctle = -2*pi*f_z*10^(kacdc_dB/20);
k_ctle = -p2_ctle;
bilinear_fs = 2*fb*oversampling;
p2d = (1+p2_ctle/bilinear_fs)./(1-p2_ctle/bilinear_fs);
p1d = (1+p1_ctle/bilinear_fs)./(1-p1_ctle/bilinear_fs);
zd = (1+z_ctle/bilinear_fs)./(1-z_ctle/bilinear_fs);
kd = (bilinear_fs-z_ctle)/((bilinear_fs-p1_ctle)*(bilinear_fs-p2_ctle));
B_filt =k_ctle*kd*poly([zd, -1]);
A_filt=poly([p1d, p2d]);
impulse_response=filter(B_filt,A_filt,ir_in);

function [ pdf ] = get_pdf_from_sampled_signal( input_vector, L, BinSize )
% Create PDF from interference vector using successive delta-set convolutions.
%   input_vector = list of values of samples
%   return
%   pdf.x
%   pdf.y
%   pdf.vec
%   pdf.bin

if max(input_vector) > BinSize
    input_vector=input_vector(abs(input_vector)>BinSize);
end
% for i = 1:length(input_vector)
%    if abs(input_vector(i)) < BinSize , input_vector(i)=0; end
%end
input_vector(abs(input_vector)<BinSize) = 0;

%% Equation 93A-39 %%
values = 2*(0:L-1)/(L-1)-1;
prob = ones(1,L)/L;

%% Initialize pdf to delta at 0
pdf=d_cpdf(BinSize, 0, 1);

for k = 1:length(input_vector)
    pdfn=d_cpdf(BinSize, abs(input_vector(k))*values, prob);
    pdf=conv_fct(pdf, pdfn);
end

function pdf=d_cpdf( binsize, values, probs)
%  p=cpdf(type, ...)
%
% CPDF is a probability mass function for discrete distributions or an 
% approxmation of a PDF for continuous distributions.
% 
% cpdf is internally normalized so that the sum of probabilities is 1
% (regardless of bin size).

% Internal fields:
% Min: *bin number* of minimum value.
% BinSize: size of PDF bins. Bin center is the representative value.
% Vec: vector of probabilities per bin.

values=binsize*round(values/binsize);
t=(min(values):binsize:max(values));
pdf.Min=min(values)/binsize;
pdf.y=zeros(size(t));
for k=1:length(values)
    [UNUSED_OUTPUT, bin]=min(abs(t-values(k))); %#ok<ASGLU>
    pdf.y(bin) = pdf.y(bin)+probs(k);
end

pdf.BinSize=binsize;
pdf.y=pdf.y/sum(pdf.y);
if any(~isreal(pdf.y)) || any(pdf.y<0)
    error('PDF must be real and nonnegative');
end
support=find(pdf.y);
pdf.y=pdf.y(support(1):support(end));
pdf.Min=pdf.Min+(support(1)-1);
pdf.x=(pdf.Min:-pdf.Min)*binsize;

function p=conv_fct(p1, p2)
if p1.BinSize ~= p2.BinSize
    error('bin size must be equal')
end

p=p1;
p.BinSize=p1.BinSize;
p.Min=p1.Min+p2.Min;
p.y=conv(p1.y, p2.y);
p.x =p.Min*p.BinSize:p.BinSize:-p.Min*p.BinSize;

function pdf = normal_dist(sigma,nsigma,binsize)
pdf.BinSize=binsize;
pdf.Min=-round(nsigma*sigma/binsize);
pdf.x=(pdf.Min:-pdf.Min)*binsize;
pdf.y=exp(-pdf.x.^2/(2*sigma^2+eps));
pdf.y=pdf.y/sum(pdf.y);

function [ param OP ]= read_ParamConfigFile(paramFile)
%warning('off','MATLAB:xlsread:Mode'); % suppress warning messages for reading the settings file from XLS
try
    [status,sheets] = xlsfinfo(paramFile); %#ok<ASGLU>
    if isempty(strfind(cell2mat(sheets), 'COM_Settings'))
        error('COM_Settings tab not in file')
    end
    [na1, na2, parameter] = xlsread(paramFile,'COM_Settings','',''); %#ok<ASGLU> % Import data from the settings file (imports the entire sheet)
catch ME %#ok<NASGU>
    warning('off','MATLAB:xlsread:Mode'); % suppress warning messages for reading the settings file from XLS
    [na1, na2, parameter] = xlsread(paramFile,'','','basic'); %#ok<ASGLU> % Import data from the settings file (imports the entire sheet)
end

% Default values are given for parameters when they are common to all clauses in 802.3bj and 803.2bm.

param.fb = xls_parameter(parameter, 'f_b')*1e9;
param.max_start_freq = xls_parameter(parameter, 'f_min')*1e9;
param.max_freq_step = xls_parameter(parameter, 'Delta_f')*1e9;
param.tx_ffe_c0_min = xls_parameter(parameter, 'c(0)', false);
param.tx_ffe_cm1_values = xls_parameter(parameter, 'c(-1)', true); % eval if string
param.tx_ffe_cp1_values = xls_parameter(parameter, 'c(1)', true); % eval if string
param.ndfe = xls_parameter(parameter, 'N_b');

param.ctle_gdc_values = xls_parameter(parameter, 'g_DC', true); % eval if string
param.CTLE_fp1 = 1e9*xls_parameter(parameter, 'f_p1', true, param.fb/4); % fp1 is in GHz
param.CTLE_fp2 = 1e9*xls_parameter(parameter, 'f_p2', true, param.fb); % fp2 is in GHz
param.CTLE_fz = 1e9*xls_parameter(parameter, 'f_z', true, param.fb/4); % fz is in GHz

param.a_thru = xls_parameter(parameter, 'A_v');
param.a_fext = xls_parameter(parameter, 'A_fe');
param.a_next = xls_parameter(parameter, 'A_ne');
param.levels = xls_parameter(parameter, 'L');
param.specBER = xls_parameter(parameter, 'DER_0');
param.pass_threshold = xls_parameter(parameter, 'COM Pass threshold');
param.sigma_RJ = xls_parameter(parameter, 'sigma_RJ');
param.A_DD = xls_parameter(parameter, 'A_DD');
param.eta_0 = xls_parameter(parameter, 'eta_0');

param.SNDR = xls_parameter(parameter, 'SNR_TX');
param.R_LM = xls_parameter(parameter, 'R_LM');

param.samples_per_ui = xls_parameter(parameter, 'M', 32);
% This will keep bmax length 0 if Nb=0
param.bmax(1:param.ndfe) = xls_parameter(parameter, 'b_max(1)');
param.bmax(2:param.ndfe) = xls_parameter(parameter, 'b_max(2..N_b)', true);

% eval if string for all three - can use different for TX and RX
param.C_pkg_board = xls_parameter(parameter, 'C_p', true)*1e-9; % C_p in nF
param.C_diepad = xls_parameter(parameter, 'C_d', true)*1e-9; % C_d in nF
param.R_diepad = xls_parameter(parameter, 'R_d', true);

param.Z0 = xls_parameter(parameter, 'R_0', 50);
param.z_p_tx_cases = xls_parameter(parameter, 'z_p (TX)', true); % eval if string
param.z_p_next_cases = xls_parameter(parameter, 'z_p (NEXT)', true); % eval if string
param.z_p_fext_cases = xls_parameter(parameter, 'z_p (FEXT)', true); % eval if string
param.z_p_rx_cases = xls_parameter(parameter, 'z_p (RX)', true); % eval if string

% Table 93A-3 parameters
param.pkg_gamma0_a1_a2 = xls_parameter(parameter, 'package_tl_gamma0_a1_a2', true, [0 1.734e-3 1.455e-4]);
param.pkg_tau = xls_parameter(parameter, 'package_tl_tau', false, 6.141e-3);
param.pkg_Z_c = xls_parameter(parameter, 'package_Z_c', false, 78.2);

% Table 92-12 parameters
param.brd_gamma0_a1_a2 = xls_parameter(parameter, 'board_tl_gamma0_a1_a2', true, [0 4.114e-4 2.547e-4]); % eval string, default
param.brd_tau = xls_parameter(parameter, 'board_tl_tau', false, 6.191e-3);
param.brd_Z_c = xls_parameter(parameter, 'board_Z_c', false, 109.8);

param.z_bp_tx = xls_parameter(parameter, 'z_bp (TX)', false, 151);
param.z_bp_next = xls_parameter(parameter, 'z_bp (NEXT)', false, 72);
param.z_bp_fext = xls_parameter(parameter, 'z_bp (FEXT)', false, 72);
param.z_bp_rx = xls_parameter(parameter, 'z_bp (RX)', false, 151);

% Unofficial parameters
param.snpPortsOrder = xls_parameter(parameter, 'Port Order', true, [1 3 2 4]);

% Deprecated parameters - affect only frequency domain analysis. better not change.
param.f_v = xls_parameter(parameter, 'f_v', false, 4);
param.f_f = xls_parameter(parameter, 'f_f', false, 4);
param.f_n = xls_parameter(parameter, 'f_n', false, 4);
param.f_r = xls_parameter(parameter, 'f_r', false, 4);

% Operational control variables
OP.include_pcb = xls_parameter(parameter, 'Include PCB (table 92-13)', false, 0);
OP.INCLUDE_CTLE = xls_parameter(parameter, 'INCLUDE_CTLE', false, 1);
OP.INCLUDE_FILTER = xls_parameter(parameter, 'INCLUDE_TX_RX_FILTER', false, 1);
OP.force_pdf_bin_size = xls_parameter(parameter, 'Force PDF bin size', false, 0);
OP.BinSize = xls_parameter(parameter, 'PDF bin size', false, 1e-5);
OP.DEBUG = xls_parameter(parameter, 'DIAGNOSTICS', false, false);
OP.DISPLAY_WINDOW = xls_parameter(parameter, 'DISPLAY_WINDOW', false, true);
OP.CSV_REPORT = xls_parameter(parameter, 'CSV_REPORT', false, true);
OP.SAVE_RESP = xls_parameter(parameter, 'SAVE_RESP', false, false);
OP.SAVE_FIGURES=xls_parameter(parameter, 'SAVE_FIGURES', false, false);
OP.SAVE_FIGURE_to_CSV=xls_parameter(parameter, 'SAVE_FIGURE_to_CSV', false, false);
OP.GET_FD = xls_parameter(parameter, 'Display frequency domain', false, false);
OP.INC_PACKAGE = xls_parameter(parameter, 'INC_PACKAGE', false, true);
OP.IDEAL_RX_TERM = xls_parameter(parameter, 'IDEAL_RX_TERM', false, false);
OP.IDEAL_TX_TERM = xls_parameter(parameter, 'IDEAL_TX_TERM', false, false);
OP.EXTERNAL = xls_parameter(parameter, 'USE_EXTERNAL_PARAM',false, false);
OP.RESULT_DIR = regexprep(xls_parameter(parameter, 'RESULT_DIR'), '\\', filesep);
OP.BREAD_CRUMBS = xls_parameter(parameter, 'BREAD_CRUMBS',false, false);
OP.ENFORCE_CAUSALITY = xls_parameter(parameter, 'Enforce Causality', false, 0);
OP.EC_REL_TOL = xls_parameter(parameter, 'Enforce Causality REL_TOL', false, 1e-2);
OP.EC_DIFF_TOL = xls_parameter(parameter, 'Enforce Causality DIFF_TOL', false, 1e-3);
OP.EC_PULSE_TOL = xls_parameter(parameter, 'Enforce Causality pulse start tolerance', false, 0.01);
OP.pkg_len_select = xls_parameter(parameter, 'z_p select', true, 1);  % eval if string
OP.RX_CALIBRATION = xls_parameter(parameter, 'RX_CALIBRATION', false, false);
OP.sigma_bn_STEP = xls_parameter(parameter, 'Sigma BBN step', false, 5e-3);
OP.BBN_Q_factor = xls_parameter(parameter, 'BBN Q factor', false, 5);
OP.force_BBN_Q_factor = xls_parameter(parameter, 'Force BBN Q factor', false, false);
OP.transmitter_transition_time = xls_parameter(parameter, 'T_r', false, 8e-3);
OP.LIMIT_JITTER_CONTRIB_TO_DFE_SPAN = xls_parameter(parameter, 'LIMIT_JITTER_CONTRIB_TO_DFE_SPAN', false, false);
OP.impulse_response_truncation_threshold = xls_parameter(parameter, 'Impulse response truncatio threshold', false, 1e-3);
OP.interp_sparam_mag = xls_parameter(parameter, 'S-parameter magnitude extrapolation policy', false, 'linear_trend_to_DC');
OP.interp_sparam_phase = xls_parameter(parameter, 'S-parameter phase extrapolation policy', false, 'extrap_cubic_to_dc_linear_to_inf');

% Parameters for error burst probability calculation. Not officially used.
OP.use_simple_EP_model = xls_parameter(parameter, 'Use simple error propagation model', false, false);
OP.nburst = xls_parameter(parameter, 'Max burst length calculated', false, 0);
OP.COM_EP_margin = xls_parameter(parameter, 'Error propagation COM margin', false, 0);


function p=xls_parameter(param_sheet, param_name, eval_if_string, default_value)
% helper function to read parameter values from XLS file. Uses names to find values.
if nargin<3, eval_if_string=0; end
[row, col]=find(strcmp(param_sheet, param_name));
if numel(row)*numel(col)==0
    if nargin<4
        missingParameter(param_name);
    else
        p = default_value;
    end
elseif numel(row)*numel(col)>1
    % if there are several occurrences, use the first, but warn
    warning('COM:XLS_parameter:MultipleOccurrence', ...
        '%d occurrences of "%s" found. Using the first', numel(row), param_name);
    p=param_sheet{row(1), col(1)+1};
else
    p=param_sheet{row, col+1};
end
if ischar(p) && eval_if_string
    p=eval(p);
end

function missingParameter (parameterName)
error( 'error:badParameterInformation', ...
    'The data for mandatory parameter %s is missing or incorrect' , parameterName);

function [sdd21p] = s21_pkg(sdd21, sdd11, sdd22, faxis, param, OP, channel_type, channel_number)
% concatenates package reflections with s21,s11,and s22 with spec return loss (gammas)
% faxis is the frequency array
% sdd21, sdd11, sdd22 are the corresponding array of differential parameters
Z0=param.Z0;

% The following three parameters have possibly different values for TX and
% RX (so can be 2-element vectors).
R_diepad = param.R_diepad;
C_diepad = param.C_diepad;
C_pkg_board = param.C_pkg_board;

% generate TX package according to channel type.
switch upper(channel_type)
    case 'THRU'
        pkg_len = param.Pkg_len_TX;
    case 'NEXT'
        pkg_len = param.Pkg_len_NEXT;
    case 'FEXT'
        pkg_len = param.Pkg_len_FEXT;
end
[ s11out, s12out, s21out, s22out ]= make_pkg(faxis, pkg_len, C_diepad(1), C_pkg_board(1), param);

% RX package length is assumed to be the same for all channel types.
[ s11in, s12in, s21in, s22in ]= make_pkg(faxis, param.Pkg_len_RX, C_diepad(2), C_pkg_board(2), param);

debug_plot_package=0; % set breakpoint and change manually if desired
if debug_plot_package
    figure; %#ok<UNRCH>
    subplot(2,1,1);
    plot(faxis/1e9, 20*log10(abs(s12out)), 'disp', 'TX package')
    hold on
    plot(faxis/1e9, 20*log10(abs(s21in)), 'disp', 'RX package');
    xlabel('f [GHz]'); ylabel('IL [dB]'); recolor_plots_in_axis(gca);
    subplot(2,1,2);
    plot(faxis/1e9, unwrap(angle(s12out)), 'disp', 'TX package')
    hold on
    plot(faxis/1e9, unwrap(angle(s21in)), 'disp', 'RX package');
    xlabel('f [GHz]'); ylabel('Phase [rad]'); recolor_plots_in_axis(gca);
end

if OP.IDEAL_TX_TERM || (OP.RX_CALIBRATION == 1 && channel_number == 2)
    gamma_tx=0;
else
    gamma_tx=(R_diepad(1)-Z0)/(R_diepad(1)+Z0);% equation 93A-17
end
if OP.IDEAL_RX_TERM
    gamma_rx=0;
else
    gamma_rx=(R_diepad(2)-Z0)/(R_diepad(2)+Z0);% equation 93A-17
end

if OP.INC_PACKAGE==0
    sdd21p= sdd21;
else
    if OP.RX_CALIBRATION == 1 && channel_number == 2
        %   for calibration do not include the transmitter package
        [s11out_rx, UNUSED_OUTPUT, s21out_rx, s22out_rx ] = combines4p( sdd11, sdd21, sdd21, sdd22, s22out, s12out, s21out, s11out ); %#ok<ASGLU> % s22 is ball side of package
        %% Equation 93A-18
        sdd21p= s21out_rx.*(1-gamma_tx).*(1+gamma_rx)./(1.- s11out_rx.*gamma_tx - s22out_rx.*gamma_rx  -s21out_rx.^2.*gamma_tx.*gamma_rx +s11out_rx.*s22out_rx.*gamma_tx.*gamma_rx);
    else
        %% Equations 93A-4 to 93A-7
        if ~OP.IDEAL_TX_TERM
            [sdd11, UNUSED_OUTPUT, sdd21, sdd22] = combines4p( s11out, s12out, s21out, s22out, sdd11, sdd21, sdd21, sdd22 ); %#ok<ASGLU>
        else % limited transition time RX may be used in RX compliance
            H_t = exp(-(pi*faxis/1e9*OP.transmitter_transition_time/1.6832).^2); %% Equation 93A-46 %%
            sdd21 = sdd21 .* H_t;
        end
        if ~OP.IDEAL_RX_TERM
            [sdd11, UNUSED_OUTPUT, sdd21, sdd22] = combines4p( sdd11, sdd21, sdd21, sdd22, s22in, s12in, s21in, s11in ); %#ok<ASGLU> % s22 is ball side of package
        end
        %% Equation 93A-18
        sdd21p= sdd21.*(1-gamma_tx).*(1+gamma_rx)./(1.- sdd11.*gamma_tx - sdd22.*gamma_rx  -sdd21.^2.*gamma_tx.*gamma_rx +sdd11.*sdd22.*gamma_tx.*gamma_rx);
    end
end

function [ s11out, s12out, s21out, s22out ] = make_pkg(f, pkg_len, cpad, cball, param)
f(f<eps)=eps;
tau = param.pkg_tau; gamma0_a1_a2=param.pkg_gamma0_a1_a2; Lenscale=pkg_len; zref=param.Z0;
%% Equation 93A-8
s11pad= -1i*2*pi.*f*cpad*zref./(2+1i*2*pi.*f*cpad*zref);
s21pad= 2./(2+1i*2*pi.*f*cpad*zref);

[ S11, S12, S21, S22 ] = synth_tline(f, param.pkg_Z_c, zref, gamma0_a1_a2, tau, Lenscale); %#ok<NASGU,ASGLU>
[ s11out1, s12out1, s21out1, s22out1 ]= ...
    combines4p(  s11pad, s21pad, s21pad, s11pad, S11, S21, S21, S11 ); % first part of equation 93A-15

%% Equation 93A-8
s11ball= -1i*2*pi.*f*cball*zref./(2+1i*2*pi.*f*cball*zref);
s21ball= 2./(2+1i*2*pi.*f*cball*zref);
[ s11out, s12out, s21out, s22out ]= ...
    combines4p( s11out1, s12out1, s21out1, s22out1, s11ball, s21ball, s21ball, s11ball );% second part of equation 93A-15

function [ s11out, s12out, s21out, s22out ] = add_brd(chdata, param)
% Used in Clause 92 for adding board trace between TP0 and TP2
switch chdata.type
    case 'THRU'
        z_bp_tx = param.z_bp_tx;
    case 'NEXT'
        z_bp_tx = param.z_bp_next;
    case 'FEXT'
        z_bp_tx = param.z_bp_fext;
end

[ s11tx, s12tx, s21tx, s22tx ] = synth_tline(chdata.faxis, param.brd_Z_c, param.Z0, param.brd_gamma0_a1_a2, param.brd_tau, z_bp_tx);
[ s11rx, s12rx, s21rx, s22rx ] = synth_tline(chdata.faxis, param.brd_Z_c, param.Z0, param.brd_gamma0_a1_a2, param.brd_tau, param.z_bp_rx);

[ s11out1, s12out1, s21out1, s22out1 ]=combines4p(  s11tx, s12tx, s21tx, s22tx, chdata.sdd11_raw, chdata.sdd21_raw, chdata.sdd12_raw, chdata.sdd22_raw );
[ s11out, s12out, s21out, s22out ]=combines4p(  s11out1, s12out1, s21out1, s22out1,  s11rx, s12rx, s21rx, s22rx);

function [s11, s12, s21, s22] = synth_tline(f, Z_c, Z_0, gamma_coeff, tau, d)
f_GHz=f/1e9;
%% Equation 93A-10 %%
gamma_1 = gamma_coeff(2)*(1+1i);
%% Equation 93A-11 %%
gamma_2 = gamma_coeff(3)*(1-2i/pi*log(f_GHz)) + 2i*pi*tau;
%% Equation 93A-9 %%
gamma = gamma_coeff(1)+gamma_1.*sqrt(f_GHz)+gamma_2.*f_GHz;
gamma(f_GHz==0) = gamma_coeff(1);

%% Equation 93A-12 %%
rho_rl=(Z_c-2*Z_0)/(Z_c+2*Z_0);

exp_gamma_d = exp(-d*gamma);
%% Equations 93A-13 and 93A-14 %%
s11 = rho_rl*(1-exp_gamma_d.^2)./(1-rho_rl^2*exp_gamma_d.^2);
s21 = (1-rho_rl^2)*exp_gamma_d./(1-rho_rl^2*exp_gamma_d.^2);
s12 = s21;
s22 = s11;

function [ s11out, s12out, s21out, s22out ] = combines4p( s11in1, s12in1, s21in1, s22in1, s11in2, s12in2, s21in2, s22in2)

s1=zeros(2,2,length(s11in1)); s2=s1; t3=s1;
for i=1:length(s11in1)
    s1(:,:,i)=[s11in1(i) s12in1(i); s21in1(i) s22in1(i) ];
    s2(:,:,i)=[s11in2(i) s12in2(i); s21in2(i) s22in2(i) ];
end
t1=stot(s1);
t2=stot(s2);
for i=1:length(s11in1)
    t3(:,:,i)=t1(:,:,i)*t2(:,:,i);
end
s3=ttos(t3);
s11out=s3(1,1,:);
s11out=transpose(s11out(:));
s12out=s3(1,2,:);
s12out=transpose(s12out(:));
s21out=s3(2,1,:);
s21out=transpose(s21out(:));
s22out=s3(2,2,:);
s22out=transpose(s22out(:));

function t_params = stot(s_params)
% p 67 R. Mavaddat. (1996). Network scattering parameter. Singapore: World Scientific.
% ISBN 978-981-02-2305-2. http://books.google.com/?id=287g2NkRYxUC&lpg=PA65&dq=T-parameters+&pg=PA67.
[s11, s12, s21, s22] = deal(s_params(1,1,:), s_params(1,2,:), s_params(2,1,:), s_params(2,2,:));
delta = (s11.*s22-s12.*s21);
s21(s21==0)=eps;
t_params = [1./s21, -s22./s21; s11./s21, -delta./s21];

function s_params = ttos(t_params)
% p 67 R. Mavaddat. (1996). Network scattering parameter. Singapore: World Scientific.
% ISBN 978-981-02-2305-2. http://books.google.com/?id=287g2NkRYxUC&lpg=PA65&dq=T-parameters+&pg=PA67.
[t11, t12, t21, t22] = deal(t_params(1,1,:), t_params(1,2,:), t_params(2,1,:), t_params(2,2,:));
delta = t11.*t22-t21.*t12;
t11(t11==0)=eps;
s_params = [t21./t11, delta./t11; 1./t11, -t12./t11];

function [ sigma_NE ] = get_sigma_noise( H_ctf, param, chdata, sigma_bn )
% the FEXT channel for calibration basically a DC connection unlike normal
% FEXT channels which are nearly open at DC channels
H_r = 1./polyval([1 2.613126 3.414214 2.613126 1], 1i*chdata(2).faxis/(param.f_r*param.fb));
idxfbby2=find( chdata(2).faxis(:) >= param.fb/2, 1);
if size(chdata,2) >= 2
    Hnoise_channel=chdata(2).sdd21;% rx package is already included tx is not
else
    Hnoise_channel=1;
end
H_np=Hnoise_channel.*H_ctf.*H_r;
%% Equation 93A-48 %%
sigma_NE = sigma_bn*sqrt(mean(abs(H_np(1:idxfbby2).^2)));

function sigma_XT = get_xtlk_noise( upsampled_txffe, type, param, chdata )
% calculate crosstalk sigma at worst phase, per equation 93A-33, for
% channels of chosen type (FEXT and NEXT are treated separately).

sigma_XT_sqd=0;
for k=1:param.num_s4p_files
    if isequal(chdata(k).type, type)
        sigma_i_sqd=0;
        effective_channel = filter(upsampled_txffe, 1, chdata(k).ctle_imp_response);
        N=round(length(effective_channel)/param.samples_per_ui)-3;
        sbr=filter(ones(param.samples_per_ui, 1), 1, effective_channel);
        for m=1:param.samples_per_ui
            h=sbr(param.samples_per_ui*(1:N)+m);
            sigma_i_sqd=max(param.sigma_X^2*sum(h.^2),sigma_i_sqd);
        end
        sigma_XT_sqd = sigma_XT_sqd + sigma_i_sqd;
    end
end
sigma_XT=sqrt(sigma_XT_sqd);

function [chdata, param] = get_s4p_files(param, OP, num_fext, num_next, file_list)
% filename parsing and acquisition
%------------------------------------------------------------------
%----------put files names into chdata structure ---------
% The thru file has the index of 1
% crosstalk file are indexed from 2
% nxi is incremented each time a file is read in  so that nxi will end
filepath=[]; % path name for file
nxi=0; % file index
% get the THRU file
if size(file_list,2) ~= 0
    file_list(1)=strrep(file_list(1),'\', filesep); % OS file convention conversion
    [filepath, basename, fileext]=fileparts(file_list{1});

else
    if OP.RX_CALIBRATION == 1
        h = msgbox('enter test channel s-parameter file'); set(h,'Color',[1 .85 0]);
        movegui(h,'northeast')
    end
    dir=fullfile(filepath, '*.s4p');
    [basename,filepath]=uigetfile(dir,'input thru channel response .s4p ');
    if filepath == 0
        error('No Thru file')
    end
    [UNUSED_OUTPUT, basename, fileext]=fileparts(basename); %#ok<ASGLU>
end
nxi=nxi+1;
chdata(nxi).filename = fullfile(filepath,  [basename fileext]);
chdata(nxi).ext = fileext;
[UNUSED_OUTPUT, dirname]=fileparts(filepath); %#ok<ASGLU>
%    chdata(nxi).base=[pth(max(strfind(pth,filesep))+1:end) '--' basename]; % add 1 directory back to basename
chdata(nxi).base=[dirname '--' basename]; % add 1 directory back to basename
chdata(nxi).A=param.a_thru; % pam encoding amplitude reduction is do in reporting but is comprehending in crosstalk PDF
chdata(nxi).type='THRU';
chdata(nxi).ftr=param.fb*param.f_v;
param.base=chdata(nxi).base; %for print out function that don't pass chdata

% now get FEXT file names into chdata structure
kxi=nxi;
for nxi=kxi+1:num_fext+kxi
    lastfilepath=filepath;
    if size(file_list,2) ~= 0
        [filepath, basename, fileext]=fileparts(file_list{nxi});
    else
        if OP.RX_CALIBRATION == 1
            h = msgbox('enter noise channel s-parameter file'); set(h,'Color',[1 .85 0]);
            movegui(h,'northeast')
        end
        dir=fullfile(filepath, '*.s4p');
        if OP.RX_CALIBRATION == 1
            [basename,filepath]=uigetfile(dir,'input noise channel response .s4p');
        else
            [basename,filepath]=uigetfile(dir,['input fext channel response .s4p #', num2str(nxi-kxi)]);
        end
        if filepath==0
            error('Not enough NEXT files')
        end
        [UNUSED_OUTPUT, basename, fileext]=fileparts(basename); %#ok<ASGLU>
    end
    if isempty( filepath), filepath=lastfilepath; end
    chdata(nxi).filename = fullfile(filepath,  [basename fileext]);
    chdata(nxi).ext = fileext;
    [UNUSED_OUTPUT, dirname]=fileparts(filepath); %#ok<ASGLU>
    chdata(nxi).base=[dirname '--' basename];
    chdata(nxi).A=param.a_fext;
    chdata(nxi).ftr=param.fb*param.f_f;
    chdata(nxi).type='FEXT';
end
% now get NEXT file names into chdata structure
kxi=num_fext+kxi;
for nxi=kxi+1:num_next+kxi
    lastfilepath=filepath;
    if size(file_list,2) ~= 0
        [filepath, basename, fileext]=fileparts(file_list{nxi});
    else
        dir=fullfile(filepath, '*.s4p');
        [basename,filepath]=uigetfile(dir,['input next channel response .s4p ', num2str(nxi-kxi)]);
        if filepath==0
            error('Not enough NEXT files')
        end
        [UNUSED_OUTPUT, basename, fileext]=fileparts(basename); %#ok<ASGLU>
    end
    if isempty( filepath), filepath=lastfilepath; end
    chdata(nxi).filename = fullfile(filepath,  [basename fileext]);
    chdata(nxi).ext = fileext;
    [UNUSED_OUTPUT, dirname]=fileparts(filepath); %#ok<ASGLU>
    chdata(nxi).base=[dirname '--' basename];
    chdata(nxi).A=param.a_next;
    chdata(nxi).ftr=param.fb*param.f_n;
    chdata(nxi).type='NEXT';
end

function [chdata, param] = read_s4p_files(param, OP, chdata)
%% extract s-parameter and convert to differential mode
% extract s-parameter data from files and apply tx and rx filters as well as package filters
num_files=length(chdata);
if ~OP.DISPLAY_WINDOW, fprintf('reading file '); end
for i=1:num_files
    if OP.DISPLAY_WINDOW; hwaitbar=waitbar(0);end
    progress = i/num_files;
    if OP.DISPLAY_WINDOW
        waitbar(progress, hwaitbar, ['Processing ' chdata(i).base]); figure(hwaitbar); drawnow;
    else
        fprintf('%i ',i);
    end

    % Skip reading file if it was already read (multiple test cases)
    if (~isfield(chdata(i), 'faxis')) || isempty(chdata(i).faxis)

        switch chdata(i).ext
            case '.s4p'
                if ~isfield(param,'snpPortsOrder') || isempty(param.snpPortsOrder)
                    param.snpPortsOrder = [1 3 2 4]; % default order normally used.
                elseif length(param.snpPortsOrder) ~= 4
                    error( 'warning:sNpFilePortMismatch', ...
                        '\n\t The number of ports defined (%G) does not match the sNp file type (%s)', ...
                        length(param.snpPortsOrder), ...
                        chdata(i).ext ...
                        );
                end
                % read function returns differnetial mode parameters
                [Sch SDDch] = read_p4_s4params(chdata(i).filename,  0, 0, param.snpPortsOrder, OP);
                chdata(i).fmaxi = length(Sch.freq);

                if Sch.freq(chdata(i).fmaxi) < param.fb
                    warning('COM:read_s4p:MaxFreqTooLow', ...
                        'In %s: the maximum frequency provided, %g, is less than the signaling rate: %g', ...
                        chdata(i).filename, Sch.freq(end), param.fb);
                end
                if Sch.freq(1) > param.max_start_freq
                    warning('COM:read_s4p:StartFreqTooHigh', ...
                        'In %s: minimum frequency, %.2g GHz, is larger than the recommended %.2g GHz', ...
                        chdata(i).filename, Sch.freq(1)/1e9, param.max_start_freq/1e9);
                end
                freqstep=diff(Sch.freq);
                % ignore frequency differences up to 1 Hz - possible numerical artifacts
                if max(freqstep)-min(freqstep) > 1
                    warning('COM:read_s4p:NonUniformFreqSpacing', 'In %s: non-uniform frequency steps: min=%.3g GHz, max=%.3g GHz', ...
                        chdata(i).filename, min(freqstep)/1e9, max(freqstep)/1e9);
                end
                if max(freqstep) - param.max_freq_step > 1
                    warning('COM:read_s4p:FreqStepTooHigh', 'In %s: frequency step, %.2g GHz, is larger than the recommended %.2g GHz', ...
                        chdata(i).filename, max(freqstep)/1e9, param.max_freq_step/1e9);
                end

                chdata(i).faxis = Sch.freq;
                chdata(i).sdd12_raw = transpose(SDDch(1:chdata(i).fmaxi,1,2));
                chdata(i).sdd21_raw = transpose(SDDch(1:chdata(i).fmaxi,2,1));
                chdata(i).sdd22_raw = transpose(SDDch(1:chdata(i).fmaxi,2,2));
                chdata(i).sdd11_raw = transpose(SDDch(1:chdata(i).fmaxi,1,1));

            case '.s12p' % not fully inplemented, this is for a 3 pair model only
                if ~isfield(param,'snpPortsOrder') || isempty(param.snpPortsOrder)
                    param.snpPortsOrder = [1 2 3 4 6 7 8 9 10 11 12]; % default order normally used.
                elseif length(param.snpPortsOrder) ~= 12
                    error( 'warning:sNpFilePortMismatch', ...
                        '\n\t The number of ports defined (%G) does not match the sNp file type (%s)', ...
                        length(param.snpPortsOrder), ...
                        chdata(i).ext ...
                        );
                end
                Sch = readdataSnPx(chdata(i).filename, 12);
                chdata(i).fmaxi = length(Sch.freq);
                chdata(i).faxis = Sch.freq;
                if isequal([1 3 5 7 9 11 2 4 6 8 10 12], param.snpPortsOrder)
                    chdata(i).sdd21_raw = reduce((Sch.cs(6,5,1:chdata(i).fmaxi)-Sch.cs(8,5,1:chdata(i).fmaxi)-Sch.cs(6,7,1:chdata(i).fmaxi)+Sch.cs(8,7,1:chdata(i).fmaxi))./2); %dd21=s21-s41-s23+s43
                    chdata(i).sdd11_raw = reduce((Sch.cs(5,5,1:chdata(i).fmaxi)-Sch.cs(7,5,1:chdata(i).fmaxi)-Sch.cs(5,7,1:chdata(i).fmaxi)+Sch.cs(7,7,1:chdata(i).fmaxi))./2); %dd11=s11-s31-s13+s33
                    chdata(i).sdd22_raw = reduce((Sch.cs(6,6,1:chdata(i).fmaxi)-Sch.cs(6,8,1:chdata(i).fmaxi)-Sch.cs(8,6,1:chdata(i).fmaxi)+Sch.cs(8,8,1:chdata(i).fmaxi))./2); %dd22=s22-s41-s14+s44
                elseif isequal([1 2 3 4 5 6 7 8 9 10 11 12], param.snpPortsOrder)
                    chdata(i).sdd21_raw = reduce((Sch.cs(9,3,1:chdata(i).fmaxi)-Sch.cs(10,3,1:chdata(i).fmaxi)-Sch.cs(9,4,1:chdata(i).fmaxi)+Sch.cs(10,4,1:chdata(i).fmaxi))./2);  %dd21=s21-s41-s23+s43
                    chdata(i).sdd11_raw = reduce((Sch.cs(3,3,1:chdata(i).fmaxi)-Sch.cs(4,9,1:chdata(i).fmaxi)-Sch.cs(3,4,1:chdata(i).fmaxi)+Sch.cs(4,4,1:chdata(i).fmaxi) )/2);    %dd11=s11-s31-s13+s33
                    chdata(i).sdd22_raw = reduce((Sch.cs(9,9,1:chdata(i).fmaxi)-Sch.cs(10,3,1:chdata(i).fmaxi)-Sch.cs(9,10,1:chdata(i).fmaxi)+Sch.cs(10,10,1:chdata(i).fmaxi) )/2);%dd22=s22-s41-s14+s44
                else
                    display('INVALID ENTRY FOR TX PORTS.  Use [1 2 3 4 5 6] or [1 3 5 7 9 11]')
                    return
                end
            otherwise
        end
        %
        % differential response is sdd21 <------------------------------
        % differential return losses are sdd11 and sdd22
        if OP.include_pcb
            % add boards to sdd
            [chdata(i).sdd11_raw, chdata(i).sdd12_raw, chdata(i).sdd21_raw, chdata(i).sdd22_raw] = add_brd(chdata(i), param);
        end
        chdata(i).sdd11=chdata(i).sdd11_raw;
        chdata(i).sdd22=chdata(i).sdd22_raw;
    end
    if  OP.INC_PACKAGE ~= 0,
        % the new sdd
        % changes to s21_pkg for calibration is to not include Tx package in the sp21p
        chdata(i).sdd21p= s21_pkg(chdata(i).sdd21_raw, chdata(i).sdd11_raw, chdata(i).sdd22_raw, chdata(i).faxis, param, OP, chdata(i).type, i);
        chdata(i).sdd21=chdata(i).sdd21p;
    else
        chdata(i).sdd21=chdata(i).sdd21_raw;
    end
    chdata(i).sdd21f=chdata(i).sdd21_raw; % used for FD analysis i.e. not filtered
end
if ~OP.DISPLAY_WINDOW, fprintf('\n'); end

function plot_bathtub_curves(hax, max_signal, sci_pdf, cci_pdf, isi_and_xtalk_pdf, noise_pdf, combined_interference_and_noise_pdf, bin_size)
cursors = d_cpdf(bin_size,max_signal*[-1 1], [1 1]/2);
signal_and_isi_pdf = conv_fct(cursors, sci_pdf);
signal_and_xtalk_pdf = conv_fct(cursors, cci_pdf);
signal_and_channel_noise_pdf = conv_fct(cursors, isi_and_xtalk_pdf);
signal_and_system_noise_pdf = conv_fct(cursors, noise_pdf);
signal_and_total_noise_pdf = conv_fct(cursors, combined_interference_and_noise_pdf);

semilogy(signal_and_isi_pdf.x, abs(cumsum(signal_and_isi_pdf.y)-0.5) ,'r','Disp','ISI', 'parent', hax)
hold on
semilogy(signal_and_xtalk_pdf.x, abs(cumsum(signal_and_xtalk_pdf.y)-0.5) ,'b','Disp','Xtalk', 'parent', hax)
semilogy(signal_and_channel_noise_pdf.x, abs(cumsum(signal_and_channel_noise_pdf.y)-0.5) ,'c','Disp','ISI+Xtalk', 'parent', hax)
semilogy(signal_and_system_noise_pdf.x, abs(cumsum(signal_and_system_noise_pdf.y)-0.5) ,'m','Disp','Jitter, TX and system noise', 'parent', hax)
semilogy(signal_and_total_noise_pdf.x, abs(cumsum(signal_and_total_noise_pdf.y)-0.5) ,'k','Disp','total noise PDF', 'parent', hax)
hc=semilogy(max_signal*[-1 -1 1 1], [0.5 1e-20 1e-20 0.5], '--ok');
set(get(get(hc,'Annotation'),'LegendInformation'), 'IconDisplayStyle','off');

ylabel(hax, 'Probability')
xlabel(hax, 'volts')
legend(hax, 'show')

function recolor_plots(ax)
colors='brgcmk';
ch=flipud(get(ax, 'children'));

for k=1:length(ch)
    set(ch(k), 'Color', colors(mod(k-1, length(colors))+1));
    set(ch(k), 'LineWidth', 2*floor((k-1)/length(colors))+1);
end
legend (ax, 'off');
warning('off', 'MATLAB:legend:PlotEmpty');
set(legend (ax, 'show'), 'interp', 'none');

function csv_string = str2csv(c)
% convert a cell array of strings to a csv string
cell_tmp = cell(2, length(c));
cell_tmp(1,:)=c;
cell_tmp(2,:) = {','};
cell_tmp{2,end} = '';
csv_string=strcat(cell_tmp{:});

function [ h ] = savefigs( param, OP )

%% find the figures
hw = waitbar(0,'Saving figures...');
h = findobj(0, 'Type', 'figure');
for ii=1:length(h)
    figname= get(h(ii), 'Name'); % use the figure name as file name
    if isempty(strfind(figname,param.base))
        figname = [figname '_' param.base]; %#ok<AGROW>
    end
    figname = ['f_' num2str(h(ii)) '_' figname]; %#ok<AGROW>
    figname = strrep(figname,':','-');
    figname = strrep(figname,' ','_');
    if OP.SAVE_FIGURES==1
        saveas(h(ii), fullfile(OP.RESULT_DIR, [figname '.fig']));
    end
    %% get x y data
    if OP.SAVE_FIGURE_to_CSV==1
        h_L = findobj(h(ii),'Type','line'); % find handles to all the lines
        M=[]; %ncol=1;
        for nk=1:length(h_L)
            % get x and data for a line.
            x_data=get(h_L(nk),'xdata')';
            y_data=get(h_L(nk),'ydata')';
            % .........>> need to get data in the line structure (legend or label) for headers
            M=[M; x_data; y_data]; %#ok<AGROW>
        end
        csvwrite([OP.RESULT_DIR figname '.csv'],M);
        %      clear M y x header h_L
    end
    waitbar(ii/length(h),hw)

end

close(hw)
