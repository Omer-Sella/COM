function result=optimize_fom(OP, param, chdata, sigma_bn,do_C2M)
%% input
% chdata(1).uneq_imp_response is the impulse response input expected to be normalized to At, peak drive voltage
% baud_rate - baud rate in seconds
% param.samples_per_ui = samples per UI of IR
% param.max_ctle - maximum ac to dc gain in dB
% param.tx_ffe(1) - maximum pre cursor (positive value)
% param.tx_ffe(2) - maximum post cursor (positive value)
% param.tx_ffe_step - sweep step size for tx pre and post taps
% param.ndfe - number of reference dfe taps
% do_C2M.  set to 0 for standard optimize_fom.  set to 1 for optimize_fom_for_C2M
% output
% result.eq.txle - [ precusor curosr postcursor]: pre and post are negative
% result.eq.ctle - index of CTLE parameters in table
% result.IR - impulse response
% result.avail_signal - maximum signal after equalization
% result.avail_sig_index - index in result.IR of max signal
% result.best_FOM - best raw ISI


min_number_of_UI_in_response=40;
baud_rate=1/param.ui;
% H_r = 1./polyval([1 2.613126 3.414214 2.613126 1], 1i*chdata(1).faxis/(param.f_r*param.fb));
f=chdata(1).faxis;

%Read user input of ts_sample_adj_range
%if one value was entered, go from 0 to that value
%if 2 values were entered, go from the 1st value to the 2nd value
if length(param.ts_sample_adj_range)==1
    param.ts_sample_adj_range(2)=param.ts_sample_adj_range(1);
    param.ts_sample_adj_range(1)=0;
end
full_sample_range=param.ts_sample_adj_range(1):param.ts_sample_adj_range(2);

H_bt=Bessel_Thomson_Filter(param,f,OP.Bessel_Thomson);
H_bw=Butterworth_Filter(param,f,OP.Butterworth);
H_RCos=Raised_Cosine_Filter(param,f,OP.Raised_Cosine);% experiment with RCos
% need to include H_RCos in noise and when computing the system ir for thru
% and crosstalk
H_r=H_bw.*H_bt.*H_RCos;
%% Bill Kirkland, need to get auto correlation of H_r.*HCTLE
% Get f vector from 0 to Fs/2-delta_f.
N_fft_by2 = 512;
f_xc = ((1:N_fft_by2)-1)/N_fft_by2*baud_rate*param.samples_per_ui/2;
H_bt_xc=Bessel_Thomson_Filter(param,f_xc,OP.Bessel_Thomson);
H_bw_xc=Butterworth_Filter(param,f_xc,OP.Butterworth);
H_RCos_xc=Raised_Cosine_Filter(param,f_xc,OP.Raised_Cosine);
H_r_xc=H_bw_xc.*H_bt_xc.*H_RCos_xc;
%%

% system noise H_sy PSD
if OP.USE_ETA0_PSD
    fspike=1e9;
    % requires communication tool box if used
    H_sy=sinc(sqrt(2)*(chdata(1).faxis-fspike)/(fspike)).^2;
else
    H_sy=ones(1,length(chdata(1).faxis));
end

%Build txffe values dynamically
%any param field that is "tx_ffe_cm<X>_values" is a precursor
%any param field that is "tx_ffe_cp<X>_values" is a postcursor
%where <X> is any integer
param_fields=fieldnames(param);
num_pre=length(find(~cellfun('isempty',regexp(param_fields,'tx_ffe_cm\d+_values'))));
num_post=length(find(~cellfun('isempty',regexp(param_fields,'tx_ffe_cp\d+_values'))));
num_taps=num_pre+num_post;
cur=num_pre+1;
%txffe_cell combines all the txffe values into a single cell array
%It is ordered: [precursorN precursorN-1 ... precursor1 postcursor1 ... postcursorN-1 postcursorN]
txffe_cell=cell(1,num_taps);
for k=num_pre:-1:1
    idx=num_pre-k+1;
    this_tx_field=sprintf('tx_ffe_cm%d_values',k);
    txffe_cell{idx}=param.(this_tx_field);
end
for k=1:num_post
    idx=k+num_pre;
    this_tx_field=sprintf('tx_ffe_cp%d_values',k);
    txffe_cell{idx}=param.(this_tx_field);
end
%total number of txffe runs is the product of the lengths of each tap
txffe_lengths=cellfun('length',txffe_cell);
if isempty(txffe_cell)
    num_txffe_runs=1;
else
    num_txffe_runs=prod(txffe_lengths);
end
%txffe_sweep_indices are used in the LOCAL_SEARCH block
%any tap with length=1 can be ignored
%Also is statistically likely that taps with greater number of values
%will exceed the LOCAL SEARCH criteria, so searching those first is faster
txffe_sweep_indices=find(txffe_lengths>1);
[~,length_sort]=sort(txffe_lengths(txffe_sweep_indices),'descend');
txffe_sweep_indices=txffe_sweep_indices(length_sort);
num_txffe_sweep_indices=length(txffe_sweep_indices);

gdc_values = param.ctle_gdc_values;
Gffe_values = param.cursor_gain;
switch param.CTLE_type
    case 'CL93'
    case 'CL120d'
        g_DC_HP_values =param.g_DC_HP_values;
    case 'CL120e'
        f_HP_Z=param.f_HP_Z;
        f_HP_P=param.f_HP_P;
        
end
best_ctle = [];
best_FOM = -inf;
best_txffe = [];
delta_sbr = [];
PSD_results=[];
MMSE_results=[];
best_bmax=param.bmax;
%AJG021820
best_bmin=param.bmin;
h_J=[];
pxi=0;
if OP.DISPLAY_WINDOW
    hwaitbar=waitbar(0);
else
    fprintf('FOM search ');
end
FOM=0;
if ~OP.RxFFE
    Gffe_values=0;
end
param.ndfe_passed=param.ndfe;
old_loops=0;
new_loops=0;

%GDC Qual construction
gqual= param.gqual;
g2qual=param.g2qual;
if ~strcmp(param.CTLE_type,'CL120d')
    qual=ones(1,length(gdc_values));
else
    if isempty(gqual) && isempty(g2qual)
        qual=ones(length(g_DC_HP_values),length(gdc_values));
    else
        qual=zeros(length(g_DC_HP_values),length(gdc_values));
        
        %prepare gqual and g2qual
        [g2qual,si]=sort(g2qual,'descend');
        gqual=gqual(si,:);
        tmp=g2qual;
        g2qual=zeros(length(tmp),2);
        for kk=1:length(tmp)
            if kk==1
                g2qual(kk,:)=[tmp(kk)+eps tmp(kk)];
            else
                g2qual(kk,:)=[tmp(kk-1) tmp(kk)];
            end
            gqual(kk,:)=sort(gqual(kk,:),'descend');
        end
        
        %Qual Construction
        for jj=1:length(g_DC_HP_values)
            for ii=1:length(gdc_values)
                for kk=1:size(gqual,1)
                    if g_DC_HP_values(jj) >= g2qual(kk,2) && g_DC_HP_values(jj) < g2qual(kk,1)
                        if gdc_values(ii) >= gqual(kk,2) && gdc_values(ii) < gqual(kk,1)
                            qual(jj,ii)=1;
                            break;
                        end
                    end
                end
            end
        end
    end
end

progress_interval=0.025;
if do_C2M
    loop_count=[1 2];
    T_O=floor((param.T_O/1000)*param.samples_per_ui);
    T_O=max(0,T_O);
else
    loop_count=1;
    T_O=0;
end
switch param.CTLE_type
    case 'CL93'
        lf_indx=1;
    case 'CL120d'
        lf_indx=length(g_DC_HP_values);
    case 'CL120e'
        lf_indx=1;
end
runs=length(gdc_values)*lf_indx*length(Gffe_values)*num_txffe_runs;
if  OP.Optimize_loop_speed_up == 1
    OP.BinSize = 1e-4;
    OP.impulse_response_truncation_threshold = 1e-3;
end

%Used to speed up FFE by only performing circshift when necessary
pulse_ctle_circshift=[];
ctle_response_updated=1;

%Used to speed up get_xtlk_noise by pre-calculating all the phase shift exponentials
calc_exp_phase=0;

%calculate cur index and pre/post indices outside of the loop
cur_start=cur;
precursor_indices=[];
postcursor_indices=[];
auto_count_trigger=0;
for kv=1:num_taps
    if ~auto_count_trigger && length(txffe_cell{kv})==1 && txffe_cell{kv}==0
        %precursor values fill the beginning of the vector.  Any empty precursor means
        %cursor position must be subtracted by 1
        if kv<cur_start
            cur=cur-1;
        end
    else
        %non empty value:  add to precursor or postcursor indices depending on position
        %in the vector
        if kv<cur_start
            auto_count_trigger=1;
            precursor_indices=[precursor_indices kv];
        else
            auto_count_trigger=0;
            postcursor_indices=[postcursor_indices kv];
        end
    end
end
if ~isempty(postcursor_indices)
    postcursor_indices=postcursor_indices(1):postcursor_indices(end);
end

%Calculate the full grid matrix of all txffe combinations
if isempty(txffe_cell)
    TXFFE_grid=0;
    FULL_tx_index_vector=1;
else
    TXFFE_grid=Full_Grid_Matrix(txffe_cell);
    %Also calculate the full grid matrix for the index used in each txffe combination
    %(the index is used in the LOCAL SEARCH block)
    for k=1:num_taps
        txffe_index_cell{k}=1:txffe_lengths(k);
    end
    FULL_tx_index_vector=Full_Grid_Matrix(txffe_index_cell);
end

%pre-calculate cursor to save time
txffe_cursor_vector=1-sum(abs(TXFFE_grid),2);

%pre-calculate full txffe for each iteration to save time
precursor_matrix=TXFFE_grid(:,precursor_indices);
postcursor_matrix=TXFFE_grid(:,postcursor_indices);
txffe_matrix = [precursor_matrix txffe_cursor_vector  postcursor_matrix];

if OP.TDMODE
    uneq_field='uneq_pulse_response';
    ctle_field='ctle_pulse_response';
else
    uneq_field='uneq_imp_response';
    ctle_field='ctle_imp_response';
end

%Speed up search for max(sbr)
if OP.TDMODE
    [~,init_max]=max(chdata(1).uneq_pulse_response);
else
    [~,init_max]=max(filter(ones(1,param.samples_per_ui),1,chdata(1).uneq_imp_response));
end
UI_max_window=20;
start_max_idx=init_max-UI_max_window*param.samples_per_ui;
if start_max_idx<1
    start_max_idx=1;
end
end_max_idx=init_max+UI_max_window*param.samples_per_ui;
if end_max_idx>length(chdata(1).(uneq_field))
    end_max_idx=length(chdata(1).(uneq_field));
end

itick_skips=0;
itick_cases=0;
FOM_TRACKER(1:length(Gffe_values),1:length(gdc_values),1:lf_indx,1:size(TXFFE_grid,1),1:length(full_sample_range))=0;
for i=loop_count
    
    for Gffe_index=1:length(Gffe_values)
        param.current_ffegain=Gffe_values(Gffe_index);
        for ctle_index=1:length(gdc_values)
            g_dc = gdc_values(ctle_index);
            kacdc = 10^(g_dc/20);
            CTLE_fp1 = param.CTLE_fp1(ctle_index);
            CTLE_fp2 = param.CTLE_fp2(ctle_index);
            CTLE_fz = param.CTLE_fz(ctle_index);
            switch param.CTLE_type
                case 'CL93'
%
                case 'CL120d'
%
                case 'CL120e'
                    HP_Z = param.f_HP_Z(ctle_index);
                    HP_P = param.f_HP_P(ctle_index);
            end
              %% HF Boost
            ctle_gain = (kacdc + 1i*chdata(1).faxis/CTLE_fz) ./ ...
                ((1+1i*chdata(1).faxis/CTLE_fp1).*(1+1i*chdata(1).faxis/CTLE_fp2));
%% Mid Frequency Boost
            ctle_gain_xc = (kacdc + 1i*f_xc/CTLE_fz) ./ ...
                ((1+1i*f_xc/CTLE_fp1).*(1+1i*f_xc/CTLE_fp2)); % Bill Kirkland          
            for  g_LP_index=1:lf_indx
                
                %GDC Qual Check
                if qual(g_LP_index,ctle_index)==0
                    pxi=pxi+num_txffe_runs;
                    continue;
                end
                
                switch param.CTLE_type
                    case 'CL93'
                        H_low=1;
                        kacde_DC_low=1;
                    case 'CL120d'
                        g_DC_low = g_DC_HP_values(g_LP_index);
                        f_HP=param.f_HP(g_LP_index);
                        kacde_DC_low = 10^(g_DC_low/20);
                        H_low=(kacde_DC_low +  1i*chdata(1).faxis/f_HP)./(1 + 1i*chdata(1).faxis/f_HP);
                        H_low_xc = (kacde_DC_low +  1i*f_xc/f_HP)./(1 + 1i*f_xc/f_HP);% Bill Kirkland
                    case 'CL120e' % z1 has been adusted on read in
                        H_low=(1 +  1i*chdata(1).faxis/HP_Z)./(1 + 1i*chdata(1).faxis/HP_P);
                        H_low_xc=(1 +  1i*f_xc/HP_Z)./(1 + 1i*f_xc/HP_P); % Bill Kirkland
                end
                H_ctf=H_low.*ctle_gain;
                switch upper(OP.FFE_OPT_METHOD)
                    case 'WIENER-HOPF'
                                        %% Bill Kirkland
                        H_ctf_xc = H_low_xc.*ctle_gain_xc;
                        H_rx_ctle_xc = H_r_xc.*H_ctf_xc;
                        % use Fourier Transform pair for correlation as we have to
                        % take ifft of H_r anyways.
                        % onesided and two sided responses - tricky, tricky, tricky
                        Var_eta0 =  param.eta_0*f_xc(end)/1e9;
                        XC_rx_ctle = ifft (H_rx_ctle_xc.*conj(H_rx_ctle_xc),2*length(H_rx_ctle_xc),'symmetric');
                        Noise_XC = Var_eta0.*XC_rx_ctle(1:param.samples_per_ui:N_fft_by2);

                        if OP.Do_White_Noise
                            Noise_XC = Noise_XC(1);
                        end
                    otherwise
                        Noise_XC=[];            
                end


                
                if OP.INCLUDE_CTLE==1
                    for k=1:param.num_s4p_files
                        ir_peak = max(abs(chdata(k).(uneq_field)));
                        ir_last  = find(abs(chdata(k).(uneq_field))>ir_peak*OP.impulse_response_truncation_threshold, 1, 'last');
                        chdata(k).(uneq_field) = chdata(k).(uneq_field)(1:ir_last);
                        chdata(k).(ctle_field) = TD_CTLE(chdata(k).(uneq_field), baud_rate ...
                            , CTLE_fz, CTLE_fp1, CTLE_fp2, g_dc, param.samples_per_ui);
                        switch param.CTLE_type
                            case 'CL93'
                            case 'CL120d'
                                chdata(k).(ctle_field) = TD_CTLE(chdata(k).(ctle_field), baud_rate, f_HP, f_HP,100e100 , g_DC_low , param.samples_per_ui);
                            case 'CL120e' % z1 has been adusted on read in
                                chdata(k).(ctle_field) = TD_CTLE(chdata(k).(ctle_field), baud_rate, HP_Z,HP_P,100e100 , 0 , param.samples_per_ui);
                        end
                    end
                    %set the flag to show ctle response was updated
                    ctle_response_updated=1;
                else
                    for k=1:param.num_s4p_files
                        chdata(k).(ctle_field) = chdata(k).(uneq_field);
                    end
                end
                for k=1:param.num_s4p_files
                    chdata(k).sdd21ctf=chdata(k).sdd21.*H_ctf; % sdd21 is a VTF, includes H_t, H_f, and package
                end
                %% Equation 93A-22 %%
                %         figure(1000)
                %         semilogx(chdata(1).faxis/1e9,db(H_ctf))
                %         hold on
                if OP.RX_CALIBRATION
                    ctle_gain2 = (kacdc + 1i*chdata(2).faxis/CTLE_fz) ./ ...
                        ((1+1i*chdata(2).faxis/CTLE_fp1).*(1+1i*chdata(2).faxis/CTLE_fp2));
                    switch param.CTLE_type
                        case 'CL93'
                            H_low2=1;
                        case 'CL120d'
                            g_DC_low = g_DC_HP_values(g_LP_index);
                            f_HP=param.f_HP(g_LP_index);
                            kacde_DC_low = 10^(g_DC_low/20);
                            H_low2=(kacde_DC_low +  1i*chdata(1).faxis/f_HP)./(1 + 1i*chdata(1).faxis/f_HP);
                        case 'CL120e' % z1 has been adusted on read in
                            H_low2=(1 +  1i*chdata(1).faxis/HP_Z)./(1 + 1i*chdata(1).faxis/HP_P);
                    end
                    H_ctf2=H_low2.*ctle_gain2;
                end
                % RIM 11-30-2020 moved to a subfunction
                [sigma_N] = get_sigma_eta_ACCM_noise(chdata,param,H_sy,H_r,H_ctf);
                if OP.RX_CALIBRATION
                    sigma_ne = get_sigma_noise( H_ctf2, param, chdata, sigma_bn); %% Equation 93A-48 %%
                    sigma_NEXT=sqrt(param.eta_0*sum( abs(H_sy(2:end).^2 .* H_r(2:end).*2 .* H_ctf(2:end).^2 ) .* diff(chdata(1).faxis)/1e9));% changed from /chdata(1).faxis(end) B. Kirkland S. Elnagar 11/6/2021
                else
                    %% Equations 93A-33 and 93A-34  for NEXT - independent of TXFFE setting %%
                    sigma_NEXT =  get_xtlk_noise( [0 1 0], 'NEXT', param, chdata );
                    sigma_ne=0;
                end
                
                if param.GDC_MIN ~= 0 && gdc_values(ctle_index) + g_DC_HP_values(g_LP_index) > param.GDC_MIN
                    pxi=pxi+num_txffe_runs;
                    continue; % change per 0.3k draft 2.3
                end
                %%
                PSD_results=[];
                if strcmp(OP.FFE_OPT_METHOD,'MMSE') && OP.RxFFE
                    OP.WO_TXFFE=1;
                    PSD_results=get_PSDs(PSD_results,[],[],[],gdc_values(ctle_index),g_DC_low,param,chdata,OP);
                end
                %TXFFE Loop
                %Originally this was a separate for loop for each tap, but it is now all contained in the TXFFE_grid matrix to use a single modular loop
                for TK=1:size(TXFFE_grid,1)

                    pxi=pxi+1;
                    progress = pxi/runs;
                    if OP.DISPLAY_WINDOW
                        if ~mod(pxi,floor(runs*progress_interval))
                            waitbar(progress, hwaitbar, 'Linear equalization tuning'); figure(hwaitbar); drawnow;
                        end
                    else
                        if ~mod(pxi,floor(runs*progress_interval)), fprintf('%i%% ', round(progress*100) );end
                    end
                    
                    %get the cursor for this iteration
                    txffe_cur=txffe_cursor_vector(TK);
                    
                    % Skip combinations with small values of c(0), not guaranteed to be supported by all transmitters.
                    if txffe_cur<param.tx_ffe_c0_min
                        continue;
                    end
                    old_loops=old_loops+1;
                    
                    %get the index used for each tap on this iteration
                    %this is needed for the LOCAL SEARCH block
                    tx_index_vector=FULL_tx_index_vector(TK,:);
                    
                    %Original LOCAL SEARCH Block:
                    %Keeping this one as commented code because it is a bit more readable than the Modular Block below
                    %But unlike the Modular Block, this one does not work if additional TXFFE taps are added
%                     % speedup "local search" heuristic - Adee Ran 03-17-2020
%                     % skip configurations more than
%                     % 2 steps away from current "best" point on any grid direction
%                     %  Matt Brown 11/19/2021 for cp2 and cp3
%                     if param.LOCAL_SEARCH>0 && ~isinf(best_FOM) && ...
%                                ((k_cp2>1 && length(cp3_values)>1 && abs(k_cp3-find(cp3_values==best_txffe(cur+3)))>param.LOCAL_SEARCH) ...
%                             || (k_cp1>1 && length(cp2_values)>1 && abs(k_cp2-find(cp2_values==best_txffe(cur+2)))>param.LOCAL_SEARCH) ...
%                             || (k_cm1>1 && length(cp1_values)>1 && abs(k_cp1-find(cp1_values==best_txffe(cur+1)))>param.LOCAL_SEARCH) ...
%                             || (k_cm2>1 && length(cm1_values)>1 && abs(k_cm1-find(cm1_values==best_txffe(cur-1)))>param.LOCAL_SEARCH) ...
%                             || (k_cm3>1 && length(cm2_values)>1 && abs(k_cm2-find(cm2_values==best_txffe(cur-2)))>param.LOCAL_SEARCH) ...
%                             || (k_cm4>1 && length(cm3_values)>1 && abs(k_cm3-find(cm3_values==best_txffe(cur-3)))>param.LOCAL_SEARCH) ...
%                             || (g_LP_index>1 && length(cm4_values)>1 && abs(k_cm4-find(cm4_values==best_txffe(cur-4)))>param.LOCAL_SEARCH) ...
%                             || (ctle_index>1 && abs(g_LP_index-best_G_high_pass)>param.LOCAL_SEARCH))
% 
%                         continue;
%                     end

                    %Modular LOCAL_SEARCH block:
                    % speedup "local search" heuristic - Adee Ran 03-17-2020
                    % skip configurations more than 2 steps away from current "best" point on any grid direction
                    skip_it=0;
                    if param.LOCAL_SEARCH>0 && ~isinf(best_FOM)
                        %instead of looping across all taps, only loop across
                        %those with length>1 (txffe_sweep_indices).
                        %It saves time since this block is encountered so often
                        for kj=1:num_txffe_sweep_indices
                            kv=txffe_sweep_indices(kj);
                            if kv==1
                                previous_loop_val=g_LP_index;
                            else
                                previous_loop_val=tx_index_vector(kv-1);
                            end
                            if previous_loop_val>1
                                best_index_this_tap=best_txffe_index(kv);
                                if abs(tx_index_vector(kv)-best_index_this_tap)>param.LOCAL_SEARCH
                                    skip_it=1;
                                    break;
                                end
                            end
                        end
                        
                        if ~skip_it && ctle_index>1 && abs(g_LP_index-best_G_high_pass)>param.LOCAL_SEARCH
                            skip_it=1;
                        end
                    end
                    if skip_it
                        continue;
                    end
                    %End Modular LOCAL SEARCH block
                    
                    new_loops=new_loops+1;
                    
                    %fetch txffe for this iteration
                    txffe=txffe_matrix(TK,:);
                    
                    %The phase shift exponentials used in get_xtlk_noise are independent of
                    %everything except number of taps and cursor position
                    %So it can be calculated 1 time here to avoid thousands of re-calcs
                    if ~calc_exp_phase
                        calc_exp_phase=1;
                        for k=1:length(txffe)
                            phase_memory(:,k)=exp(-1j*2*pi*(k-cur).*f/param.fb);
                        end
                        if OP.RxFFE
                            for k=-1*param.RxFFE_cmx:param.RxFFE_cpx
                                phase_memoryRXFFE(:,k+param.RxFFE_cmx+1)=exp(-1j*2*pi*(k+1).*f/param.fb);
                            end
                            phase_memory=[phase_memory phase_memoryRXFFE];
                        end
                    end
                    
                    %% Unequalized Pulse Reponse & circshift for FFE
                    %Perform circshift for FFE only when CTLE is updated.  The FFE_Fast function takes
                    %in the pre-shifted matrix
                    if ctle_response_updated
                        ctle_response_updated=0;
                        num_pre=cur-1;
                        %Another speed up:  the unequalized pulse is also only unique for each CTLE update
                        %Calculating here reduces number of convolutions by thousands
                        if OP.TDMODE
                            pulse_ctle=chdata(1).(ctle_field)(:);
                        else
                            %uneq_pulse=filter(ones(param.samples_per_ui, 1), 1, chdata(1).(ctle_field)(:));
                            %"conv2" is faster than filter. Just need to chop off extra points at the end
                            pulse_ctle=conv2(chdata(1).(ctle_field)(:),ones(param.samples_per_ui, 1));
                            pulse_ctle=pulse_ctle(1:end-param.samples_per_ui+1);
                        end
                        for k=1:length(txffe)
                            pulse_ctle_circshift(:,k)=circshift(pulse_ctle,[(k-1-num_pre)*param.samples_per_ui 0]);
                        end
                    end
                    
                    %% Apply TXFFE to pre-shifted pulse response
                    %[sbr] = FFE( txffe, cur-1, param.samples_per_ui, uneq_pulse);
                    sbr=FFE_Fast(txffe,pulse_ctle_circshift);
                    sbr_from_txffe=sbr;
                    sbr1=sbr;
                    
                    
                    %% Find Sample Location
                    % If RXFFE is included, the sample location will be found again below
                    [cursor_i,no_zero_crossing,sbr_peak_i,zxi]=cursor_sample_index(sbr,param,OP,start_max_idx:end_max_idx);
                    if param.ts_anchor==0
                        %keep MM
                    elseif param.ts_anchor==1
                        %peak sample
                        cursor_i=sbr_peak_i;
                        no_zero_crossing=0;
                    elseif param.ts_anchor==2
                        %max DV
                        possible_cursor=sbr(sbr_peak_i-param.samples_per_ui:sbr_peak_i+param.samples_per_ui);
                        possible_precursor=sbr(sbr_peak_i-2*param.samples_per_ui:sbr_peak_i);
                        [max_diff,d_idx]=max(possible_cursor-possible_precursor);
                        cursor_i=sbr_peak_i-param.samples_per_ui+d_idx-1;
                        no_zero_crossing=0;
                    else
                        error('ts_anchor parameter must be 0, 1, or 2');
                    end
                    if no_zero_crossing
                        continue;
                    end
                    raw_cursor_i=cursor_i;
                    
                    %%%%%%%%%%
                    %%%%%%%%%%
                    %%%%%%%%%%
                    %NEW ITICK LOOP (not indenting everything yet)
                    [~,si]=sort(abs(full_sample_range));
                    best_positive_itick_FOM=-inf;
                    best_negative_itick_FOM=-inf;
                    best_positive_itick_in_loop=[];
                    best_negative_itick_in_loop=[];
                    best_itick_FOM=-inf;
                    best_itick_in_cluster=[];
                    best_cluster=[];
                    
                    %box_search:  take the middle of 5 point windows, then use the best of those to search the rest of that window
                    %middle_search:  start from itick=0 and work outwards in negative and positive direction. stop searching when FOM starts to go down (using LOCAL SEARCH) 
                    box_search=0;
                    middle_search=1;
                    
                    if box_search
                        box_size=5;
                        box_mid=floor(box_size/2);
                        cluster=full_sample_range(1)+box_mid:box_size:full_sample_range(end);
                        CL=length(cluster);
                        loop_range=1:CL+box_mid*2;
                    elseif middle_search
                        loop_range=si;
                    else
                        loop_range=1:length(full_sample_range);
                    end
                    
                    for itickn=loop_range
                        if box_search
                            if itickn<=CL
                                itick=cluster(itickn);
                            else
                                if itickn==CL+1
                                    best_cluster=setdiff([best_itick_in_cluster-box_mid:best_itick_in_cluster+box_mid],best_itick_in_cluster);
                                end
                                if isempty(best_cluster)
                                    continue;
                                end
                                itick=best_cluster(itickn-CL);
                            end
                        else
                            itick=full_sample_range(itickn); 
                        end
                    
                    itick_cases=itick_cases+1;
                    
                    sbr=sbr_from_txffe;
                    cursor_i = raw_cursor_i+itick;
                    
                    %Local Search for +/- itick sweep
                    if middle_search && param.LOCAL_SEARCH>0
                        if itick>=0 && ~isinf(best_positive_itick_FOM) && abs(best_positive_itick_in_loop-itick)>=param.LOCAL_SEARCH
                            itick_skips=itick_skips+1;
                            continue;
                        end
                        if itick<=0 && ~isinf(best_negative_itick_FOM) && abs(best_negative_itick_in_loop-itick)>=param.LOCAL_SEARCH
                            itick_skips=itick_skips+1;
                            continue;
                        end
                    end
                    
                    triple_transit_time = round(sbr_peak_i*2/param.samples_per_ui)+20;
                    if min_number_of_UI_in_response < triple_transit_time
                        min_number_of_UI_in_response = triple_transit_time;
                    end
                    
                    cursor = sbr(cursor_i);
                    
                    %% RXFFE
                    if OP.RxFFE
                        % [ sbr, C]=force(sbr,param,OP,cursor_i);
                        %[ sbr, C]=force(sbr,param,OP,zxi+param.samples_per_ui);
                        %if isrow(sbr), sbr=sbr';end
                        
                        %AJG:  do not return sbr here (run time improvement)
                        %UPDATE:  use cursor_i in RXFFE instead of zero crossing + 1 UI
                        %[ ~, C]=force(sbr,param,OP,zxi+param.samples_per_ui,[],0);
                        % [ ~, C, floating_tap_locations]=force(sbr,param,OP,cursor_i,[],0);
                        % sbr at this point include the current setting
                        %                under consideration of txffe h21 ctf and fr
                        switch upper(OP.FFE_OPT_METHOD)
                            case 'MMSE'
                                OP.WO_TXFFE=0;
                                PSD_results=get_PSDs(PSD_results,sbr,cursor_i,txffe,gdc_values(ctle_index),g_DC_low,param,chdata,OP);
                                S_n=PSD_results.S_rn+PSD_results.S_tn+PSD_results.S_xn+PSD_results.S_jn;
                                if 0 % for debug
                                    plot(PSD_results.fvec(1:param.num_ui_RXFF_noise)/param.fb,10*log10(PSD_results.S_rn*1000/100) ,'disp','Srn')
                                    hold on
                                    plot(PSD_results.fvec(1:param.num_ui_RXFF_noise)/param.fb,10*log10(PSD_results.S_xn*1000/100) ,'disp','Sxn')
                                    plot(PSD_results.fvec(1:param.num_ui_RXFF_noise)/param.fb,10*log10(PSD_results.S_tn*1000/100) ,'disp','Stn')
                                    plot(PSD_results.fvec(1:param.num_ui_RXFF_noise)/param.fb,10*log10(PSD_results.S_jn*1000/100) ,'disp','Sjn')
                                    plot(PSD_results.fvec(1:param.num_ui_RXFF_noise)/param.fb,10*log10(PSD_results.S_n*1000/100) ,'disp','Sn')
                                     xlim([0 0.5])
                                    % ylim([-190 -160])
                                    set(gcf,'defaulttextinterpreter','none')
                                    xlabel('Normalized Frequency')
                                    ylabel('PSD dBm/Hz')
                                    hold on
                                    grid on
                                    legend show
                                    title('PSD')
                                end
                                MMSE_results = MMSE(PSD_results,sbr,cursor_i, param, OP ) ;
                                % floating_tap_locations=MMSE_results.MLSE_results;
                                C=MMSE_results.C;
                                FOM=MMSE_results.FOM;
                                floating_tap_locations=MMSE_results.floating_tap_locations;
                            otherwise
                                [ ~, C, floating_tap_locations]=force(sbr,param,OP,cursor_i,[],0, chdata, txffe, Noise_XC);
                        end
                        %Now there is a stand alone function for determining if RXFFE taps are illegal
                        %This is because the "force" function will also do a legality check when "backoff" is enabled
                        if ~(strcmp(OP.FFE_OPT_METHOD,'MMSE') && OP.RxFFE)
                            if RXFFE_Illegal(C,param)
                                continue;
                            end
                        end
                        %AJG:  speed up:  calculate sbr after checks for illegal taps (many tap combinations are illegal and "FFE" is time consuming)
                        sbr=FFE(C,param.RxFFE_cmx,param.samples_per_ui,sbr);
                        if isrow(sbr), sbr=sbr';end
                        
                        %% second guess at cursor location (t_s)  - based on approximate zero crossing
                        %This entire block is now inside the "if OP.RxFFE" to avoid unnecessary re-calculation
                        %UPDATE:  NOT RESAMPLING AFTER RXFFE
%                         [cursor_i,no_zero_crossing,sbr_peak_i]=cursor_sample_index(sbr,param,OP,start_max_idx:end_max_idx);
%                         if no_zero_crossing
%                             continue;
%                         end
                        
                        cursor = sbr(cursor_i);
                    end
                    A_p=sbr(sbr_peak_i);
                    %% 93A.1.6 step c defines A_s %%
                    A_s = param.R_LM*cursor/(param.levels-1);
                    if isempty(delta_sbr)
                        delta_sbr = sbr;
                    end
                    sbr=sbr(:);
                    %% Equation 93A-27 "otherwise" case %% param.N_bmax is param.ndfe if groups are not used
                    
                    if(param.Floating_DFE), param.ndfe=param.N_bmax; end
                    far_cursors = sbr(cursor_i-T_O+param.samples_per_ui*(param.ndfe+1):param.samples_per_ui:end);
                    t=((cursor_i+param.samples_per_ui*(param.ndfe+1):param.samples_per_ui:length(sbr))-(cursor_i+param.samples_per_ui*(param.ndfe+1)))*...
                        param.ui/param.samples_per_ui;
                    precursors = sbr(cursor_i-param.samples_per_ui:-param.samples_per_ui:1);
                    precursors = precursors(end:-1:1);
                    
                    %                                 % Error message if the sbr is not long enough for the specified range of Nb
                    %                                 if length(sbr) < cursor_i+param.samples_per_ui*(param.ndfe+1)
                    %                                     close(hwaitbar);
                    %                                     error('Pulse Response contains %d samples after the cursor. Specified Nb requires %d samples after the cursor.' ...
                    %                                         , length(sbr)-cursor_i, param.samples_per_ui*(param.ndfe+1));
                    %                                 end
                    
                    
                    
                    %% skip this case if FOM has no chance of beating old FOM
                    %this is also done below but with excess_dfe_cursors included.
                    %excess_dfe_cursors requires the floating DFE computation which is
                    %time consuming, so checking here can have significant run time improvements
                    sigma_ISI_ignoreDFE = param.sigma_X*norm([precursors;  far_cursors]);
                    if ~(strcmp(OP.FFE_OPT_METHOD,'MMSE') && OP.RxFFE)
                        if (20*log10(A_s/sigma_ISI_ignoreDFE) < best_FOM)
                            continue
                        end
                    end
                    
                    %% Equation 93A-27, when 1<=n<=N_b
                    %required length = cursor + all DFE UI + 1 additional UI
                    sbr_required_length=cursor_i+param.samples_per_ui*(param.ndfe+1);
                    if length(sbr)<sbr_required_length
                        sbr(end+1:sbr_required_length)=0;
                    end
                    dfecursors=sbr(cursor_i+param.samples_per_ui*(1):param.samples_per_ui:cursor_i+param.samples_per_ui*(param.ndfe));
                    if param.dfe_delta ~= 0
                        dfecursors_q=floor(abs(dfecursors/sbr(cursor_i))./param.dfe_delta).*param.dfe_delta.*sign(dfecursors)*sbr(cursor_i);
                        
                    else
                        dfecursors_q=dfecursors;
                    end
                    if param.Floating_DFE
                        %% floating taps
                        postcurors= sbr(cursor_i+param.samples_per_ui:param.samples_per_ui:end);
                        
                        [floating_tap_locations, floating_tap_coef, hisi, bmax]= floatingDFE( postcurors ,param.ndfe_passed,param.N_bf,param.N_bg,param.N_bmax, param.bmaxg, sbr(cursor_i), param.dfe_delta  );
                        
                        newbmax= [ param.bmax bmax(param.ndfe_passed+1:param.N_bmax)].';
                        param.use_bmax=newbmax;
                        %AJG021820
                        param.use_bmin=[param.bmin bmax(param.ndfe_passed+1:param.N_bmax)*-1].';
                    else
                        param.use_bmax=param.bmax;
                        %AJG021820
                        param.use_bmin=param.bmin;
                    end
                    
                    %AJG021820
                    actual_dfecursors=dfe_clipper(dfecursors_q,sbr(cursor_i)*param.use_bmax(:),sbr(cursor_i)*param.use_bmin(:));
                    if do_C2M
                        dfecursors_windowed=sbr(cursor_i-T_O+param.samples_per_ui*(1):param.samples_per_ui:cursor_i+param.samples_per_ui*(param.ndfe)-T_O);
                        % readjust SBR
                        if 0
                            %PR_DFE_center not currently used, so this is in "if 0" statement
                            PR_DFE_center=sbr;
                            for n=1:param.ndfe
                                %                                                 for ix=-param.samples_per_ui/2: param.samples_per_ui/2
                                %                                                     i_sample=ix+n*param.samples_per_ui+cursor_i;
                                %                                                     dper=sbr(i_sample)- actual_dfecursors(n);
                                %                                                     PR_DFE_center(i_sample)=dper;
                                %                                                 end
                                i_sample=(-param.samples_per_ui/2: param.samples_per_ui/2)+n*param.samples_per_ui+cursor_i;
                                PR_DFE_center(i_sample)=sbr(i_sample)-actual_dfecursors(n);
                            end
                        end
                        excess_dfe_cursors=dfecursors_windowed-actual_dfecursors;
                    else
                        excess_dfe_cursors=dfecursors-actual_dfecursors;
                    end
                    dfetaps=actual_dfecursors/sbr(cursor_i);
                    
                    if length(dfetaps) >= param.N_tail_start && param.N_tail_start ~=0
                        tail_RSS=norm(dfetaps(param.N_tail_start:end));
                        if  tail_RSS ~= 0
                            if tail_RSS >= param.B_float_RSS_MAX
                                param.use_bmax(param.N_tail_start:end)= ...
                                    min(tail_RSS, param.B_float_RSS_MAX) * sign(dfetaps(param.N_tail_start:end)).*dfetaps(param.N_tail_start:end) /tail_RSS;
                                %AJG021820
                                param.use_bmin(param.N_tail_start:end)= ...
                                    min(tail_RSS, param.B_float_RSS_MAX) * -1 * sign(dfetaps(param.N_tail_start:end)).*dfetaps(param.N_tail_start:end) /tail_RSS;
                            end
                        end
                        
                        %AJG021820
                        actual_dfecursors=dfe_clipper(dfecursors_q,sbr(cursor_i)*param.use_bmax(:),sbr(cursor_i)*param.use_bmin(:));
                        if do_C2M
                            excess_dfe_cursors=dfecursors_windowed-actual_dfecursors;
                        else
                            excess_dfe_cursors=dfecursors-actual_dfecursors;
                        end
                        dfetaps=actual_dfecursors/sbr(cursor_i);
                        
                    else
                        tail_RSS=0;
                    end
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
                    if ~OP.SNR_TXwC0
                        %% Equation 93A-30 %%
                        % since A_s = param.R_LM*cursor/(param.levels-1), cursor=(param.levels-1)*A_s/param.R_LM
                        sigma_TX = (param.levels-1)*A_s/param.R_LM*10^(-param.SNR_TX/20);
                    else
                        sigma_TX = (param.levels-1)*A_s/txffe(cur)/param.R_LM*10^(-param.SNR_TX/20);% SNER_TX mod from Adee
                    end
                    %% Equation 93A-31 %%
                    sigma_ISI = param.sigma_X*norm([precursors; excess_dfe_cursors; far_cursors]);
                    ISI_N=param.sigma_X*norm( far_cursors);
                    %% break if FOM has no chance of beating old e
                    OP.exe_mode=1;
                    if ~(strcmp(OP.FFE_OPT_METHOD,'MMSE') && OP.RxFFE)
                        switch OP.EXE_MODE
                            case 0
                            case 1
                                if (20*log10(A_s/sigma_ISI) < best_FOM)
                                    continue
                                end
                            case 2
                                if (20*log10(A_s/sigma_ISI) < best_FOM)
                                    break
                                end
                        end
                    end
                    %% Equation 93A-32 %%
                    sigma_J = norm([param.A_DD param.sigma_RJ])*param.sigma_X*norm(h_J);
                    
                    %% Equations 93A-33 and 93A-34 for FEXT (depends on TXFFE setting) %%
                    if OP.RX_CALIBRATION
                        sigma_XT=0;
                    else
                        if ~OP.RxFFE
                            % sigma_FEXT =  get_xtlk_noise( 0, 'FEXT', param ,chdata);
                            % sigma_XT = norm([sigma_NEXT sigma_FEXT]);
                            [sigma_XT,~,~] =  get_xtlk_noise( txffe, 'FEXT', param ,chdata,phase_memory); %with three outputs, the sigma_XT includes both FEXT and NEXT zhilei huang 01/11/2019
                            %% Equation 93A-36 denominator (actually its sqrt)
                        else % John Ewen: 13/12/20018
                            [sigma_XT,~,~] =  get_xtlk_noise( txffe, 'FEXT', param ,chdata, phase_memory,C); %with three outputs, the sigma_XT includes both FEXT and NEXT zhilei huang 01/11/2019
                            %                                         sigma_FEXT_ffe=  get_xtlk_noise( 0, 'FEXT', param, chdata, C );
                            %                                         if ~OP.RxFFE
                            %                                             sigma_NEXT_ffe =  get_xtlk_noise(0, 'NEXT', param, chdata);
                            %                                         else
                            %                                             sigma_NEXT_ffe =  get_xtlk_noise(0, 'NEXT', param, chdata , C);
                            %                                         end
                            %                                         sigma_XT = norm([sigma_NEXT_ffe sigma_FEXT_ffe])*sqrt((param.levels^2-1)/(3*(param.levels-1)^2));
                        end
                    end
                    if ~(strcmp(OP.FFE_OPT_METHOD,'MMSE') && OP.RxFFE)
                        if OP.RxFFE % modify sigma_N with rx noise from the rx ffe
                            index_f2=find(chdata(1).faxis(:)>param.fb,1,'first');
                            if isempty(index_f2), index_f2=length(chdata(1).faxis);end
                            f=chdata(1).faxis;
                            H_Rx_FFE=zeros(1,length(f));
                            for ii=-param.RxFFE_cmx:param.RxFFE_cpx
                                %H_Rx_FFE=C(ii+param.RxFFE_cmx+1).*exp(-1j*2*pi*(ii+1).*f/param.fb)+H_Rx_FFE;
                                if C(ii+param.RxFFE_cmx+1)==0
                                    %speed up:  skip cases when rxffe=0
                                    continue;
                                end
                                if ii+1==0
                                    %speed up:  ii+1=0, so just scalar addition and avoid exp calc
                                    H_Rx_FFE = H_Rx_FFE + C(ii+param.RxFFE_cmx+1);
                                else
                                    H_Rx_FFE=C(ii+param.RxFFE_cmx+1).*transpose(phase_memory(:,ii+param.RxFFE_cmx+1+length(txffe)))+H_Rx_FFE;
                                end
                            end
                            sigma_N =sqrt(param.eta_0*sum(H_sy(2:end) .* abs(H_r(2:end) .* H_ctf(2:end) .*H_Rx_FFE(2:end)).^2 .* diff(chdata(1).faxis)/1e9)); % changed from /chdata(1).faxis(end) B. Kirkland S. Elnagar 11/6/2021
                        end
                    end
                    %% Equation 93A-36 (note log argument is voltage rather than power ratio)
                    total_noise_rms  = norm([sigma_ISI sigma_J sigma_XT sigma_N sigma_TX sigma_ne]);
                    if do_C2M
                        if param.Noise_Crest_Factor == 0
                            ber_q = sqrt(2)*erfcinv(2*param.specBER);
                        else
                            ber_q=param.Noise_Crest_Factor;
                        end
                        if OP.force_pdf_bin_size
                            delta_y = OP.BinSize;
                        else
                            delta_y = min(A_s/1000, OP.BinSize);
                        end
                        ne_noise_pdf = normal_dist(0, ber_q, delta_y);
                        cci_pdf = normal_dist(0, ber_q, delta_y);
                        chdata(1).eq_pulse_response=sbr;
                        tmp_result.t_s= cursor_i;
                        tmp_result.A_s=A_s;
                        EH_1st= 2*(A_s-erfcinv(param.specBER*2)*2/sqrt(2)*total_noise_rms);
                        if EH_1st <=  param.Min_VEO_Test/1000 -.001
                            %                                         sprintf( '%g.1 As ..  %g.1 EH\n',A_s*1000,EH_1st*1000)
                            continue
                        else
                            %                                          sprintf( ' OK before %g As ..  %g EH\n',A_s*1000,EH_1st*1000)
                        end
                        Struct_Noise.sigma_N=sigma_N;
                        Struct_Noise.sigma_TX=sigma_TX;
                        Struct_Noise.cci_pdf=cci_pdf;
                        Struct_Noise.ber_q=ber_q;
                        Struct_Noise.ne_noise_pdf=ne_noise_pdf;
                        [Left_EW,Right_EW,eye_contour,EH_T_C2M,EH_B_C2M]=COM_eye_width(chdata,delta_y,tmp_result,param,OP,Struct_Noise,1);
                        EH=EH_T_C2M-EH_B_C2M;
                        N_i=(A_s*2-EH)/2;
                        if EH <= param.Min_VEO_Test/1000
                            %                                         sprintf( 'After  As=%.1f ..  EH=%.1f  EH_1st=%.1f  \n',A_s*1000,EH*1000, EH_1st*1000)
                            continue
                        else
                            %                                         sprintf( '<strong> After  As=%.1f ..  EH=%.1f  EH_1st=%.1f </strong> \n',A_s*1000,EH*1000, EH_1st*1000)
                        end
                        if ~(strcmp(OP.FFE_OPT_METHOD,'MMSE') && OP.RxFFE) % MMSE defines its own FOM
                            FOM =20*log10(A_s/N_i);
                        end
                    else
                        if ~(strcmp(OP.FFE_OPT_METHOD,'MMSE') && OP.RxFFE) % MMSE defines its own FOM
                            FOM = 20*log10(A_s/total_noise_rms);
                        end
                        %                                             if strfind(param.CTLE_type,'CL120e')
                        %                                                 FOM = A_s*C(param.RxFFE_cmx+1)-total_noise_rms_ffe;
                    end
                    if 0 % for loop analysis
                        result.FOM_array(new_loops)=FOM;
                    end
                    
                    if FOM>best_itick_FOM
                        best_itick_FOM=FOM;
                        best_itick_in_cluster=itick;
                    end
                    
                    if itick>=0 && FOM>best_positive_itick_FOM
                        best_positive_itick_FOM=FOM;
                        best_positive_itick_in_loop=itick;
                    end
                    if itick<=0 && FOM>best_negative_itick_FOM
                        best_negative_itick_FOM=FOM;
                        best_negative_itick_in_loop=itick;
                    end
                    
                    itick_index=find(itick==full_sample_range);
                    FOM_TRACKER(Gffe_index,ctle_index,g_LP_index,TK,itick_index)=FOM;
                    
                    if (FOM > best_FOM)
                        best_current_ffegain=param.current_ffegain;
                        best_txffe = txffe;
                        %along with best_txffe, save the indices of the best_txffe
                        %(saves time in LOCAL SEARCH block)
                        best_txffe_index=tx_index_vector;
                        best_sbr = sbr;
                        best_ctle = ctle_index;
                        best_G_high_pass =g_LP_index;
                        best_FOM = FOM;
                        best_cursor_i = cursor_i;
                        best_itick = itick;                       
                        if ~OP.TDMODE
                            [ effective_channel ] = FFE( txffe , cur-1, param.samples_per_ui, chdata(1).ctle_imp_response );
                            best_IR=effective_channel;
                        end
                        best_sigma_N = sigma_N;
                        best_h_J = h_J;
                        best_A_s=A_s;
                        best_A_p=A_p;
                        best_ISI=ISI_N;
                        best_bmax=param.use_bmax;
                        %AJG021820
                        best_bmin=param.use_bmin;
                        best_tail_RSS=tail_RSS;
                        best_dfetaps=dfetaps;
                        if param.Floating_DFE
                            best_floating_tap_locations=floating_tap_locations;
                            best_floating_tap_coef=floating_tap_coef;
                        end
                        if param.Floating_RXFFE
                            best_floating_tap_locations=floating_tap_locations;
                            % best_floating_tap_coef=floating_tap_coef;
                        end                        
                        if OP.RxFFE
                            best_RxFFE=C;
                            best_PSD_results=PSD_results;
                            best_MMSE_results=MMSE_results;
                        end
                    end
                    end
                end
                
            end
        end
    end
    if do_C2M
        if  best_FOM == -inf
            param.Min_VEO_Test=0;
        else
            break
        end
    end
end
if 0
    fprintf('old loops = %d\n',old_loops);
    fprintf('new loops = %d\n',new_loops);
    display(sprintf('\n :loops = %g',pxi))
end

%turn this on to review if FOM changes sign more than once in an itick loop
if 0
    DIR_CHANGE={};
    for m=1:length(Gffe_values)
        for n=1:length(gdc_values)
            for k=1:lf_indx
                FOM_this_mat=squeeze(FOM_TRACKER(m,n,k,:,:));
                %x reveals if FOM on a particular row (locked txffe, moving itick) goes up or down
                %1 = goes up, -1=goes down
                x=sign(diff(FOM_this_mat')');
                %y = change in sign on x.  the location of a "2" is where FOM changes direction
                y=abs(diff(x'))';
                %the goal is the FOM only changes direction once. so count the occurences of the 2
                for j=1:size(FOM_this_mat,1)
                    z{j}=find(y(j,:)==2);
                end
                zL=cellfun('length',z);
                %return any row where FOM changed direction more than once
                DIR_CHANGE{j,k}=find(zL>1);
            end
        end
    end
    multi_direction_change=find(~cellfun('isempty',DIR_CHANGE))
end

if ~exist('best_cursor_i', 'var')% take last setting
    result.eq_failed=true;
    display('equalization failed')
    best_bmax=param.bmax;
    %AJG021820
    best_bmin=param.bmin;
    best_tail_RSS=0;
    best_current_ffegain=0;
    best_txffe = txffe;
    best_sbr = sbr;
    best_ctle = ctle_index;
    if OP.RxFFE
        best_PSD_results=PSD_results;
        best_MMSE_results=MMSE_results;
        best_RxFFE=C;
    end
    best_G_high_pass =g_LP_index;
    best_FOM = FOM;
    %if this block is reached, the last encountered EQ parameters are used
    %if it so happened that there was no zero crossing in the last encountered EQ set, then cursor_i will be empty
    %EQ search has failed, so it is not important to give an exact sample location, so just use the peak of the pulse
    if isempty(cursor_i)
        [~,cursor_i]=max(sbr);
    end
    best_cursor_i = cursor_i;
    best_itick = itick;
    if ~OP.TDMODE
        [ effective_channel ] = FFE( txffe , cur-1, param.samples_per_ui, chdata(1).ctle_imp_response );
        best_IR=effective_channel;
    end
    best_sigma_N = sigma_N;
    best_h_J = h_J;
    best_A_p=max(sbr);
    best_ISI=1;
    best_dfetaps= sbr( cursor_i+param.samples_per_ui*(1):param.samples_per_ui:cursor_i+param.samples_per_ui*(param.ndfe))/sbr(cursor_i) ;
    best_A_s= sbr( cursor_i);
    if param.Floating_DFE
        best_floating_tap_locations=[];
        best_floating_tap_coef=[];
    end
    if do_C2M
        return
    end
    %     return
else
    result.eq_failed=false; % RIM 12/30/2023
end

best_cursor = best_sbr(best_cursor_i);
% report during debug
PRin=filter(ones(param.samples_per_ui, 1),1, chdata(1).uneq_imp_response);
%If sbr was zero padded, then PRin needs to do so as well)
if length(PRin)<length(best_sbr)
    PRin(end+1:length(best_sbr))=0;
end
f=1e8:1e8:100e9;

H_bt=Bessel_Thomson_Filter(param,f,OP.Bessel_Thomson);
H_bw=Butterworth_Filter(param,f,OP.Butterworth);
H_RCos=Raised_Cosine_Filter(param,f,OP.Raised_Cosine);% experiment with RCos
% need to include H_RCos in noise and when computing the system ir for thru
% and crosstalk
H_r=H_bw.*H_bt.*H_RCos;

ctle_gain1 = (10^(gdc_values(best_ctle)/20) + 1i*f/param.CTLE_fz(best_ctle)) ./ ...
    ((1+1i*f/param.CTLE_fp1(best_ctle)).*(1+1i*f/param.CTLE_fp2(best_ctle)));

switch param.CTLE_type
    case 'CL93'
        H_low=1;
    case 'CL120d'
        H_low=(10^(param.g_DC_HP_values(best_G_high_pass)/20) +  1i*f/param.f_HP(best_G_high_pass))./(1 + 1i*f/param.f_HP(best_G_high_pass));
    case 'CL120e'
        H_low=(1 +  1i*f/param.f_HP_Z(best_ctle))./ (1 + 1i*f/param.f_HP_P(best_ctle));
end
ctle_gain=H_low.*ctle_gain1.*H_r;



%lsbr=length(sbr);
%use length of best_sbr in case zero padding was performed
%check "sbr_required_length" variable
lsbr=length(best_sbr);
t=0:param.ui/param.samples_per_ui:(lsbr-1)*param.ui/param.samples_per_ui;

sampled_best_sbr_precursors_t = (best_cursor_i/param.samples_per_ui:-1:1/param.samples_per_ui)*param.ui;
sampled_best_sbr_precursors_t = sampled_best_sbr_precursors_t(end:-1:2); % exclude cursor
sampled_best_sbr_precursors = best_sbr(round(sampled_best_sbr_precursors_t/param.ui*param.samples_per_ui));
sampled_best_sbr_postcursors_t = (best_cursor_i:param.samples_per_ui:lsbr)/param.samples_per_ui*param.ui;
sampled_best_sbr_postcursors_t = sampled_best_sbr_postcursors_t(2:end); % exclude cursor
sampled_best_sbr_postcursors = best_sbr(round(sampled_best_sbr_postcursors_t/param.ui*param.samples_per_ui));
sampled_best_sbr_dfecursors_t = (best_cursor_i/param.samples_per_ui+(1:param.ndfe_passed))*param.ui;
if param.Floating_DFE
    sampled_best_sbr_fdfecursors_t = (best_cursor_i/param.samples_per_ui+(best_floating_tap_locations))*param.ui;
end
% apply max tap value constraint
dfe_cursors =    sampled_best_sbr_postcursors(1:param.ndfe);
dfe_SBRcursors = sampled_best_sbr_postcursors(1:param.ndfe);
if  isrow(best_bmax) == 1, best_bmax=best_bmax.';end

%AJG021820
if  isrow(best_bmin) == 1, best_bmin=best_bmin.';end
DFE_taps_mV=dfe_clipper(dfe_cursors,best_cursor*best_bmax(1:param.ndfe),best_cursor*best_bmin(1:param.ndfe));
if param.Floating_DFE
    FDFE_taps_mV=DFE_taps_mV(best_floating_tap_locations);
end

sampled_best_sbr_postcursors(1:param.ndfe) = dfe_SBRcursors-DFE_taps_mV;
Symbol_Adj = (param.levels-1);% 3A.1.6
if OP.DEBUG ~=0
    if OP.DISPLAY_WINDOW && ~OP.RX_CALIBRATION
        % display pulse responses in one axis per test case.
        switch upper(OP.TIME_AXIS)
            case 'S' % RIM 11-13-2023 added user selectable xaxis
                xnorm=1;
                xaxis_label='seconds';
                offset=0;
            case 'UI'
                xnorm=param.ui;
                xaxis_label='UI';
                offset=t(best_cursor_i)/xnorm;
            otherwise
                xnorm=1;
                xaxis_label='seconds';
                offset=0;      
        end
        figure_name = sprintf('PKG %d: Equalization effect: %s : ',  OP.pkg_len_select( param.package_testcase_i), param.base);
        fig=findobj('Name', figure_name);
        if isempty(fig), fig=figure('Name', figure_name); end
        figure(fig);set(gcf,'Tag','COM');
        movegui(fig,'north')
        %figure(fig.Number);
        % hax = subplot(length(OP.pkg_len_select), 1, param.package_testcase_i);
        if OP.RxFFE
            ax1=subplot(2,1,1);
        end
        plot(t/xnorm-offset,best_sbr/Symbol_Adj,'disp', '1/2 Symbol 0-3 Equalized PR');
        hold on
        
        PRplt(1:param.samples_per_ui+5)=PRin(1); % line up with the ffe introduced delay
        PRplt(param.samples_per_ui+6:length(t))=PRin(1:length(t)-param.samples_per_ui-5);
        plot((t-param.ui-t(6))/xnorm-offset,PRplt/Symbol_Adj,'r','disp', '1/2 Symbol 0-3 Unequalized (tp5d) PR');
        stem(t(best_cursor_i)/xnorm-offset,best_sbr(best_cursor_i)/Symbol_Adj,'g','disp','Cursor (sample point)');
        title(sprintf('PKG Case %d',  OP.pkg_len_select( param.package_testcase_i)));
        ylabel('volts')
        xlabel(xaxis_label)
        grid on
        legend show
        legend( 'Location', 'best')
        plot((sampled_best_sbr_precursors_t-param.ui/param.samples_per_ui)/xnorm-offset, sampled_best_sbr_precursors'/Symbol_Adj, 'kx', 'disp','Pre cursors');
        plot((sampled_best_sbr_postcursors_t-param.ui/param.samples_per_ui)/xnorm-offset, sampled_best_sbr_postcursors'/Symbol_Adj, 'ko', 'disp','Post cursors');
        if param.ndfe_passed ~=0
            stem((sampled_best_sbr_dfecursors_t-param.ui/param.samples_per_ui)/xnorm-offset,DFE_taps_mV(1:param.ndfe_passed)/Symbol_Adj','m', 'LineWidth',2,'disp','DFE-canceled cursors');
        end
        if param.Floating_DFE
            stem((sampled_best_sbr_fdfecursors_t-param.ui/param.samples_per_ui)/xnorm-offset,FDFE_taps_mV/Symbol_Adj, 'MarkerFaceColor','red','MarkerEdgeColor','m','LineWidth',1,'disp','FDFE-canceled cursors');
        end
        if OP.RxFFE
            ax2=subplot(2,1,2);
            if param.Floating_RXFFE
                stem((t(best_cursor_i+(best_floating_tap_locations)*param.samples_per_ui))/xnorm-offset,best_RxFFE(best_floating_tap_locations)...
                    ,'filled','disp','RxFFE floating FFE taps')
                hold on
            end
            stem((t(best_cursor_i+param.samples_per_ui*(-param.RxFFE_cmx:param.RxFFE_cpx)))/xnorm-offset,best_RxFFE(1:param.RxFFE_cmx+param.RxFFE_cpx+1)...
                ,'filled','disp','RxFFE fixted FFE taps')
            legend show
            zoom xon
            linkaxes([ax1 ax2],'x')
        end


        grid on
        legend show
        legend( 'Location', 'best')
        zoom xon
        % set(hax, 'tag', 'EQE');
        %
        figure(110);set(gcf,'Tag','COM');
        set(gcf, 'Name', 'CTLE selection');
        movegui(gcf, 'southeast');
        semilogx(f,20*log10(abs(ctle_gain)), 'disp', sprintf('Case %d', param.package_testcase_i));
        hold on
        semilogx(f,20*log10(abs(H_r)), 'disp', sprintf('Rx filter Case %d', param.package_testcase_i));
        semilogx(f,20*log10(abs(H_low.*ctle_gain1)), 'disp', sprintf('CTF Case %d', param.package_testcase_i));
        fbaud_tick=find(f >= baud_rate, 1);
        fnq_tick=find(f >= baud_rate/2, 1);
        stem(f(fnq_tick),20*log10(abs(ctle_gain(fnq_tick))),'g', 'handlevisibility', 'off');
        stem(f(fbaud_tick),20*log10(abs(ctle_gain(fbaud_tick))),'g', 'handlevisibility', 'off');
        recolor_plots(gca);
        title('CTF/w Rx Filter Response')
        ylabel('dB')
        xlabel('Hz')
        legend show
    end
    display(['FOM:                ' ,num2str(best_FOM, 2),' dB']);
    display(['TXFFE coefficients: ' ,mat2str(best_txffe) ] );
    display(['SNR ISI:                ' ,num2str(20*log10(best_A_p/best_ISI), 2),' dB']);
    display(['CTLE DC gain:       ' ,num2str(gdc_values(best_ctle)), ' dB']);
    display(['CTF peaking gain:  ' ,num2str(20*log10(max(abs(ctle_gain))), 2), ' dB']);
    display(['Symbol Available signal:   ' ,num2str(best_cursor/Symbol_Adj)]);
end
if OP.DEBUG && OP.DISPLAY_WINDOW && OP.RX_CALIBRATION==0
    eqe_axes = findobj('tag', 'EQE');
    if ~isempty(eqe_axes), linkaxes(eqe_axes, 'xy'); end
end
if OP.DISPLAY_WINDOW
    close(hwaitbar);
else
    fprintf('\n');
end

% % eq_data
result.cur=cur;
result.txffe = best_txffe;
result.ctle = best_ctle;
result.best_G_high_pass=best_G_high_pass;
result.DFE_taps = best_dfetaps; %relative
result.DFE_taps_i = best_cursor_i+(1:param.ndfe)*param.samples_per_ui;
if param.Floating_DFE
    result.floating_tap_locations=best_floating_tap_locations;
    result.floating_tap_coef=best_floating_tap_coef;
end
if param.Floating_RXFFE
    result.floating_tap_locations=best_floating_tap_locations;
end
result.A_s = best_A_s;
result.t_s = best_cursor_i;
result.itick = best_itick;
result.sigma_N = best_sigma_N;
result.h_J = best_h_J;
result.FOM = best_FOM;
if ~OP.TDMODE
    %If sbr was zero padded, then best_IR needs to do so as well)
    if length(best_IR)<length(best_sbr)
        best_IR(end+1:length(best_sbr))=0;
    end
    result.IR = best_IR;
end
result.t=t;
result.sbr=best_sbr;
if OP.RxFFE
    result.RxFFE=best_RxFFE;
    result.PSD_results=best_PSD_results;
    result.MMSE_results=best_MMSE_results;
end



% changed RIM 4/17/2019 use sum(IR) for Vf of originl IR at N_b UI
% updated RIM 12/17/2021 
result.A_p = max(chdata(1).uneq_pulse_response);
its=find(chdata(1).uneq_pulse_response>=max(chdata(1).uneq_pulse_response),1,'first');
PR=chdata(1).uneq_pulse_response;
iend = its+param.N_v*param.samples_per_ui-param.samples_per_ui/2; % from eq: 163A-3
ibeg= its-param.D_p*param.samples_per_ui-param.samples_per_ui/2;% from eq: 163A-3
if iend >= length(PR)
    iend = length (PR);
end
if ibeg < 1
    ibeg = 1;
end
PR=PR(ibeg:iend);
result.A_f = sum(PR/param.samples_per_ui); %% eq 163A-3
SRn=PR;
for ik=1:floor(length(PR)/param.samples_per_ui)
    SPR=circshift(PR,param.samples_per_ui*ik);
    SPR(1:ik*param.samples_per_ui)=0;
    SRn=SRn+ SPR;
end
codedebug=0;
if codedebug
    fig=figure('Name', 'step and pulse response for code debug');
    figure(fig);set(gcf,'Tag','COM');
    UI=(1:length(SRn))/param.samples_per_ui-param.D_p;
    plot(UI,SRn)
    hold on
    plot(UI,PR)
    xlim([-param.D_p param.N_v])
    grid on;hold off;
    result.step=SRn;
end
i20=find(SRn>=0.20*result.A_f,1,'first');
i80=find(SRn>=0.80*result.A_f,1,'first');
result.Tr_measured_from_step=(i80-i20)/(param.fb*param.samples_per_ui);
result.Pmax_by_Vf=result.A_p/result.A_f;
result.ISI =best_ISI;
result.SNR_ISI=20*log10(best_A_p/best_ISI);
result.best_current_ffegain=best_current_ffegain;
result.best_bmax=best_bmax;
%AJG021820
result.best_bmin=best_bmin;
result.tail_RSS=best_tail_RSS;
function param=parameter_size_adjustment(param,OP)

make_length2={'C_pkg_board' 'C_diepad' 'L_comp' 'C_bump' 'tfx' 'C_v' 'C_0' 'C_1' 'pkg_Z_c' 'brd_Z_c' 'R_diepad'};
make_length_WCPORTZ={'a_thru' 'a_fext' 'a_next' 'SNDR'};
make_length_GDC={'CTLE_fp1' 'CTLE_fp2' 'CTLE_fz' 'f_HP_Z' 'f_HP_P'};
make_length_DCHP={'f_HP'};
make_length_ncases={'AC_CM_RMS'};

%ncases used by make_length_ncases fields
[ncases, mele]=size(param.z_p_tx_cases); % need find the number of test cases RIM 01-08-20

%PORTZ_mult used by make_length_WCPORTZ fields
pkg_sel_vec=ones(1,max(OP.pkg_len_select));
if OP.WC_PORTZ
    PORTZ_mult=[1 1];
else
    PORTZ_mult=pkg_sel_vec;
end

%Parameters that have length = 2
for j=1:length(make_length2)
    if numel(param.(make_length2{j}))==1
        param.(make_length2{j}) = param.(make_length2{j})*[1 1];
    end
end

%Parameters that have length = ncases
for j=1:length(make_length_ncases)
    if numel(param.(make_length_ncases{j}))==1
        param.(make_length_ncases{j}) = param.(make_length_ncases{j})*ones(1,ncases);
    end
end

%Parameters that have length = length(ctle_gdc_values)
for j=1:length(make_length_GDC)
    if numel(param.(make_length_GDC{j}))==1
        param.(make_length_GDC{j}) = param.(make_length_GDC{j})*ones(size(param.ctle_gdc_values));
    end
end

%Parameters that have length = length(g_DC_HP_values)
for j=1:length(make_length_DCHP)
    if numel(param.(make_length_DCHP{j}))==1
        param.(make_length_DCHP{j}) = param.(make_length_DCHP{j})*ones(size(param.g_DC_HP_values));
    end
end

%Parameters that have length associated with PORTZ_mult
for j=1:length(make_length_WCPORTZ)
    if numel(param.(make_length_WCPORTZ{j}))==1
        param.(make_length_WCPORTZ{j}) = param.(make_length_WCPORTZ{j})*PORTZ_mult;
    end
end
function sgm = pdf2sgm(pdf)
avg = sum(pdf.x .* pdf.y);
sgm = sqrt(sum((pdf.x - avg).^2 .* pdf.y));
% end yasuo patch


%% adding tx packgage
function cdf=pdf_to_cdf(pdf)

%Transform PDF to CDF
%The CDF is natively calculated from negative-to-positive voltage.
%This only gives BER calculation for bottom eye.  Need to also
%calculate a CDF of reversed PDF to get top eye.  The final CDF is the
%min of top and bottom CDF values.
%If only interested in one side, a simple cumsum on y is all that is needed.

cdf.yB=cumsum(pdf.y);
cdf.yT=fliplr(cumsum(fliplr(pdf.y)));
cdf.y=min([cdf.yB(:) cdf.yT(:)],[],2);
cdf.x=pdf.x;
function plot_bathtub_curves(hax, max_signal, sci_pdf, cci_pdf, isi_and_xtalk_pdf, noise_pdf,jitt_pdf, combined_interference_and_noise_pdf, bin_size)
cursors = d_cpdf(bin_size,max_signal*[-1 1], [1 1]/2);
signal_and_isi_pdf = conv_fct(cursors, sci_pdf);
signal_and_xtalk_pdf = conv_fct(cursors, cci_pdf);
signal_and_channel_noise_pdf = conv_fct(cursors, isi_and_xtalk_pdf);
signal_and_system_noise_pdf = conv_fct(cursors, noise_pdf);
signal_and_system_jitt_pdf = conv_fct(cursors, jitt_pdf);
signal_and_total_noise_pdf = conv_fct(cursors, combined_interference_and_noise_pdf);
%% Added by Bill Kirkland, June 14, 2017
cursors_l = cursors; cursors_l.y(cursors_l.x>0) = 0;
cursors_r = cursors; cursors_r.y(cursors_r.x<0) = 0;
signal_and_total_noise_pdf_l = conv_fct(cursors_l, combined_interference_and_noise_pdf);
signal_and_total_noise_pdf_r = conv_fct(cursors_r, combined_interference_and_noise_pdf);

semilogy(signal_and_isi_pdf.x, abs(cumsum(signal_and_isi_pdf.y)-0.5) ,'r','Disp','ISI', 'parent', hax)
hold on
semilogy(signal_and_xtalk_pdf.x, abs(cumsum(signal_and_xtalk_pdf.y)-0.5) ,'b','Disp','Xtalk', 'parent', hax)
semilogy(signal_and_channel_noise_pdf.x, abs(cumsum(signal_and_channel_noise_pdf.y)-0.5) ,'c','Disp','ISI+Xtalk', 'parent', hax)
semilogy(signal_and_system_noise_pdf.x, abs(cumsum(signal_and_system_noise_pdf.y)-0.5) ,'m','Disp','Jitter, SNR_TX,RL_M, eta_0 noise', 'parent', hax)
semilogy(signal_and_system_jitt_pdf.x, abs(cumsum(signal_and_system_jitt_pdf.y)-0.5) ,'g','Disp','Jitter noise', 'parent', hax)

%% Added by Bill Kirkland, June 14, 2017
% modification allows bathtub curves to cross over and hence one can
% directly read the noise component.
%semilogy(signal_and_total_noise_pdf.x, abs(cumsum(signal_and_total_noise_pdf.y)-0.5) ,'k','Disp','total noise PDF', 'parent', hax)
vbt_l = abs(0.5-cumsum(signal_and_total_noise_pdf_l.y));
vbt_r = fliplr(0.5-(cumsum(fliplr(signal_and_total_noise_pdf_r.y))));
semilogy(signal_and_total_noise_pdf_l.x, vbt_l ,'k','Disp','total noise PDF left', 'parent', hax)
semilogy(signal_and_total_noise_pdf_r.x, vbt_r ,'k','Disp','total noise PDF right', 'parent', hax)

hc=semilogy(max_signal*[-1 -1 1 1], [0.5 1e-20 1e-20 0.5], '--ok');
set(get(get(hc,'Annotation'),'LegendInformation'), 'IconDisplayStyle','off');

ylabel(hax, 'Probability')
xlabel(hax, 'volts')
legend(hax, 'show')
% testing code
if 0
    figure_name =  'COM curves';
    fig=findobj('Name', figure_name);
    if isempty(fig), fig=figure('Name', figure_name); end
    figure(fig);set(gcf,'Tag','COM');
    grid on
    semilogy(db(max_signal./(signal_and_total_noise_pdf_r.x)), 2*vbt_r ,'Disp','SNR CDF')
    hold on
    semilogy(db(max_signal./(max_signal-signal_and_total_noise_pdf_r.x)), 2*vbt_r ,'Disp','COM CDF')
    ylim([ 1e-6 0.25])
    xlim([0 30])
    grid on
end