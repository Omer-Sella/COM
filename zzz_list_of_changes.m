function zzz_list_of_changes
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
% added output_args.peak_uneq_pulse_mV
% added output_args.cable_loss when "Include PCB" is not 0 in the config file
% added: tap c(-2) c(2) and c(3)
% added: g_DC_HP and f_HP_PZ
% added: new value for "Include PCB" = 2 for cable Rx compliance test only Rx host added
% added: BREAD_CRUMBS is 1 then a mat file with the structures params and OP is created in the results directory
% added: T_r_filter_type for RITT testing when IDEAL_TX_TERM is 1:
% added T_r_meas_point for RITT, if 0, measurement was at tp0, if 1 measurement was tp0a
% 0 is for is for Gaussian filter and 1 is for a 4th order Bessel-Thomson filter
% fixed INCLUDE_CTLE=0 to really remove from computation
% r161a fixed matlab version issue when for OP.INCLUDE_CTLE=0
% r162 adjusting RITT rise time to Mike Dudek's recommendations also always enable risetime filter if T_r_filter_type=1
% r162 tx and rx package impedance {Zc)
% r162a Gaussian equation corrected
% r163 cast snr_tx with package test case
% r164 fix pdf for very low noise and lo pass filter enhancements
% r164 add zero gain at nqyist CTLE as in CL12e
% r165 add simpler congfig command called FORCE_TR (force risetime)
% r200 cm3 and cm4 added cp3 removed
% r200 fixed problem in s21_to_impulse_DC when s parameter have a DC entry
% r200 ILD_FOM updated to EQ93A-55  ERL adde
% r200 improved phase interpolation for return loss time conversion
% r200a db = @(x) 20*log10(abs(x)) added so sig processing toolbox won't be required
% r200b Fixed error in bifurcation of Tx/Rx Rd%
% r200c missed on fix for interpolation
% r210 new ERL with time gating function
% r224 update ERL with from D3.1
% r226  fix s2p reading problem
%    change SNR_ISI_XTK_normalized_1_sigma to SNR_ISI (with nulled noise)
%    Fix Rx calibration issue
%    added ERL limit and Nd
% r227 adding Pmax/Vf and peak of eq pulse, fixed issue with rx testing
%     INC_PACKAGE=0 not fully supported message
%     if N=0 use TDR_duration
%     red display text for fail ERL and COM
% r228 fixed ERL pass fail report, default Grr_limit to 1
% r230 add rx ffe
% r231 change crosstalk noise to icn like to speed things up
% r231 change default OP.impulse_response_truncation_threshold to 1e-3 from
% 1e-5mof-
% r232 fix default for Rx eq so old spead sheets work
% r234 fix inadvertent typo for clause 120e ctle and problem with TXffe loop time reduction
% r235 adding dfe quantization changed to normalized DFE taps reported
% r236 adding ffe gain loop and resample after RxFFE
% r240 added output for C2M and setting defaults for some FFE eq
% r241 force FFE main cursor to 1 and remove sum of taps = 1
% r250 adding more complex package
% r251 post cursor fix for DFE in force() and ffe backoff
% r251 remove TDR threshold noise filter
% r252 add rx FFE filter to receiver noise filter
% r252 change ICN in the xtk noise calculation to end at fb rather than fb/2
% r253 a few bug fixes in force from i indexing and for no ffe postcursors
% r254 precursor check fix in optimize_fom % mod fix in force
% r254 help to align columns in csv file
% r254 accept syntax for 2 tline flex package model
% r256 speed up optimize FOM
% r256 fix problem reading in config file from q/a
% r256 added code from Yasou Hidaka for reading in parameter an and printing out noise
% r257 fixed extrapolation of channel with lower bandwidths in s21_to_impulse_DC
% r257 in get_xtlk_noise in optimize_FOM: reomove crosstalk double counting and apply TXFFE is  FEXT
% r258 EXE_MODE switch 12/21 0:legacy 1:fast 2:superfast
% r258 CDR switch 'MM' or 'mod-MM'
% r258 correction for asymentirc tx/Rx packages
% r258 revamped display results display window
% r259 fix problem if Min_VEO is set in spreadsheet.
% r259 fix problem in optimize_FOM. get_xtlk_noise need to have 3 output
%      parameter else only FEXT is considered for FOM.zhilei huang 01/11/2019
% r259 putting COM_db and IL last in output to terminal
% r259 msgtext change to msg for C2C case other cases not vetted but not presently used
% r259 use N_bx for ERL rather than Nb (ndfe))
% r259 added TDR_W_TXPKG which performs TDR and ERL with the Tx package added
% r260 r259 used rd for the reciever to terminate the package. It was changed to the rd of the transmitter
% r260 used eta_0 PSD equation for sigma_n
% r260 fix IL graph legend to w/pkg and Tr
% r260 define tfx for each port
% r262 fx parameter passing parsing for mod_string revert to 2.57 no COM computational impact
% r262 Report estimate for DER for channel (Yasuo 2/30/19)
% r262 reset on exit default text interpreter to tex
% r262 localize run timer (John Buck 1/17/19)
% r262 set db as internal function in force to avoid tool box
% r262 changed loop for Grr and Gloss in get_tdr so that nbx and tfx works when beta x = 0
% r263 added to output_args RL structure and report "struct" in csv file
% r264 added EW estimate
% r266 using unequalized IR for Vf and Vf to compute ratio of Vp/Vf
% r267 added floating taps with param.N_bf, param.N_bg, param.N_bmax, param.bmaxg.  groups not used for ERL
% r268 added sequential/co-optimization switch for floating tap banks OP.FT_COOP default 0 i.e. sequential
% r269 changed  param.N_bmax to  param.N_f
% r270 implement JingBo Li's and Howard Heck's floating tap method
% r270 modification by Adam Healey for Ls and Cb termination (aka t-coil emulation)
% r270 added c_0 and c_1 for CA in add_brd
% r272 fixed version syntax problem in output_args RL report
% r272 fixed eye width computation problem crosstalk was missed in pervious versions
% r272 removed eye width report if doing a Rx calibration
% r273 better alignment and control for ICN reporting
% r273 fixed PSXTK graph
% r275 fixed delay adjustment for ERL/TDR in get_TDR (Adam Healey 09/06/2019)
% r276 go back to reporting channel IL results (output_args.IL_dB_channel_only_at_Fnq) with board added read_s4p_files (as in r270)
% r276 chdata(i).Aicn=param.a_icn_fext should have been chdata(i).Aicn=param.a_icn_next for the next selection. Since in most spec's they are the same there is little no impact in results
% r276 test for output_args for isfield(chdata(1),'sdd22_raw')
% r276 change divisor for ICN and FOM_ILD to param.f2 from param.fb, may raise ICN and ILD value reported in r275
% r276 C_1 was instantiated as C_0. This was fixed
% r276 fixed rounding problem in reporting of loss at f_nq
% r276 power limit (RSS) for tail DFE taps (B_float_RSS_MAX, N_tail_start)
% r277 added nv for deterining steady state voltage for fitting compatibility
% r278 added b_min to support asymmetric bmax
% r278 added kappa1 and kappa2 to scale package to channel reflection for ERL experiments
% r278 added keyword OP.SHOW_BRD which includes added board in TDR and ERL
% r292 speed up for FOM search (Adee Ran) implemented by Adam Gregory.
% r292 param.LOCAL_SEARCH set to is the heuristic step distance keyword is 'Local Search'
% r292 fixing TDR for different impedance references in get_TDR and s2p file compatibility
% r292 eq. 93A-19 and 93-20 code implementation  bug when include .3by change% to fix edge rate equation 93A-46 (h_T). no effect if Rd=50 or IL > 5 dB
% r292 H_t implemented  in s21_pkg
% r292 plot and report for die to die IL remove the Tr effect "IL with pkgs & Tr filter" goes to "IL with pkgs"
% r292 add GDC_MIN to optimize_FOM
% r293 fix if ndfe-0 and ERL only and s2p issue
% r293a investigate the Tukey filtering
% r293a if fix if bmin is missing
% r294 fix problems reading s2p files for ERL computation
% r294 align Tukey_Window with .3ck definition for ERL and TDR computations
% r294 add parameter param.Noise_Crest_Factor. Default is not to use
% r294 add gdc and gdc2 range limitations
% r295 add VEC Pass threshold
% r295 removed close force all. Tagged all figures with "COM"
% r295 consolidated print in new function "end_display_control"
% r295 report pre/pmax for Txffe
% r295 speed up test cases by not re-reading in s4p files
% r297 add  provisions for AC_CM_RMS for through CM (experimental)
% r299 add keyword T_O (param.T_O) and  samples_for_C2M (parsm.samples_for_C2M) for new C2M VEC and EH computations
% r310 refine VEC and EH for C2M from Adam Gregory in
% r315 added keyword for Bessel_Thomson and Butterworth(default) filter. 
%      cdf_to_ber_contour,COM_eye_width,combine_pdf_same_voltage_axis,
%      optimize_fom_for_C2M. pdf_to_cdf, conv_fct_MeanNotZero, and get_pdf_full
% r311 added RILN
% r314 when T_O is not zero 3 eyes are used to compute VEC and VEO
% r315 Bessel_Thomson keyword is added mostly for measuring  Pmax, Vf, and SNDR
% r316 remove DC computation for RX Calibration loops
% r317 for SAVE_TD to include EQ and unEQ FIR
% r317 clean up bessel thomson and butterworth filter logic for ERL and normal COM 
% r318 if min_VEO_test fails to find a solution the loop is restarted with min_VEO_test to near zero. Makes sure COM returns results 
% r320 fixed RX_CALIBRATION which was broken in r310
% r320 speed up for C2M by moving managing optimize loop distribution of computations
% r320 for C2M added Gaussian window keyword, Gaussian_histogram_window, for T_O and keyword, QL which is  at Q limit at +/-T_O
% r320 removed external feature and replace with TDMODE
% r320 added TDMODE which allows for the use of pulse resonance files (CSV) instead of s4p files
% r330 changed FOM ILN to use a complex fit and compute FOM_ILN in the time domain330 added tfx to N for ERL
% r335 fixed typo in when processing the bessel thompson filter option
% r335 process in CD mode instead of DC mode to get CM noise at Rx
% r335 compute and report CD_CM_RMS
% r335 fixed where output_arg is save i.e. move to end
% r335 refine interp_Sparam to do zero fill instead of extrapolation
% r335 change raw IL plot to not include boards
% r335 set T_0 to zero if not C2M
% r335 change for s parameter interp: check fit sigma, if not OK zero fill 
% r335 added actual sdd12 (instead fo mirroring sdd21) to s21_pkg and 21dc_pkg and read_s4p_files
% r335 TD_RILN changes from Hansel Dsilva
% r335 Fixed sigma_N for RxFFE
% r335 added more to self documenting keyword capability from read_ParamConfigFile and xls_parameter routines
% r335 added c(2) and C(3) back to read_ParamConfigFile
% r335 Optimize_loop_speed_up keyword option added.  Mostly speeds up c2m(vsr)
% r335 corrected GDC_MIN per 0.3ck D2.3
% r335 sigma_r replaces Qr which replaced QL for Gaussian histogram window
% r340 fix for when post cursor taps 2 and 3 are used (from Matt Brown)
% r370 speed up
% r370 fix for floating tap missing locations
% r370 variable Tx FFE taps 
% r370 package die load with ladder circuit
% r370 mods for SNDR_tx exporation using keyword SNR_TXwC0
% r380 fix for Rx Calibaration (error introduced going from 3.4 to 3.7)
% r380 added capabablity to enable a raised cosine Rx filter0
% r380 keyword added: RC_Start, RC_end, Raised_Cosine
% r380 added plot for VTF 
% r385 added capability for additional Tx FFE per package
% r385 keyword added: PKG_Tx_FFE_preset default is 0 i.e. noop
% r385 SAVE_CONFIG2MAT set to 0 as default i.e. don't create a conifig mat file(
% r388 Adjusted Rx caliberation for CL 162 i.e. adding Tx noise (sigma hn) instead of Rx noise line
% r389 Improvement by A. Ran for reporting loss at Nq
% r389 Fixed typo: changed VIM to VMP 
% r400 fixed PR with zero pad extension
% r400 keyword MLSE and SNRADJ_EQUA for future work
% r400 replaced function db with instances of 20*log10(abs(...))
% r410 widen voltage distriution for normal_dist doubled max Q
% r410 improve reading in of config files
% r410 renormalize s-parameter if not 50 ohm ref
% r410 reference for RXFFE changed to MM from UI+zero first precursor
% r410 remove RL from output_args bc not need and too much storage allocation
% r410 s21^2 changed to s12*s21  in s21_pkg. Corrected VTF needed for non-passive sparameters
% r420 updated equalization figures. if Rxffe is use a subplot of Rx FFE taps is graphed
% r420 updated equalization figures. Now separate per pkg case in optimize_fom
% r420 updade force to account for pulse responces with short delays in force
% r420 added Tx/Rx p/n skew with keywords Txpskew, Txnskew, Rxpskew, Pxnskew
% r420 add common mode outputs: VMC_H_mV and SCMR_dB from CDF of CD PR and DD PR
% r420 fixed and added control for RXFFE_TAP_CONSTRAINT and RXFFE_FLOAT_CTL 
% r420 Wiener-Kofp MMSE optimization for RxFFE 
