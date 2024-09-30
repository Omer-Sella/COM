function [chdata, param] = process_sxp(param, OP, chdata, SDDch)
num_files=length(chdata);
if ~OP.DISPLAY_WINDOW, fprintf('reading file '); end
for i=1:num_files
    if param.package_testcase_i==1  && i==1
        if OP.TDR && i==1
            S.Frequencies=chdata(i).faxis;
            S.Impedance=100;
            if ~OP.SHOW_BRD
                Sfield='_orig';
            else
                Sfield='_raw';
            end
            S.Parameters(1,1,:)=chdata(i).(['sdd11' Sfield]) ;
            if ~param.FLAG.S2P
                S.Parameters(1,2,:)=chdata(i).(['sdd12' Sfield]) ;
                S.Parameters(2,1,:)=chdata(i).(['sdd21' Sfield]) ;
                S.Parameters(2,2,:)=chdata(i).(['sdd22' Sfield]) ;
                S.NumPorts=2; % rim 2/26/2019 correct from S.NumPorts=4;
            else
                S.NumPorts=1;
            end
            if OP.TDR_W_TXPKG
                if OP.ERL == 2
                    error('Cannot add pacakge to s2p files. ERL==2 not supportted if TDR_W_TXPKG = 1')
                end
                R_diepad = param.R_diepad;
                % RX package length is assumed to be the same for all
                % channel types.  swap Cp and Cd for Tx. TDR_W_TXPKG only
                % for  Rx pkg
                [ s11in, s12in, s21in, s22in]=make_full_pkg('TX',S.Frequencies,param,'THRU');
                [ S.Parameters(1,1,:), S.Parameters(1,2,:), S.Parameters(2,1,:), S.Parameters(2,2,:)]  = ...
                    combines4p(  s11in, s12in, s21in, s22in, ...
                    S.Parameters(1,1,:), S.Parameters(1,2,:), S.Parameters(2,1,:), S.Parameters(2,2,:) );
                %                     S=sparameters(S.Parameters,S.Frequencies,100);
                S=SL(S,S.Frequencies,R_diepad(1)*2);
                chdata(i).TX_RL=S.Parameters(2,2,:);
                S.Parameters(1,1,:)=chdata(i).sdd11_orig; % when looking at Tx don't include package
            end
            
            % need to combine S wiht is page and channel
            if param.FLAG.S2P
                port_sel=1;
            else
                port_sel=[1 2];
                if OP.AUTO_TFX
                    [ fir4del, tu] =get_RAW_FIR(squeeze(chdata(i).sdd12_orig),S.Frequencies,OP,param);
                    pix=find(fir4del==max(fir4del),1);
                     param.tfx(2)=2*tu(pix);
                end
            end
            OP.impulse_response_truncation_threshold=1e-5; %Only for TDR not returned out of "process_sxp" function
            for ipsl=1:length(port_sel) % do for both port if s4p
                for izt=1:length(param.Z_t) % do for all tdr impedances
                    param.RL_sel=port_sel(ipsl); % this used in get_TDR
                    %                     OP.interp_sparam_phase='interp_to_DC'; % better for return loss
                    %                     OP.interp_sparam_mag='trend_to_DC';
                    OP.interp_sparam_mag='linear_trend_to_DC';
                    % OP.interp_sparam_mag='extrap_to_DC_or_zero';
                    OP.interp_sparam_phase='extrap_cubic_to_dc_linear_to_inf';
                    TDR_results(izt,ipsl) = get_TDR(S, OP, param,param.Z_t(izt),ipsl);
                    if ipsl ==1
                        chdata(i).TDR11(izt).ZSR=[TDR_results(izt,1).tdr];
                        chdata(i).TDR11(izt).t=TDR_results(izt,1).t;
                        chdata(i).TDR11(izt).avgZport=[TDR_results(izt,1).avgZport];
                        if OP.PTDR, chdata(i).PDTR11(izt).ptdr=TDR_results(izt,1).ptdr_RL;end
                    else
                        chdata(i).TDR22(izt).ZSR=[TDR_results(izt,2).tdr];
                        chdata(i).TDR22(izt).t=TDR_results(izt,2).t;
                        chdata(i).TDR22(izt).avgZport=[TDR_results(izt,2).avgZport];
                        if OP.PTDR, chdata(i).PDTR22(izt).ptdr=TDR_results(izt,2).ptdr_RL;end
                    end
                    if OP.PTDR && i==1
                        if ipsl ==1
                            chdata(i).TDR11(izt).ERL=[TDR_results(izt,1).ERL];
                            chdata(i).TDR11(izt).ERLRMS=[TDR_results(izt,1).ERLRMS];
                        else
                            if ~param.FLAG.S2P
                                chdata(i).TDR22(izt).ERL=[TDR_results(izt,2).ERL];
                                chdata(i).TDR22(izt).ERLRMS=[TDR_results(izt,2).ERLRMS];
                            else
                                chdata(i).TDR22(izt).ERL=[];
                                chdata(i).TDR22(izt).ERLRMS=[];
                            end
                        end
                    else
                        chdata(i).TDR11(izt).ERL=[];
                        chdata(i).TDR22(izt).ERL=[];
                        chdata(i).TDR11(izt).ERLRMS=[];
                        chdata(i).TDR22(izt).ERLRMS=[];
                    end
                end
            end
        end
        
    end
    if OP.DISPLAY_WINDOW && OP.DEBUG && OP.TDR
        h=figure(180);set(gcf,'Tag','COM');
        if param.package_testcase_i==1 && i == 1
            if i==1
                htabgroup = uitabgroup(h);
                htab1 = uitab(htabgroup, 'Title', 'TDR TX');
                htab3 = uitab(htabgroup, 'Title', 'PTDR TX');
                hax1 = axes('Parent', htab1);
                hax3 = axes('Parent', htab3);
                if ~param.FLAG.S2P
                    htab2 = uitab(htabgroup, 'Title', 'TDR RX');
                    htab4 = uitab(htabgroup, 'Title', 'PTDR RX');
                    hax2 = axes('Parent', htab2);
                    hax4 = axes('Parent', htab4);
                end
            end
            set(h,'CurrentAxes',hax1)
            hold on
            plot(chdata(i).TDR11(izt).t(:),chdata(i).TDR11(izt).ZSR,'disp',[ chdata(i).base ' Tx port']);
            hold off
            legend (hax1, 'off');grid on;zoom xon;
            set(legend (hax1, 'show'), 'interp', 'none');
            
            if ~param.FLAG.S2P
                set(h,'CurrentAxes',hax2)
                hold on
                plot(chdata(i).TDR22(izt).t(:),chdata(i).TDR22(izt).ZSR,'disp',[ chdata(i).base ' Rx port']);
                hold off
                legend (hax2, 'off');grid on;zoom xon;
                set(legend (hax2, 'show'), 'interp', 'none');
            end
            
            set(h,'CurrentAxes',hax3)
            hold on
            if OP.PTDR
                for izt=1:length(param.Z_t)
                    msg=[ chdata(i).base ' Tx port Zt='  num2str(param.Z_t(izt),3) '  ERL=', num2str(chdata(i).TDR11(izt).ERL, 3) 'db'];
                    plot(chdata(i).TDR11(izt).t/param.ui,chdata(i).PDTR11(izt).ptdr,'disp',msg);
                    msg=['PTDR Zt=' num2str(param.Z_t(izt)/param.ui,3) ' worst sampled noise cursors ' 'Tx port Zt='  num2str(param.Z_t(izt),3) ];
                    stem(TDR_results(izt,1).WC_ptdr_samples_t/param.ui,TDR_results(izt,1).WC_ptdr_samples,'disp',msg);
                end
            end
            hold off
            legend (hax3, 'off');grid on;zoom xon;
            set(legend (hax3, 'show'), 'interp', 'none');
            if ~param.FLAG.S2P
                set(h,'CurrentAxes',hax4)
                hold on
                if OP.PTDR
                    for izt=1:length(param.Z_t)
                        msg=[ chdata(i).base ' Tx port Zt='  num2str(param.Z_t(izt),3) '  ERL=', num2str(chdata(i).TDR22(izt).ERL, 3) 'db'];
                        plot(chdata(i).TDR22(izt).t/param.ui,chdata(i).PDTR22(izt).ptdr,'disp',msg);
                        msg=['PTDR Zt=' num2str(param.Z_t(izt),3)/param.ui ' worst sampled noise cursors ' 'Rx port Zt='  num2str(param.Z_t(izt),3) ];
                        stem(TDR_results(izt,2).WC_ptdr_samples_t/param.ui,TDR_results(izt,2).WC_ptdr_samples,'disp',msg);
                    end
                end
                hold off
                legend (hax4, 'off');grid on;zoom xon;
                set(legend (hax4, 'show'), 'interp', 'none');
            end
        end
    end
    if param.FLAG.S2P, return; end
end