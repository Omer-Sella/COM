function results= get_cm_noise(M,PR,L,BER,OP)

if ~exist('OP')
    OP.DC_norm_test=0;
    OP.DISPLAY_WINDOW=1;
end
param.BinSize=1e-5;
PR_test=-inf;
PR_fom_best=-inf;
% hwaitbar=waitbar(0);
for ki=1:M
    progress = ki/M;
    % if OP.DISPLAY_WINDOW
    %     waitbar(progress, hwaitbar, 'DM to CM computing'); figure(hwaitbar); drawnow;
    % else
    %     if ~mod(progress*100,1), fprintf('%i%% ', progress*100 );end
    % end
    tps=PR(ki:M:end);
    if OP.DC_norm_test
        PR_fom=(norm(tps));
    else
        testpdf=get_pdf_from_sampled_signal( tps,L, param.BinSize*10 );
        cdf_test=cumsum(testpdf.y);
        PRn_test=(-testpdf.x(find(cdf_test>=BER,1,'first')));
        PR_fom=PRn_test;
    end
    if PR_fom > PR_fom_best
        PR_fom_best=PR_fom;
        best_ki=ki;
    end
    if ~OP.DC_norm_test
        results.DCn=PR_fom_best;
        results.DCn_pdf=testpdf;
        results.DCn_cdf=cdf_test;
    else
        results.DCn=PR_fom_best;
    end
    results.DCn_p2p=max(PR)-min(PR);
end

