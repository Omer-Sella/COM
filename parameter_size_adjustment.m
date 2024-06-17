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
