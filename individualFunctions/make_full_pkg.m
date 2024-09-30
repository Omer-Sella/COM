function [ s11out, s12out, s21out, s22out]=make_full_pkg(type,faxis,param,channel_type,mode,include_die)

%This function makes the TX or RX package.  The type input must be
%'TX' or 'RX'
%If the mode argument is omitted, mode='dd' is assumed.  Currently
%mode='dc' is only used when making the TX package for AC CM noise
%inclusion.  The Rx package for 'dc' mode is still generated using
%the same parameters as 'dd' mode
%channel_type should be 'THRU' 'FEXT' or 'NEXT'
%
%One instance of package block looks like this (if no elements are set to 0):
%-------------Lcomp----------Tline---------------
%   |                   |               |
%   Cpad                Cbump           Cball
%   |                   |               |
%------------------------------------------------

if nargin<6
    %optional input "include_die"=0 allows die parameters to be forced to 0
    %this includes Cpad, Lcomp, and Cbump
    include_die=1;
end
if nargin<5
    mode='dd';
end


if ~isempty(param.PKG_NAME)
    %The gamma and tau parameters do not currently have a separate Tx and Rx home to live (they were locked for both sides originally)
    %so they are swapped in depending on if Tx or Rx is set for type
    %Note that param is not returned from this function, so the swap does not persist
    swap_fields = {'pkg_gamma0_a1_a2' 'pkg_tau'};
    if strcmpi(type,'tx')
        pkg_name = param.PKG_NAME{1};
    elseif strcmpi(type,'rx')
        pkg_name = param.PKG_NAME{2};
    else
        error('Pkg type must be Tx or Rx');
    end
    pkg_parameter_struct = param.PKG.(pkg_name);

    
    for j=1:length(swap_fields)
        param.(swap_fields{j}) = pkg_parameter_struct.(swap_fields{j});
    end
    
end

C_diepad = param.C_diepad;
C_pkg_board = param.C_pkg_board;
% [ahealey] Unpack optional compensating L and "bump" C model parameters.
L_comp = param.L_comp;
C_bump = param.C_bump;
if ~include_die
    %best to multiply by 0.  that way vectors maintain original size
    C_diepad=C_diepad*0;
    L_comp=L_comp*0;
    C_bump=C_bump*0;
end
% [ahealey] End of modifications.
% generate TX package according to channel type.
[ncases, mele]=size(param.z_p_next_cases);

%Syntax update for C_diepad and L_comp
%Allow a chain of values to be entered as a matrix:
%[L_Tx1 L_Tx2 L_Tx3 ; L_Rx1 L_Rx2 L_Rx3]
if isvector(C_diepad)
    Cd_Tx=C_diepad(1);
    Cd_Rx=C_diepad(2);
    L_comp_Tx=L_comp(1);
    L_comp_Rx=L_comp(2);
    num_blocks=mele;
else
    Cd_Tx=C_diepad(1,:);
    Cd_Rx=C_diepad(2,:);
    L_comp_Tx=L_comp(1,:);
    L_comp_Rx=L_comp(2,:);
    num_blocks=mele+length(Cd_Tx)-1;
end
extra_LC=length(Cd_Tx)-1;
%note:  "insert_zeros" is empty if length(Cd_Tx) = 1
insert_zeros=zeros([1 extra_LC]);

%Updated technique of building Tx/Rx packages
%each index corresponds to the package segment
switch type
    case 'TX'
        switch mele
            case 1
                Cpad=Cd_Tx;
                Lcomp=L_comp_Tx;
                Cbump=C_bump(1);
                Cball=C_pkg_board(1);
                Zpkg=param.pkg_Z_c(1);
            case 4
                Cpad=[Cd_Tx 0 0 0];
                Lcomp=[L_comp_Tx 0 0 0];
                Cbump=[C_bump(1) 0 0 0];
                Cball=[0 0 param.C_v(1) C_pkg_board(1)];
                Zpkg=param.pkg_Z_c(1,:);
            otherwise
                error('package syntax error')
        end
        switch upper(channel_type)
            case 'THRU'
                Len=param.Pkg_len_TX;
            case 'NEXT'
                Len=param.Pkg_len_NEXT;
            case 'FEXT'
                Len=param.Pkg_len_FEXT;
        end
    case 'RX'
        switch mele
            case 1
                Cpad=Cd_Rx;
                Lcomp=L_comp_Rx;
                Cbump=C_bump(2);
                Cball=C_pkg_board(2);
                Zpkg=param.pkg_Z_c(2);
            case 4
                Cpad=[Cd_Rx 0 0 0];
                Lcomp=[L_comp_Rx 0 0 0];
                Cbump=[C_bump(2) 0 0 0];
                Cball=[0 0 param.C_v(2) C_pkg_board(2)];
                Zpkg=param.pkg_Z_c(2,:);
            otherwise
                error('package syntax error')
        end
        switch upper(channel_type)
            case 'THRU'
                Len=param.Pkg_len_RX;
            case 'NEXT'
                Len=param.Pkg_len_RX;
            case 'FEXT'
                Len=param.Pkg_len_RX;
        end
end

%Insert the extra 0 at the front end of Cball, Cbump, Len, and Zpkg
Cball=[insert_zeros Cball];
Cbump=[insert_zeros Cbump];
Len=[insert_zeros Len];
Zpkg=[insert_zeros Zpkg];

% debug_string='';
% for j=1:length(Zpkg)
%     if Cpad(j)~=0
%         debug_string=[debug_string sprintf(', Cd=%0.4g',Cpad(j))];
%     end
%     if Lcomp(j)~=0
%         debug_string=[debug_string sprintf(', Ls=%0.4g',Lcomp(j))];
%     end
%     if Cbump(j)~=0
%         debug_string=[debug_string sprintf(', Cb=%0.4g',Cbump(j))];
%     end
%     if Len(j)~=0
%         debug_string=[debug_string sprintf(', Len=%0.4g Zc=%0.3g',Len(j),Zpkg(j))];
%     end
%     if Cball(j)~=0
%         debug_string=[debug_string sprintf(', Cp=%0.4g',Cball(j))];
%     end
% end
% if length(debug_string)>2
%     debug_string=debug_string(3:end);
% end

% tx package
pkg_param=param;
if strcmpi(mode,'dc')
    % change tx package to CC mode
    pkg_param.Z0=pkg_param.Z0/2;
    Cpad=Cpad*2;
    Cball=Cball*2;
    Zpkg=Zpkg*2;
    Lcomp=Lcomp/2;
    Cbump=Cbump*2;
end
switch num_blocks
    case 1
        [ s11out, s12out, s21out, s22out ]= make_pkg(faxis, Len(1), Cpad(1), Cball(1),Zpkg(1), pkg_param, Lcomp(1), Cbump(1));
    otherwise
        for j=1:num_blocks
            [spkg11,spkg12,spkg21,spkg22]=make_pkg(faxis, Len(j),  Cpad(j),Cball(j) ,Zpkg(j), pkg_param, Lcomp(j),Cbump(j));
            if j==1
                s11out=spkg11; s12out=spkg12; s21out=spkg21; s22out=spkg22;
            else
                [ s11out, s12out, s21out, s22out ]=combines4p(  s11out, s12out, s21out, s22out, spkg11,spkg12,spkg21,spkg22   );
            end
        end
end
function [ s11out, s12out, s21out, s22out ] = make_pkg(f, pkg_len, cpad, cball, pkg_z, param, varargin)
f(f<eps)=eps;
tau = param.pkg_tau; gamma0_a1_a2=param.pkg_gamma0_a1_a2; Lenscale=pkg_len; zref=param.Z0;
%% Equation 93A-8
s11pad= -1i*2*pi.*f*cpad*zref./(2+1i*2*pi.*f*cpad*zref);
s21pad= 2./(2+1i*2*pi.*f*cpad*zref);

% [ahealey] Add compensating L and shunt C (bump) when requested.
s12pad = s21pad;
s22pad = s11pad;
if nargin > 6
    lcomp = varargin{1};
    if lcomp>0
        s11comp = (1i*2*pi*f*lcomp/zref)./(2+1i*2*pi*f*lcomp/zref);
        s21comp = 2./(2+1i*2*pi*f*lcomp/zref);
        [s11pad, s12pad, s21pad, s22pad] = combines4p( ...
            s11pad, s12pad, s21pad, s22pad, ...
            s11comp, s21comp, s21comp, s11comp);
    end
end
if nargin > 7
    cbump = varargin{2};
    if cbump>0
        s11bump = -1i*2*pi.*f*cbump*zref./(2+1i*2*pi.*f*cbump*zref);
        s21bump = 2./(2+1i*2*pi.*f*cbump*zref);
        [s11pad, s12pad, s21pad, s22pad] = combines4p( ...
            s11pad, s12pad, s21pad, s22pad, ...
            s11bump, s21bump, s21bump, s11bump);
    end
end
% [ahealey] End of modifications.

[ S11, S12, S21, S22 ] = synth_tline(f, pkg_z, zref, gamma0_a1_a2, tau, Lenscale); %#ok<NASGU,ASGLU>
% [ahealey] Symmetry cannot be assumed with more complex termination models.
% [ s11out1, s12out1, s21out1, s22out1 ]= ...
%     combines4p(  s11pad, s21pad, s21pad, s11pad, S11, S21, S21, S11 ); % first part of equation 93A-15
[s11out1, s12out1, s21out1, s22out1] = combines4p( ...
    s11pad, s12pad, s21pad, s22pad, ...
    S11, S21, S21, S11);
% [ahealey] End of modifications.

%% Equation 93A-8
s11ball= -1i*2*pi.*f*cball*zref./(2+1i*2*pi.*f*cball*zref);
s21ball= 2./(2+1i*2*pi.*f*cball*zref);
[ s11out, s12out, s21out, s22out ]= ...
    combines4p( s11out1, s12out1, s21out1, s22out1, s11ball, s21ball, s21ball, s11ball );% second part of equation 93A-15

function missingParameter (parameterName)
error( 'error:badParameterInformation', ...
'The data for mandatory parameter %s is missing or incorrect' , parameterName);

function pdf = normal_dist(sigma,nsigma,binsize)
pdf.BinSize=binsize;
pdf.Min=-round(2*nsigma*sigma/binsize); % RIM 03/03/2023 capture more of the tails
pdf.x=(pdf.Min:-pdf.Min)*binsize;
pdf.y=exp(-pdf.x.^2/(2*sigma^2+eps));
pdf.y=pdf.y/sum(pdf.y);
