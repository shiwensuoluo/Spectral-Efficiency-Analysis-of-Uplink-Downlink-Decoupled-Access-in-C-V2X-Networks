clear all;  %  important %nakagami-m   coverage probability analitical expression
            %jlf            

n=80
for i = 1:n
% % ============ SPATIAL MODEL ========= %nagani-m channeel
lambda_l = 10/pi; % (LINE_DENSITY/PI); Line density - km^-1
lambda_r = 15; %(NUMBER OF RECEIVING VEHICULAR NODES PER KM)
lambda_2 = i*0.5; %NUMBER OF TIER 2 NODES PER KM
lambda_1 = 5; %NUMBER OF TIER 1 NODES PER KM SQ. 

% % ============ TRANSMIT PARAMETERS ========= %
p1_dbm = 46; % Transmit power of tier 1 node in dBm
p1 = 10^(.1*(p1_dbm-30));
p2_dbm = 20; % Transmit power of tier 2 node in dBm
p2 = 10^(.1*(p2_dbm-30));

pv_dbm = 20;
pv = 10^(.1*(pv_dbm-30));

%%gain增益
mlg_1 = 1; % Mainlobe gain TIER 1 node  1
%slg_1 = .01; % Sidelobe gain TIER 1 node

mlg_2 = 1; % Mainlobe gain TIER 2 node
slg_2 = .01; % Sidelobe gain TIER 2 node

v_m = 0.01;  %vehicle gain to typical
v_s0 = 1;  %vehicle gain to other
v_s1 = 0.01;%vehicle gain to other



b1_db = 0; % Selection bias of TIER 1 node in dB  in dl
b1 = 10^(.1*b1_db);  
b2_db = 0; % Selectio  bias of TIER 2 node in dB in  dl
b2 = 10^(.1*b2_db);


A_MS = p1*b1*mlg_1/(p2*b2*mlg_2); %dl
B_MS = pv*v_m/(pv*v_s0);          %ul
%qc = .05; % Probability with which the mainlobe of interfering TIER 1 node is directed towards the typical receiver


% % ============ PROPAGATION PARAMETERS =================================================== %
alpha1 = 4; % PATH_LOSS EXPONENT 
alpha2 = 4;  %tier 2

% ~~~ NAKAGAMI-M FADING PARAMETERS ~~~ %
m1 = 1; % Tier 1 本来是1
m20 = 2; % Tier 2 typical line
m21 = 1; % Tier 2 other lines
% m = 1;

% ~~~ SHADOWING PARAMETERS ~~~~ %
ln_mu1 = 0;% Mean of log-normal shadowing gain in dB for TIER 1
ln_sig1 = 4; % Std deviation of log-normal shadowing gain in dB for TIER 1

ln_mu20 = 0; % Mean of log-normal shadowing gain in dB for TIER 2 TYP LINE
ln_sig20 = 2; % Std deviation of log-normal shadowing gain in dB for TIER 2 TYP LINE

ln_mu21 = 0;% Mean of log-normal shadowing gain in dB for TIER 2 OTHER LINES
ln_sig21 = 4;% Std deviation  of log-normal shadowing gain in dB for TIER 2 OTHER LINES

% % % ======================================================================================

sir_threshold_dB = 0;%-50:5:50;
sir_threshold = 10.^(.1*sir_threshold_dB);

%---------------END of PARAMETERS ------------------------------------------------------
alpha = 4 ;

sf_c = exp( 2/alpha1*log(10)/10*ln_mu1 + (2/alpha1*log(10)/10)^2*ln_sig1^2/2 );   %tier 1

sf_u0 = exp( 1/alpha2*log(10)/10*ln_mu20 + (1/alpha2*log(10)/10)^2*ln_sig20^2/2);%tier 2 line 
sf_a = exp( 2/alpha2*log(10)/10*ln_mu21 + (2/alpha2*log(10)/10)^2*ln_sig21^2/2 ); % tier 2 other
%替代理论后的密度
sh_lambda_1 = sf_c*lambda_1; %1 
sh_lambda_2 = sf_u0*lambda_2; %20
sh_lambda_a = sf_a*(pi*lambda_l*lambda_2);%21



%%原文的计算覆盖
% peu = sh_lambda_2*(1/(sh_lambda_1*((p2*b2)*mlg_2/(p1*b1))^(-2/alpha)))^.5*...
%     exp(sh_lambda_2^2/(sh_lambda_1*pi*((p2*b2)*mlg_2/(p1*b1))^(-2/alpha)))*...
%     erfc((sh_lambda_2^2/(sh_lambda_1*pi*((p2*b2)*mlg_2/(p1*b1))^(-2/alpha)))^.5 );
% pec = 1- peu;


cdfrm = @(xM) 1 - exp(-sh_lambda_1*pi*xM.^2);    %MBS 
pdfrm = @(xM) 2*pi*sh_lambda_1*xM.*exp(-sh_lambda_1*pi*xM.^2);

cdfrs = @(xS) 1 - exp(-2*sh_lambda_2*xS);  %sbs
pdfrs = @(xS) 2*sh_lambda_2*exp(-2*sh_lambda_2*xS);


%case2的概率
pr_case2 = sh_lambda_2*(1/(sh_lambda_1*(B_MS)^(2/alpha)))^.5*...
     exp(sh_lambda_2^2/(sh_lambda_1*pi*(B_MS)^(2/alpha)))*...
     erfc((sh_lambda_2^2/(sh_lambda_1*pi*(B_MS)^(2/alpha)))^.5 )-...
     sh_lambda_2*(1/(sh_lambda_1*(A_MS)^(2/alpha)))^.5*...    
     exp(sh_lambda_2^2/(sh_lambda_1*pi*(A_MS)^(2/alpha)))*...
     erfc((sh_lambda_2^2/(sh_lambda_1*pi*(A_MS)^(2/alpha)))^.5 );

 %distance distribution  of case2  DL
 pdf_case2_DL =@(x) pdfrm(x)/pr_case2.*(exp(-sh_lambda_2*pi*2*(1/A_MS)^(1/alpha2)*x^(alpha1/alpha2))-...
                                        exp(-sh_lambda_2*pi*2*(1/B_MS)^(1/alpha2)*x^(alpha1/alpha2)));
   
% pdf_case2_UL =@(x) pdfrs(x)/pr_case2.*(exp(-sh_lambda_1*pi*(B_MS)^(2/alpha1)*x^(alpha2/alpha1))-...
%                                         exp(-sh_lambda_1*pi*(A_MS)^(2/alpha1)*x^(alpha2/alpha1)));     
%pdfr_g_u = @(r) pdfru(r)/peu.*(1-cdfrc( ((p2*b2)*mlg_2/(p1*b1))^(-1/alpha)*r)) ; %+ 1/peu*( cdfru( ((p1*b1)/(p2*b2)/mlg_2)^(-1/alpha)*r) - cdfru(r) ).*pdfrc(r);

%pdfr_g_c = @(r) pdfrc(r)/pec.*(1- cdfru( ((p1*b1)/(p2*b2)/mlg_2)^(-1/alpha)*r)) ;





lapI_m = @(j,x_m) exp (-2*pi*sh_lambda_1*integral(@(x) ( 1 - (1+ j*(p1)*mlg_1*x.^(-alpha1)/m1).^(-m1)).*x , x_m, Inf,'ArrayValued', true));
lapI_s0 = @(j,x_m) exp (-2*sh_lambda_2*integral(@(x) (1 - ((1+ j*p2*mlg_2*x.^(-alpha2)/m20).^(-m20))),(1/A_MS)^(1/alpha2)*x_m^(alpha1/alpha2), Inf,'ArrayValued', true));
lapI_s1 = @(j,x_m) exp( -2*pi*sh_lambda_a*integral(@(x) ( 1 - (1+ j*pv*slg_2*x.^(-alpha2)/m21).^(-m21)).*x,0, Inf,'ArrayValued', true));

lapI_case2_DL = @(j,x_m) lapI_m(j,x_m) .*lapI_s0(j,x_m) .*lapI_s1(j,x_m) ;

% pc_e1 = @(r) 0;
 se_case2_DL = 0;
pc_e1 = @(t,r) 0;
f1 = @(s,x) -2*pi*sh_lambda_1*(1 - (1+ s*(p1)*mlg_1*x.^(-alpha1)/m1).^(-m1)).*x;

f2 = @(s,x) -2*sh_lambda_2*( 1 - (1+ s*(p2*mlg_2)*x.^(-alpha2)/m20).^(-m20));
f3 = @(s,x) -2*pi*sh_lambda_a*( 1 - (1+ s*pv*slg_2*x.^(-alpha2)/m21).^(-m21)).*x;

fg_mf = @(r,s) integral(@(x) f1(s,x), r,Inf,'ArrayValued', true) + ...
    integral(@(x) f2(s,x),(1/A_MS)^(1/alpha2)*r^(alpha1/alpha2), Inf,'ArrayValued', true) +...
    integral(@(x) f3(s,x),0,Inf,'ArrayValued', true);




if m1==1
   % pc_e1 = @(b) (integral(@(r) lapI_g_c(m1*b*r.^(alpha)/(p1*mlg_1),r).*pdfr_g_c(r), 0, Inf, 'ArrayValued', true));
   pc_e1 = @(x_m) (integral( @(t)  lapI_case2_DL(m1*x_m.^(alpha1)*(exp(t)-1)/(p1*mlg_1),x_m),0, Inf, 'ArrayValued', true));
   se_case2_DL= integral(@(r) pc_e1(r).*pdf_case2_DL(r),0,inf,'ArrayValued', true)
else  
    syms s r x t;
    assume(s,'real');
    assume(r,'positive');
%     assume(t,'positive');
    g1 = int( f1(s,x),x,r,Inf);%ints是对x求定积分
    g2 = int( f2(s,x),x,(1/A_MS)^(1/alpha2)*r^(alpha1/alpha2), Inf);
    g3 = int( f3(s,x),x,0,Inf);
    fg = g1+g2+g3;
    
     for n=0:m1-1
       derk1_I_e1 = @(r,s) 0;
        for k=0:n
            for j=0:k
                a1_str = char( diff( (fg)^(k-j), 's', n)); % To convert sym expr to function handle
                a1_str = strrep(a1_str, 'int(', 'integral(@(x)');
                a1_str = strrep(a1_str, 'x,','');
                a1_str = strrep(a1_str, '*','.*');   
                a1_str = strrep(a1_str, '/','./');
                a1_str = strrep(a1_str, '^','.^');
                a1_mf = eval(['@(r,s) ' a1_str]); % Converting Sym expr to Matlab Function handle
                derk1_I_e1 = @(r,s) derk1_I_e1(r,s) + (-1)^j/factorial(k)*nchoosek(k,j)*...
                                                                     a1_mf(r,s).*(fg_mf(r,s)).^j;
            end     
        end

        dern_I_e1 = @(s,r) lapI_case2_DL(s,r).*derk1_I_e1(r,s);
        
        j_m =@(t,r) m1*(exp(t)-1)/(p1*mlg_1)*r.^(alpha1);%(exp(t)-1)/(p1*mlg_1);%jlf加的
        

        pc_e11 = @(t,r)(-1*j_m(t,r)).^n/factorial(n).*dern_I_e1(j_m(t,r),r);
        
       
        pc_e1 = @(t,r) pc_e1(t,r)+ pc_e11(t,r);       
     end
      pc_e = @(r) (integral(@(t) pc_e1(t,r), 0,19, 'ArrayValued', true));
      se_case2_DL = se_case2_DL + integral(@(r) pc_e(r).*pdf_case2_DL(r),0,inf, 'ArrayValued', true); %
        

    
end

se_lambda_case2_DL(i) = se_case2_DL;
aa(i) = i*0.5/5.0;
i = i + 1;
end
  



