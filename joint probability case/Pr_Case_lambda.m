clear all;  %  important %nakagami-m   coverage probability analitical expression
            %jlf 
n=30
for i = 1:n
% % ============ SPATIAL MODEL ========= %nagani-m channeel
lambda_l = 10/pi; % (LINE_DENSITY/PI); Line density - km^-1
lambda_r = 15; %(NUMBER OF RECEIVING VEHICULAR NODES PER KM)
lambda_2 = i*0.5; %NUMBER OF TIER 2 NODES PER KM
lambda_1 = 5; %NUMBER OF TIER 1 NODES PER KM SQ. 原来是0.5

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

gv = 1 ;



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
alpha1 = 2; % PATH_LOSS EXPONENT   NLOS是非视距4   LOS是2
alpha2 = 2;  %tier 2

alpha = 2;  %NLOS 是4   los是2

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


sf_c = exp( 2/alpha1*log(10)/10*ln_mu1 + (2/alpha1*log(10)/10)^2*ln_sig1^2/2 );   %tier 1

sf_u0 = exp( 1/alpha2*log(10)/10*ln_mu20 + (1/alpha2*log(10)/10)^2*ln_sig20^2/2);%tier 2 line 
sf_a = exp( 2/alpha2*log(10)/10*ln_mu21 + (2/alpha2*log(10)/10)^2*ln_sig21^2/2 ); % tier 2 other
%替代理论后的密度
sh_lambda_1 = sf_c*lambda_1; %1 
sh_lambda_2 = sf_u0*lambda_2; %20
sh_lambda_a = sf_a*(pi*lambda_l*lambda_2);%21





cdfrm = @(xM) 1 - exp(-sh_lambda_1*pi*xM.^2);    %MBS 
pdfrm = @(xM) 2*pi*sh_lambda_1*xM.*exp(-sh_lambda_1*pi*xM.^2);

cdfrs = @(xS) 1 - exp(-2*sh_lambda_2*xS);  %sbs
pdfrs = @(xS) 2*sh_lambda_2*exp(-2*sh_lambda_2*xS);





%case1的概率
pr_case1 = 1 - sh_lambda_2*(1/(sh_lambda_1*(B_MS)^(2/alpha)))^.5*...
     exp(sh_lambda_2^2/(sh_lambda_1*pi*(B_MS)^(2/alpha)))*...
     erfc((sh_lambda_2^2/(sh_lambda_1*pi*(B_MS)^(2/alpha)))^.5 );
 
%case 2 概率
pr_case2 = sh_lambda_2*(1/(sh_lambda_1*(B_MS)^(2/alpha)))^.5*...
     exp(sh_lambda_2^2/(sh_lambda_1*pi*(B_MS)^(2/alpha)))*...
     erfc((sh_lambda_2^2/(sh_lambda_1*pi*(B_MS)^(2/alpha)))^.5 )-...
     sh_lambda_2*(1/(sh_lambda_1*(A_MS)^(2/alpha)))^.5*...    
     exp(sh_lambda_2^2/(sh_lambda_1*pi*(A_MS)^(2/alpha)))*...
     erfc((sh_lambda_2^2/(sh_lambda_1*pi*(A_MS)^(2/alpha)))^.5 );

 
 %case4的概率
pr_case4 = sh_lambda_2*(1/(sh_lambda_1*(A_MS)^(2/alpha)))^.5*...
     exp(sh_lambda_2^2/(sh_lambda_1*pi*(A_MS)^(2/alpha)))*...
     erfc((sh_lambda_2^2/(sh_lambda_1*pi*(A_MS)^(2/alpha)))^.5 );
 
 
Pr_case_lambda(1,i) = pr_case1;
Pr_case_lambda(2,i) = pr_case2;
Pr_case_lambda(3,i) = pr_case4;
x(i) = i*0.5/5;
i = i + 1;
end
 
%% 
% Pr_case_lambda=[1,Pr_case_lambda(1,:);0,Pr_case_lambda(2,:);0,Pr_case_lambda(3,:)];

% x = 1:1:n;
%box on
plot(x,Pr_case_lambda(1,:));
hold on;
plot(x,Pr_case_lambda(2,:));
hold on;
plot(x,Pr_case_lambda(3,:));
hold on;
plot(x,zeros(1,n));
grid on;
axis([0 2.5 0 1]);