clear;
%%%use monte carlo to campute the coverage and rate
%参数    计算额了频谱效率   覆盖概率
d =5;  % OBSERVATION WINDOW RADIUS  original——40

n=61
for iiii = 1:4:n
    iiii
% % ============ SPATIAL MODEL ========= %nagani-m channeel
lambda_l = 10/pi; % (LINE_DENSITY/PI); Line density - km^-1
lambda_r = 15; %(NUMBER OF RECEIVING VEHICULAR NODES PER KM)
lambda_2 = iiii*0.5; %NUMBER OF TIER 2 NODES PER KM 原来是4
lambda_1 = 5; %NUMBER OF TIER 1 NODES PER KM SQ. 原来是0.5



% % ============ TRANSMIT PARAMETERS ========= %
p1_dbm = 46; % Transmit power of tier 1 node in dBm
p1 = 10^(.1*(p1_dbm-30));
p2_dbm = 20; % Transmit power of tier 2 node in dBm
p2 = 10^(.1*(p2_dbm-30));

pv_dbm = 20;
pv = 10^(.1*(pv_dbm-30));


mlg_2 = 1; % Mainlobe gain TIER 2 node
slg_2 = .01; % Sidelobe gain TIER 2 node

mlg_1 = 1; % Mainlobe gain TIER 1 node
slg_1 = .01; % Sidelobe gain TIER 1 node

gv  = 1;

v_m = 0.01;

v_s0 = 1;
v_s1 = 0.01;

b1_db = 0; % Selection bias of TIER 1 node in dB
b1 = 10^(.1*b1_db);
b2_db = 0; % Selection bias of TIER 2 node in dB
b2 = 10^(.1*b2_db);

qc = .05; % Probability with which the mainlobe of interfering TIER 1 node is directed towards the typical receiver


% % ============ PROPAGATION PARAMETERS =================================================== %
alpha1 = 4;
alpha2 =4;
alpha = 4; % PATH_LOSS EXPONENT

% ~~~ NAKAGAMI-M FADING PARAMETERS ~~~ %
m1 = 1; % Tier 1
m20 = 2; % Tier 2 typical line
m21 = 1; % Tier 2 other lines

% ~~~ SHADOWING PARAMETERS ~~~~ %
ln_mu1 = 0;% Mean of log-normal shadowing gain in dB for TIER 1
ln_sig1 = 4; % Std deviation of log-normal shadowing gain in dB for TIER 1

ln_mu20 = 0; % Mean of log-normal shadowing gain in dB for TIER 2 TYP LINE
ln_sig20 = 2; % Std deviation of log-normal shadowing gain in dB for TIER 2 TYP LINE

ln_mu21 = 0;% Mean of log-normal shadowing gain in dB for TIER 2 OTHER LINES
ln_sig21 = 4;% Std deviation  of log-normal shadowing gain in dB for TIER 2 OTHER LINES

% % % ======================================================================================

num_iterations = 1e1;%原文是1e4

sir_threshold_dB = 0;
sir_threshold = 10^(.1*sir_threshold_dB);

EVAL_RATE = 1;%%原来是0
bw = 10; % BW in MHz  这个是每个基站的带宽  
rate_threshold = 10; % in Mbps
%-------------------------------------------------------------------
%======================Simulation=============================

sir_measured = zeros(1, num_iterations);
load_measured = zeros(1, num_iterations); 




%%仿真布置
numc = poissrnd(lambda_1*pi*(d^2),1,num_iterations); % Number of TIER 1 nodes in observation window for all iterations   服从lambade分布的  
numl = poissrnd(lambda_l*(2*d*pi),1,num_iterations); % Number of LINES in observation window for all iterations



for j1=1:num_iterations
    Nl = numl(j1); % Number of lines 从泊松点的数字里随机选出一个值
    Nc = numc(j1); % Number of TIER 1 nodes
    
    
    % GENERATING THE RANDOM LINES
    rho1 = [0 -d+2*d*rand(1,Nl) ]; %[ (rv+(d-rv)*urv).*sign(urv-.5) 0];  这里多了一个0保证第一部分都是在典型线上   
    [~,ind2] = sort(abs(rho1));% - rho1(end);  ind2是排序对应的位置
    rho = rho1(ind2) - 2*d*(rho1(ind2)>d) + 2*d*(rho1(ind2)<-d);%这步感觉就多此一举
    theta =pi*rand(1,length(rho)); %%极坐标 的值
    len_chord = 2*(d^2-(rho.^2)).^.5; % Length of chord of each line
    
    % ============ vehicular user generation ===================
    numv = poissrnd(lambda_r*len_chord);
    numv = numv + (numv==0);%将为零的道路上补上一辆车
    total_v = sum(numv);
    
    clens = cumsum(numv);%%从0到1逐步累加
    idx=zeros(1,clens(end));idx2=idx; idx3=idx;
    idx([1 clens(1:end-1)+1]) = diff([0 len_chord]); %%mei条路之间长度差值
    len_c = cumsum(idx);    %路宽
    idx2([1 clens(1:end-1)+1]) = diff([0 rho]);
    rho_vec = cumsum(idx2);%到原点距离
    idx3([1 clens(1:end-1)+1]) = diff([0 theta]);
    theta_vec = cumsum(idx3);%角度
    
    x1 = (-len_c/2 + len_c.*rand(1,total_v));%位置  将所有车的位置都表示在一个矩阵里
    x = zeros(size(x1));
    for k1=1:length(rho) %rho是线的数目 对每一根线上的车的位置排一下
        x( sum(numv(1:k1-1))+1: sum(numv(1:k1))) = sort(x1(sum(numv(1:k1-1))+1: sum(numv(1:k1))));%%这是对应的每条路的车进行位置排序
    end
    beta = atan(x./rho_vec);%求每个车所在的位置
    gamma = (theta_vec+beta);
    u = (rho_vec.^2+x.^2).^.5;%车到原点的距离
    signm = sign(rho_vec) + 1*(rho_vec==0);%sign正为一，负为-1，零为0    距离原点为0的设为1,避免直接成了0
    pts_v = ([u.*cos(gamma).*signm] + 1j*[u.*sin(gamma).*signm]).';  %表示成极坐标形式
    
    
    % ============ TIER 2 generation ===================
    numu = poissrnd(lambda_2*len_chord);
    numu = numu + (numu==0);
    total_u = sum(numu);
    
    clens = cumsum(numu);
    idx=zeros(1,clens(end));idx2=idx; idx3=idx;
    idx([1 clens(1:end-1)+1]) = diff([0 len_chord]);
    len_c = cumsum(idx);  %%路长
    idx2([1 clens(1:end-1)+1]) = diff([0 rho]);
    rho_vec = cumsum(idx2);  %线到原点的距离
    idx3([1 clens(1:end-1)+1]) = diff([0 theta]);
    theta_vec = cumsum(idx3);  %线的 角度
    
    x1 = [ -len_c/2 + len_c.*rand(1,total_u)];
    x = zeros(size(x1));
    for k1=1:length(rho)
        x( sum(numu(1:k1-1))+1: sum(numu(1:k1))) = sort(x1(sum(numu(1:k1-1))+1: sum(numu(1:k1))));%%sort是从小到大排列
    end
    beta = atan(x./rho_vec);
    gamma = (theta_vec+beta);
    u = (rho_vec.^2+x.^2).^.5; %点到原点的距离
    signm = sign(rho_vec) + 1*(rho_vec==0);   %表示点在零轴的左边还是右边
    pts_u = ([ u.*cos(gamma).*signm] + 1j*[ u.*sin(gamma).*signm]).';%极坐标表示  signm排序后把从小到大弄反了又
    dist_u = abs(pts_u);  
    %信道H
    chnl_u = [gamrnd(m20, 1/m20, length(dist_u(1:numu(1))), 1); gamrnd(m21, 1/m21, length(dist_u(numu(1)+1:end)),1) ];
    
    % ================== Cellular BSs generation =======================
    rad_c = sqrt(rand(Nc,1))*d; %Nc层1节点数目
    phi_c = rand(Nc,1)*(2*pi);  %角度
    xc = [ rad_c.*cos(phi_c)];
    yc = [ rad_c.*sin(phi_c)]; % x =[ -d + 2*d*rand(1,Nc)];% y =[ -d + 2*d*rand(1,Nc)];
    pts_c = xc+1j*yc;  %极坐标
    dist_c = abs(pts_c); %距离
    % Rc(j1) = min(dist_c);
    %  if Nc>0
    [vx,vy] = voronoi(xc,yc);%泰森多边形 vx: Voronoi 边角的 x 坐标，以列向量形式返回。vy:Voronoi 边角的 y 坐标，以列向量形式返回。
    [verts, cs] = voronoin([xc(:) yc(:)]);%基于由矩阵 P 表示的 N 维点确定一个 Voronoi 图并返回其 Voronoi 顶点 v 和 Voronoi 元胞 c。
                                             %c：是索引，以元胞数组形式返回。c 的每个元素都包含构成 Voronoi 元胞的 Voronoi 顶点 v 的行索引。
    %蜂窝基站的信道
    chnl_c = gamrnd(m1, 1/m1, length(dist_c), 1);%%gamma分布函数  信道参数H
    %  end

    
    
    % ==================== ============== ============= ======== ========= ========== =
    %     %     [val, ind] = min(dist_c);
    typ_verts_tmp = verts(cs{1},:);
    [ax,ay] =poly2cw(typ_verts_tmp(:,1) , typ_verts_tmp(:,2) );  %对点逆时针排序
    typ_verts =  [ax ay];
    typ_perm(j1) = sum(abs(diff( [typ_verts(:,1)+1j*typ_verts(:,2); typ_verts(1,1)+1j*typ_verts(1,2)])));%%%后面这个是做了个闭环，把第一项也补上了
    typ_area(j1) = polyarea(typ_verts(:,1), typ_verts(:,2));%%计算多边形的面积
    % ==================== ============== ============= ======== ========= ========== =
    %     %     [val, ind] = min(dist_c);
    typ_verts_tmp = verts(cs{1},:);
    [ax,ay] =poly2cw(typ_verts_tmp(:,1) , typ_verts_tmp(:,2) );
    typ_verts =  [ax ay];
    typ_perm(j1) = sum(abs(diff( [typ_verts(:,1)+1j*typ_verts(:,2); typ_verts(1,1)+1j*typ_verts(1,2)])));
    typ_area(j1) = polyarea(typ_verts(:,1), typ_verts(:,2));
    
    Ru(j1) = min(dist_u(1:numu(1)));%dist_u 是层2节点
    % ================== Load & SIR Computation =======================
  
   %val是值 ind是索引位置
   

    for iL = 1: numl(j1)+1
        
%           iL
     
                  
        for iv = 1:numv(iL)
            
            
            shad_u0 = 10.^(.1*(ln_mu20 + ln_sig20*randn(numu(iL),1)));   %典型线上的
            shad_c = 10.^(.1*(ln_mu1 + ln_sig1*randn(length(dist_c), 1))); %原点的
            
       %%%%shadow     
            shad_u0u = [10.^(.1*(ln_mu21 + ln_sig21*randn(sum(numu(1:iL-1)),1)));...
                      10.^(.1*(ln_mu20 + ln_sig20*randn(numu(iL),1)));...
                      10.^(.1*(ln_mu21 + ln_sig21*randn(sum(numu(iL+1:end)),1)))];%shad_u0指在典型线上
            shad_c = 10.^(.1*(ln_mu1 + ln_sig1*randn(length(dist_c), 1))); %原点的
            
      %%%mainlobe 
           bf_gains = [slg_2*ones(sum(numu(1:iL-1)),1);
                       mlg_2*ones(numu(iL),1);
                       slg_2*ones(sum(numu(iL+1:end)),1)];%%numu1默认是在典型线上的
                  
           mlg_1 = 1; % Mainlobe gain TIER 1 node
      %%%%%%%%%%channel parameter  
          chnl_u = [gamrnd(m20, 1/m20, numu(iL), 1); 
                      gamrnd(m21, 1/m21, total_u-numu(iL),1) ];
          chnl_c = gamrnd(m1, 1/m1, length(dist_c), 1);%%gamma分布函数  信道参数H  
                  
           
            [DLval,DLind] = max( ...
                          [p2*mlg_2*b2*shad_u0.*(abs(pts_v(sum(numv(1:iL-1))+iv)-...
                                                     pts_u(sum(numu(1:iL-1))+1:sum(numu(1:iL-1))+numu(iL)))).^(-alpha2);...
                          p1*mlg_1*b1*shad_c.*(abs(pts_v(sum(numv(1:iL-1))+iv)-pts_c)).^(-alpha1)] ...
                             );%%找这个一大列最大值。   下行  %最大功率的SBS和MBS
                         
                         

            if DLind <= length(dist_u(sum(numu(1:iL-1))+1:sum(numu(1:iL-1))+numu(iL)))  
                pr_decouple(1,sum(numv(1:iL-1))+iv)= 1;%第一行下行，第2行是上行。1指sbs；  2 指MBS
                
   
                sig_pwr = p2*bf_gains(DLind+sum(numu(1:iL-1)))*shad_u0u(DLind+sum(numu(1:iL-1)))*chnl_u(DLind).*...
                          abs(pts_v(sum(numv(1:iL-1))+iv)-pts_u(sum(numu(1:iL-1))+DLind)).^(-alpha2); %abs(tagged_rsu).^(-alpha);

                total_pwr = [p2*bf_gains.*shad_u0u.*...
                                 [chnl_u(numu(iL)+1:numu(iL)+sum(numu(1:iL-1)));chnl_u(1:numu(iL));chnl_u(sum(numu(1:iL))+1:end)]...  
                                 .*(abs(pts_v(sum(numv(1:iL-1))+iv)-pts_u)).^(-alpha2);...
                                 p1*mlg_1.*chnl_c.*shad_c.*(abs(pts_v(sum(numv(1:iL-1))+iv)-pts_c)).^(-alpha1)] ;
%                 sir = sig_pwr/(total_pwr - sig_pwr);
                sir = sig_pwr/sum( [total_pwr(1:find(total_pwr==sig_pwr)-1);total_pwr(find(total_pwr==sig_pwr)+1:end)]);
                sir_decouple(1,sum(numv(1:iL-1))+iv)= sir;
                
            else
                pr_decouple(1,sum(numv(1:iL-1))+iv)= 2;
                
                sig_pwr = p1*mlg_1*shad_c(DLind-numu(iL)).*chnl_c(DLind-numu(iL)).*(abs(pts_v(sum(numv(1:iL-1))+iv)-pts_c(DLind-numu(iL)))).^(-alpha1) ;   
                total_pwr = sum( [p2*bf_gains.*shad_u0u.*...
                                 [chnl_u(numu(iL)+1:numu(iL)+sum(numu(1:iL-1)));chnl_u(1:numu(iL));chnl_u(sum(numu(1:iL))+1:end)]...
                                 .*(abs(pts_v(sum(numv(1:iL-1))+iv)-pts_u)).^(-alpha2);...
                                 p1*mlg_1.*chnl_c.*shad_c.*abs(pts_v(sum(numv(1:iL-1))+iv)-pts_c).^(-alpha1)] );
                sir = sig_pwr/(total_pwr - sig_pwr);
                sir_decouple(1,sum(numv(1:iL-1))+iv)= sir;  
                
            end
            

          UL_shad_u0u = [10.^(.1*(ln_mu21 + ln_sig21*randn(sum(numv(1:iL-1)),1)));...
                      10.^(.1*(ln_mu20 + ln_sig20*randn(numv(iL),1)));...
                      10.^(.1*(ln_mu21 + ln_sig21*randn(sum(numv(iL+1:end)),1)))];%shad_u0指在典型线上
          UL_shad_c = 10.^(.1*(ln_mu1 + ln_sig1*randn(length(pts_v), 1))); %原点的
            
      %%%mainlobe 
          bf_gains = [v_s1*ones(sum(numv(1:iL-1)),1);
                       v_s0*ones(numv(iL),1);
                       v_s1*ones(sum(numv(iL+1:end)),1)];%%numu1默认是在典型线上的
                   
         v_m = 0.01; % Mainlobe gain TIER 1 node
      %%%%%%%%%%channel parameter  
         chnl_u = [gamrnd(m20, 1/m20, numv(iL), 1); 
                      gamrnd(m21, 1/m21, total_v-numv(iL),1) ];
         chnl_c = gamrnd(m1, 1/m1, length(pts_v), 1);%%gamma分布函数  信道参数H  
                  
            
            
            %%%%%%%%%%%%%%% UP Link

            
            [ULval,ULind] = max( ...
                          [pv*v_s0*shad_u0.*(abs(pts_v(sum(numv(1:iL-1))+iv)-pts_u(sum(numu(1:iL-1))+1:sum(numu(1:iL-1))+numu(iL)))).^(-alpha2);...
                          pv*v_m*shad_c.*(abs(pts_v(sum(numv(1:iL-1))+iv)-pts_c)).^(-alpha1)] ...
                          );%%找这个一大列最大值。  上行
            if ULind <= length(pts_u(sum(numu(1:iL-1))+1:sum(numu(1:iL-1))+numu(iL)))  
                pr_decouple(2,sum(numv(1:iL-1))+iv)= 1;%第一行下行，第2行是上行。1指sbs；  2 指MBS
                
                
                
                sig_pwr = pv*bf_gains(iv+sum(numv(1:iL-1))).*UL_shad_u0u(iv+sum(numv(1:iL-1))).*chnl_u(iv).*...
                             abs(pts_v(sum(numv(1:iL-1))+iv)-pts_u(sum(numu(1:iL-1))+ULind)).^(-alpha2); %abs(tagged_rsu).^(-alpha);
                
                     
                total_pwr = [pv*bf_gains.*UL_shad_u0u.*...
                                 [chnl_u(numv(iL)+1:numv(iL)+sum(numv(1:iL-1)));chnl_u(1:numv(iL));chnl_u(sum(numv(1:iL))+1:end)]...
                                 .*(abs(pts_v-pts_u(sum(numu(1:iL-1))+ULind))).^(-alpha2)];
%                 sir = sig_pwr/(total_pwr - sig_pwr);
                sir = sig_pwr/sum( [total_pwr(1:find(total_pwr==sig_pwr)-1);total_pwr(find(total_pwr==sig_pwr)+1:end)]);
                
                sir_decouple(2,sum(numv(1:iL-1))+iv)= sir;
                
            else
                pr_decouple(2,sum(numv(1:iL-1))+iv)= 2;
                
             
                
                sig_pwr = pv*v_m*UL_shad_c(sum(numv(1:iL-1))+iv)*chnl_c(sum(numv(1:iL-1))+iv)*...  
                         (abs(  pts_v(sum(numv(1:iL-1))+iv) - pts_c(ULind-numu(iL)))  ).^(-alpha1);

                total_pwr = sum( [pv*v_m.*UL_shad_c.*...
                                 chnl_c.*...
                                 (abs(pts_v-pts_c(ULind-numu(iL)))).^(-alpha1)]);
                             
                sir = sig_pwr/(total_pwr - sig_pwr);
                sir_decouple(2,sum(numv(1:iL-1))+iv)= sir;   
            end
            
           
       
        end  
        
    end  
    

    
    
    [case1 posti_case1] = find_case(pr_decouple , 2, 2); % d m   u-M
    [case2 posti_case2] = find_case(pr_decouple , 2, 1); % d-M   u-S
%    [case3 posti_case3] = find_case(pr_decouple , 1, 2); %d-S    u-M
    [case4 posti_case4] = find_case(pr_decouple , 1, 1); %d-S   u-S
    
    
    
        D_se(j1) = sum(log(1+sir_decouple(1,:))) / total_v;
    U_se(j1) = sum(log(1+sir_decouple(2,:))) / total_v;
    
    
    pr_case(1:3,j1) = [case1;case2;case4]./ (length(posti_case1)+length(posti_case2)+length(posti_case4));
    
    sir_case(1:6,j1) = [sum(log(1+sir_decouple(1,posti_case1)))/length(posti_case1);...%1 下行
                        sum(log(1+sir_decouple(2,posti_case1)))/length(posti_case1);...  %1 上行
                        sum(log(1+sir_decouple(1,posti_case2)))/length(posti_case2);...
                        sum(log(1+sir_decouple(2,posti_case2)))/length(posti_case2);...
%                         sum(log2(1+sir_decouple(posti_case3)))/length(posti_case3);...
                        sum(log(1+sir_decouple(1,posti_case4)))/length(posti_case4);
                        sum(log(1+sir_decouple(2,posti_case4)))/length(posti_case4)];
    cov_NLOS_decouple(1:4,j1)= [ (  length(find(sir_decouple(1,posti_case1)>=3.98))+...
                                                  length(find(sir_decouple(1,posti_case2)>=3.98)) )/total_v;...%宏基站 下行
                                                  length(find(sir_decouple(2,posti_case1)>=3.98)) /total_v;...%宏基站 行
                                              length(find(sir_decouple(1,posti_case4)>=3.98)) /total_v;...%小基站 下行
                                             ( length(find(sir_decouple(2,posti_case2)>=3.98))+...
                                                  length(find(sir_decouple(2,posti_case4)>=3.98)) )/total_v;   ]%小基站 上行
                                              
                                              
    
end

se_NLOS_decouple(1:2,iiii) =[mean(D_se);mean(U_se)]

pr_NLOS_decouple(1:3,iiii) = mean(pr_case,2) 
sir_NLOS_decouple(1:6,iiii) = mean(sir_case,2)%频谱效率

cov_pr_3_98_NLOS_decouple(1:4,iiii)=mean(cov_NLOS_decouple,2)

pr_decouple=0;
sir_decouple=0;

end

% se_case_NLOS_decouple = log2(1+sir_NLOS_decouple);
figure
se_decouple_NLOS = ((sir_NLOS_decouple(1,1:4:61)+sir_NLOS_decouple(2,1:4:61)).*pr_NLOS_decouple(1,1:4:61)...
                +(sir_NLOS_decouple(3,1:4:61)+sir_NLOS_decouple(4,1:4:61)).*pr_NLOS_decouple(2,1:4:61)+...
                (sir_NLOS_decouple(5,1:4:61)+sir_NLOS_decouple(6,1:4:61)).*pr_NLOS_decouple(3,1:4:61));%*0.8;
plot(0.1:0.4:6.1,se_decouple_NLOS);

figure 
plot(0.1:0.4:6.1,cov_pr_3_98_NLOS_decouple(:,1:4:61));

function [num ind] = find_case(A,a,b)
      ii = 1;
  for i=1:length(A)
      if A(1,i) == a && A(2,i)==b
         ind(ii)=i;
         ii= ii+1;
         
      end
  end
  num =  ii -1;

end






% prob_E1 = sum(assocc)/num_iterations;%事件1的发生概率
% prob_E2 = sum(assocv)/num_iterations;%事件2的发生概率




