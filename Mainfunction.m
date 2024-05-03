%% 文章：计及风光不确定性含碳排放和碳惩罚的虚拟电厂优化调度策略
%_________________________________________________________________________%
%                                                                         %
%  Yalmip/GUROBI                                                          %
%                                                                         %
%_________________________________________________________________________%
%                                                                         %
%  Main Function                                                          %
%  Developed in MATLAB R2021b                                             %
%                                                                         %
%  Author and programmer: Jijun Shui (E-mail: shuijijun@mail.shiep.edu.cn)%
%                                                                         %
%  Main paper:  Shui Jijun, Peng Daogang, Song Yankan, et al              %
%               Optimal Scheduling Strategy of Virtual Power Plant with Carbon Emission and Carbon Penalty Considering Uncertainty of Wind Power and Photovoltaic Power,                                     %
%               Journal of System Simulation                              %
%               DOI: 10.16182/j.issn1004731x.joss.23-0840                 %
%                                                                         %
%_________________________________________________________________________%
%% 格式化
clear
close all
clc
%% 仅考虑电负荷，不考虑热/冷负荷，即不是多能的系统
%% 场景生成与削减后的最大概率场景作为算例场景，场景数据导入
% P_tpvwtload 第一行为时间，第二行为光伏出力；第三行为风电出力；第四行为用户负荷
[P_tpvwtload] = ScenarioMethod;
time = 1:24;                       % 时间轴24h，步长为1h
%% 典型日数据拆分
%% 可以根据自己需求选择是否拆分
P_pvScenario = P_tpvwtload(2,:);   % 光伏出力数据
P_wtScenario = P_tpvwtload(3,:);   % 风电出力数据
P_loadScenario = P_tpvwtload(4,:); % 用户负荷数据
%%
T = 24;
%% 光伏模型及其约束
% 决策变量
P_pv = sdpvar(1,T);                % 光伏的实际出力
% 约束条件
C = [];
C = [C,P_pv <= P_pvScenario;       % 光伏实际出力约束
       P_pv >= 0;                  % 光伏实际出力约束
       ];
%% 风电模型及其约束
% 决策变量
P_wt = sdpvar(1,T);                % 风电的实际出力
% 约束条件
C = [C,P_wt <= P_wtScenario;       % 光伏实际出力约束
       P_wt >= 0;                  % 光伏实际出力约束
       ];
%% 蓄电池模型及其约束
% 决策变量
P_essd = sdpvar(1,T);              % 蓄电池的实际出力 kW
W_essd = sdpvar(1,T);              % 蓄电池存储能量 kWh
SOC_essd = sdpvar(1,T);            % 蓄电池荷电状态
% 初始值及模型参数设置
miu_C = 0.9;                       % 蓄电池储能效率
miu_D = 0.9;                       % 蓄电池放能效率
SOC_min = 0.2;                     % 蓄电池最小荷电状态
SOC_max = 0.95;                    % 蓄电池最大荷电状态
P_essd_CD_max = 300;               % 蓄电池最大充放功率
SOC_essd0t = 0.5;                  % 蓄电池初始荷电状态：0.5
W_essd_0t = 500;                   % 蓄电池初始储能：500 kWh
W_essd_Max = 1000;                 % 蓄电池容量 1000 kWh
theta_s = 0.05;                    % 蓄电池的自损耗率
% 约束条件
C = [C,W_essd(1) == W_essd_0t;
       SOC_essd(1) == SOC_essd0t;
       P_essd(T) == 0;
       SOC_essd == W_essd / W_essd_Max;W_essd(2:T) - W_essd(1:T-1)*(1-theta_s) == -P_essd(1:T-1) * miu_D;
       % 充电为负，放电为正
       % 不等式约束
       % 充放功率约束
       P_essd >= (-P_essd_CD_max);
       P_essd <= P_essd_CD_max;
       % 储能约束
       W_essd >= 0;
       W_essd <= W_essd_Max;
       % 荷电状态约束
       SOC_essd >= SOC_min;
       SOC_essd <= SOC_max;
       ];
%% 燃气轮机模型及其约束
% 决策变量
P_chp = sdpvar(1,T);              % 燃气轮机的实际出力 kW
Fuel_chp = sdpvar(1,T);           % 燃气轮机消耗的燃料所提供的能量 kWh
Carbon_chp = sdpvar(1,T);         % 燃气轮机的碳排放量，单位：g
P_chp_heat = sdpvar(1,T);         % 燃气轮机的供热功率
% 初始值及模型参数设置
miu_e = 0.8;                      % 发电效率
miu_h = 0.5;                      % 制热效率
CarbonVal = 500;                  % 注入碳强度，即每消耗kWh燃料所产生的碳排放
% 约束条件
C = [C,P_chp == Fuel_chp * miu_e;
       P_chp_heat == Fuel_chp * miu_h;
       Carbon_chp == Fuel_chp * CarbonVal;
       P_chp <= 1000;
       P_chp >= 100;
       P_chp_heat >= 0;
       (P_chp(2:T) - P_chp(1:T-1)) <= 500;
       (P_chp(2:T) - P_chp(1:T-1)) >= -500;
       Carbon_chp >= 0;
       Fuel_chp >= 0;
       ];
%% 向主网买/购电及其约束
% 决策变量
P_buyfromGrid = sdpvar(1,T);
P_saleforGrid = sdpvar(1,T);
% 约束条件
C = [C,P_buyfromGrid .* P_saleforGrid == 0;   % 无法在某一个时刻同时向主网买电和卖电
       P_buyfromGrid >= 0;
       P_saleforGrid >= 0;
       abs(P_buyfromGrid) <= 200;
       abs(P_saleforGrid) <= 200;
       ];
%% 
%% 用户负荷约束
% 约束条件
C = [C,P_loadScenario + P_saleforGrid == P_pv + P_wt + P_buyfromGrid + P_essd + P_chp; % 虚拟电厂功率平衡
       ];
%% 计价规则
price_load = [0.6,0.6,0.6,0.5,0.5,0.6,0.8,0.8,1,1,0.8,0.8,1.2,1.2,1.2,0.8,0.8,0.8,1,1,1,0.8,0.8,0.8];   % 用户电负荷分时电价
price_Grid = [0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.6,0.6,0.6,0.6,1,1,1,1,0.6,0.6,0.6,0.6,0.6,0.6,0.3,0.3]; % 主网上网购电分时电价
price_heat = 0.157 * ones(length(price_load),1);                            % 供热价格
% 天然气价格
price_Gas = 4.3;                                                            % 元/立方
%% 
Gas2kWh = 5;                                                                % 单位：kWh
% 碳排放成本
CarbonPerkg = 1.25;                                                         % 碳排放成本取1.25元/kg
CarbonPerkgDeC = 0.5;                                                       % 碳排放惩罚成本1.5
%% 目标函数
F1 = sdpvar(1,1);                                                           % 日前优化调度目标函数
Pi_buyfromGrid = sdpvar(1,1);                                               % 购电费用
Pi_saleforGrid = sdpvar(1,1);                                               % 卖电收益
Pi_chp = sdpvar(1,1);                                                       % 燃气轮机的燃料费用
Pi_loadScenario = sdpvar(1,1);                                              % 用户负荷的用电收益
Pi_carboncost = sdpvar(1,1);                                                % 碳排放成本
Pi_heat = sdpvar(1,1);                                                      % 供热收益
Pi_yunwei = sdpvar(1,1);% 运维成本
Pi_yunwei_pv = sdpvar(1,1);% PV运维成本
price_yunwei_pv = 0.025;
Pi_yunwei_wt = sdpvar(1,1);% WT运维成本
price_yunwei_wt = 0.015;
Pi_yunwei_eesd = sdpvar(1,1);% EESD运维成本
price_yunwei_eesd = 0.01;
Pi_carbonDeC = sdpvar(1,1);% 碳排放惩罚成本
price_yunwei_chp = 0.021;
Pi_yunwei_chp = sdpvar(1,1);% CHP运维成本
C = [C,F1 == -(Pi_saleforGrid + Pi_loadScenario - Pi_buyfromGrid - Pi_chp - Pi_carboncost + Pi_heat - Pi_yunwei - Pi_carbonDeC); % 收益约束
       Pi_buyfromGrid == sum(P_buyfromGrid .* price_Grid);                  % 购电成本
       Pi_saleforGrid == sum(P_saleforGrid .* price_Grid);                  % 上网收益
       Pi_loadScenario == sum(P_loadScenario .* price_load);                % 售电收益
       Pi_chp == sum(Fuel_chp) / Gas2kWh * price_Gas;                       % 天然气费用
       Pi_carboncost == sum(Carbon_chp) / 1000 * CarbonPerkg;               % 碳排放成本
       Pi_heat == sum(P_chp_heat) * price_heat;                             % 供热收益
       Pi_yunwei == Pi_yunwei_pv + Pi_yunwei_wt + Pi_yunwei_eesd + Pi_yunwei_chp;
       Pi_yunwei_pv == sum(P_pv) * price_yunwei_pv;
       Pi_yunwei_wt == sum(P_wt) * price_yunwei_wt;
       Pi_yunwei_eesd == sum(abs(P_essd)) * price_yunwei_eesd;
       Pi_carbonDeC == sum(Carbon_chp) / 1000 * CarbonPerkgDeC;             % 碳排放惩罚成本
       Pi_yunwei_chp == sum(sum(P_chp)+sum(P_chp_heat)) * price_yunwei_chp;
       ];
%% 求解器配置与模型求解
%% 求解器配置
ops = sdpsettings('verbose',2,'solver','gurobi','showprogress',1);
ops.gurobi.LogFile ='path_to_log_file.log';
%% 进行求解
result = optimize(C,F1,ops);
if result.problem == 0
    disp('求解成功！')
else
    disp('求解错误！')
end
%% value 值转化
% 发用功率数值转化
P_pv = value(P_pv);
P_wt = value(P_wt);
P_buyfromGrid = value(P_buyfromGrid);
P_saleforGrid = value(P_saleforGrid);
P_essd = value(P_essd);
P_chp = value(P_chp);
% 费用数值转化
F1 = value(F1);
Pi_buyfromGrid = value(Pi_buyfromGrid);
Pi_saleforGrid = value(Pi_saleforGrid);
Pi_loadScenario = value(Pi_loadScenario);
Pi_chp = value(Pi_chp);
Pi_carboncost = value(Pi_carboncost);
Pi_heat = value(Pi_heat);
Pi_yunwei = value(Pi_yunwei);
Pi_yunwei_pv = value(Pi_yunwei_pv);
Pi_yunwei_wt = value(Pi_yunwei_wt);
Pi_yunwei_eesd = value(Pi_yunwei_eesd);
Pi_carbonDeC = value(Pi_carbonDeC);
Pi_yunwei_chp = value(Pi_yunwei_chp);
%
Fuel_chp = value(Fuel_chp);
Carbon_chp = value(Carbon_chp);
W_essd = value(W_essd);
SOC_essd = value(SOC_essd);
%% 输出及绘图
fprintf('\n\n总收益（不包含碳惩罚成本）为：%.3f元\n\n',-(F1-Pi_carbonDeC))
fprintf('总收益（包含碳惩罚成本）为：%.3f元\n',-F1)
fprintf('购电成本为：%.3f元\n',Pi_buyfromGrid)
fprintf('天然气成本为：%.3f元\n',Pi_chp)
fprintf('上网收益为：%.3f元\n',Pi_saleforGrid)
fprintf('售电收益为：%.3f元\n',Pi_loadScenario)
fprintf('供热收益为：%.3f元\n',Pi_heat)
fprintf('碳排放成本为：%.3f元\n',Pi_carboncost)
fprintf('碳惩罚成本为：%.3f元\n',Pi_carbonDeC)
fprintf('运维总成本为：%.3f元\n',Pi_yunwei)
fprintf('PV运维成本为：%.3f元\n',Pi_yunwei_pv)
fprintf('WT运维成本为：%.3f元\n',Pi_yunwei_wt)
fprintf('EESD运维成本为：%.3f元\n',Pi_yunwei_eesd)
fprintf('CHP运维成本为：%.3f元\n\n',Pi_yunwei_chp)
%% 绘图
%% 根据自己的需求绘图
figure(5)
hold on
plot(time,P_pv,'LineWidth',2)
plot(time,P_wt,'LineWidth',2)
plot(time,P_chp,'LineWidth',2)
plot(time,P_essd,'LineWidth',2)
plot(time,P_loadScenario,'LineWidth',2)
plot(time,P_buyfromGrid,'LineWidth',2)
plot(time,P_saleforGrid,'LineWidth',2)
legend('光伏出力','风电出力','CHP出力','储能出力','用户电负荷','主网出力','上网')
xlabel('时间/h')
ylabel('功率/kW')
grid on
%% 所有设备的出力和用电情况
figure(6)
hold on
p1 = bar(P_loadScenario,'stacked');
p2 = bar(P_chp,'stacked');
p3 = bar(P_wt,'stacked');
p4 = bar(P_buyfromGrid,'stacked');
p5 = bar(P_pv,'stacked');
p6 = bar(P_saleforGrid,'stacked');
p7 = bar(P_essd,'stacked');
xlabel('时间/h')
ylabel('功率/kW')
lgd1=legend([p1,p2,p3,p4,p5,p6,p7],'用户电负荷','CHP出力','风电出力','主网出力','光伏出力','上网','储能出力','orientation','horizontal','location','north');
set(lgd1,'FontName','宋体','FontSize',15);
legend boxoff
%% 蓄电池的状态
figure(7)
subplot(3,1,1)
bar(P_essd,'stacked')
subtitle('蓄电池充放功率变化')
xlabel('时间/h')
ylabel('功率/kW')
grid on
subplot(3,1,2)
bar(W_essd,'stacked')
subtitle('蓄电池储能变化')
xlabel('时间/h')
ylabel('储能/kWh')
grid on
subplot(3,1,3)
bar(SOC_essd,'stacked')
subtitle('蓄电池荷电状态变化')
xlabel('时间/h')
ylabel('荷电状态/%')
grid on
%% 燃气轮机
figure(8)
subplot(2,1,1)
bar(Fuel_chp,'stacked')
subtitle('CHP燃料耗能')
xlabel('时间/h')
ylabel('能量/kWh')
grid on
subplot(2,1,2)
bar(Carbon_chp,'stacked')
subtitle('CHP碳排放量')
xlabel('时间/h')
ylabel('碳排放量/g')
grid on
%% 弃风弃光
figure(9)
hold on
plot(time,P_pvScenario)
plot(time,P_wtScenario)
plot(time,P_pv)
plot(time,P_wt)


