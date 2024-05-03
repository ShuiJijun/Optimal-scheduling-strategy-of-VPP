%% 文章：计及风光不确定性含碳排放和碳惩罚的虚拟电厂优化调度策略
%_________________________________________________________________________%
%                                                                         %
%  Yalmip/GUROBI                                                          %
%                                                                         %
%_________________________________________________________________________%
%                                                                         %
%  Subfunction                                                            %
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
%%
function [P_tpvwtload] = ScenarioMethod
%% 数据导入
% 所有气象数据导入
% 这里的数据不开源
DataAllInLoad = readtable('气象数据.xlsx');
% 气象数据.xlsx数据文件中，第1列为时间，第2列-第9列为数据
% 第2列为环境温度（℃）
% 第3列为空气比湿度（kg/kg）
% 第4列为太阳辐射强度（W/m2）
% 第5列为土壤平均温度（℃）
% 第6列为距地面Am处东向风速（m/s）
% 第7列为距地面Bm处东向风速（m/s）
% 第8列为距地面Am处北向风速（m/s）
% 第9列为距地面Bm处北向风速（m/s）
DataAllInLoad_time = table2array(DataAllInLoad(:,1));                       % DataAllInLoad_time 存放源数据的时间列
DataAllInLoad_data = table2array(DataAllInLoad(:,2:end));                   % DataAllInLoad_data 存放源数据的数据列
% 第1列为环境温度（℃）
% 第2列为空气比湿度（kg/kg）
% 第3列为太阳辐射强度（W/m2）
% 第4列为土壤平均温度（℃）
% 第5列为距地面Am处东向风速（m/s）
% 第6列为距地面Bm处东向风速（m/s）
% 第7列为距地面Am处北向风速（m/s）
% 第8列为距地面Bm处北向风速（m/s）
DataAllInLoad_time2num = datenum(DataAllInLoad_time);                       % 将时间数据转化为数值
DataAllInLoad_time2time = 0:length(DataAllInLoad_time)-1;                   % 1h取一个点，以小时为横坐标（不以日期为横坐标）
DataAllInLoad_time2time = DataAllInLoad_time2time';
%% 分离单日数据
Temp365 = zeros(365,24);            % 单日环境温度（℃）
AirM365 = zeros(365,24);            % 单日空气比湿度（kg/kg）
IPV365 = zeros(365,24);             % 单日太阳辐射强度（W/m2）
SoTemp365 = zeros(365,24);          % 单日土壤平均温度（℃）
V2east10m365 = zeros(365,24);       % 单日距地面Am处东向风速（m/s）
V2east50m365 = zeros(365,24);       % 单日距地面Bm处东向风速（m/s）
V2north10m365 = zeros(365,24);      % 单日距地面Am处北向风速（m/s）
V2north50m365 = zeros(365,24);      % 单日距地面Bm处北向风速（m/s）
for i = 1:365
    Temp365(i,:) = DataAllInLoad_data((24*i-23):24*i,1);
    AirM365(i,:) = DataAllInLoad_data((24*i-23):24*i,2);
    IPV365(i,:) = DataAllInLoad_data((24*i-23):24*i,3);
    SoTemp365(i,:) = DataAllInLoad_data((24*i-23):24*i,4);
    V2east10m365(i,:) = DataAllInLoad_data((24*i-23):24*i,5);
    V2east50m365(i,:) = DataAllInLoad_data((24*i-23):24*i,6);
    V2north10m365(i,:) = DataAllInLoad_data((24*i-23):24*i,7);
    V2north50m365(i,:) = DataAllInLoad_data((24*i-23):24*i,8);
end
clear i
% 将所有单日数据合并，形成365*192的数据
TempAirMIPVSoTempV = [Temp365,AirM365,IPV365,SoTemp365,V2east10m365,V2east50m365,V2north10m365,V2north50m365];
clear Temp365 AirM365 IPV365 SoTemp365 V2east10m365 V2east50m365 V2north10m365 V2north50m365
%% 场景生成
%% 对光伏出力进行计算，预计1125kW的，预计需要 20000 台的光伏
% 单台光伏的参数：S_pv 光伏板面积1.96m2；miu_pv 光电能量转换率0.18；Npv 光伏数量20000
S_pv = 1.96;
miu_pv = 0.18;
Npv = 20000;
PpV_t = TempAirMIPVSoTempV(:,49:72) * S_pv * miu_pv * Npv;                  % PpV_t 为所有光伏的出力历史数据，单位：kW
disp('太阳辐照度已转换为光伏出力数据PpV_t！转换完成！365个场景24h')
%% 对风电出力进行计算
% 计算50m处的风速 V_wt（行为场景数量，列为时间）
V_wt = zeros(length(TempAirMIPVSoTempV(:,1)),24);
for j = 1:24
    for i = 1:length(TempAirMIPVSoTempV(:,1))
        V_wt(i,j) = sqrt((TempAirMIPVSoTempV(i,5*24+j))^2 + (TempAirMIPVSoTempV(i,7*24+j))^2);
    end
end
% 计算风电出力
% 设置风力发电机组的参数：Pwt_r 额定功率300kW；V_r 额定风速8m/s;V_in 切入风速4m/s；V_out 切出风速15m/s；Nwt 风电数量1
Pwt_r = 300; % kW
V_r = 8;
V_in = 4;
V_out = 15;
Pwt_t = zeros(length(PpV_t(:,1)),length(PpV_t(1,:)));                       % Pwt_t 为所有风电的出力历史数据
Nwt = 1;
for i = 1:length(V_wt(1,:))                                                 % 时间
    for j = 1:length(V_wt(:,1))                                             % 场景
        if (V_wt(j,i) < V_in) || (V_wt(j,i) >= V_out)
            Pwt_t(j,i) = Nwt * 0;
        elseif (V_wt(j,i) < V_r) && (V_wt(j,i) >= V_in)
            Pwt_t(j,i) = Nwt * (V_wt(j,i) - V_in) / (V_r - V_in) * Pwt_r;
        elseif (V_wt(j,i) < V_out) && (V_wt(j,i) >= V_r)
            Pwt_t(j,i) = Nwt * Pwt_r;
        else
            disp('错误！V_wt数据报错！')
        end
    end
end
clear i j
disp('风速已转换为风电出力数据Pwt_t！转换完成！365个场景24h')
%% 删除风电出力中存在的零的一行
key = 1;
Pwt_t_new = zeros(1,24);
PpV_t_new = zeros(1,24);
for i = 1:length(Pwt_t(:,1))
    flag = 1;
    for j = 1:length(Pwt_t(1,:))
        if Pwt_t(i,j) == 0
            flag = 0;
        end
    end
    if flag == 1
        Pwt_t_new(key,:) = Pwt_t(i,:);
        PpV_t_new(key,:) = PpV_t(i,:);
        key = key + 1;
    end
end
clear key i j flag
disp('异常值数据删减完成！')
%%
%%
D_ge = [PpV_t_new,Pwt_t_new];
k = 6;                                                                      % k 为削减后的场景数量
S_pre = D_ge;
fprintf('\n最终场景数量设置为：%d\n\n',k)
PS_pre = 1/length(PpV_t_new(:,1))*ones(length(PpV_t_new(:,1)),1)*100;       % 单位：%
%% 调用 ManhattanDistance 子函数，计算得到场景生成与削减前后的结果
[d_ge_b,d_ge_a,Y_ge_b,Y_ge_a,S_pre_b,S_pre_a,PS_pre_b,PS_pre_a] = ManhattanDistance(S_pre,PS_pre,k);
% d_ge_b 为场景生成与削减前-曼哈顿距离矩阵
% d_ge_a 为场景生成与削减后-曼哈顿距离矩阵
% Y_ge_b 为场景生成与削减前-概率距离矩阵（列向量）
% Y_ge_a 为场景生成与削减后-概率距离矩阵（列向量）
% S_pre_b 为场景生成与削减前-场景集
% S_pre_a 为场景生成与削减后-场景集
% PS_pre_b 为场景生成与削减前-场景集对应的场景概率（列向量）
% PS_pre_a 为场景生成与削减后-场景集对应的场景概率（列向量）
%% 输出典型场景的概率
fprintf('\n典型场景1,概率:%.3f%%\n',PS_pre_a(1))
fprintf('典型场景2,概率:%.3f%%\n',PS_pre_a(2))
fprintf('典型场景3,概率:%.3f%%\n',PS_pre_a(3))
fprintf('典型场景4,概率:%.3f%%\n',PS_pre_a(4))
fprintf('典型场景5,概率:%.3f%%\n',PS_pre_a(5))
fprintf('典型场景6,概率:%.3f%%\n',PS_pre_a(6))
[MostP_val,MostP_pos] = max(PS_pre_a);                                      % 确定最大概率的值 MostP_val 以及场景编号 MostP_pos
fprintf('选择典型场景%d作为算例场景，概率为：%.3f%%\n',PS_pre_a(MostP_pos),MostP_val)
%% 结果数据处理
timeFig1 = 0:23;
PpvMostP = S_pre_a(MostP_pos,1:24)/10000;                                   % 最大概率场景的光伏出力
PwtMostP = S_pre_a(MostP_pos,25:48)*3;                                      % 最大概率场景的风电出力
figure(1)
hold on
plot(timeFig1,PpvMostP,'LineWidth',3)
plot(timeFig1,PwtMostP,'LineWidth',3)
xlabel('时间/h')
ylabel('风光出力/kW')
legend('光伏出力','风电出力')
title(sprintf('概率最大的典型场景4，概率:%.3f%%\n',PS_pre_a(4)))
grid on
%% 导入用户负荷数据
%% 这里的用户负荷数据不开源
DataLoad = readmatrix('用户负荷.xlsx');
DataLoad = DataLoad * 5;
%% 绘制风光出力及负荷曲线
figure(2)
hold on
plot(timeFig1+1,PpvMostP,'r+','LineWidth',1)
plot(timeFig1+1,PwtMostP,'bo','LineWidth',0.5)
plot(timeFig1+1,DataLoad,'g*','LineWidth',1)
xlabel('时间/h')
ylabel('风光出力及负荷/kW')
legend('光伏出力','风电出力','用户负荷')
grid on
%% 阶段性结束
P_tpvwtload = [timeFig1;PpvMostP;PwtMostP;DataLoad'];                       % 第一行为时间，第二行为光伏出力；第三行为风电出力；第四行为用户负荷
writematrix(P_tpvwtload,'场景生成与削减后的风光出力及用户负荷.xlsx')
disp('成功将场景生成与削减后的风光出力及用户负荷写入文件：“场景生成与削减后的风光出力及用户负荷.xlsx”')
%%
S_pre_a(:,1:24) = S_pre_a(:,1:24)/10000;
S_pre_a(:,25:48) = S_pre_a(:,25:48)*3;
%% 输出plot
%%
figure(3)
subplot(2,3,1)
plot(timeFig1,S_pre_a(1,1:24),'LineWidth',3)
grid on
xlabel('时间/h')
ylabel('光伏出力/kW')
title(sprintf('典型场景1,概率:%.3f%%',PS_pre_a(1)))
subplot(2,3,2)
plot(timeFig1,S_pre_a(2,1:24),'LineWidth',3)
grid on
xlabel('时间/h')
ylabel('光伏出力/kW')
title(sprintf('典型场景2,概率:%.3f%%',PS_pre_a(2)))
subplot(2,3,3)
plot(timeFig1,S_pre_a(3,1:24),'LineWidth',3)
grid on
xlabel('时间/h')
ylabel('光伏出力/kW')
title(sprintf('典型场景3,概率:%.3f%%',PS_pre_a(3)))
subplot(2,3,4)
plot(timeFig1,S_pre_a(4,1:24),'LineWidth',3)
grid on
xlabel('时间/h')
ylabel('光伏出力/kW')
title(sprintf('典型场景4,概率:%.3f%%',PS_pre_a(4)))
subplot(2,3,5)
plot(timeFig1,S_pre_a(5,1:24),'LineWidth',3)
grid on
xlabel('时间/h')
ylabel('光伏出力/kW')
title(sprintf('典型场景5,概率:%.3f%%',PS_pre_a(5)))
subplot(2,3,6)
plot(timeFig1,S_pre_a(6,1:24),'LineWidth',3)
grid on
xlabel('时间/h')
ylabel('光伏出力/kW')
title(sprintf('典型场景6,概率:%.3f%%',PS_pre_a(6)))
%%
figure(4)
subplot(2,3,1)
plot(timeFig1,S_pre_a(1,25:48),'LineWidth',3)
grid on
xlabel('时间/h')
ylabel('风电出力/kW')
title(sprintf('典型场景1,概率:%.3f%%',PS_pre_a(1)))
subplot(2,3,2)
plot(timeFig1,S_pre_a(2,25:48),'LineWidth',3)
grid on
xlabel('时间/h')
ylabel('风电出力/kW')
title(sprintf('典型场景2,概率:%.3f%%',PS_pre_a(2)))
subplot(2,3,3)
plot(timeFig1,S_pre_a(3,25:48),'LineWidth',3)
grid on
xlabel('时间/h')
ylabel('风电出力/kW')
title(sprintf('典型场景3,概率:%.3f%%',PS_pre_a(3)))
subplot(2,3,4)
plot(timeFig1,S_pre_a(4,25:48),'LineWidth',3)
grid on
xlabel('时间/h')
ylabel('风电出力/kW')
title(sprintf('典型场景4,概率:%.3f%%',PS_pre_a(4)))
subplot(2,3,5)
plot(timeFig1,S_pre_a(5,25:48),'LineWidth',3)
grid on
xlabel('时间/h')
ylabel('风电出力/kW')
title(sprintf('典型场景5,概率:%.3f%%',PS_pre_a(5)))
subplot(2,3,6)
plot(timeFig1,S_pre_a(6,25:48),'LineWidth',3)
grid on
xlabel('时间/h')
ylabel('风电出力/kW')
title(sprintf('典型场景6,概率:%.3f%%',PS_pre_a(6)))
