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
function [d_ge_b,d_ge_a,Y_ge_b,Y_ge_a,S_pre_b,S_pre_a,PS_pre_b,PS_pre_a] = ManhattanDistance(S_pre,PS_pre,num)
% 输入
% S_pre 为场景集 行为场景集个数，列为某一个场景的维数
% num 为最终剩余的场景个数，进行 length(S_pre(:,1)) - num 次削减
% 输出
% d_ge_b 为初始曼哈顿距离矩阵
% d_ge_a 为输出的曼哈顿距离矩阵
% Y_ge_b 为初始概率距离矩阵（列向量）
% Y_ge_a 为输出的概率距离矩阵（列向量）
% S_pre_b 为初始场景集
% S_pre_a 为削减后的场景集
% PS_pre_b 为初始概率矩阵
% PS_pre_a 为输出概率矩阵
S_pre_b = S_pre;
S_pre_a = S_pre_b;
PS_pre_b = PS_pre;
PS_pre_a = PS_pre_b;
d_ge_b = zeros(length(S_pre_b(:,1)),length(S_pre_b(:,1)));                  % 方阵
d_ge_a = d_ge_b;
Y_ge_b = zeros(length(S_pre(:,1)),1);
Y_ge_a = Y_ge_b;
%% 全局缩减
% 削减次数times
times = length(S_pre(:,1)) - num;
%% 场景削减
for epochs = 1:times
    %% 需要更新的参数：d_ge_a、Y_ge_a、S_pre_a、PS_pre_a                   %这个非常重要
    %% 只计算上半对角线距离
    for i = 1:length(S_pre_a(:,1))
        for j = i+1:length(S_pre_a(:,1))
            scence_1 = S_pre_a(i,:);                                        % 场景1
            scence_2 = S_pre_a(j,:);                                        % 场景2
            d_ge_a(i,j) = sum(abs(scence_1-scence_2));                      % 曼哈顿距离计算
            d_ge_a(j,i) = d_ge_a(i,j);
        end
    end
    for k = 1:length(S_pre_a(:,1))
        Y_ge_a(k) = PS_pre_a(k)*sum(d_ge_a(k,:));                           % 这里写d_ge(i,:)或者d_ge(:,i)都可以
    end
    %% 查找待削减场景
    [~,~,pomin1,pomin2] = min12(Y_ge_a);                                    % 查找最小的待削减场景g；查找次要小的待削减场景h
    %% 已确定待削减场景的位置编号
    %% 场景生成，只是场景对应的概率叠加了
    %% 删减 pomin2 位置的场景，保留 pomin1 位置的场景
    if (pomin1 < pomin2) && (pomin2 < length(PS_pre_a))
        PS_pre_a = [PS_pre_a(1:pomin1-1);(PS_pre_a(pomin1)+PS_pre_a(pomin2));PS_pre_a(pomin1+1:pomin2-1);PS_pre_a(pomin2+1:end)];
        d_ge_a = [d_ge_a(1:pomin2-1,1:pomin2-1),d_ge_a(1:pomin2-1,pomin2+1:end);d_ge_a(pomin2+1:end,1:pomin2-1),d_ge_a(pomin2+1:end,pomin2+1:end)];
        S_pre_a = [S_pre_a(1:pomin2-1,:);S_pre_a(pomin2+1:end,:)];
    elseif (pomin1 < pomin2) && (pomin2 == length(PS_pre_a))
        PS_pre_a = [PS_pre_a(1:pomin1-1);(PS_pre_a(pomin1)+PS_pre_a(pomin2));PS_pre_a(pomin1+1:pomin2-1)];
        d_ge_a = [d_ge_a(1:pomin2-1,1:pomin2-1)];
        S_pre_a = [S_pre_a(1:pomin2-1,:);S_pre_a(pomin2+1:end,:)];
    elseif (pomin1 > pomin2) && (pomin1 < length(PS_pre_a))
        PS_pre_a = [PS_pre_a(1:pomin2-1);PS_pre_a(pomin2+1:pomin1-1);(PS_pre_a(pomin1)+PS_pre_a(pomin2));PS_pre_a(pomin1+1:end)];
        d_ge_a = [d_ge_a(1:pomin2-1,1:pomin2-1),d_ge_a(1:pomin2-1,pomin2+1:end);d_ge_a(pomin2+1:end,1:pomin2-1),d_ge_a(pomin2+1:end,pomin2+1:end)];
        S_pre_a = [S_pre_a(1:pomin2-1,:);S_pre_a(pomin2+1:end,:)];
    elseif (pomin1 > pomin2) && (pomin1 == length(PS_pre_a))
        PS_pre_a = [PS_pre_a(1:pomin2-1);PS_pre_a(pomin2+1:pomin1-1);(PS_pre_a(pomin1)+PS_pre_a(pomin2))];
        d_ge_a = [d_ge_a(1:pomin2-1,1:pomin2-1),d_ge_a(1:pomin2-1,pomin2+1:end);d_ge_a(pomin2+1:end,1:pomin2-1),d_ge_a(pomin2+1:end,pomin2+1:end)];
        S_pre_a = [S_pre_a(1:pomin2-1,:);S_pre_a(pomin2+1:end,:)];
    else
        disp('Wrong!')
    end
    Y_ge_a = zeros(length(PS_pre_a),1);                                     % 重置 Y_ge_a
    for k = 1:length(PS_pre_a)
        Y_ge_a(k) = PS_pre_a(k)*sum(d_ge_a(k,:));
    end
end
disp('场景生成与削减完成！结果已返回！')
end
