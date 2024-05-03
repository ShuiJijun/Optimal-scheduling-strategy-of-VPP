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
function [val1,val2,pos1,pos2] = min12(S)
% min12 函数用于求解列向量或者行向量中最小的两个值以及坐标
% val1 为最小值
% val2 为次小值
[val1,pos1] = min(S);
S_except1 = S;
S_except1(pos1) = inf;
[val2,pos2] = min(S_except1);
end