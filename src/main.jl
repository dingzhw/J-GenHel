# TO RUN:
# include("main.jl")
#
# AUTHOR
# Dr. Umberto Saetti, Assistant Professor, Department of Aerospace Engineering,
# Auburn university
# saetti@auburn.edu
#
# LAST UPDATED
# 1/3/2022
#
# DESCRIPTION
# Complex helicopter model based on the General Helicopter (GenHel) flight
# dynamics simulation model (Ref. 1). This model is representative of a utility
# helicopter similar to a UH-60. The model contains a 6-DoF rigid-body dynamic
# model of the fuselage, nonlinear aerodynamic lookup tables for the fuselage,
# rotor blades, and empennage, rigid flap and lead-lag rotor blade dynamics,
# a three-state Pitt-Peters inflow model, and a one-state tail rotor model.
# The trim, linearization routines, and overall architecture of the simulation
# follows the teachings of Dr. Joe Horn at Penn State. If using the code,
# please cite Ref. 2.
#
# REFERENCES
# 1) Howlett, J. J., “UH-60A Black Hawk Engineering Simulation Program.
#    Volume 1: Mathematical Model,” Tech. rep. NASA-CR-166309, 1980.
# 2) Saetti, U., and Horn, J. F., "Flight Simulation and Control using the Julia
#    Language", AIAA Scitech Forum, San Diego, CA, Jan 3-7, 2022.
#    DOI: https://arc.aiaa.org/doi/10.2514/6.2022-2354.
#
# STATES             | NINDICES    | DESCRIPTION
# ______________________________________________________________________________
# u v w              |  1  2  3    | body velocities [ft/s]
# p q r              |  4  5  6    | angula rates [rad]
# phi theta psi      |  7  8  9    | Euler angles [rad]
# x  y  z            | 10 11 12    | position [ft]
# b0 b1c b1s b0d     | 13 14 15 16 | flapping angles [rad]
# db0 db1c db1s db0d | 17 18 19 20 | flapping angles derivatives [rad/s]
# z0 z1c z1s z0d     | 21 22 23 24 | lead-lag angles [rad]
# dz0 dz1c dz1s dz0d | 25 26 27 28 | lead-lag angles derivative [rad/s]
# lam0 lam1c lam1s   | 29 30 31    | inflow angles [rad]
# psi                | 32          | rotor azimuth [rad]
# dynamic twist      | 33          | dynamic twist force? [?]
# lam0T              | 34          | tail rotor inflow [rad]
#
# CONTROL INPUTS     | INDICED     | DESCRIPTION
# ______________________________________________________________________________
# delta_lat          | 1           | lateral stick
# delta_lon          | 2           | longitudinal stick
# delta_col          | 3           | collective stick
# delta_ped          | 4           | pedals
#
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
#
# ## 脚本设置与包导入
#
# ------------------------------------------------------------------------------

# 导入所需的 Julia 包
# using DifferentialEquations # 微分方程求解器，当前未启用
using LinearAlgebra     # 线性代数工具
using SparseArrays      # 稀疏数组支持
using ControlSystems    # 控制系统库，用于传递函数和状态空间模型
using MAT               # 读取 MATLAB 的 .mat 文件
using Interpolations    # 插值计算
using ApproxFun         # 函数逼近
using Printf            # 格式化输出
# ENV["MPLBACKEND"]="tkagg" # 设置绘图后端，当前未启用
using Plots             # 强大的绘图库
using JLD               # Julia 数据格式，用于保存和加载变量
# using PyPlot            # 基于 Matplotlib 的绘图库，当前未启用
# using GR                # 另一个绘图后端，当前未启用

# 打印程序标题和作者信息
@printf("\n                  J-GenHel           \n")
@printf("\n          Author: Dr. Umberto Saetti\n")
@printf("\n-----------------------------------------------\n")

# ------------------------------------------------------------------------------
#
# ## 导入自定义函数与模型常数
#
# ------------------------------------------------------------------------------

# 导入项目中的各个功能模块
include("control_mixing.jl")    # 控制输入混合逻辑
include("table_lookup.jl")      # 空气动力学查表功能
include("stabSched.jl")         # 安定面调度逻辑
include("atan2.jl")             # 自定义 atan2 函数
include("eqnmot.jl")            # 运动方程
include("interpp1.jl")          # 一维插值
include("linearize.jl")         # 线性化模块
include("interpp2.jl")          # 二维插值
include("wrapper.jl")           # 频率响应相位角包装函数
include("trimmer.jl")           # 配平计算模块
include("ModRed.jl")            # 模型降阶（通用）
include("ModRed8.jl")           # 8 状态模型降阶
include("ModRed10.jl")          # 10 状态模型降阶
include("GenHel.jl")            # GenHel 直升机核心动力学模型
include("GenHel_DI.jl")         # GenHel 动态反转闭环动力学模型
include("rk4.jl")               # 四阶龙格-库塔积分器
include("finp.jl")              # 控制输入函数
include("simulate.jl")          # 开环仿真模块
include("simulate_DI.jl")       # 闭环仿真模块

# 加载 H60 (UH-60) 直升机模型的物理常数
include("H60_constants.jl")

# ------------------------------------------------------------------------------
#
# ## 模型配平与线性化
#
# ------------------------------------------------------------------------------

# 设置配平飞行条件
VXTRIM=80*1.688   # 前飞速度 [ft/s] (80 kts)
VYTRIM=0.0        # 侧飞速度 [ft/s]
VZTRIM=0.0        # 垂直速度 [ft/s]
PSIDTRIM=0.0      # 偏航速率 [rad/s]

# 设置初始姿态和位置
PSITRIM = 0.0*D2R # 初始偏航角 [rad]
ALTTRIM = 0.0     # 初始高度 [ft]
XNTRIM = 0.0      # 初始北向位置 [ft]
YETRIM = 0.0      # 初始东向位置 [ft]

# 初始化配平计算所需的变量
include("trimInit.jl")

# 调用配平函数计算配平状态(x0)和控制量(u0)
x0, u0 = trimmer(GenHel!,x0,u0,xdot_targ,t)

# 在配平点进行线性化，得到状态矩阵 A 和控制矩阵 B
A, B = linearize(GenHel!,x0,u0,xdot0,t)

# ------------------------------------------------------------------------------
#
# ## 光谱分析 (特征值分析)
#
# ------------------------------------------------------------------------------

# 计算全阶模型的特征值
eigs=eigvals(A)

# 将模型降至 8 状态
sysRed = ModRed8(A,B)
A8=sysRed[1]
B8=sysRed[2]
eigs8=eigvals(A8)

# 将模型降至 10 状态
sysRed = ModRed10(A,B)
A10=sysRed[1]
B10=sysRed[2]
eigs10=eigvals(A10)

# (可选) 保存状态矩阵和特征值到 JLD 文件
#save("output/AB_J-GenHel.jld","A",A,"B",B,"A8",A8,"B8",B8)
#save("output/eigs_J-GenHel.jld","eigs",eigs,"eigs8",eigs8)

# 绘制不同模型的特征值以进行比较
gr(size=(800,600)) # 设置绘图窗口大小
plt_eigs=plot(real(eigs),imag(eigs),seriestype=:scatter,
      label="全阶模型", legendfontsize=12, legend=:topleft,
      marker = (:circle, :royalblue, 6), markerstrokecolor = :royalblue)
plot!(real(eigs10),imag(eigs10),seriestype=:scatter,
      label="降阶模型 (10-状态)", legendfontsize=12,
      marker = (:utriangle, :brown3, 6), markerstrokecolor = :brown3,
      reuse = false, color=:purple)
plot!(real(eigs8),imag(eigs8),seriestype=:scatter,
      label="降阶模型 (8-状态)", legendfontsize=12,
      marker = (:star5, :forestgreen, 6), markerstrokecolor = :forestgreen)
display(plt_eigs)
xaxis!("实部",xguidefontsize=14)
yaxis!("虚部",yguidefontsize=14)

# ------------------------------------------------------------------------------
#
# ## 频率响应分析
#
# ------------------------------------------------------------------------------

# 创建状态空间模型 (全阶)
sys=ss(A,B,Matrix{Float64}(I, 37, 37),zeros(37,4))
pdlat=tf(sys[4,1]) # 提取横向杆输入到滚转角速度的传递函数

# 创建状态空间模型 (10-状态)
sys10=ss(A10,B10,Matrix{Float64}(I, 10, 10),zeros(10,4))
pdlat10=tf(sys10[4,1])

# 创建状态空间模型 (8-状态)
sys8=ss(A8,B8,Matrix{Float64}(I, 8, 8),zeros(8,4))
pdlat8=tf(sys8[4,1])

# 计算伯德图数据 (幅频和相频响应)
mag, phase, w = bode(pdlat, 10.0.^(range(-1,stop=2,length=1000)))
mag10, phase10, w10 = bode(pdlat10, 10.0.^(range(-1,stop=2,length=1000)))
mag8, phase8, w8 = bode(pdlat8, 10.0.^(range(-1,stop=2,length=1000)))

# 绘制伯德图进行比较
gr(size=(800,600))
p1_mag=plot(vec(w), 20*log10.(vec(mag)/perc2in[1,1]),
      label="全阶模型", reuse = false, line=(3, :royalblue, :solid),
      legendfontsize=12)
plot!(vec(w10), 20*log10.(vec(mag10)/perc2in[1,1]), label="降阶模型 (10-状态)",
                          reuse = false, line=(3, :brown3, :dash))
plot!(vec(w8), 20*log10.(vec(mag8)/perc2in[1,1]), label="降阶模型 (8-状态)",
                          reuse = false, line=(3, :forestgreen, :dashdot))
xaxis!(:log10)
yaxis!("幅值 [dB]",yguidefontsize=14)
xlims!((0.1,30))
ylims!((-50,0))

p1_phase=plot(vec(w), wrapper(vec(phase)).-360, label="", line=(3, :royalblue, :solid))
plot!(vec(w10), wrapper(vec(phase10)).-360, label="", line=(3, :brown3, :dash))
plot!(vec(w8), wrapper(vec(phase8)).-360, label="", line=(3, :forestgreen, :dashdot))
yaxis!("相位 [deg]",yguidefontsize=14)
xaxis!(:log10)
xlims!((0.1,30))

p1=plot(p1_mag, p1_phase, layout=(2,1), leftmargin=3Plots.mm)
display(p1)

# ------------------------------------------------------------------------------
#
# ## 时域仿真
#
# ------------------------------------------------------------------------------

# 设置仿真参数
dt=0.01   # 时间步长 [s]
Tsim=10   # 仿真总时长 [s]

# 运行开环仿真
state_OL, time_OL = simulate(GenHel!,finp,Tsim,dt,x0,xdot0,u0)

# 设置动态反转（闭环）控制器参数
NDISTATES=9 # 动态反转控制器的状态数量

# 构建闭环仿真的初始状态向量
x0_DI=[x0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; u0[1]; u0[2]; u0[4]]

# 构建闭环仿真的初始状态导数向量
xdot0_DI=[xdot0; zeros(NDISTATES,1)]

# 运行闭环仿真
state_CL, time_CL = simulate_DI(GenHel_DI!,finp,Tsim,dt,x0,xdot0,u0)

# 将开环和闭环的仿真结果保存到 JLD 文件
save("output/resp_GenHel.jld","time_OL",time_OL,"time_CL",time_CL,
    "state_OL",state_OL,"state_CL",state_CL)

# ------------------------------------------------------------------------------
#
# ## 仿真结果绘图
#
# ------------------------------------------------------------------------------

# 绘制姿态角（滚转、俯仰、偏航）响应
gr(size=(800,600))
i=7 # 状态索引：滚转角 phi
p1=plot(time_OL,state_OL[i,:]*180/pi,label="开环",legend=:topright,
    legendfontsize=12,line=(3, :royalblue, :solid))
plot!(time_OL,state_CL[i,:]*180/pi,label="闭环",legend=:topright,
    legendfontsize=12,line=(3, :brown3, :dash))
yaxis!("φ (滚转角) [deg]",yguidefontsize=14)
xlims!((0,10))

i=8 # 状态索引：俯仰角 theta
p2=plot(time_OL,state_OL[i,:]*180/pi,label="",line=(3, :royalblue, :solid))
plot!(time_OL,state_CL[i,:]*180/pi,label="",line=(3, :brown3, :dash))
yaxis!("θ (俯仰角) [deg]",yguidefontsize=14)
xlims!((0,10))

i=9 # 状态索引：偏航角 psi
p3=plot(time_OL,state_OL[i,:]*180/pi,label="",line=(3, :royalblue, :solid))
plot!(time_OL,state_CL[i,:]*180/pi,label="",line=(3, :brown3, :dash))
xaxis!("时间 [s]",xguidefontsize=14)
yaxis!("ψ (偏航角) [deg]",yguidefontsize=14)
xlims!((0,10))

p4=plot(p1, p2, p3, layout=(3,1), leftmargin=3Plots.mm)
display(p4)
# (可选) 保存图像
#savefig(".\\Plots\\att_J-GenHel.svg")

# 绘制角速率（p, q, r）响应
gr(size=(800,600))
i=4 # 状态索引：滚转角速率 p
p5=plot(time_OL,state_OL[i,:],label="开环",legend=:bottomleft,
    legendfontsize=12,line=(3, :royalblue, :solid))
plot!(time_OL,state_CL[i,:],label="闭环",legend=:bottomleft,
    legendfontsize=12,line=(3, :brown3, :dash))
yaxis!("p (滚转角速率) [rad/s]",yguidefontsize=14)
xlims!((0,10))

i=5 # 状态索引：俯仰角速率 q
p6=plot(time_OL,state_OL[i,:],label="",line=(3, :royalblue, :solid))
plot!(time_OL,state_CL[i,:],label="",line=(3, :brown3, :dash))
yaxis!("q (俯仰角速率) [rad/s]",yguidefontsize=14)
xlims!((0,10))

i=6 # 状态索引：偏航角速率 r
p7=plot(time_OL,state_OL[i,:],label="",line=(3, :royalblue, :solid))
plot!(time_OL,state_CL[i,:],label="",line=(3, :brown3, :dash))
xaxis!("时间 [s]",xguidefontsize=14)
yaxis!("r (偏航角速率) [rad/s]",yguidefontsize=14)
xlims!((0,10))

p8=plot(p5, p6, p7, layout=(3,1), leftmargin=3Plots.mm)
display(p8)
# (可选) 保存图像
#savefig(".\\Plots\\ang_J-GenHel.svg")
