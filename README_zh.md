# J-GenHel：一个基于 Julia 的直升机飞行动力学仿真

J-GenHel 是一个使用 Julia 语言编写的综合性、非线性直升机飞行动力学仿真模型。它基于 NASA 的通用直升机 (GenHel) 飞行动力学模型，具体实现以类似于 UH-60 黑鹰的通用直升机为蓝本。

该代码由奥本大学的 Umberto Saetti 博士开发，并基于论文 "Flight Simulation and Control using the Julia Language" 的研究成果 (详见 [如何引用](#如何引用))。

## 主要特性

*   **高保真度模型:** 实现了机身的六自由度 (6-DoF) 刚体动力学模型。
*   **非线性空气动力学:** 为机身、旋翼桨叶和尾翼使用了非线性的空气动力学查表。
*   **精细的旋翼动力学:** 包括刚性挥舞和摆振的桨叶动力学。
*   **先进的入流模型:** 主旋翼采用三状态的 Pitt-Peters 入流模型，尾桨采用单状态模型。
*   **完整的分析套件:** 提供了一系列工具用于：
    *   在不同飞行条件下对飞机进行配平。
    *   对模型进行线性化以获得状态空间表达式。
    *   模型降阶。
    *   光谱分析 (特征值绘图)。
    *   频率响应分析。
    *   时域仿真 (包括开环和带动态逆控制器的闭环仿真)。

## 快速上手

### 环境要求

您的系统需要安装 [Julia](https://julialang.org/downloads/)。

此外，还需要以下 Julia 包：
*   `LinearAlgebra`
*   `SparseArrays`
*   `ControlSystems`
*   `MAT`
*   `Interpolations`
*   `ApproxFun`
*   `Plots`
*   `JLD`

### 安装

1.  克隆此仓库。
2.  打开 Julia REPL (交互式命令行) 并进入项目目录。
3.  安装所需的包。您可以在 REPL 中输入 `]` 进入包管理器模式，然后运行：
    ```julia
    add LinearAlgebra SparseArrays ControlSystems MAT Interpolations ApproxFun Plots JLD
    ```

### 运行仿真

主脚本 `src/main.jl` 运行一个完整的配平、分析和仿真流程。要运行它，请启动一个 Julia 会话，进入 `src` 目录，然后执行：

```julia
include("main.jl")
```

脚本将会生成多个图表，显示系统的特征值、频率响应和时间仿真结果，同时也会将输出数据保存到 `output/` 目录中。

## 项目结构

*   `src/`: 包含所有 Julia 源代码文件 (`.jl`)。
    *   `main.jl`: 运行仿真和分析的主入口脚本。
    *   `GenHel.jl`: 定义直升机运动方程的核心函数。
    *   `H60_constants.jl`: 定义了类 UH-60 模型的物理常数。
    *   其他文件包含各种辅助函数 (例如，用于控制混合、配平、线性化等)。
*   `data/`: 包含仿真所需的数据文件 (`.mat`, `.jld`)，如气动数据表、旋翼数据和其他模型参数。
*   `output/`: 用于保存仿真结果和生成数据的默认目录。

## 如何引用

如果您在研究中使用了此代码，请引用以下出版物：

> Saetti, U., and Horn, J. F., "Flight Simulation and Control using the Julia Language", AIAA Scitech Forum, San Diego, CA, Jan 3-7, 2022.
> DOI: [https://arc.aiaa.org/doi/10.2514/6.2022-2354](https://arc.aiaa.org/doi/10.2514/6.2022-2354)

原始的 UH-60 数学模型在以下文献中描述：

> Howlett, J. J., “UH-60A Black Hawk Engineering Simulation Program. Volume 1: Mathematical Model,” Tech. rep. NASA-CR-166309, 1980.

## 许可证

本项目采用 MIT 许可证。详情请见 [LICENSE](src/LICENSE) 文件。

## 作者

*   **Dr. Umberto Saetti**
    *   奥本大学航空航天工程系助理教授
    *   邮箱: saetti@auburn.edu
