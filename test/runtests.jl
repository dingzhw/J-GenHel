# test/runtests.jl

using Test

# ==============================================================================
# ## 测试套件主入口
#
# 本文件是项目测试套件的入口点。
# 它会按顺序执行所有的单元测试和集成测试。
#
# 要运行测试，请在项目根目录下打开 Julia REPL，然后执行：
#
# ```julia
# using Pkg
# Pkg.activate(".")
# Pkg.test()
# ```
#
# 或者，直接运行此文件：
#
# ```julia
# include("test/runtests.jl")
# ```
# ==============================================================================

@testset "J-GenHel 测试套件" begin
    println("开始运行工具函数测试...")
    @testset "工具函数测试 (Utility Functions)" begin
        include("test_utils.jl")
    end

    println("\n开始运行集成测试...")
    @testset "集成/冒烟测试 (Integration/Smoke Test)" begin
        include("test_integration.jl")
    end
end
