# test/test_integration.jl

using Test

# ==============================================================================
# ## 集成测试 (Smoke Test)
#
# 这个文件包含一个高层级的“冒烟测试”。
# 它会尝试运行整个 `src/main.jl` 脚本，以确保没有发生崩溃或致命错误。
# 这有助于捕捉可能由代码更改引入的重大问题。
#
# 注意：此测试会执行完整的仿真和绘图流程，可能需要一些时间。
# ==============================================================================

@testset "主脚本 `main.jl` 冒烟测试" begin
    # 我们期望 `main.jl` 脚本能够从头到尾无错误地运行。
    # `@test_nowarn` 会捕捉任何错误或异常。
    # 由于 `main.jl` 在 `src` 目录中，我们需要进入该目录来运行它，
    # 以确保所有的 `include` 路径都正确。

    println("正在执行 `src/main.jl` ... (这可能需要一分钟左右)")

    # 切换到 src 目录执行脚本，然后切回
    original_dir = pwd()
    try
        cd("../src")
        # 使用 @test_nowarn 宏来断言代码块不会抛出任何异常
        @test_nowarn include("main.jl")
        println("`main.jl` 脚本成功运行完毕。")
    catch e
        # 如果出现错误，打印错误信息以便调试
        println("运行 `main.jl` 时发生错误:")
        showerror(stdout, e, catch_backtrace())
        # 强制测试失败
        @test false
    finally
        # 确保我们总是能切回原始目录
        cd(original_dir)
    end
end
