# test/test_utils.jl

# 引入测试库和需要测试的函数
using Test
include("../src/control_mixing.jl")
include("../src/atan2.jl")
include("../src/wrapper.jl")

# ==============================================================================
# ## 单元测试：工具函数
#
# 这个文件包含了对项目中一些小型、独立的工具函数的单元测试。
# ==============================================================================

@testset "control_mixing.jl" begin
    # 测试控制混合逻辑
    # 输入一个 4x1 的向量 p，代表 [lat, lon, col, ped] 的百分比 (0-100)
    # 输出一个 4x1 的向量 controls，代表 [th1c, th1s, th0, thtr] 的实际控制量

    # 测试用例 1: 零输入
    p_zero = [0.0, 0.0, 0.0, 0.0]
    controls_zero = control_mixing(p_zero)
    @test controls_zero[1] ≈ 0.0
    @test controls_zero[2] ≈ 15.0
    @test controls_zero[3] ≈ 2.0
    @test controls_zero[4] ≈ -10.0

    # 测试用例 2: 满输入 (100%)
    p_full = [100.0, 100.0, 100.0, 100.0]
    controls_full = control_mixing(p_full)
    @test controls_full[1] ≈ -15.0
    @test controls_full[2] ≈ 0.0
    @test controls_full[3] ≈ 21.0
    @test controls_full[4] ≈ 20.0

    # 测试用例 3: 50% 输入 (中间值)
    p_mid = [50.0, 50.0, 50.0, 50.0]
    controls_mid = control_mixing(p_mid)
    @test controls_mid[1] ≈ -7.5
    @test controls_mid[2] ≈ 7.5
    @test controls_mid[3] ≈ 11.5
    @test controls_mid[4] ≈ 5.0
end

@testset "atan2.jl" begin
    # 测试自定义的 atan2 函数
    # 它应该与标准库的 atan(y, x) 行为一致
    @test atan2(1, 1) ≈ atan(1, 1)
    @test atan2(1, -1) ≈ atan(1, -1)
    @test atan2(-1, -1) ≈ atan(-1, -1)
    @test atan2(-1, 1) ≈ atan(-1, 1)
    @test atan2(0, 1) ≈ atan(0, 1)
    @test atan2(0, -1) ≈ atan(0, -1)
    @test atan2(1, 0) ≈ atan(1, 0)
    @test atan2(-1, 0) ≈ atan(-1, 0)
end

@testset "wrapper.jl" begin
    # 测试相位包装函数，确保输出在 0 到 360 度之间

    # 测试用例 1: 在范围内的值
    @test wrapper(90.0) ≈ 90.0

    # 测试用例 2: 负值
    @test wrapper(-90.0) ≈ 270.0

    # 测试用例 3: 大于 360 的值
    @test wrapper(450.0) ≈ 90.0

    # 测试用例 4: 数组输入
    input_angles = [90.0, -90.0, 450.0]
    expected_output = [90.0, 270.0, 90.0]
    @test wrapper(input_angles) ≈ expected_output
end
