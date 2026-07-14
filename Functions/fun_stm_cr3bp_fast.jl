
# 【使い方】
# このファイルを include するだけでよい（fun_stm_cr3bp.jl は include 不要）。
#   include("Orbital_Design_Library/Functions/fun_stm_cr3bp_fast.jl")
#   prob = ODEProblem(fun_stm_cr3bp!, X0, tspan, (mu,))       # 元ライブラリ互換の名前で呼んでもOK
#   prob = ODEProblem(fun_stm_cr3bp_fast!, X0, tspan, (mu,))  # fast の名前で呼んでもOK
# どちらの名前で呼んでも、実体はこの高速実装（fun_stm_cr3bp_fast!）が動く。
#
# 【注意】
# 同じスコールで元の fun_stm_cr3bp.jl も include すると、fun_stm_cr3bp! の
# 定義が二重になり、後から include した方で上書きされる。
#   fun_stm_cr3bp_fast.jl → fun_stm_cr3bp.jl の順で include: 最終的にスロー版が有効
#   fun_stm_cr3bp.jl → fun_stm_cr3bp_fast.jl の順で include: 最終的に高速版が有効
# 高速化を確実に効かせたい場合は、元の fun_stm_cr3bp.jl は include しないこと。

function fun_stm_cr3bp_fast!(dx, x, parameter, t)
    mu = parameter[1]
    x1, x2, x3, x4, x5, x6 = x[1], x[2], x[3], x[4], x[5], x[6]
    om = 1 - mu
    d1 = x1 + mu            # = X[1] + mu
    d2 = x1 - om            # = X[1] - 1 + mu
    r1 = sqrt(d1^2 + x2^2 + x3^2)
    r2 = sqrt(d2^2 + x2^2 + x3^2)
    r1_3 = r1^3; r2_3 = r2^3; r1_5 = r1^5; r2_5 = r2^5

    # 状態微分（fun_cr3bp! と同一式）
    dx[1] = x4; dx[2] = x5; dx[3] = x6
    dx[4] =  2x5 + x1 - om*d1/r1_3 - mu*d2/r2_3
    dx[5] = -2x4 + x2 - om*x2/r1_3 - mu*x2/r2_3
    dx[6] =           - om*x3/r1_3 - mu*x3/r2_3

    # 擬ポテンシャルのヘッセ（fun_A_cr3bp と同一式）
    Uxx = 1 - om/r1_3 - mu/r2_3 + 3*om*d1^2/r1_5 + 3*mu*d2^2/r2_5
    Uyy = 1 - om/r1_3 - mu/r2_3 + 3*om*x2^2/r1_5 + 3*mu*x2^2/r2_5
    Uzz =   - om/r1_3 - mu/r2_3 + 3*om*x3^2/r1_5 + 3*mu*x3^2/r2_5
    Uxy = 3*om*d1*x2/r1_5 + 3*mu*d2*x2/r2_5
    Uxz = 3*om*d1*x3/r1_5 + 3*mu*d2*x3/r2_5
    Uyz = 3*om*x2*x3/r1_5 + 3*mu*x2*x3/r2_5

    # STM: dPhi = A * Phi（割当ゼロ）
    A = @SMatrix [ 0.0  0.0  0.0  1.0  0.0  0.0;
                   0.0  0.0  0.0  0.0  1.0  0.0;
                   0.0  0.0  0.0  0.0  0.0  1.0;
                   Uxx  Uxy  Uxz  0.0  2.0  0.0;
                   Uxy  Uyy  Uyz -2.0  0.0  0.0;
                   Uxz  Uyz  Uzz  0.0  0.0  0.0 ]
    Phi  = reshape(view(x,  7:42), 6, 6)
    dPhi = reshape(view(dx, 7:42), 6, 6)
    mul!(dPhi, A, Phi)
    return nothing
end

# 元ライブラリ互換エイリアス:
# 補正系関数（fun_differential_correction_cr3pb!, fun_manifold_cr3bp,
# fun_multiple_shooting_cr3bp, fun_null_cr3bp 等）が呼ぶ fun_stm_cr3bp! を、
# この高速実装に委譲する。入力 (dx, x, parameter, t)・出力 dx[1:42] は
# 元の fun_stm_cr3bp! と完全に同一。
fun_stm_cr3bp!(dx, x, parameter, t) = fun_stm_cr3bp_fast!(dx, x, parameter, t)