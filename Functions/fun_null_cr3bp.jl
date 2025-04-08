function fun_null_cr3bp(x0, t0, mu)
    # 初期状態にSTM (6x6 単位行列) を埋め込む
    Φ0 = Matrix{Float64}(I, 6, 6)
    X0 = vcat(x0, reshape(Φ0, 36))  # 6 + 36 = 42次元

    # 時間スパン
    tspan = (0.0, t0)

    # 解く
    prob = ODEProblem(fun_stm_cr3bp!, X0, tspan, mu)
    sol = solve(prob, Vern9(), abstol=1e-14, reltol=1e-14)

    # 最終時刻の値取得
    final_state = sol.u[end]
    X = final_state[1:6]
    Φ = reshape(final_state[7:end], 6, 6)
    dx = zeros(Float64,6)
    fun_cr3bp!(dx, X, mu, 0.0)
    f_x = dx

    DF = [
        Φ[2,1] Φ[2,3] Φ[2,5] f_x[2];
        Φ[4,1] Φ[4,3] Φ[4,5] f_x[4];
        Φ[6,1] Φ[6,3] Φ[6,5] f_x[6];
    ]

    R = rank(DF)
    # Find a non-trivial vector x that satisfies DF * x = 0.
    N = vec(nullspace(DF))
    return N, R
end
