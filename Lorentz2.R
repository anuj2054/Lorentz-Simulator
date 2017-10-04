library("deSolve")
# library('rootSolve')

parameters <- span=""> c(Kt = 0.35, Kz = 0.052, alpha = 9.622e-05, beta = 0.0007755, 
    mu = 4 * 10^10, gamma = 0.06364, c1 = 0.00166, c2 = 0.00022, F0 = 6.65, 
    F1 = 2, F2 = 47.9, G0 = -3.6, G1 = 1.24, G2 = 3.81, w = 2 * 3.14 * 1, a = 0.25, 
    b = 4, Ta2 = 1, T0 = 298.15, V1 = 0.832 * 10^16, V2 = 2.592 * 10^16, V3 = 10.3 * 
        10^16)

state <- span=""> c(X = 1, Y = -1, Z = 0.2, T1 = 0.95, T2 = 0.93, T3 = 0.94, S1 = 35, 
    S2 = 36, S3 = 37)

times <- span=""> seq(0, 10000, by = 1)
# Fm=vector('numeric',length(times)) Gm=vector('numeric',length(times))
Lorenz <- span=""> function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
        # rate of change
        q <- span=""> mu * (alpha * (T2 - T1) - beta * (S2 - S1))
        Fm <- span=""> F0 + F1 * cos(w * t) + F2 * (T2 - T1)/T0
        Gm <- span=""> G0 + G1 * cos(w * t) + G2 * (T1/T0)
        Qs <- span=""> c1 + c2 * (Y^2 + Z^2)
        Ta1 <- span=""> Ta2 - gamma * X
        dX <- span=""> -Y^2 - Z^2 - a * X + a * Fm
        dY <- span=""> X * Y - b * X * Z - Y + Gm
        dZ <- span=""> b * X * Y + X * Z - Z
        dT1 <- span=""> (0.5 * q * (T2 - T3) + Kt * (Ta1 - T1) - Kz * (T1 - T3))/V1
        dT2 <- span=""> (0.5 * q * (T3 - T1) + Kt * (Ta2 - T2) - Kz * (T2 - T3))/V2
        dT3 <- span=""> (0.5 * q * (T1 - T2) + Kt * (T1 - T3) + Kz * (T2 - T3))/V3
        dS1 <- span=""> (0.5 * q * (S2 - S3) - Kz * (S1 - S3) - Qs)/V1
        dS2 <- span=""> (0.5 * q * (S3 - S1) - Kz * (S2 - S3) + Qs)/V2
        dS3 <- span=""> (0.5 * q * (S1 - S2) + Kz * (S1 - S3) + Kz * (S2 - S3))/V3

        # return the rate of change
        list(c(dX, dY, dZ, dT1, dT2, dT3, dS1, dS2, dS3))
    })  # end with(as.list ...
}

################################################ 3
out <- span=""> ode(y = state, times = times, func = Lorenz, parms = parameters)
# head(out)

layout(matrix(c(1, 2, 3), nrow = 3))
plot(out[, "X"][1:500], xlab = "time", ylab = "X")
plot(out[, "Y"][1:500], xlab = "time", ylab = "Y")
plot(out[, "Z"][1:500], xlab = "time", ylab = "Z")

