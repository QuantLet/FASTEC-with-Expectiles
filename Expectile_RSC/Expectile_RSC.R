delf = function(u, tau, delta) {
    abs(tau - as.numeric(u + delta < 0)) * (u + delta)^2 - abs(tau - as.numeric(u < 
        0)) * u^2 - 2 * (tau - as.numeric(u < 0)) * abs(u) * delta
}

del = seq(-0.5, 0.5, length = 50)
tau = 0.9
u   = 0.1
plot(del, delf(u = u, tau = tau1, delta = del), type = "l", lty = 2, lwd = 2, xlab = "delta", 
    ylab = "")
lines(del, delf(u = -u, tau = tau1, delta = del), lty = 3, lwd = 2)
lines(del, min(tau1, 1 - tau1) * del^2, col = "red")
