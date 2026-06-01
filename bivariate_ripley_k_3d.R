################################################################################
# Bivariate Ripley's K-function in 3D — rectangular box domain
# Implementation of Statistical spatial analysis for cryo-electron tomography
# Antonio Martinez-Sanchez, et al.
# Written by Paolo Lambre, 2026
################################################################################

cap_vol <- function(r, d) {
  if (d >= r) return(0.0)
  h <- r - d
  pi * h^2 * (3*r - h) / 3
}

wedge_vol <- function(r, a, b) {
  if (a >= r || b >= r || (a^2 + b^2) >= r^2) return(0.0)
  h_u <- sqrt(r^2 - b^2)
  tryCatch({
    integrate(function(u) {
      sapply(u, function(ui) {
        rem_r2 <- r^2 - ui^2
        if (rem_r2 <= b^2) return(0.0)
        rem_r <- sqrt(rem_r2)
        theta <- 2 * acos(b / rem_r)
        0.5 * rem_r^2 * (theta - sin(theta))
      })
    }, a, h_u, rel.tol=1e-4, abs.tol=1e-6)$value
  }, error = function(e) 0.0)
}

semi_wedge_vol <- function(r, a, b, c) {
  if (a >= r || b >= r || c >= r || (a^2 + b^2 + c^2) >= r^2) return(0.0)
  u_hi <- sqrt(max(r^2 - b^2 - c^2, 0))
  if (u_hi <= a) return(0.0)
  tryCatch({
    integrate(function(u) {
      sapply(u, function(ui) {
        rem_r2 <- r^2 - ui^2
        if (rem_r2 <= (b^2 + c^2)) return(0.0)
        v_hi <- sqrt(rem_r2 - c^2)
        if (v_hi <= b) return(0.0)
        integrate(function(v) {
          sapply(v, function(vi) {
            arg <- rem_r2 - vi^2
            if (arg < c^2) return(0.0)
            sqrt(arg) - c
          })
        }, b, v_hi, rel.tol=1e-4, abs.tol=1e-6)$value
      })
    }, a, u_hi, rel.tol=1e-4, abs.tol=1e-6)$value
  }, error = function(e) 0.0)
}

clipped_sphere_volume <- function(r, point, cube_bounds) {
  dx0 <- point[1] - cube_bounds[1]; dx1 <- cube_bounds[2] - point[1]
  dy0 <- point[2] - cube_bounds[3]; dy1 <- cube_bounds[4] - point[2]
  dz0 <- point[3] - cube_bounds[5]; dz1 <- cube_bounds[6] - point[3]
  
  V <- (4/3) * pi * r^3
  
  caps <- cap_vol(r,dx0) + cap_vol(r,dx1) + cap_vol(r,dy0) + cap_vol(r,dy1) + cap_vol(r,dz0) + cap_vol(r,dz1)
  wedges <- wedge_vol(r,dx0,dy0) + wedge_vol(r,dx0,dy1) + wedge_vol(r,dx1,dy0) + wedge_vol(r,dx1,dy1) +
    wedge_vol(r,dx0,dz0) + wedge_vol(r,dx0,dz1) + wedge_vol(r,dx1,dz0) + wedge_vol(r,dx1,dz1) +
    wedge_vol(r,dy0,dz0) + wedge_vol(r,dy0,dz1) + wedge_vol(r,dy1,dz0) + wedge_vol(r,dy1,dz1)
  semi_wedges <- semi_wedge_vol(r,dx0,dy0,dz0) + semi_wedge_vol(r,dx0,dy0,dz1) +
    semi_wedge_vol(r,dx0,dy1,dz0) + semi_wedge_vol(r,dx0,dy1,dz1) +
    semi_wedge_vol(r,dx1,dy0,dz0) + semi_wedge_vol(r,dx1,dy0,dz1) +
    semi_wedge_vol(r,dx1,dy1,dz0) + semi_wedge_vol(r,dx1,dy1,dz1)
  
  max(V - caps + wedges - semi_wedges, 0.0)
}

bivariate_k_3d <- function(ref_points, eval_points, r_values, cube_bounds, precomputed_vols = NULL) {
  domain_vol <- (cube_bounds[2]-cube_bounds[1]) * (cube_bounds[4]-cube_bounds[3]) * (cube_bounds[6]-cube_bounds[5])
  lambda_e   <- nrow(eval_points) / domain_vol
  
  K <- numeric(length(r_values))
  L <- numeric(length(r_values))
  H <- numeric(length(r_values))
  
  dist_mat <- sqrt(
    outer(ref_points[,1], eval_points[,1], "-")^2 +
      outer(ref_points[,2], eval_points[,2], "-")^2 +
      outer(ref_points[,3], eval_points[,3], "-")^2
  )
  
  for (ri in seq_along(r_values)) {
    r <- r_values[ri]
    total_count <- sum(dist_mat <= r)
    
    if (!is.null(precomputed_vols)) {
      total_vol <- precomputed_vols[ri]
    } else {
      total_vol <- sum(apply(ref_points, 1, function(x) clipped_sphere_volume(r, x, cube_bounds)))
    }
    
    if (total_vol > 0) {
      K[ri] <- ((4 * pi * r^3) / 3) * total_count / (lambda_e * total_vol)
    } else {
      K[ri] <- NA_real_
    }
    L[ri] <- (3 * K[ri] / (4 * pi))^(1/3)
  }
  H <- L - r_values
  list(r = r_values, K = K, L = L, H = H)
}

csr_envelope <- function(ref_points, n_eval, r_values, cube_bounds, n_sims = 19, conf = 0.95, seed = 42) {
  set.seed(seed)
  K_sims <- matrix(NA_real_, nrow = n_sims, ncol = length(r_values))
  L_sims <- matrix(NA_real_, nrow = n_sims, ncol = length(r_values))
  H_sims <- matrix(NA_real_, nrow = n_sims, ncol = length(r_values))
  
  precomputed_vols <- numeric(length(r_values))
  for (ri in seq_along(r_values)) {
    precomputed_vols[ri] <- sum(apply(ref_points, 1, function(x) clipped_sphere_volume(r_values[ri], x, cube_bounds)))
  }
  
  for (s in seq_len(n_sims)) {
    sim_eval <- cbind(
      runif(n_eval, cube_bounds[1], cube_bounds[2]),
      runif(n_eval, cube_bounds[3], cube_bounds[4]),
      runif(n_eval, cube_bounds[5], cube_bounds[6])
    )
    res <- bivariate_k_3d(ref_points, sim_eval, r_values, cube_bounds, precomputed_vols = precomputed_vols)
    K_sims[s, ] <- res$K
    L_sims[s, ] <- res$L
    H_sims[s, ] <- res$H
  }
  
  lo <- (1 - conf) / 2; hi <- 1 - lo
  list(
    r = r_values,
    K_low  = apply(K_sims, 2, quantile, probs = lo, na.rm = TRUE),
    K_high = apply(K_sims, 2, quantile, probs = hi, na.rm = TRUE),
    K_mean = colMeans(K_sims, na.rm = TRUE),
    L_low  = apply(L_sims, 2, quantile, probs = lo, na.rm = TRUE),
    L_high = apply(L_sims, 2, quantile, probs = hi, na.rm = TRUE),
    L_mean = colMeans(L_sims, na.rm = TRUE),
    H_low  = apply(L_sims, 2, quantile, probs = lo, na.rm = TRUE),
    H_high = apply(L_sims, 2, quantile, probs = hi, na.rm = TRUE),
    H_mean = colMeans(L_sims, na.rm = TRUE),
    
    K_sims = K_sims,
    L_sims = L_sims,
    H_sims = H_sims
  )
}

global_rank_envelope <- function(res, env, alpha = 0.05) {
  rank_envelope_test <- function(Tlist, alpha = 0.05) {
    s_plus_1 <- length(Tlist)
    s <- s_plus_1 - 1
    m <- length(Tlist[[1]])
    
    R_lower <- matrix(NA_real_, nrow = s_plus_1, ncol = m)
    R_upper <- matrix(NA_real_, nrow = s_plus_1, ncol = m)
    
    for (j in seq_len(m)) {
      vals <- sapply(Tlist, function(x) x[j])
      R_lower[, j] <- rank(vals, ties.method = "average")
      R_upper[, j] <- rank(-vals, ties.method = "average")
    }
    
    R_star <- pmin(R_lower, R_upper)
    
    R_extreme <- apply(R_star, 1, min)
    
    R1 <- R_extreme[1]
    p_minus <- sum(R_extreme < R1) / (s + 1)
    p_plus  <- sum(R_extreme <= R1) / (s + 1)
    
    a_k <- sapply(1:max(R_extreme), function(k) sum(R_extreme < k) / (s + 1))
    k_alpha <- max(which(a_k < alpha))
    
    Tmat <- do.call(rbind, Tlist)
    
    T_low <- apply(Tmat, 2, function(col) sort(col)[k_alpha])
    T_upp <- apply(Tmat, 2, function(col) sort(col, decreasing = TRUE)[k_alpha])
    
    T1 <- Tlist[[1]]
    outside <- any(T1 < T_low | T1 > T_upp)
    
    decision <- if (p_plus <= alpha) {
      "Reject H0"
    } else if (p_minus > alpha) {
      "Do not reject H0"
    } else {
      "Borderline (touching envelope)"
    }
    
    list(
      p_interval = c(p_minus = p_minus, p_plus = p_plus),
      k_alpha = k_alpha,
      envelope = list(lower = T_low, upper = T_upp),
      R_extreme = R_extreme,
      decision = decision
    )
  }

  mat_to_list <- function(mat) split(mat, row(mat))

  Tlist_K <- c(list(res$K), mat_to_list(env$K_sims))
  Tlist_L <- c(list(res$L), mat_to_list(env$L_sims))
  Tlist_H <- c(list(res$H), mat_to_list(env$H_sims))

  out_K <- rank_envelope_test(Tlist_K, alpha)
  out_L <- rank_envelope_test(Tlist_L, alpha)
  out_H <- rank_envelope_test(Tlist_H, alpha)

  list(
    alpha = alpha,
    r = env$r,
    
    K = list(
      p_interval = out_K$p_interval,
      decision   = out_K$decision,
      envelope   = out_K$envelope,
      R_extreme  = out_K$R_extreme
    ),
    
    L = list(
      p_interval = out_L$p_interval,
      decision   = out_L$decision,
      envelope   = out_L$envelope,
      R_extreme  = out_L$R_extreme
    ),
    
    H = list(
      p_interval = out_H$p_interval,
      decision   = out_H$decision,
      envelope   = out_H$envelope,
      R_extreme  = out_H$R_extreme
    )
  )
}

plot_bivariate_k <- function(result, envelope = NULL, type = "L", ...) {
  r <- result$r
  y <- switch(type,
              K = result$K,
              L = result$L,
              H = result$H
  )
  y0 <- switch(type,
               K = (4 * pi * r^3 / 3),
               L = r,
               H = rep(0, length(r))
  )
  
  y_all <- c(y, y0)
  if (!is.null(envelope)) {
    elo <- switch(type, K = envelope$K_low, L = envelope$L_low, H = envelope$H_low)
    ehi <- switch(type, K = envelope$K_high, L = envelope$L_high, H = envelope$H_high)
    y_all <- c(y_all, elo, ehi)
  }
  
  ylab <- switch(type,
                 K = expression(K^re*(r)),
                 L = expression(L^re*(r)),
                 H = expression(H^re*(r))
  )
  title <- switch(type,
                  K = "Bivariate K(r)",
                  L = "Bivariate L(r) [uncentered]",
                  H = "Bivariate H(r) [centered L]"
  )
  
  plot(r, y, type = "l", col = "blue", lwd = 2,
       ylim = range(y_all, na.rm = TRUE),
       xlab = "r", ylab = ylab, main = title, ...)
  
  if (!is.null(envelope)) {
    polygon(c(r, rev(r)), c(elo, rev(ehi)),
            col = adjustcolor("grey70", alpha.f = 0.4), border = NA)
    lines(r, elo, col = "grey50", lty = 2)
    lines(r, ehi, col = "grey50", lty = 2)
  }
  
  lines(r, y0, col = "black", lty = 3, lwd = 1.5)
  lines(r, y,  col = "blue",  lwd = 2)
  
  legend("topleft",
         legend = c("Observed", "CSR expected", "95% envelope"),
         col    = c("blue", "black", "grey50"),
         lty    = c(1, 3, 2), lwd = c(2, 1.5, 1),
         bty    = "n")
}
