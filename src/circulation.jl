"""
    convergence!(T1, fields::ClimateFields, timestate, ws::CirculationWorkspace)

Moisture flux convergence from `T1` (specific humidity, `[kg/kg]`) and the
current `fields.omegaclim` (vertical velocity), writing the tendency into
`ws.dX_conv`. Implements Eq. 18 from Stassen et al. (2019).
"""
function convergence!(T1, fields::ClimateFields, timestate, ws::CirculationWorkspace)
    omega = @view fields.omegaclim[:, :, timestate.ityr]

    @. ws.dX_conv = -T1 * omega * const_factor
    return nothing
end

"""
    diffusion!(T1, h_scl, fields::ClimateFields, ws::CirculationWorkspace, timestate)

Meridional + zonal diffusion of `T1` (temperature or humidity), writing the
tendency into `ws.dX_diff`. `h_scl` (`z_air` or `z_vapor`) selects the
topographic weighting field.
"""
function diffusion!(T1, h_scl, fields::ClimateFields, ws::CirculationWorkspace, timestate)
    # Zero output buffer (we will accumulate into it)
    fill!(ws.dX_diff, 0.0)

    # Topographic scaling (choose based on scale height)
    wz = if h_scl == z_air
        fields.wz_air
    elseif h_scl == z_vapor
        fields.wz_vapor
    else
        error("Invalid h_scl = $h_scl (must be z_air or z_vapor)")
    end

    # Precomputed geometry/coefficients
    dxlat = dxlat_grid
    ccy = ccy_diff
    ccx = ccx_diff

    # Pre‑cached neighbour indices
    jm1, jp1 = lon_jm1, lon_jp1
    jm2, jp2 = lon_jm2, lon_jp2
    jm3, jp3 = lon_jm3, lon_jp3

    # ----- Precompute k‑independent terms for the poles -----
    # For k == 1 (North Pole)
    @views @. ws.term_north = ccy * wz[:, 2] * (T1[:, 2] - T1[:, 1])
    # For k == ydim (South Pole)
    @views @. ws.term_south = ccy * wz[:, ydim-1] * (T1[:, ydim-1] - T1[:, ydim])

    for k in 1:ydim
        # ----- Meridional diffusion -----
        if k == 1
            @turbo for i in 1:xdim
                ws.dX_diff[i, k] += wz[i, k] * ws.term_north[i]
            end
        elseif k == ydim
            @turbo for i in 1:xdim
                ws.dX_diff[i, k] += wz[i, k] * ws.term_south[i]
            end
        else
            # Mid‑latitudes: no precomputation possible (depends on k‑1, k+1)
            @turbo for i in 1:xdim
                ws.dX_diff[i, k] += wz[i, k] * ccy * (
                    wz[i, k-1] * (T1[i, k-1] - T1[i, k]) +
                    wz[i, k+1] * (T1[i, k+1] - T1[i, k])
                )
            end
        end

        # ----- Zonal diffusion -----
        if dxlat[k] > 2.5e5   # mid‑latitudes, normal time step
            @turbo for j in 1:xdim
                jm1v = jm1[j];
                jp1v = jp1[j]
                jm2v = jm2[j];
                jp2v = jp2[j]
                jm3v = jm3[j];
                jp3v = jp3[j]

                dTx = ccx[k] * 0.05 * (
                    10.0 * (wz[jm1v, k] * (T1[jm1v, k] - T1[j, k]) +
                            wz[jp1v, k] * (T1[jp1v, k] - T1[j, k])) +
                    4.0 * (wz[jm2v, k] * (T1[jm2v, k] - T1[jm1v, k]) +
                           wz[jm1v, k] * (T1[j, k] - T1[jm1v, k])) +
                    4.0 * (wz[jp1v, k] * (T1[j, k] - T1[jp1v, k]) +
                           wz[jp2v, k] * (T1[jp2v, k] - T1[jp1v, k])) +
                    1.0 * (wz[jm3v, k] * (T1[jm3v, k] - T1[jm2v, k]) +
                           wz[jm2v, k] * (T1[jm1v, k] - T1[jm2v, k])) +
                    1.0 * (wz[jp2v, k] * (T1[jp1v, k] - T1[jp2v, k]) +
                           wz[jp3v, k] * (T1[jp3v, k] - T1[jp2v, k]))
                )
                ws.dX_diff[j, k] += wz[j, k] * dTx
            end
        else   # polar regions – sub‑timestepping
            # Number of sub‑steps for stability (precomputed, depends only on k)
            time2 = POLAR_DIFF_TIME2[k]
            ccx2 = POLAR_DIFF_CCX2[k]

            # Copy current row into temporary buffer
            ws.T1h .= @view T1[:, k]

            for _ in 1:time2
                # Jacobi
                @turbo for j in 1:xdim
                    jm1v = jm1[j];
                    jp1v = jp1[j]
                    jm2v = jm2[j];
                    jp2v = jp2[j]
                    jm3v = jm3[j];
                    jp3v = jp3[j]

                    ws.dTxh[j] = ccx2 * 0.05 * (
                        10.0 * (wz[jm1v, k] * (ws.T1h[jm1v] - ws.T1h[j]) +
                                wz[jp1v, k] * (ws.T1h[jp1v] - ws.T1h[j])) +
                        4.0 * (wz[jm2v, k] * (ws.T1h[jm2v] - ws.T1h[jm1v]) +
                               wz[jm1v, k] * (ws.T1h[j] - ws.T1h[jm1v])) +
                        4.0 * (wz[jp1v, k] * (ws.T1h[j] - ws.T1h[jp1v]) +
                               wz[jp2v, k] * (ws.T1h[jp2v] - ws.T1h[jp1v])) +
                        1.0 * (wz[jm3v, k] * (ws.T1h[jm3v] - ws.T1h[jm2v]) +
                               wz[jm2v, k] * (ws.T1h[jm1v] - ws.T1h[jm2v])) +
                        1.0 * (wz[jp2v, k] * (ws.T1h[jp1v] - ws.T1h[jp2v]) +
                               wz[jp3v, k] * (ws.T1h[jp3v] - ws.T1h[jp2v]))
                    )
                end
                @turbo for j in 1:xdim
                    dq = ifelse(ws.dTxh[j] <= -ws.T1h[j], -0.9 * ws.T1h[j], ws.dTxh[j])
                    ws.T1h[j] += dq
                end
            end

            # Add total change (scaled by outer wz) to output buffer
            @views @. ws.dX_diff[:, k] += wz[:, k] * (ws.T1h - T1[:, k])
        end
    end

    return nothing
end

"""
    advection!(T1, h_scl, fields::ClimateFields, ws::CirculationWorkspace, timestate, cfg::PhysicsConfig)

Meridional + zonal advection of `T1` (temperature or humidity), writing the
tendency into `ws.dX_adv`. Gated by `cfg.log_hadv`/`cfg.log_vadv` depending on
`h_scl`.
"""
function advection!(T1, h_scl, fields::ClimateFields, ws::CirculationWorkspace, timestate, cfg::PhysicsConfig)
    # Disable advection for water vapour or heat according to switches
    if (h_scl == z_vapor && !cfg.log_vadv) || (h_scl == z_air && !cfg.log_hadv)
        fill!(ws.dX_adv, 0.0)
        return nothing
    end

    # Pre-zero the output buffer (we will accumulate into it)
    fill!(ws.dX_adv, 0.0)

    # Extract 2D views for current time step
    vclim_p_t = @view fields.vclim_p[:, :, timestate.ityr]
    vclim_m_t = @view fields.vclim_m[:, :, timestate.ityr]
    uclim_p_t = @view fields.uclim_p[:, :, timestate.ityr]
    uclim_m_t = @view fields.uclim_m[:, :, timestate.ityr]

    # Topographic scaling (choose based on scale height)
    wz = if h_scl == z_air
        fields.wz_air
    elseif h_scl == z_vapor
        fields.wz_vapor
    else
        error("Invalid h_scl = $h_scl (must be z_air or z_vapor)")
    end

    # Precomputed constants
    dxlat = dxlat_grid
    ccy = ccy_adv
    ccx = ccx_adv
    is_polar = IS_POLAR

    @inbounds for k in 1:ydim
        # ----- Meridional (v) advection -----
        if k == 1          # North Pole
            @turbo for j in 1:xdim
                v_p = vclim_p_t[j, k]
                ws.dX_adv[j, k] += ccy * v_p * (
                    wz[j, 2] * (T1[j, 1] - T1[j, 2]) +
                    wz[j, 3] * (T1[j, 1] - T1[j, 3])
                ) / 3.0
            end
        elseif k == 2
            @turbo for j in 1:xdim
                v_m = vclim_m_t[j, k]
                v_p = vclim_p_t[j, k]
                ws.dX_adv[j, k] += ccy * (
                    -v_m * wz[j, 1] * (T1[j, 2] - T1[j, 1]) +
                    v_p * (wz[j, 3] * (T1[j, 2] - T1[j, 3]) +
                           wz[j, 4] * (T1[j, 2] - T1[j, 4])) / 3.0
                )
            end
        elseif k >= 3 && k <= ydim-2
            @turbo for j in 1:xdim
                km1, km2 = k-1, k-2
                kp1, kp2 = k+1, k+2
                v_m = vclim_m_t[j, k]
                v_p = vclim_p_t[j, k]
                ws.dX_adv[j, k] += ccy * (
                    -v_m * (wz[j, km1] * (T1[j, k] - T1[j, km1]) +
                            wz[j, km2] * (T1[j, k] - T1[j, km2])) +
                    v_p * (wz[j, kp1] * (T1[j, k] - T1[j, kp1]) +
                           wz[j, kp2] * (T1[j, k] - T1[j, kp2]))
                ) / 3.0
            end
        elseif k == ydim-1
            @turbo for j in 1:xdim
                km1, km2 = k-1, k-2
                kp1 = k+1
                v_m = vclim_m_t[j, k]
                v_p = vclim_p_t[j, k]
                ws.dX_adv[j, k] += ccy * (
                    -v_m * (wz[j, km1] * (T1[j, k] - T1[j, km1]) +
                            wz[j, km2] * (T1[j, k] - T1[j, km2])) / 3.0 +
                    v_p * wz[j, kp1] * (T1[j, k] - T1[j, kp1])
                )
            end
        else               # k == ydim (South Pole)
            @turbo for j in 1:xdim
                km1, km2 = k-1, k-2
                v_m = vclim_m_t[j, k]
                ws.dX_adv[j, k] += ccy * (
                    -v_m * (wz[j, km1] * (T1[j, k] - T1[j, km1]) +
                            wz[j, km2] * (T1[j, k] - T1[j, km2]))
                ) / 3.0
            end
        end

        # ----- Zonal (u) advection -----
        if !is_polar[k]   # mid‑latitudes, normal timestep
            @turbo for j in 1:xdim
                jm1, jp1 = lon_jm1[j], lon_jp1[j]
                jm2, jp2 = lon_jm2[j], lon_jp2[j]
                u_m = uclim_m_t[j, k]
                u_p = uclim_p_t[j, k]
                ws.dX_adv[j, k] += ccx[k] * (
                    -u_m * (wz[jm1, k] * (T1[j, k] - T1[jm1, k]) +
                            wz[jm2, k] * (T1[j, k] - T1[jm2, k])) +
                    u_p * (wz[jp1, k] * (T1[j, k] - T1[jp1, k]) +
                           wz[jp2, k] * (T1[j, k] - T1[jp2, k]))
                ) / 3.0
            end
        else # polar regions – sub‑timestepping
            # Number of sub‑steps (CFL stability. Precomputed, depends only on k)
            time2 = POLAR_ADV_TIME2[k]
            ccx2 = POLAR_ADV_CCX2[k]

            # Copy current row into temporary buffer
            ws.T1h .= @view T1[:, k]

            for _ in 1:time2
                # Jacobi
                @turbo for j in 1:xdim
                    jm1, jp1 = lon_jm1[j], lon_jp1[j]
                    jm2, jp2 = lon_jm2[j], lon_jp2[j]
                    jm3, jp3 = lon_jm3[j], lon_jp3[j]
                    u_m = uclim_m_t[j, k]
                    u_p = uclim_p_t[j, k]

                    ws.dTxh[j] = ccx2 * (
                        -u_m * (10.0 * wz[jm1, k] * (ws.T1h[j] - ws.T1h[jm1]) +
                                4.0 * wz[jm2, k] * (ws.T1h[jm1] - ws.T1h[jm2]) +
                                1.0 * wz[jm3, k] * (ws.T1h[jm2] - ws.T1h[jm3])) +
                        u_p * (10.0 * wz[jp1, k] * (ws.T1h[j] - ws.T1h[jp1]) +
                               4.0 * wz[jp2, k] * (ws.T1h[jp1] - ws.T1h[jp2]) +
                               1.0 * wz[jp3, k] * (ws.T1h[jp2] - ws.T1h[jp3]))
                    ) / 20.0
                end
                @turbo for j in 1:xdim
                    # Stability clamp (avoid negative water vapour)
                    dq = ifelse(ws.dTxh[j] <= -ws.T1h[j], -0.9 * ws.T1h[j], ws.dTxh[j])
                    ws.T1h[j] += dq
                end
            end

            # Add total change to the output buffer
            @views @. ws.dX_adv[:, k] += ws.T1h - T1[:, k]
        end
    end

    return nothing
end

"""
    circulation!(X_in, h_scl, dX_out, fields::ClimateFields, ws::CirculationWorkspace, timestate, cfg::PhysicsConfig)

Sub-steps `X_in` through `ntime` iterations of [`diffusion!`](@ref),
[`advection!`](@ref), and [`convergence!`](@ref) (each gated by the relevant
`cfg.log_*` switch), writing the total change into `dX_out`. The sub-step
loop is a genuine sequential recurrence and is not parallelized.
"""
function circulation!(X_in, h_scl, dX_out, fields::ClimateFields, ws::CirculationWorkspace, timestate, cfg::PhysicsConfig)
    # Early exit if atmospheric processes disabled
    if (!cfg.log_atmos_dmc || !cfg.log_crcl_dmc || !cfg.log_crcl_drsp)
        fill!(dX_out, 0.0)
        return nothing
    end

    # Precompute flags
    do_diff_v = cfg.log_vdif && h_scl == z_vapor
    do_diff_h = cfg.log_hdif && h_scl == z_air
    do_adv_v = cfg.log_vadv && h_scl == z_vapor
    do_adv_h = cfg.log_hadv && h_scl == z_air
    do_conv = cfg.log_conv && h_scl == z_vapor

    copyto!(ws.X_work, X_in)

    fill!(ws.dX_diff, 0.0)
    fill!(ws.dX_adv, 0.0)
    fill!(ws.dX_conv, 0.0)

    for _tt in 1:ntime
        do_diff_v && diffusion!(ws.X_work, h_scl, fields, ws, timestate)
        do_diff_h && diffusion!(ws.X_work, h_scl, fields, ws, timestate)
        do_adv_v && advection!(ws.X_work, h_scl, fields, ws, timestate, cfg)
        do_adv_h && advection!(ws.X_work, h_scl, fields, ws, timestate, cfg)
        do_conv && convergence!(ws.X_work, fields, timestate, ws)

        @turbo for j in 1:ydim
            for i in 1:xdim
                ws.X_work[i, j] += ws.dX_diff[i, j] + ws.dX_adv[i, j] + ws.dX_conv[i, j]
            end
        end
    end

    # Final difference
    @. dX_out = ws.X_work - X_in

    return nothing
end
