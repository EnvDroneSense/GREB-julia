"""
    seaice!(Ts0, fields::ClimateFields, timestate, cfg::PhysicsConfig)

Updates `fields.cap_surf` (surface heat capacity) for ocean points based on
`Ts0`-derived ice fraction, blending land/open-ocean/ice capacities. No-op if
`cfg.log_ocean_dmc` is false; skips the ice-albedo blend if `cfg.log_ice` is
false.
"""
# ── notebook cell 9bff59c5-4631-4091-8230-989a835788e5  (orig lines 1650-1688) ──
function seaice!(Ts0, fields::ClimateFields, timestate, cfg::PhysicsConfig)
    mld = @view fields.mldclim[:, :, timestate.ityr]
    z_topo = fields.z_topo
    glacier = fields.glacier
    cap_surf = fields.cap_surf

    if !cfg.log_ocean_dmc
        return     # No ice feedback: skip sea ice calculation
    end

    # Compute ice‑dependent heat capacity for ocean points
    @turbo for i in 1:xdim, j in 1:ydim
        is_ocean = z_topo[i, j] < 0.0
        T = Ts0[i, j]
        mld_val = mld[i, j]
        cap_open = cap_ocean * mld_val

        # Ice fraction (0 = no ice, 1 = full ice)
        ice_frac = ifelse(T <= To_ice1, 1.0,
            ifelse(T >= To_ice2, 0.0,
                1.0 - (T - To_ice1) * inv_To_ice_range))

        # Blend between land (ice) and open ocean capacities
        cap_with_ice = cap_land * ice_frac + cap_open * (1.0 - ice_frac)

        # Apply only to ocean points; keep land points unchanged
        cap_surf[i, j] = ifelse(is_ocean, cap_with_ice, cap_surf[i, j])
    end

    # Override for experiments without ice‑albedo feedback
    if !cfg.log_ice
        @turbo for i in 1:xdim, j in 1:ydim
            cap_surf[i, j] = ifelse(z_topo[i, j] > 0.0, cap_land, cap_ocean * mld[i, j])
        end
        return
    end

    # Glacier override: ice sheets have land heat capacity
    @. cap_surf = ifelse(glacier > 0.5, cap_land, cap_surf)
end

"""
    deep_ocean!(Ts, To, fields::ClimateFields, timestate, cfg::PhysicsConfig, ws::CirculationWorkspace)

Computes surface/deep-ocean coupling tendencies (`dT_ocean`, `dTo`) from
mixed-layer-depth entrainment/detrainment and turbulent mixing, active only
where the point is ocean and above the sea-ice threshold. Returns zeros if
`cfg.log_ocean_dmc`/`cfg.log_ocean_drsp` disable ocean coupling.
"""
# ── notebook cell 625089e2-ef77-4821-a6d6-d0a0f88207f2  (orig lines 1703-1750) ──
function deep_ocean!(Ts, To, fields::ClimateFields, timestate, cfg::PhysicsConfig, ws::CirculationWorkspace)
    # Use pre-allocated zero buffers
    dT_ocean = ws.dT_ocean_buf
    dTo = ws.dTo_buf

    # no deep-ocean coupling
    if !cfg.log_ocean_dmc || !cfg.log_ocean_drsp
        fill!(dT_ocean, 0.0)
        fill!(dTo, 0.0)
        return (dT_ocean=dT_ocean, dTo=dTo)
    end

    z_topo = fields.z_topo
    z_ocean = fields.z_ocean

    # ── Change in mixed-layer depth ─────────────────────────
    mld_now = @view fields.mldclim[:, :, timestate.ityr]
    mld_prev = timestate.ityr > 1 ? @view(fields.mldclim[:, :, timestate.ityr-1]) : @view(fields.mldclim[:, :, nstep_yr])

    # Zero buffers first (one pass, cheap)
    fill!(dT_ocean, 0.0)
    fill!(dTo, 0.0)

    # ── Entrainment & detrainment & turbulent mixing ──────
    @turbo for i in 1:xdim, j in 1:ydim
        active = (z_topo[i, j] < 0.0) & (Ts[i, j] >= To_ice2)
        h_now = mld_now[i, j]
        h_prev = mld_prev[i, j]
        dh = h_now - h_prev
        z_deep = z_ocean[i, j]
        z_rem = z_deep - h_now

        # Entrainment/detrainment contributions (only when active)
        dTo_entr = ifelse(active & (dh < 0.0), c_effmix * (-dh / z_rem) *
                                               (Ts[i, j] - To[i, j]), 0.0)
        dT_ocean_entr = ifelse(active & (dh > 0.0), c_effmix * (dh / h_now) *
                                                    (To[i, j] - Ts[i, j]), 0.0)

        # Turbulent mixing (only when active)
        Tx = ifelse(Ts[i, j] > To_ice2, Ts[i, j], To_ice2)
        dTo_turb = ifelse(active, turb_coeff * (Tx - To[i, j]) / z_rem, 0.0)
        dT_ocean_turb = ifelse(active, turb_coeff * (To[i, j] - Tx) / h_now, 0.0)

        # Combine (buffer was zeroed before loop)
        dTo[i, j] = dTo_entr + dTo_turb
        dT_ocean[i, j] = dT_ocean_entr + dT_ocean_turb
    end
    return (dT_ocean=dT_ocean, dTo=dTo)
end
