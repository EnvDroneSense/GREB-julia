# Hoisted from hydro!'s body (IMPROVEMENTS.md §4.5): these depend only on
# module-level constants, not on any per-call argument, so recomputing them
# as locals every timestep was pure waste.
const _HYDRO_CONST_FACTOR1 = 3.75e-3
const _HYDRO_CONST_FACTOR2 = 17.08085
const _HYDRO_CONST_FACTOR3 = 234.175
const _HYDRO_GUST_LAND = 4.0
const _HYDRO_GUST_OCEAN = 9.0
const _HYDRO_CE_LAND = 0.25 * ce
const _HYDRO_CE_OCEAN = 0.58 * ce
const _HYDRO_CONST_LATENT = cq_latent * ρ_air * ce

"""
    hydro!(Ts, q, fields::ClimateFields, timestate, cfg::PhysicsConfig, ws::CirculationWorkspace)

Computes latent heat flux and evaporation/rain tendencies. `cfg.log_eva`
(-1/0/1/2) selects the wind-gust/coefficient parameterization used for
evaporation; `cfg.log_rain` (via `cfg.c_q`/`c_rq`/`c_omega`/`c_omegastd`, set
by [`set_hydrology_parameters!`](@ref)) controls the rain regression.
Returns `(Q_lat, Q_lat_air, dq_eva, dq_rain)`.
"""
function hydro!(Ts, q, fields::ClimateFields, timestate, cfg::PhysicsConfig, ws::CirculationWorkspace)
    c_q = cfg.c_q
    c_rq = cfg.c_rq
    c_omega = cfg.c_omega
    c_omegastd = cfg.c_omegastd

    fill!(ws.Q_lat_buf, 0.0)
    fill!(ws.Q_lat_air_buf, 0.0)
    fill!(ws.dq_eva_buf, 0.0)
    fill!(ws.dq_rain_buf, 0.0)

    if !cfg.log_atmos_dmc || !cfg.log_hydro_dmc || !cfg.log_hydro_drsp
        return (Q_lat=ws.Q_lat_buf, Q_lat_air=ws.Q_lat_air_buf,
            dq_eva=ws.dq_eva_buf, dq_rain=ws.dq_rain_buf)
    end

    z_topo = fields.z_topo
    wz_air = fields.wz_air
    wz_vapor = fields.wz_vapor
    u = @view fields.uclim[:, :, timestate.ityr]
    v = @view fields.vclim[:, :, timestate.ityr]
    swet = @view fields.swetclim[:, :, timestate.ityr]
    omega = @view fields.omegaclim[:, :, timestate.ityr]
    omegastd = @view fields.omegastdclim[:, :, timestate.ityr]
    rain_limit = fields.rain_limit  
    apply_rain_limit = cfg.log_rain == 1

    const_factor1 = _HYDRO_CONST_FACTOR1
    const_factor2 = _HYDRO_CONST_FACTOR2
    const_factor3 = _HYDRO_CONST_FACTOR3
    gust_land = _HYDRO_GUST_LAND
    gust_ocean = _HYDRO_GUST_OCEAN
    cE_land = _HYDRO_CE_LAND
    cE_ocean = _HYDRO_CE_OCEAN
    const_latent = _HYDRO_CONST_LATENT

    Q_lat = ws.Q_lat_buf
    Q_lat_air = ws.Q_lat_air_buf
    dq_eva = ws.dq_eva_buf
    dq_rain = ws.dq_rain_buf

    if cfg.log_eva == -1
        @turbo for j in 1:ydim
            for i in 1:xdim
                # Saturation humidity
                T = Ts[i, j] - 273.15
                qs = max(const_factor1 * exp(const_factor2 * T / (T + const_factor3)) * wz_air[i, j], 1e-8)
                ws.qs[i, j] = qs
                rq = q[i, j] / qs
                # Relative humidity
                ws.rq[i, j] = rq
                
                # Evaporation and rain tendencies
                u_val = u[i, j]; v_val = v[i, j]
                wind = sqrt(u_val*u_val + v_val*v_val)
                wind = sqrt(wind*wind + ifelse(z_topo[i, j] > 0.0, gust_land, gust_ocean))
                qlat = (q[i, j] - qs) * wind * const_latent * swet[i, j]
                Q_lat[i, j] = qlat

                drain = (c_q + c_rq * rq + c_omega * omega[i, j] + c_omegastd * omegastd[i, j]) * cq_rain * q[i, j]
                limit_val = rain_limit[i, j]
                drain = ifelse(apply_rain_limit & (drain >= limit_val), limit_val, drain)
                dq_rain[i, j] = drain

                dq_eva[i, j] = -qlat / cq_latent / r_qviwv
                Q_lat_air[i, j] = -drain * cq_latent * r_qviwv
            end
        end
    elseif cfg.log_eva == 0
        ws_view = @view fields.wsclim[:, :, timestate.ityr]
        @turbo for j in 1:ydim
            for i in 1:xdim
                T0 = Ts[i, j] - 273.15
                qs0 = max(const_factor1 * exp(const_factor2 * T0 / (T0 + const_factor3)) * wz_air[i, j], 1e-8)
                ws.qs[i, j] = qs0
                rq = q[i, j] / qs0
                ws.rq[i, j] = rq

                Tskin = ifelse(z_topo[i, j] > 0.0, Ts[i, j] + 5.0, Ts[i, j] + 1.0)
                Tskin = ifelse(Tskin < 200.0, 200.0, Tskin)
                ws.Tskin[i, j] = Tskin
                T = Tskin - 273.15
                qs_val = const_factor1 * exp(const_factor2 * T / (T + const_factor3)) * wz_air[i, j]

                ws_base = ws_view[i, j]
                ws.ws_base[i, j] = ws_base
                gust = ifelse(z_topo[i, j] > 0.0, 132.25, 29.16)
                wind = sqrt(ws_base*ws_base + gust)

                cE = ifelse(z_topo[i, j] > 0.0, cE_land, cE_ocean)
                ws.cE_buf[i, j] = cE
                qlat = cE * wind * ρ_air * cq_latent * (q[i, j] - qs_val) * swet[i, j]
                Q_lat[i, j] = qlat

                drain = (c_q + c_rq * rq + c_omega * omega[i, j] + c_omegastd * omegastd[i, j]) * cq_rain * q[i, j]
                limit_val = rain_limit[i, j]
                drain = ifelse(apply_rain_limit & (drain >= limit_val), limit_val, drain)
                dq_rain[i, j] = drain

                dq_eva[i, j] = -qlat / cq_latent / r_qviwv
                Q_lat_air[i, j] = -drain * cq_latent * r_qviwv
            end
        end
    elseif cfg.log_eva == 1
        gust_land_1 = gust_land + 144.0
        gust_ocean_1 = gust_ocean + 50.41  # 7.1^2
        @turbo for j in 1:ydim
            for i in 1:xdim
                T = Ts[i, j] - 273.15
                qs = max(const_factor1 * exp(const_factor2 * T / (T + const_factor3)) * wz_air[i, j], 1e-8)
                ws.qs[i, j] = qs
                rq = q[i, j] / qs
                ws.rq[i, j] = rq

                u_val = u[i, j]; v_val = v[i, j]
                wind = sqrt(u_val*u_val + v_val*v_val)
                wind = sqrt(wind*wind + ifelse(z_topo[i, j] > 0.0, gust_land_1, gust_ocean_1))
                coeff = ifelse(z_topo[i, j] > 0.0, 0.04, 0.73)
                qlat = (q[i, j] - qs) * wind * cq_latent * ρ_air * coeff * ce * swet[i, j]
                Q_lat[i, j] = qlat

                drain = (c_q + c_rq * rq + c_omega * omega[i, j] + c_omegastd * omegastd[i, j]) * cq_rain * q[i, j]
                limit_val = rain_limit[i, j]
                drain = ifelse(apply_rain_limit & (drain >= limit_val), limit_val, drain)
                dq_rain[i, j] = drain

                dq_eva[i, j] = -qlat / cq_latent / r_qviwv
                Q_lat_air[i, j] = -drain * cq_latent * r_qviwv
            end
        end
    elseif cfg.log_eva == 2
        ws_view = @view fields.wsclim[:, :, timestate.ityr]
        gust_land_2 = 81.0  # 9.0^2
        gust_ocean_2 = 16.0  # 4.0^2
        @turbo for j in 1:ydim
            for i in 1:xdim
                T = Ts[i, j] - 273.15
                qs = max(const_factor1 * exp(const_factor2 * T / (T + const_factor3)) * wz_air[i, j], 1e-8)
                ws.qs[i, j] = qs
                rq = q[i, j] / qs
                ws.rq[i, j] = rq

                wind = ws_view[i, j]
                wind = sqrt(wind*wind + ifelse(z_topo[i, j] > 0.0, gust_land_2, gust_ocean_2))
                coeff = ifelse(z_topo[i, j] > 0.0, 0.56, 0.79)
                qlat = (q[i, j] - qs) * wind * cq_latent * ρ_air * coeff * ce * swet[i, j]
                Q_lat[i, j] = qlat

                drain = (c_q + c_rq * rq + c_omega * omega[i, j] + c_omegastd * omegastd[i, j]) * cq_rain * q[i, j]
                limit_val = rain_limit[i, j]
                drain = ifelse(apply_rain_limit & (drain >= limit_val), limit_val, drain)
                dq_rain[i, j] = drain

                dq_eva[i, j] = -qlat / cq_latent / r_qviwv
                Q_lat_air[i, j] = -drain * cq_latent * r_qviwv
            end
        end
    else
        error("Unknown log_eva value: $(cfg.log_eva). Valid values: -1, 0, 1, 2")
    end

    return (Q_lat=Q_lat,
        Q_lat_air=Q_lat_air,
        dq_eva=dq_eva,
        dq_rain=dq_rain)
end
