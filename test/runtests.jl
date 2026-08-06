using GREB
using Test

# Smoke tests that do NOT require the (large, external) JLD2 input data.
# They check that the package loads, its types build, and the grid/constants
# are intact after the notebook -> package extraction. Full integration runs
# (which need `greb_dataset_jld2/`) are demonstrated in examples/run_greb.jl.

const DATA_DIR = joinpath(@__DIR__, "..", "greb_dataset_jld2")

@testset "GREB.jl" begin

    @testset "grid constants" begin
        @test GREB.xdim == 96
        @test GREB.ydim == 48
        @test GREB.nstep_yr == 730
    end

    @testset "min_T_K is -40°C, a real cold-extreme floor" begin
        # Fortran's Tmin_limit=40 is raw Kelvin (greb.model.mscm.f90:470-477)
        # -- colder than anywhere on Earth ever gets, so it can never
        # physically bind. Kept at 233.15 K (-40°C) intentionally: Fortran
        # isn't always the right reference to match literally.
        @test GREB.min_T_K ≈ 233.15
        Ts = fill(220.0, GREB.xdim, GREB.ydim)
        @test all(≈(233.15), max.(Ts, GREB.min_T_K))
    end

    @testset "PhysicsConfig" begin
        cfg = PhysicsConfig()
        @test cfg isa PhysicsConfig

        for exp in (:full_model, :constant_topo, :co2_double, :co2_quadruple,
                    :elnino, :lanina, :rcp85)
            c = create_experiment_config(exp)
            @test c isa PhysicsConfig
            @test c.experiment == exp
        end

        @test create_experiment_config(:co2_double).co2_concentration == 680.0
    end

    @testset "set_hydrology_parameters! writes into cfg, not a disconnected global" begin
        # Regression: this used to declare `global c_q, c_rq, c_omega, c_omegastd`
        # and assign those instead of cfg.c_q etc. hydro! reads cfg.c_q, so
        # log_rain silently had zero effect on every run regardless of its
        # value. Check cfg's own fields match the HYDRO_PARAMS table directly.
        for log_rain in (-1, 0, 1, 2, 3)
            cfg = create_experiment_config(:full_model)
            cfg.log_rain = log_rain
            set_hydrology_parameters!(cfg)
            expected = GREB.HYDRO_PARAMS[log_rain]
            @test (cfg.c_q, cfg.c_rq, cfg.c_omega, cfg.c_omegastd) == expected
        end

        # NCEP adjustment branch (log_rain == 0 && log_clim == 1) overrides the table
        cfg = create_experiment_config(:full_model)
        cfg.log_rain = 0
        cfg.log_clim = 1
        set_hydrology_parameters!(cfg)
        @test (cfg.c_q, cfg.c_rq, cfg.c_omega, cfg.c_omegastd) == (-1.27, 1.99, -16.54, 21.15)
    end

    @testset "set_hydrology_parameters! errors on invalid log_rain" begin
        # Regression: an out-of-range log_rain (e.g. a typo) used to silently
        # fall back to the "Original GREB" parameters via get(dict, key, default)
        # with no warning. Now it should fail loudly instead.
        cfg = create_experiment_config(:full_model)
        cfg.log_rain = 4
        @test_throws ErrorException set_hydrology_parameters!(cfg)
    end

    @testset "circulation! do_conv gate matches cfg.log_conv (not inverted)" begin
        # Regression: do_conv used to read `cfg.log_conv == 0` while every
        # sibling flag (do_diff_v, do_diff_h, do_adv_v, do_adv_h) reads `== 1`.
        # Since log_conv::Bool defaults to true, moisture-flux convergence was
        # silently OFF for every default-configured run. Give omegaclim a
        # nonzero pattern so convergence! has something to contribute, then
        # confirm log_conv actually gates it (and in the right direction).
        fields = ClimateFields()
        fields.omegaclim .= 0.01
        ts = TimeState(1, 1)
        q_in = fill(0.01, GREB.xdim, GREB.ydim)
        dq_on = similar(q_in)
        dq_off = similar(q_in)

        # Separate workspaces: convergence! only overwrites ws.dX_conv
        # when do_conv is true, so a shared workspace would leak the
        # "on" call's nonzero dX_conv into the "off" call's accumulation.
        cfg_on = create_experiment_config(:full_model)
        cfg_on.log_vdif = false
        cfg_on.log_vadv = false
        cfg_on.log_conv = true
        circulation!(q_in, GREB.z_vapor, dq_on, fields, CirculationWorkspace(), ts, cfg_on)

        cfg_off = create_experiment_config(:full_model)
        cfg_off.log_vdif = false
        cfg_off.log_vadv = false
        cfg_off.log_conv = false
        circulation!(q_in, GREB.z_vapor, dq_off, fields, CirculationWorkspace(), ts, cfg_off)

        @test !all(iszero, dq_on)
        @test all(iszero, dq_off)
        @test dq_on != dq_off
    end

    @testset "circulation! doesn't leak stale ws.dX_conv into a later h_scl==z_air call" begin
        # Regression: circulation! relied on each kernel to zero its own
        # ws.dX_* output, but convergence! only writes ws.dX_conv when
        # do_conv is true — which is structurally always false for
        # h_scl==z_air. A shared workspace reused across a vapor call (where
        # do_conv is true) and a later air call would carry the vapor call's
        # nonzero ws.dX_conv straight into the air call's sub-step
        # accumulation. Fortran zeroes dx_diffuse/dx_advec/dx_conv once per
        # circulation() call (greb.model.mscm.f90:856-858); disable diffusion
        # and advection for the air call so its only possible nonzero
        # contribution would be a leaked dX_conv.
        fields = ClimateFields()
        fields.omegaclim .= 0.01
        ts = TimeState(1, 1)
        ws = CirculationWorkspace()

        cfg_vapor = create_experiment_config(:full_model)
        cfg_vapor.log_conv = true
        q_in = fill(0.01, GREB.xdim, GREB.ydim)
        dq_out = similar(q_in)
        circulation!(q_in, GREB.z_vapor, dq_out, fields, ws, ts, cfg_vapor)
        @test !all(iszero, dq_out)  # confirm convergence! actually ran and left ws.dX_conv nonzero

        cfg_air = create_experiment_config(:full_model)
        cfg_air.log_hdif = false
        cfg_air.log_hadv = false
        Ta_in = fill(280.0, GREB.xdim, GREB.ydim)
        dTa_out = similar(Ta_in)
        circulation!(Ta_in, GREB.z_air, dTa_out, fields, ws, ts, cfg_air)
        @test all(iszero, dTa_out)
    end

    @testset "diffusion!/advection! polar sub-stepping uses Jacobi, not Gauss-Seidel" begin
        # Regression: Fortran computes the whole row's increment (dTxh) from
        # the OLD T1h first, then applies the clamp and adds it to the whole
        # row at once (Jacobi: greb.model.mscm.f90:983-1036 for diffusion,
        # :1227-1228 for advection). The prior Julia code wrote
        # `ws.T1h[j] += dq` inside the same @turbo loop that computed dq, so
        # later j read T1h values already updated by earlier j in the same
        # sweep (Gauss-Seidel) — a different numerical scheme, and unsound
        # under @turbo (which assumes no cross-iteration dependency).
        #
        # This test independently re-implements the polar branch's exact
        # stencil as a plain (non-@turbo, provably order-safe) Jacobi
        # reference and checks diffusion!'s actual output matches it — a
        # Gauss-Seidel reimplementation would diverge from this reference
        # for the nonuniform input below.
        k = 5  # dxlat_grid[5] <= 2.5e5: triggers the polar sub-stepping branch
        fields = ClimateFields()
        for j in 1:GREB.ydim, i in 1:GREB.xdim
            fields.wz_air[i, j] = 0.6 + 0.3 * sin(i / 5.0) * cos(j / 3.0)
        end
        T1 = [280.0 + 8.0 * sin(i / 4.0) + 3.0 * (i % 5) for i in 1:GREB.xdim, j in 1:GREB.ydim]
        ts = TimeState(1, 1)
        ws = CirculationWorkspace()

        GREB.diffusion!(T1, GREB.z_air, fields, ws, ts)
        actual = copy(ws.dX_diff[:, k])

        # Independent Jacobi reference for row k, mirroring the production
        # stencil/coefficients exactly but with plain loops (no @turbo) so
        # order-safety is guaranteed by construction.
        wz = fields.wz_air
        dxlat = GREB.dxlat_grid
        dd = max(1, round(Int, GREB.Δt_crcl / (dxlat[k]^2 / GREB.κ)))
        dtdff2 = GREB.Δt_crcl / dd
        time2 = max(1, round(Int, GREB.Δt_crcl / dtdff2))
        ccx2 = GREB.κ * dtdff2 / dxlat[k]^2
        jm1, jp1 = GREB.lon_jm1, GREB.lon_jp1
        jm2, jp2 = GREB.lon_jm2, GREB.lon_jp2
        jm3, jp3 = GREB.lon_jm3, GREB.lon_jp3

        T1h = copy(T1[:, k])
        dTxh = similar(T1h)
        for _ in 1:time2
            for j in 1:GREB.xdim
                jm1v, jp1v = jm1[j], jp1[j]
                jm2v, jp2v = jm2[j], jp2[j]
                jm3v, jp3v = jm3[j], jp3[j]
                dTxh[j] = ccx2 * 0.05 * (
                    10.0 * (wz[jm1v, k] * (T1h[jm1v] - T1h[j]) + wz[jp1v, k] * (T1h[jp1v] - T1h[j])) +
                    4.0 * (wz[jm2v, k] * (T1h[jm2v] - T1h[jm1v]) + wz[jm1v, k] * (T1h[j] - T1h[jm1v])) +
                    4.0 * (wz[jp1v, k] * (T1h[j] - T1h[jp1v]) + wz[jp2v, k] * (T1h[jp2v] - T1h[jp1v])) +
                    1.0 * (wz[jm3v, k] * (T1h[jm3v] - T1h[jm2v]) + wz[jm2v, k] * (T1h[jm1v] - T1h[jm2v])) +
                    1.0 * (wz[jp2v, k] * (T1h[jp1v] - T1h[jp2v]) + wz[jp3v, k] * (T1h[jp3v] - T1h[jp2v]))
                )
            end
            for j in 1:GREB.xdim
                dq = ifelse(dTxh[j] <= -T1h[j], -0.9 * T1h[j], dTxh[j])
                T1h[j] += dq
            end
        end
        expected = wz[:, k] .* (T1h .- T1[:, k])

        @test actual ≈ expected
    end

    @testset "co2_part regional CO2 mask resets between runs (no leak)" begin
        # Regression: co2_part used to be a module global only ever mutated by
        # forcing()'s regional_co2_* branches, with nothing resetting it at the
        # start of a run, so a regional-CO2 experiment's masked values leaked
        # into any later, unrelated run in the same session. Now co2_part lives
        # on ClimateFields (per IMPROVEMENTS.md §1.1's state-struct refactor);
        # by default greb_model! builds a fresh instance per call so this can't
        # happen at all, but init_model! still resets it defensively for the
        # case where a caller explicitly reuses one `fields` across calls
        # (e.g. to avoid reloading real climatology) — exercise that case here.
        fields = ClimateFields()
        cfg_regional = PhysicsConfig(experiment=:regional_co2_nh)
        redirect_stdout(devnull) do
            # scnr=0: the mask is set once by init_model! before the ctrl
            # loop even starts, so a scenario year adds runtime with no
            # extra signal for this assertion.
            greb_model!(RunSpec(scnr = 0), cfg_regional; jld2_dir = "", fields = fields)
        end
        @test any(!=(1.0), fields.co2_part)  # regional run actually changed the mask

        cfg_plain = create_experiment_config(:full_model)
        redirect_stdout(devnull) do
            greb_model!(RunSpec(scnr = 0), cfg_plain; jld2_dir = "", fields = fields)
        end
        @test all(==(1.0), fields.co2_part)  # init_model! resets it back to full CO2
    end

    @testset "workspaces & accumulators" begin
        ws = CirculationWorkspace()
        @test ws isa CirculationWorkspace

        acc = MonthlyAccumulator()
        @test acc isa MonthlyAccumulator
        @test (GREB.reset!(acc); true)   # reset! runs without error

        ts = TimeState(1, 1)
        @test ts.jday == 1
        @test ts.ityr == 1
    end

    @testset "MonthlyRecord type" begin
        @test MonthlyRecord <: NamedTuple
        @test :Ts in fieldnames(MonthlyRecord)
        @test :precip in fieldnames(MonthlyRecord)
    end

    @testset "build_monthly_climatology/apply_scenario_anomalies" begin
        # Coverage gap (IMPROVEMENTS.md §3): zero direct unit tests before —
        # only exercised indirectly through 2 full-model-run tests. Hand-
        # built records with every field filled to one scalar value make the
        # averaging arithmetic trivial to check by hand.
        mkrec(v) = (Ts=fill(v, GREB.xdim, GREB.ydim), Ta=fill(v, GREB.xdim, GREB.ydim),
            To=fill(v, GREB.xdim, GREB.ydim), q=fill(v, GREB.xdim, GREB.ydim),
            albedo=fill(v, GREB.xdim, GREB.ydim), ice=fill(v, GREB.xdim, GREB.ydim),
            precip=fill(v, GREB.xdim, GREB.ydim), evap=fill(v, GREB.xdim, GREB.ydim),
            qcrcl=fill(v, GREB.xdim, GREB.ydim), sw=fill(v, GREB.xdim, GREB.ydim),
            lw=fill(v, GREB.xdim, GREB.ydim), qlat=fill(v, GREB.xdim, GREB.ydim),
            qsens=fill(v, GREB.xdim, GREB.ydim))

        @test build_monthly_climatology(MonthlyRecord[]) == MonthlyRecord[]

        # Multi-year averaging: 2 years, record idx's value == idx. Month m
        # sees values m (year 1) and m+12 (year 2) -> mean == m+6.
        two_years = MonthlyRecord[mkrec(Float64(idx)) for idx in 1:24]
        clim = build_monthly_climatology(two_years)
        @test length(clim) == 12
        for m in 1:12
            @test all(==(Float64(m + 6)), clim[m].Ts)
        end

        # Non-12-multiple count: only 5 records -> months 6..12 never occur
        # and must fall back to records[1] exactly (not e.g. NaN/zero).
        five = MonthlyRecord[mkrec(Float64(idx)) for idx in 1:5]
        clim5 = build_monthly_climatology(five)
        for m in 1:5
            @test all(==(Float64(m)), clim5[m].Ts)
        end
        for m in 6:12
            @test clim5[m] == five[1]
        end

        # apply_scenario_anomalies: 2 years of scnr records (value == idx),
        # ctrl climatology value == 100+month. Check the wraparound (idx=13
        # is month 1, same as idx=1) subtracts the correct month's climatology.
        ctrl_clim = MonthlyRecord[mkrec(100.0 + m) for m in 1:12]
        scnr = MonthlyRecord[mkrec(Float64(idx)) for idx in 1:24]
        anom = apply_scenario_anomalies(scnr, ctrl_clim)
        @test all(==(1.0 - 101.0), anom[1].Ts)
        @test all(==(13.0 - 101.0), anom[13].Ts)
        @test all(==(12.0 - 112.0), anom[12].Ts)

        # Early-return guards: empty scnr_records or empty ctrl_clim ->
        # scnr_records passed straight through, not turned into anomalies.
        @test apply_scenario_anomalies(MonthlyRecord[], ctrl_clim) == MonthlyRecord[]
        @test apply_scenario_anomalies(scnr, MonthlyRecord[]) == scnr
    end

    @testset "greb_model! runs without notebook globals" begin
        # Regression: qflux_correction!/greb_model! used to reference the Pluto
        # @bind globals `time_flux`/`jld2_dir`. They are now parameters, so the
        # model must run to completion on default (unloaded) fields. Values are
        # NaN without real JLD2 data — we only assert it runs and shapes are OK.
        cfg = create_experiment_config(:full_model)
        result = redirect_stdout(devnull) do
            greb_model!(RunSpec(scnr = 0), cfg; jld2_dir = "")
        end
        @test length(result.ctrl) == 12
        @test length(result.scnr) == 0
        @test result.ctrl[1] isa MonthlyRecord
    end

    @testset "greb_model! runs across log_eva / log_rain branches" begin
        # Regression for a `ws.cE` vs `ws.cE_buf` field-name typo in hydro!'s
        # log_eva == 0 branch: it threw `type CirculationWorkspace has no
        # field cE` at runtime, but was invisible to the test above because
        # PhysicsConfig's default (log_eva = -1) never reached that branch.
        # log_eva and log_rain gate independent branches of hydro! (evaporation
        # vs. rain), so a full cross product (20 full-model-year runs) buys no
        # extra bug-catching power over hitting every value of each axis at
        # least once — zip the two lists (cycling the shorter) instead: 5 runs
        # cover all 4 log_eva and all 5 log_rain values.
        log_evas = (-1, 0, 1, 2)
        log_rains = (-1, 0, 1, 2, 3)
        n = max(length(log_evas), length(log_rains))
        combos = collect(Iterators.take(zip(Iterators.cycle(log_evas), Iterators.cycle(log_rains)), n))
        for (log_eva, log_rain) in combos
            cfg = create_experiment_config(:full_model)
            cfg.log_eva = log_eva
            cfg.log_rain = log_rain
            result = redirect_stdout(devnull) do
                greb_model!(RunSpec(scnr = 0), cfg; jld2_dir = "")
            end
            @test length(result.ctrl) == 12
        end
    end

    @testset "hydro! log_rain==1 rain limit uses the full wz_vapor field, not a single point" begin
        # Regression: the rain limit collapsed the spatially-varying
        # wz_vapor field to a single grid point (wz_vapor[1,1], the
        # Antarctic coast) and applied that one scalar globally. Fortran
        # uses the full 2D field (greb.model.mscm.f90:755). Zero out
        # c_q/c_rq/c_omega/c_omegastd so the pre-limit dq_rain is exactly
        # 0.0 everywhere, guaranteeing the limit binds at every point
        # regardless of hydrology parameterization.
        fields = ClimateFields()
        fields.wz_vapor .= 1.0
        fields.wz_vapor[1, 1] = 1.0
        fields.wz_vapor[50, 25] = 0.1
        cfg = create_experiment_config(:full_model)
        cfg.log_rain = 1
        cfg.c_q = 0.0
        cfg.c_rq = 0.0
        cfg.c_omega = 0.0
        cfg.c_omegastd = 0.0
        ts = TimeState(1, 1)
        Ts = fill(290.0, GREB.xdim, GREB.ydim)
        q = fill(0.02, GREB.xdim, GREB.ydim)
        out = hydro!(Ts, q, fields, ts, cfg, CirculationWorkspace())

        limit_1_1 = -0.0015 / (fields.wz_vapor[1, 1] * GREB.r_qviwv * 86400.0)
        limit_50_25 = -0.0015 / (fields.wz_vapor[50, 25] * GREB.r_qviwv * 86400.0)
        @test out.dq_rain[1, 1] ≈ limit_1_1
        @test out.dq_rain[50, 25] ≈ limit_50_25
        @test out.dq_rain[1, 1] != out.dq_rain[50, 25]
    end

    @testset "hydro! log_eva branches produce distinct output (not just distinct shape)" begin
        # Regression: log_eva modes 1 and 2 used to fall through to an `else`
        # branch that silently duplicated mode -1's formula (byte-identical
        # Qlat/dq_eva output). The branch-coverage sweep above only checks
        # output shape, which is exactly why this went undetected. Give the
        # wind/wetness climatology a nonzero, non-uniform pattern so the
        # different wind-gust/coefficient parameterizations actually diverge,
        # then confirm all four log_eva values produce different dq_eva.
        fields = ClimateFields()
        fields.uclim .= 3.0
        fields.vclim .= 2.0
        fields.swetclim .= 0.5
        fields.wsclim .= 4.0

        Ts = fill(290.0, GREB.xdim, GREB.ydim)
        q = fill(0.005, GREB.xdim, GREB.ydim)
        ts = TimeState(1, 1)

        outputs = map((-1, 0, 1, 2)) do log_eva
            cfg = create_experiment_config(:full_model)
            cfg.log_eva = log_eva
            copy(hydro!(Ts, q, fields, ts, cfg, CirculationWorkspace()).dq_eva)
        end

        for i in 1:length(outputs), j in (i+1):length(outputs)
            @test outputs[i] != outputs[j]
        end
    end

    @testset "hydro! errors on invalid log_eva" begin
        Ts = fill(290.0, GREB.xdim, GREB.ydim)
        q = fill(0.005, GREB.xdim, GREB.ydim)
        cfg = create_experiment_config(:full_model)
        cfg.log_eva = 99
        @test_throws ErrorException hydro!(Ts, q, ClimateFields(), TimeState(1, 1), cfg, CirculationWorkspace())
    end

    @testset "init_model! drsp climatology overrides use the correct switches" begin
        # Regression: init_model! gated the qclim=0.0052 override on
        # log_vapor_dmc (should be log_humid_drsp) and the mldclim=d_ocean
        # override on log_ocean_dmc (should be log_ocean_drsp), and had no
        # log_clouds_drsp -> cldclim=0.7 override at all — all three verified
        # against the Fortran reference's `if(log_cloud_drsp==0) cldclim=0.7`
        # / `log_humid_drsp` / `log_ocean_drsp` block.
        fields = ClimateFields()
        cfg = create_experiment_config(:full_model)
        cfg.log_clouds_drsp = false
        cfg.log_humid_drsp = false
        cfg.log_ocean_drsp = false
        init_model!(cfg, fields)
        @test all(==(0.7), fields.cldclim)
        @test all(==(0.0052), fields.qclim)
        @test all(==(GREB.d_ocean), fields.mldclim)
    end

    @testset "init_model! hoists static regional_co2 masks (nh/sh/tropics/extratropics)" begin
        # Regression/coverage for IMPROVEMENTS.md §4.4's co2_part hoist:
        # nh/sh/tropics/extratropics masks depend only on fixed latitude-row
        # ranges, so they're now set once in init_model! instead of every
        # timestep in forcing(). Confirm the mask is right at init, and that
        # forcing() no longer touches co2_part for these four experiments.
        fields = ClimateFields()
        cfg = create_experiment_config(:full_model)
        cfg.experiment = :regional_co2_nh
        init_model!(cfg, fields)
        @test all(==(0.5), fields.co2_part[:, 1:24])
        @test all(==(1.0), fields.co2_part[:, 25:48])

        fields.co2_part[1, 30] = 0.75  # sentinel outside the nh mask; forcing() must leave it alone
        icmn_ctrl = zeros(Float64, GREB.xdim, GREB.ydim, 1)
        forcing(1, 1970, cfg, fields, icmn_ctrl)
        @test fields.co2_part[1, 30] == 0.75
        @test fields.co2_part[1, 1] == 0.5
    end

    @testset "deep_ocean! turbulent mixing stays active under sea ice (Ts < To_ice2)" begin
        # Regression: Fortran gates entrainment/detrainment on Ts>=To_ice2 but
        # turbulent mixing only on z_topo<0 (ocean), using Tx=max(To_ice2,Ts)
        # specifically so mixing stays well-defined under ice
        # (greb.model.mscm.f90:818-830). Julia previously applied the
        # combined ice-threshold mask to all four terms, silently zeroing
        # turbulent mixing under sea ice/high-latitude winter.
        fields = ClimateFields()
        fields.z_topo .= -1.0             # ocean everywhere
        fields.z_ocean .= 1000.0
        fields.mldclim .= 50.0            # uniform -> dh == 0, isolates turbulent mixing
        cfg = create_experiment_config(:full_model)
        ts = TimeState(1, 1)
        ws = CirculationWorkspace()

        Ts = fill(260.0, GREB.xdim, GREB.ydim)  # below To_ice2: entrainment/detrainment inactive
        To = fill(270.0, GREB.xdim, GREB.ydim)
        out = deep_ocean!(Ts, To, fields, ts, cfg, ws)

        @test !all(iszero, out.dTo)
        @test !all(iszero, out.dT_ocean)
    end

    @testset "seaice! glacier override applies even when log_ice is false" begin
        # Regression: seaice!'s !cfg.log_ice branch returned early, skipping
        # the glacier -> cap_land override that Fortran applies
        # unconditionally afterward (greb.model.mscm.f90:786-792).
        fields = ClimateFields()
        fields.z_topo .= -1.0     # ocean point
        fields.glacier .= 1.0     # glacier mask everywhere
        fields.mldclim .= 50.0
        fields.cap_surf .= 0.0    # distinct sentinel, easy to detect a missed override
        cfg = create_experiment_config(:full_model)
        cfg.log_ice = false
        ts = TimeState(1, 1)
        Ts0 = fill(280.0, GREB.xdim, GREB.ydim)

        seaice!(Ts0, fields, ts, cfg)
        @test all(==(GREB.cap_land), fields.cap_surf)
    end

    @testset "LWradiation! log_atmos_dmc==false only zeros LW_down, not LW_up" begin
        # Regression: LWradiation! zeroed both LW_up and LW_down when
        # log_atmos_dmc was false. Fortran only zeros LWair_down; LWair_up is
        # snapshotted before the conditional zeroing and stays at its full
        # computed value (decouples surface from atmospheric downwelling
        # feedback without touching the atmosphere's own emission term).
        fields = ClimateFields()
        ts = TimeState(1, 1)
        Ts = fill(290.0, GREB.xdim, GREB.ydim)
        Ta = fill(280.0, GREB.xdim, GREB.ydim)
        q = fill(0.005, GREB.xdim, GREB.ydim)
        CO2 = 340.0

        cfg_on = create_experiment_config(:full_model)
        out_on = LWradiation!(Ts, Ta, q, CO2, fields, ts, cfg_on, CirculationWorkspace())

        cfg_off = create_experiment_config(:full_model)
        cfg_off.log_atmos_dmc = false
        out_off = LWradiation!(Ts, Ta, q, CO2, fields, ts, cfg_off, CirculationWorkspace())

        @test all(iszero, out_off.LW_down)
        @test !all(iszero, out_off.LW_up)
        @test out_off.LW_up == out_on.LW_up
    end

    @testset "SWradiation! is allocation-free" begin
        # Regression: `@. sw[:, j] = sf * (1.0 - albedo[:, j])` looked
        # allocation-free but `albedo[:, j]` on the RHS is plain getindex
        # (not a dotview like the LHS), materializing a fresh Vector every
        # iteration — 40704 bytes/call, the entire allocation footprint of
        # tendencies! (which calls SWradiation! once per timestep). Fixed
        # with @views; every other kernel already benchmarks at 0 bytes.
        fields = ClimateFields()
        state = ModelState()
        ts = TimeState(1, 1)
        ws = CirculationWorkspace()
        cfg = create_experiment_config(:full_model)
        Ts = fill(290.0, GREB.xdim, GREB.ydim)
        SWradiation!(Ts, fields, state, ts, cfg, ws)  # warm up (compilation)
        @test @allocated(SWradiation!(Ts, fields, state, ts, cfg, ws)) == 0
    end

    @testset "qflux_correction! pulls Ts/To/q to climatology; Ta gets no correction (matches Fortran)" begin
        # Coverage gap (IMPROVEMENTS.md §3): qflux_correction!'s loop body had
        # zero test coverage — every other test uses time_flux=0. Use a
        # climatology that's uniform in time (same value at every ityr) so
        # the correction's fixed point is reached after a single timestep:
        # algebraically, ws.Ts0_buf ends up as
        # ws.Ts0_buf_uncorrected + (Tc - ws.Ts0_buf_uncorrected) == Tc exactly
        # (same for To/q), so Ts/To/q should land on Tclim/Toclim/qclim, not
        # just move closer to it.
        fields = ClimateFields()
        fields.cap_surf .= GREB.cap_ocean
        for j in 1:GREB.ydim, i in 1:GREB.xdim
            fields.Tclim[i, j, :] .= 280.0 + 5.0 * sin(i / 10.0) * cos(j / 8.0)
            fields.Toclim[i, j, :] .= 279.0
            fields.qclim[i, j, :] .= 0.006
        end
        cfg = create_experiment_config(:full_model)
        state = ModelState()
        ts = TimeState(1, 1)
        ws = CirculationWorkspace()

        Ts = fill(290.0, GREB.xdim, GREB.ydim)
        Ta = fill(290.0, GREB.xdim, GREB.ydim)
        q = fill(0.010, GREB.xdim, GREB.ydim)
        To = fill(285.0, GREB.xdim, GREB.ydim)

        GREB.qflux_correction!(340.0, Ts, Ta, q, To, fields, state, ts, cfg, ws, 1)

        @test any(!=(0.0), fields.TF_correct)
        @test any(!=(0.0), fields.ToF_correct)
        @test any(!=(0.0), fields.qF_correct)
        @test all(isfinite, fields.TF_correct)
        @test all(isfinite, fields.ToF_correct)
        @test all(isfinite, fields.qF_correct)

        @test Ts ≈ fields.Tclim[:, :, 1]
        @test To ≈ fields.Toclim[:, :, 1]
        @test q ≈ fields.qclim[:, :, 1]

        # Ta has no corresponding correction field (fields.TaF_correct
        # doesn't exist) — it just integrates freely, matching Fortran's own
        # qflux_correction exactly (greb.model.mscm.f90:566). Confirm it
        # still moved from its initial value (i.e. isn't just frozen/broken)
        # even though nothing nudges it toward a climatology target.
        @test all(isfinite, Ta)
        @test Ta != fill(290.0, GREB.xdim, GREB.ydim)
    end

    @testset "greb_model! swaps sw_solar for paleo experiments, restores after" begin
        # Regression: paleo/orbital experiments never actually loaded the
        # alternate solar-forcing table despite load_solar_forcing_jld2
        # existing for exactly this purpose (Fortran: `sw_solar =
        # sw_solar_scnr`). sw_solar now lives on ClimateFields; exercise the
        # explicit-fields-reuse case (see the co2_part test above) to confirm
        # a paleo run's swapped table doesn't leak into a later run against
        # the same `fields` instance.
        fields = ClimateFields()
        saved_sw_solar = copy(fields.sw_solar)
        tmpdir = mktempdir()
        try
            mkpath(joinpath(tmpdir, "solar_scenarios"))
            distinctive_value = 999.0
            GREB.jldopen(joinpath(tmpdir, "solar_scenarios", "solar_paleo.jld2"), "w") do file
                file["data"] = fill(distinctive_value, GREB.ydim, GREB.nstep_yr)
                file["dim_names"] = ["lat", "time"]
            end

            cfg = create_experiment_config(:paleo_231kyr)
            captured = mktemp() do path, io
                result = redirect_stdout(io) do
                    greb_model!(RunSpec(), cfg; jld2_dir = tmpdir, fields = fields)
                end
                flush(io)
                (result = result, text = read(path, String))
            end
            @test length(captured.result.scnr) == 12
            # Confirm the swap branch actually ran (fields are all-zero in
            # this test env, so asserting on Ts/anomaly magnitude would be
            # fragile — several other tests in this suite hit NaN for the
            # same reason and only check shape/completion instead).
            @test occursin("loading alternate solar forcing", captured.text)

            # sw_solar restored to its pre-run value after greb_model! returns
            @test fields.sw_solar == saved_sw_solar

            cfg_plain = create_experiment_config(:full_model)
            redirect_stdout(devnull) do
                greb_model!(RunSpec(scnr = 0), cfg_plain; jld2_dir = "", fields = fields)
            end
            @test fields.sw_solar == saved_sw_solar
        finally
            rm(tmpdir; recursive = true, force = true)
        end
    end

    @testset "greb_model! reaches every :experiment symbol forcing()/init_model! dispatch on" begin
        # Coverage gap (IMPROVEMENTS.md §3): create_experiment_config only
        # builds presets for 9 of ~29 :experiment symbols forcing()/
        # init_model! actually handle — the rest are reachable by
        # constructing PhysicsConfig(experiment=:foo) directly, just not via
        # the preset helper. forcing() only runs during the SCENARIO loop
        # (never the control loop — see src/model.jl), so scnr=1 is required
        # here: RunSpec(scnr=0) would make this whole sweep vacuous
        # (forcing() never called, nothing exercised). ctrl=0 is safe and
        # ~halves this test's runtime: compute_annual_ice_climatology on an
        # empty ctrl_output just returns zeros (confirmed directly), which
        # even :regional_co2_ocean/_land_ice's icmn_ctrl read tolerates fine
        # — nothing here asserts on control-run output, only "did it run".
        simple_symbols = (
            :a1b_scenario, :co2_10x, :co2_half, :co2_zero, :solar_cycle_11yr,
            :a1b_enhanced, :co2_sine_wave, :co2_step, :modern_solar_paleo_co2,
            :earth_sun_distance, :regional_co2_nh, :regional_co2_sh,
            :regional_co2_tropics, :regional_co2_extratropics,
            :regional_co2_ocean, :regional_co2_land_ice, :regional_co2_winter,
            :regional_co2_summer, :sst_plus1,
        )
        for sym in simple_symbols
            cfg = PhysicsConfig(experiment = sym)
            result = redirect_stdout(devnull) do
                greb_model!(RunSpec(ctrl = 0, scnr = 1), cfg; jld2_dir = "")
            end
            @test length(result.scnr) == 12
        end

        # Genuine "not yet implemented" placeholders (forcing() error()s on
        # these deliberately) — confirm they still throw, not that they
        # silently no-op.
        for sym in (:rcp26, :rcp45, :rcp60, :custom_co2)
            cfg = PhysicsConfig(experiment = sym)
            @test_throws ErrorException greb_model!(RunSpec(ctrl = 0, scnr = 1), cfg; jld2_dir = "")
        end

        # :paleo_solar_modern_co2/:obliquity/:eccentricity swap in a real
        # solar_scenarios/*.jld2 table at scenario start — synthesize a
        # minimal set (same mktempdir+jldopen pattern as the paleo-swap test
        # above) so these three are reachable without the real dataset.
        tmpdir = mktempdir()
        try
            mkpath(joinpath(tmpdir, "solar_scenarios"))
            GREB.jldopen(joinpath(tmpdir, "solar_scenarios", "solar_paleo.jld2"), "w") do file
                file["data"] = fill(999.0, GREB.ydim, GREB.nstep_yr)
                file["dim_names"] = ["lat", "time"]
            end
            GREB.jldopen(joinpath(tmpdir, "solar_scenarios", "solar_obliquity.jld2"), "w") do file
                file["data"] = fill(999.0, 1, GREB.ydim, GREB.nstep_yr)
                file["dim_names"] = ["index", "lat", "time"]
                file["coords"] = Dict(1 => [0.0])
            end
            GREB.jldopen(joinpath(tmpdir, "solar_scenarios", "solar_eccentricity.jld2"), "w") do file
                file["data"] = fill(999.0, 1, GREB.ydim, GREB.nstep_yr)
                file["dim_names"] = ["index", "lat", "time"]
                file["coords"] = Dict(1 => [0.0])
            end

            for sym in (:paleo_solar_modern_co2, :obliquity, :eccentricity)
                cfg = PhysicsConfig(experiment = sym)
                result = redirect_stdout(devnull) do
                    greb_model!(RunSpec(ctrl = 0, scnr = 1), cfg; jld2_dir = tmpdir)
                end
                @test length(result.scnr) == 12
            end
        finally
            rm(tmpdir; recursive = true, force = true)
        end
    end

    @testset "log_hydro_dmc==false freezes humidity entirely (not just eva/rain)" begin
        # Regression: Fortran zeroes the ENTIRE humidity increment when
        # log_hydro_dmc==0 (greb.model.mscm.f90:486), including dq_crcl and
        # qF_correct -- not just the eva/rain terms (already zero from
        # hydro!'s own early return). circulation!(q, ...) runs regardless
        # of log_hydro_dmc (it has its own, separate log_crcl_dmc/
        # log_atmos_dmc gates), so before this fix a nonzero dq_crcl/
        # qF_correct would still leak into q. With the whole run's
        # log_hydro_dmc off, q must never move from its initial
        # climatological value.
        if !isdir(DATA_DIR)
            @test_skip "greb_dataset_jld2/ not present"
        else
            fields = load_greb_jld2!(DATA_DIR; dataset = :ncep)
            cfg = create_experiment_config(:full_model)
            cfg.log_hydro_dmc = false
            result = redirect_stdout(devnull) do
                greb_model!(RunSpec(scnr = 0), cfg; jld2_dir = DATA_DIR, fields = fields)
            end
            q_ini = fields.qclim[:, :, GREB.nstep_yr]
            for rec in result.ctrl
                @test all(isapprox.(rec.q, q_ini; atol = 1e-9))
            end
        end
    end

    @testset "golden regression: real dataset control+scenario run matches snapshot" begin
        # Tripwire for any future refactor touching the physics kernels: a
        # real 1yr control + 1yr scenario run against the actual NCEP
        # dataset, snapshotted as monthly global-mean Ts/Ta/q. Drift beyond
        # float-reassociation noise (~1e-12, per §1.1's validation) means
        # real behavior changed, not just codegen.
        #
        # Set RUN_GOLDEN=0 to skip this one locally (it's the slowest
        # DATA_DIR-gated test); CI always runs it (default RUN_GOLDEN=1).
        if !isdir(DATA_DIR)
            @test_skip "greb_dataset_jld2/ not present"
        elseif get(ENV, "RUN_GOLDEN", "1") == "0"
            @test_skip "RUN_GOLDEN=0"
        else
            fields = load_greb_jld2!(DATA_DIR; dataset = :ncep)
            cfg = create_experiment_config(:full_model)
            result = redirect_stdout(devnull) do
                greb_model!(RunSpec(), cfg; jld2_dir = DATA_DIR, fields = fields)
            end

            gmean(x) = sum(x) / length(x)
            summarize(rec) = (Ts = gmean(rec.Ts), Ta = gmean(rec.Ta), q = gmean(rec.q))

            # Reference values recomputed after this pass's bug fixes
            # (circulation! dX_conv leak, deep_ocean! turbulent-mixing mask,
            # hydro! rain limit, seaice! glacier override, hydro! rq ordering,
            # diffusion!/advection! Jacobi fix, humidity update semantics) —
            # real behavior changes, not codegen noise. min_T_K was tried and
            # reverted (see IMPROVEMENTS.md §0.14), so these reflect the
            # original -40°C floor.
            ctrl_ref = [
                (Ts = 276.6376432475215, Ta = 279.0132752922387, q = 0.006483368265992144),
                (Ts = 276.1076273785758, Ta = 278.6101854203388, q = 0.0066878269244861795),
                (Ts = 276.3982313917197, Ta = 278.7355819769281, q = 0.006828203121749328),
                (Ts = 278.033418596411, Ta = 280.2620085698279, q = 0.006996728195722164),
                (Ts = 280.1103499917721, Ta = 282.3720691501376, q = 0.00732842859123354),
                (Ts = 281.83909739745076, Ta = 284.25313867300514, q = 0.007855042454125112),
                (Ts = 282.5491066792044, Ta = 285.12633793891365, q = 0.008282889952438756),
                (Ts = 282.33494625480506, Ta = 285.0130523882134, q = 0.008281318600014156),
                (Ts = 281.10617631470427, Ta = 283.7911113302602, q = 0.007863431467610265),
                (Ts = 279.50523603155864, Ta = 282.12404966212097, q = 0.007469852856220665),
                (Ts = 278.537108759095, Ta = 281.1654805465912, q = 0.0073083638644664516),
                (Ts = 278.2347705686773, Ta = 280.9247609592766, q = 0.0073813136191502576),
            ]
            scnr_ref = [
                (Ts = -0.007632155594618766, Ta = -0.007133529465347189, q = -4.9333976027584994e-8),
                (Ts = -0.001309167326571708, Ta = -0.0015145523125641436, q = 1.030987112661906e-8),
                (Ts = -0.00011935227909166186, Ta = -0.00015150845673007126, q = -8.811407575403715e-9),
                (Ts = 0.00011412561624893656, Ta = 0.00010169282102569695, q = 7.191046358240965e-10),
                (Ts = 0.0001557742243889431, Ta = 0.00015178409538934505, q = 1.282729300373303e-8),
                (Ts = 9.186448051926958e-5, Ta = 9.533767452687043e-5, q = 1.5921599083929904e-8),
                (Ts = 4.87049451490498e-5, Ta = 4.9187632382707847e-5, q = 1.0540113727968019e-8),
                (Ts = 3.253222138257147e-5, Ta = 3.153559790906405e-5, q = 3.6892345609650163e-9),
                (Ts = 8.00866341792994e-5, Ta = 7.24281415625589e-5, q = 4.179801343708237e-9),
                (Ts = 9.394539367011155e-5, Ta = 9.5251277080081e-5, q = 9.7252510856804e-10),
                (Ts = 7.547731804897885e-5, Ta = 7.672918591223999e-5, q = -2.129466956195276e-9),
                (Ts = 5.324451949642286e-5, Ta = 5.3642668972774134e-5, q = -3.3368785465046947e-9),
            ]

            @test length(result.ctrl) == length(ctrl_ref)
            @test length(result.scnr) == length(scnr_ref)
            for (rec, ref) in zip(result.ctrl, ctrl_ref)
                s = summarize(rec)
                @test isapprox(s.Ts, ref.Ts; atol = 1e-6)
                @test isapprox(s.Ta, ref.Ta; atol = 1e-6)
                @test isapprox(s.q, ref.q; atol = 1e-6)
            end
            for (rec, ref) in zip(result.scnr, scnr_ref)
                s = summarize(rec)
                @test isapprox(s.Ts, ref.Ts; atol = 1e-6)
                @test isapprox(s.Ta, ref.Ta; atol = 1e-6)
                @test isapprox(s.q, ref.q; atol = 1e-6)
            end
        end
    end

    @testset "read_jld2 rejects non-JLD2 input" begin
        tmp = tempname() * ".jld2"
        write(tmp, "not a jld2 file")
        @test_throws Exception read_jld2(tmp)
        rm(tmp; force = true)
    end

    @testset "load_greb_jld2!/load_flux_corrections_jld2! file-exists branches" begin
        # Coverage gap (IMPROVEMENTS.md §3): both loaders' "file present"
        # branches are only exercised when the real (gitignored)
        # greb_dataset_jld2/ happens to be present locally — never on a
        # clean checkout/CI. Synthesize a minimal dataset matching
        # load_greb_jld2!'s exact expected layout (src/io.jl), same
        # mktempdir+jldopen pattern as the paleo-swap test above.
        write2(path, v) = (mkpath(dirname(path)); GREB.jldopen(path, "w") do f
            f["data"] = fill(v, GREB.xdim, GREB.ydim); f["dim_names"] = ["lon", "lat"]
        end)
        write3(path, v) = (mkpath(dirname(path)); GREB.jldopen(path, "w") do f
            f["data"] = fill(v, GREB.xdim, GREB.ydim, GREB.nstep_yr); f["dim_names"] = ["lon", "lat", "time"]
        end)
        write_solar(path, v) = (mkpath(dirname(path)); GREB.jldopen(path, "w") do f
            f["data"] = fill(v, GREB.ydim, GREB.nstep_yr); f["dim_names"] = ["lat", "time"]
        end)

        tmpdir = mktempdir()
        try
            write2(joinpath(tmpdir, "static", "global.topography.jld2"), 1.0)
            write2(joinpath(tmpdir, "static", "greb.glaciers.jld2"), 2.0)
            write3(joinpath(tmpdir, "climatology", "ncep.tsurf.1948-2007.clim.jld2"), 3.0)
            write3(joinpath(tmpdir, "climatology", "ncep.zonal_wind.850hpa.clim.jld2"), 4.0)
            write3(joinpath(tmpdir, "climatology", "ncep.meridional_wind.850hpa.clim.jld2"), 5.0)
            write3(joinpath(tmpdir, "climatology", "ncep.atmospheric_humidity.clim.jld2"), 6.0)
            write3(joinpath(tmpdir, "climatology", "ncep.soil_moisture.clim.jld2"), 7.0)
            write3(joinpath(tmpdir, "climatology", "isccp.cloud_cover.clim.jld2"), 8.0)
            write3(joinpath(tmpdir, "climatology", "woce.ocean_mixed_layer_depth.clim.jld2"), 9.0)
            write3(joinpath(tmpdir, "climatology", "Tocean.clim.jld2"), 10.0)
            write3(joinpath(tmpdir, "climatology", "erainterim.omega.vertmean.clim.jld2"), 11.0)
            write3(joinpath(tmpdir, "climatology", "erainterim.omega_std.vertmean.clim.jld2"), 12.0)
            write3(joinpath(tmpdir, "climatology", "erainterim.windspeed.850hpa.clim.jld2"), 13.0)
            write_solar(joinpath(tmpdir, "solar", "solar_radiation.clim.jld2"), 14.0)

            # "files missing" branch: no flux-correction files present yet ->
            # load_flux_corrections_jld2! should warn and zero-fill, not error.
            fields_nocorr = load_greb_jld2!(tmpdir; dataset = :ncep)
            @test all(==(0.0), fields_nocorr.TF_correct)
            @test all(==(0.0), fields_nocorr.qF_correct)
            @test all(==(0.0), fields_nocorr.ToF_correct)
            @test all(==(3.0), fields_nocorr.Tclim)  # loader itself still worked

            # "files present" branch: add the 3 flux-correction files and reload.
            write3(joinpath(tmpdir, "climatology", "Tsurf_flux_correction.jld2"), 15.0)
            write3(joinpath(tmpdir, "climatology", "vapour_flux_correction.jld2"), 16.0)
            write3(joinpath(tmpdir, "climatology", "Tocean_flux_correction.jld2"), 17.0)

            fields = load_greb_jld2!(tmpdir; dataset = :ncep)
            @test all(==(1.0), fields.z_topo)
            @test all(==(2.0), fields.glacier)
            @test all(==(3.0), fields.Tclim)
            @test all(==(4.0), fields.uclim)
            @test all(==(14.0), fields.sw_solar)
            @test all(==(15.0), fields.TF_correct)
            @test all(==(16.0), fields.qF_correct)
            @test all(==(17.0), fields.ToF_correct)
        finally
            rm(tmpdir; recursive = true, force = true)
        end

        missing_parent = mktempdir()
        @test_throws ErrorException load_greb_jld2!(joinpath(missing_parent, "nonexistent"))
        rm(missing_parent; recursive = true, force = true)
    end

end
