using GREB
using Test

# Smoke tests that do NOT require the (large, external) JLD2 input data.
# They check that the package loads, its types build, and the grid/constants
# are intact after the notebook -> package extraction. Full integration runs
# (which need `greb_dataset_jld2/`) are demonstrated in examples/run_greb.jl.

@testset "GREB.jl" begin

    @testset "grid constants" begin
        @test GREB.xdim == 96
        @test GREB.ydim == 48
        @test GREB.nstep_yr == 730
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
            greb_model!(0, 1, 1, cfg_regional; jld2_dir = "", fields = fields)
        end
        @test any(!=(1.0), fields.co2_part)  # regional run actually changed the mask

        cfg_plain = create_experiment_config(:full_model)
        redirect_stdout(devnull) do
            greb_model!(0, 1, 0, cfg_plain; jld2_dir = "", fields = fields)
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

    @testset "greb_model! runs without notebook globals" begin
        # Regression: qflux_correction!/greb_model! used to reference the Pluto
        # @bind globals `time_flux`/`jld2_dir`. They are now parameters, so the
        # model must run to completion on default (unloaded) fields. Values are
        # NaN without real JLD2 data — we only assert it runs and shapes are OK.
        cfg = create_experiment_config(:full_model)
        result = redirect_stdout(devnull) do
            greb_model!(0, 1, 0, cfg; jld2_dir = "")
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
        for log_eva in (-1, 0, 1, 2), log_rain in (-1, 0, 1, 2, 3)
            cfg = create_experiment_config(:full_model)
            cfg.log_eva = log_eva
            cfg.log_rain = log_rain
            result = redirect_stdout(devnull) do
                greb_model!(0, 1, 0, cfg; jld2_dir = "")
            end
            @test length(result.ctrl) == 12
        end
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
                    greb_model!(0, 1, 1, cfg; jld2_dir = tmpdir, fields = fields)
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
                greb_model!(0, 1, 0, cfg_plain; jld2_dir = "", fields = fields)
            end
            @test fields.sw_solar == saved_sw_solar
        finally
            rm(tmpdir; recursive = true, force = true)
        end
    end

    @testset "read_jld2 rejects non-JLD2 input" begin
        tmp = tempname() * ".jld2"
        write(tmp, "not a jld2 file")
        @test_throws Exception read_jld2(tmp)
        rm(tmp; force = true)
    end

end
