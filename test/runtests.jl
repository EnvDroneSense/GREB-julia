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
        saved_omega = copy(GREB.omegaclim)
        try
            GREB.omegaclim .= 0.01
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
            circulation!(q_in, GREB.z_vapor, dq_on, CirculationWorkspace(), ts, cfg_on)

            cfg_off = create_experiment_config(:full_model)
            cfg_off.log_vdif = false
            cfg_off.log_vadv = false
            cfg_off.log_conv = false
            circulation!(q_in, GREB.z_vapor, dq_off, CirculationWorkspace(), ts, cfg_off)

            @test !all(iszero, dq_on)
            @test all(iszero, dq_off)
            @test dq_on != dq_off
        finally
            GREB.omegaclim .= saved_omega
        end
    end

    @testset "co2_part regional CO2 mask resets between runs (no leak)" begin
        # Regression: co2_part is a module global only ever mutated by
        # forcing()'s regional_co2_* branches; nothing reset it at the start
        # of a run, so a regional-CO2 experiment's masked values leaked into
        # any later, unrelated run in the same session.
        cfg_regional = PhysicsConfig(experiment=:regional_co2_nh)
        redirect_stdout(devnull) do
            greb_model!(0, 1, 1, cfg_regional; jld2_dir = "")
        end
        @test any(!=(1.0), GREB.co2_part)  # regional run actually changed the mask

        cfg_plain = create_experiment_config(:full_model)
        redirect_stdout(devnull) do
            greb_model!(0, 1, 0, cfg_plain; jld2_dir = "")
        end
        @test all(==(1.0), GREB.co2_part)  # init_model! resets it back to full CO2
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
        for log_eva in (-1, 0, 1), log_rain in (-1, 0, 1, 2, 3)
            cfg = create_experiment_config(:full_model)
            cfg.log_eva = log_eva
            cfg.log_rain = log_rain
            result = redirect_stdout(devnull) do
                greb_model!(0, 1, 0, cfg; jld2_dir = "")
            end
            @test length(result.ctrl) == 12
        end
    end

    @testset "read_jld2 rejects non-JLD2 input" begin
        tmp = tempname() * ".jld2"
        write(tmp, "not a jld2 file")
        @test_throws Exception read_jld2(tmp)
        rm(tmp; force = true)
    end

end
