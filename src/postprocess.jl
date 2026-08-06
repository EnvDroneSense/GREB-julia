"""
    build_monthly_climatology(records::Vector{MonthlyRecord})::Vector{MonthlyRecord}

Averages a multi-year run's monthly `records` into a 12-month climatology,
one mean per calendar month across however many years are present.
"""
function build_monthly_climatology(records::Vector{MonthlyRecord})::Vector{MonthlyRecord}
    isempty(records) && return MonthlyRecord[]

    fields = propertynames(records[1])
    nmonths = 12
    counts = zeros(Int, nmonths)
    clim_acc = [Dict{Symbol,Matrix{Float64}}() for _ in 1:nmonths]

    # Accumulate by mont
    for (idx, rec) in enumerate(records)
        mon = mod(idx - 1, nmonths) + 1
        counts[mon] += 1
        for fld in fields
            arr = getfield(rec, fld)
            if haskey(clim_acc[mon], fld)
                clim_acc[mon][fld] .+= arr
            else
                clim_acc[mon][fld] = copy(arr)
            end
        end
    end

    # Build climatology
    clim = MonthlyRecord[]
    for mon in 1:nmonths
        if counts[mon] == 0
            push!(clim, records[1])
        else
            push!(clim, NamedTuple{fields}(
                Tuple(clim_acc[mon][fld] ./ counts[mon] for fld in fields)
            ))
        end
    end
    return clim
end

"""
    apply_scenario_anomalies(scnr_records, ctrl_clim)::Vector{MonthlyRecord}

Subtracts the matching calendar month of `ctrl_clim` (from
[`build_monthly_climatology`](@ref)) from each record in `scnr_records`,
turning absolute monthly output into anomalies relative to the control run.
"""
function apply_scenario_anomalies(scnr_records::Vector{MonthlyRecord}, ctrl_clim::Vector{MonthlyRecord})::Vector{MonthlyRecord}
    isempty(scnr_records) && return scnr_records
    isempty(ctrl_clim) && return scnr_records

    fields = propertynames(scnr_records[1])
    anom = MonthlyRecord[]

    for (idx, rec) in enumerate(scnr_records)
        mon = mod(idx - 1, 12) + 1
        ref = ctrl_clim[mon]
        push!(anom, NamedTuple{fields}(
            Tuple(getfield(rec, fld) .- getfield(ref, fld) for fld in fields)))
    end
    return anom
end

"""
    compute_annual_ice_climatology(ctrl_output::Vector{MonthlyRecord})

Averages `ctrl_output`'s `ice` field into a 12-month climatology, returned as
an `(xdim, ydim, 12)` array.
"""
function compute_annual_ice_climatology(ctrl_output::Vector{MonthlyRecord})
    ice_months = zeros(Float64, xdim, ydim, 12)
    count = zeros(Int, 12)

    for (idx, rec) in enumerate(ctrl_output)
        mon = mod1(idx, 12)          # month 1..12 based on record index
        @. ice_months[:, :, mon] += rec.ice   # note: field name is `ice`, not `ice_cover`
        count[mon] += 1
    end

    for mon in 1:12
        cnt = count[mon]
        cnt > 0 || continue
        @. ice_months[:, :, mon] *= 1.0 / cnt
    end

    return ice_months
end
