module ParkinsonsDualTaskCoordination

using DatasetManager, LabDataSources, C3D, Biomechanics, Peaks, DSP, Statistics,
    StatsBase

include("generate.jl")

export analyze, getheel, centerd

function centerd(signal, events)
    mi, ma = intervalextrema(signal, events)
    mi_μ, ma_μ = circmeand(mi), circmeand(ma)
    signalcent = signal .- Ref(mi_μ + (ma_μ - mi_μ)/2)

    return signalcent
end

function asymmetry(y,x)
    return y - x
end

const IK_SETUP = joinpath(@__DIR__, "default_Setup_InverseKinematicsTool.xml")

function analyze(trial;
    genpath, logio=IOBuffer(), force_events=false, force_ik=false, verbose=true,
    start=25.0
)
    evsfn = joinpath(genpath, "Subject $(subject(trial))", "events", trial.name*".tsv")
    requiresource!(trial, "events" => V3DEventsSource(evsfn); deps=(Source{C3DFile},),
        genfunc=roerdink2008, force=force_events)
    force_ik && requiresource!(trial, TRCFile; modfun=(x -> filterpoints!(x, df)), lab_orientation)
    motfn = joinpath(genpath, "Subject $(subject(trial))", "ik", trial.name*".mot")
    requiresource!(trial, "ik" => OSimMotion(motfn);
        force=force_ik, iksetupxml=IK_SETUP, logio, compress=true)
    verbose && logio isa IOBuffer && !iszero(logio.size) && println(String(take!(logio)))

    finish = conditions(trial)[:task] == "dual" ? 105. : nothing
    events = readsegment(Segment(trial, "events"; start, finish))
    rfcidx = toindices(events["RFC"], 100)

    seg = Segment(trial, "ik"; start, finish)
    jnts = readsegment(seg; series=["hip_flexion_r", "hip_flexion_l",
        "knee_angle_r", "knee_angle_l", "arm_flex_r", "arm_flex_l", "elbow_flex_r",
        "elbow_flex_l"])

    sr = SegmentResult(seg)
    res = results(sr)

    if conditions(trial)[:ma_side] == "R"
        ma = "r"
        la = "l"
    else
        @assert conditions(trial)[:ma_side] == "L"
        ma = "l"
        la = "r"
    end

    # Proximal joint ROMs and peak flexion
    mins, maxs = intervalextrema(jnts[!, "hip_flexion_r"], rfcidx)
    res["rhip_peak_flex"] = circmeand(maxs)
    res["rhip_rom"] = circmeand(maxs) - circmeand(mins)
    res["rhip_rom_cov"] = circstdd(maxs - mins)/res["rhip_rom"]*100
    mins, maxs = intervalextrema(jnts[!, "hip_flexion_l"], rfcidx)
    res["lhip_peak_flex"] = circmeand(maxs)
    res["lhip_rom"] = circmeand(maxs) - circmeand(mins)
    res["lhip_rom_cov"] = circstdd(maxs - mins)/res["lhip_rom"]*100
    mins, maxs = intervalextrema(jnts[!, "arm_flex_r"], rfcidx)
    res["rshoulder_peak_flex"] = circmeand(maxs)
    res["rshoulder_rom"] = circmeand(maxs) - circmeand(mins)
    res["rshoulder_rom_cov"] = circstdd(maxs - mins)/res["rshoulder_rom"]*100
    mins, maxs = intervalextrema(jnts[!, "arm_flex_l"], rfcidx)
    res["lshoulder_peak_flex"] = circmeand(maxs)
    res["lshoulder_rom"] = circmeand(maxs) - circmeand(mins)
    res["lshoulder_rom_cov"] = circstdd(maxs - mins)/res["lshoulder_rom"]*100

    res["ma_hip_rom"] = res[ma*"hip_rom"]
    res["la_hip_rom"] = res[la*"hip_rom"]
    res["ma_hip_rom_cov"] = res[ma*"hip_rom_cov"]
    res["la_hip_rom_cov"] = res[la*"hip_rom_cov"]
    res["hip_rom_asym"] = asymmetry(res["ma_hip_rom"], res["la_hip_rom"])
    res["ma_shoulder_rom"] = res[ma*"shoulder_rom"]
    res["la_shoulder_rom"] = res[la*"shoulder_rom"]
    res["ma_shoulder_rom_cov"] = res[ma*"shoulder_rom_cov"]
    res["la_shoulder_rom_cov"] = res[la*"shoulder_rom_cov"]
    res["shoulder_rom_asym"] = asymmetry(res["ma_shoulder_rom"], res["la_shoulder_rom"])

    res["ma_hip_peak_flex"] = res[ma*"hip_peak_flex"]
    res["la_hip_peak_flex"] = res[la*"hip_peak_flex"]
    res["hip_peakflex_asym"] = asymmetry(res["ma_hip_peak_flex"], res["la_hip_peak_flex"])
    res["ma_shoulder_peak_flex"] = res[ma*"shoulder_peak_flex"]
    res["la_shoulder_peak_flex"] = res[la*"shoulder_peak_flex"]
    res["shoulder_peakflex_asym"] = asymmetry(res["ma_shoulder_peak_flex"], res["la_shoulder_peak_flex"])

    # Proximal interlimb coordination
    θoffset, meansd = crpensemble(jnts[!, "hip_flexion_"*ma], jnts[!, "hip_flexion_"*la], rfcidx;
                                 centerfun=centerd)
    res["hip_inter_avg_phase"] = circmeand(θoffset)
    res["hip_inter_meansd"] = circmeand(meansd)
    if (res[la*"shoulder_rom"] > 5) & (res[ma*"shoulder_rom"] > 5)
        θoffset, meansd = crpensemble(jnts[!, "arm_flex_"*ma], jnts[!, "arm_flex_"*la], rfcidx;
                                     centerfun=centerd)
        res["shoulder_inter_avg_phase"] = circmeand(θoffset)
        res["shoulder_inter_meansd"] = circmeand(meansd)
    else
        res["shoulder_inter_avg_phase"] = missing
        res["shoulder_inter_meansd"] = missing
    end

    # Proximal ipsilateral coordination
    if res["rshoulder_rom"] > 5
        θoffset, meansd = crpensemble(jnts[!, "arm_flex_r"], jnts[!, "hip_flexion_r"], rfcidx;
                                     centerfun=centerd)
        res["rprox_ipsi_avg_phase"] = circmeand(θoffset)
        res["rprox_ipsi_meansd"] = circmeand(meansd)
    else
        res["rprox_ipsi_avg_phase"] = missing
        res["rprox_ipsi_meansd"] = missing
    end
    if res["lshoulder_rom"] > 5
        θoffset, meansd = crpensemble(jnts[!, "arm_flex_l"], jnts[!, "hip_flexion_l"], rfcidx;
                                     centerfun=centerd)
        res["lprox_ipsi_avg_phase"] = circmeand(θoffset)
        res["lprox_ipsi_meansd"] = circmeand(meansd)
    else
        res["lprox_ipsi_avg_phase"] = missing
        res["lprox_ipsi_meansd"] = missing
    end

    res["ma_ipsi_avg_phase"] = res[ma*"prox_ipsi_avg_phase"]
    res["la_ipsi_avg_phase"] = res[la*"prox_ipsi_avg_phase"]
    res["ma_ipsi_meansd"] = res[ma*"prox_ipsi_meansd"]
    res["la_ipsi_meansd"] = res[la*"prox_ipsi_meansd"]
    res["ipsi_meansd_asym"] = asymmetry(res["ma_ipsi_meansd"], res["la_ipsi_meansd"])

    # Shoulder-hip contralateral coordination
    if res["rshoulder_rom"] > 5
        θoffset, meansd = crpensemble(jnts[!, "arm_flex_r"], jnts[!, "hip_flexion_l"], rfcidx;
                                     centerfun=centerd)
        res["rprox_contra_avg_phase"] = circmeand(θoffset)
        res["rprox_contra_meansd"] = circmeand(meansd)
    else
        res["rprox_contra_avg_phase"] = missing
        res["rprox_contra_meansd"] = missing
    end
    if res["lshoulder_rom"] > 5
        θoffset, meansd = crpensemble(jnts[!, "arm_flex_l"], jnts[!, "hip_flexion_r"], rfcidx;
                                     centerfun=centerd)
        res["lprox_contra_avg_phase"] = circmeand(θoffset)
        res["lprox_contra_meansd"] = circmeand(meansd)
    else
        res["lprox_contra_avg_phase"] = missing
        res["lprox_contra_meansd"] = missing
    end

    res["ma_contra_avg_phase"] = res[ma*"prox_contra_avg_phase"]
    res["la_contra_avg_phase"] = res[la*"prox_contra_avg_phase"]
    res["ma_contra_meansd"] = res[ma*"prox_contra_meansd"]
    res["la_contra_meansd"] = res[la*"prox_contra_meansd"]

    # Lower-limb intralimb coordination
    θoffset, meansd = crpensemble(jnts[!, "hip_flexion_r"], jnts[!, "knee_angle_r"], rfcidx;
                                 centerfun=centerd)
    res["rllimb_intra_avg_phase"] = circmeand(θoffset)
    res["rllimb_intra_meansd"] = circmeand(meansd)
    θoffset, meansd = crpensemble(jnts[!, "hip_flexion_l"], jnts[!, "knee_angle_l"], rfcidx;
                                 centerfun=centerd)
    res["lllimb_intra_avg_phase"] = circmeand(θoffset)
    res["lllimb_intra_meansd"] = circmeand(meansd)

    res["ma_llimb_intra_avg_phase"] = res[ma*"llimb_intra_avg_phase"]
    res["la_llimb_intra_avg_phase"] = res[la*"llimb_intra_avg_phase"]
    res["ma_llimb_intra_meansd"] = res[ma*"llimb_intra_meansd"]
    res["la_llimb_intra_meansd"] = res[la*"llimb_intra_meansd"]
    res["llimb_intra_meansd_asym"] = asymmetry(res["ma_llimb_intra_meansd"],
                                           res["la_llimb_intra_meansd"])

    # Upper-limb intralimb coordination
    if res["rshoulder_rom"] > 5
        θoffset, meansd = crpensemble(jnts[!, "arm_flex_r"], jnts[!, "elbow_flex_r"], rfcidx;
                                     centerfun=centerd)
        res["rulimb_intra_avg_phase"] = circmeand(θoffset)
        res["rulimb_intra_meansd"] = circmeand(meansd)
    else
        res["rulimb_intra_avg_phase"] = missing
        res["rulimb_intra_meansd"] = missing
    end
    if res["lshoulder_rom"] > 5
        θoffset, meansd = crpensemble(jnts[!, "arm_flex_l"], jnts[!, "elbow_flex_l"], rfcidx;
                                     centerfun=centerd)
        res["lulimb_intra_avg_phase"] = circmeand(θoffset)
        res["lulimb_intra_meansd"] = circmeand(meansd)
    else
        res["lulimb_intra_avg_phase"] = missing
        res["lulimb_intra_meansd"] = missing
    end

    res["ma_ulimb_intra_avg_phase"] = res[ma*"ulimb_intra_avg_phase"]
    res["la_ulimb_intra_avg_phase"] = res[la*"ulimb_intra_avg_phase"]
    res["ma_ulimb_intra_meansd"] = res[ma*"ulimb_intra_meansd"]
    res["la_ulimb_intra_meansd"] = res[la*"ulimb_intra_meansd"]
    res["ulimb_intra_meansd_asym"] = asymmetry(res["ma_ulimb_intra_meansd"],
                                           res["la_ulimb_intra_meansd"])

    res["pci"] = phase_coordination_index(;rfs=events["RFC"], rfo=events["RFO"],
        lfs=events["LFC"], lfo=events["LFO"])

    # Spatiotemp
    swtm = swing(events["RFC"], events["RFO"]).*100
    res["rswing_avg"] = mean(swtm)
    res["rswing_cov"] = variation(swtm)*100
    swtm = swing(events["LFC"], events["LFO"]).*100
    res["lswing_avg"] = mean(swtm)
    res["lswing_cov"] = variation(swtm)*100

    lheel, rheel = getheel(trial; start=toindices(start, 100))
    (;lsteps, rsteps) = steplength(;rfs=events["RFC"], rfo=events["RFO"], lfs=events["LFC"],
        lfo=events["LFO"], lftpos=lheel, rftpos=rheel, AP=2, VT=nothing, requiredsteps=40,
        fs=100)
    res["numlsteps"] = length(lsteps)
    res["numrsteps"] = length(rsteps)
    res["rstep_len"] = mean(rsteps)
    res["lstep_len"] = mean(lsteps)
    res["rstep_len_cov"] = variation(rsteps)*100
    res["lstep_len_cov"] = variation(lsteps)*100
    (;lsteps, rsteps) = stepwidth(;rfs=events["RFC"], lfs=events["LFC"], lftpos=lheel,
        rftpos=rheel, ML=1, fs=100)
    res["rstep_width_cov"] = variation(-rsteps)*100
    res["lstep_width_cov"] = variation(lsteps)*100

    (;lsteps, rsteps) = steptimes(;rfs=events["RFC"], lfs=events["LFC"])
    res["rstep_time"] = mean(rsteps)
    res["lstep_time"] = mean(lsteps)
    res["rstep_time_cov"] = variation(rsteps)*100
    res["lstep_time_cov"] = variation(lsteps)*100

    # Gait speed
    if hassource(trial, "dflow")
        dflow = readsegment(Segment(trial, "dflow"; start))
        res["gait_speed"] = round(mean(dflow[1:(60*100), "Right Treadmill Speed (m/s)"]); digits=3)
    end

    return sr
end

end # module
