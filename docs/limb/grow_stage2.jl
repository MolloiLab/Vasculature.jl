# Stage 2 of the leg (limb_pipeline.jl): every arterial or venous patch grown down to
# its terminal radius, each cached on its own as LIMB_CACHE_DIR/patch_<k>.jls (arteries) or vein_patch_<k>.jls (veins).
#
#     julia grow_stage2.jl [arteries|veins] [first_patch] [last_patch]      arteries if neither is given
#
# Needs XCAT_LEG_NRB and LIMB_CACHE_DIR (lab storage). The patches grow in a worker process (all threads) that includes
# the notebook, so every function and cached stage in it is reused (LIMB_NOTEBOOK overrides its path). A patch grows in
# seconds, but one whose vessels cannot attach (patch 105: an end of its feeding segment lies in a bone voxel as the
# growth indexes it) tries every step in vain for about an hour and returns its feeding segment alone. So a patch still
# growing after a minute plus a quarter millisecond per terminal (STAGE2_LIMIT_S overrides) is listed as timed out, its
# worker is stopped and a new one carries on, and a patch grown short of its terminals is listed as failed; neither is
# saved. Patches already cached are skipped, so a stopped run resumes where it left off. Each patch is written to a .part file
# and renamed, so an interrupted run leaves no broken cache; a cache that cannot be written stops the run. Progress goes
# to the console and to LIMB_CACHE_DIR/stage2_progress.log.

using Distributed, Dates, Printf

for name in ("XCAT_LEG_NRB", "LIMB_CACHE_DIR")
    haskey(ENV, name) || error("set $name first: without it the notebook falls back to another machine's paths")
end

const NOTEBOOK = get(ENV, "LIMB_NOTEBOOK", joinpath(@__DIR__, "limb_pipeline.jl"))

time_limit(nterm) = haskey(ENV, "STAGE2_LIMIT_S") ? parse(Float64, ENV["STAGE2_LIMIT_S"]) : 60 + nterm / 4000

"""
Defined on each worker: grow patch `k`, save it as `path` and report to `done`, a channel on the main process. Growth that
cannot attach its vessels ends with only the feeding segment and no error (its time limit may not catch a small patch),
so a patch short of its terminals is not saved.
"""
const ON_WORKER = quote
    """
    Windows runs a process it takes for background work on the efficiency cores alone, about 6× slower for growth on
    this machine (2026-11-03). Opt this process out (its own power-throttling state, `ProcessPowerThrottling`).
    """
    function full_speed()
        Sys.iswindows() || return true
        state = UInt32[1, 1, 0] # PROCESS_POWER_THROTTLING_STATE: version 1, execution speed controlled, not throttled
        process = ccall((:GetCurrentProcess, "kernel32"), Ptr{Cvoid}, ())
        return ccall((:SetProcessInformation, "kernel32"), Cint, (Ptr{Cvoid}, Cint, Ptr{UInt32}, UInt32), process, 4, state, 12) != 0
    end

    function stage2_patch(k, system, path, done)
        stage = :grow
        try
            result = grow_patch(k; system)
            grown = length(terminals(result.tree))
            grown >= 0.99 * result.nterm || error("grew $grown of its $(result.nterm) terminals, so it is not saved")
            stage = :save
            serialize(path * ".part", result)
            mv(path * ".part", path; force = true)
            put!(done, (ok = true, nterm = result.nterm, grown = grown, seconds = result.seconds))
        catch e
            stage === :save && rm(path * ".part"; force = true)
            put!(done, (ok = false, stage = stage, error = sprint(showerror, e)))
        end
        return nothing
    end
end

"The running worker's process id: a stuck worker never notices this process end (Ctrl+C), so it is stopped at exit."
const WORKER_PID = Ref(0)
kill_pid(pid) = run(pipeline(ignorestatus(Sys.iswindows() ? `taskkill /F /PID $pid` : `kill -9 $pid`); stdout = devnull, stderr = devnull))
atexit(() -> WORKER_PID[] > 0 && kill_pid(WORKER_PID[]))

"A worker with the notebook loaded, and its process id."
function start_worker()
    w = only(addprocs(1; exeflags = "--threads=auto"))
    WORKER_PID[] = remotecall_fetch(getpid, w)
    remotecall_fetch(ex -> (Core.eval(Main, ex); nothing), w, ON_WORKER)
    remotecall_fetch(() -> Main.full_speed(), w) || @warn "Windows may keep the worker on its efficiency cores, about 6× slower"
    remotecall_fetch(nb -> (Base.include(Main, nb); nothing), w, NOTEBOOK)
    return w, WORKER_PID[]
end

"Stop a worker that may be stuck in a loop that never yields, so it cannot be asked to exit."
function stop_worker(w, pid)
    try
        rmprocs(w; waitfor = 0)
    catch
    end
    kill_pid(pid)
    WORKER_PID[] = 0
    return nothing
end

"Grow the patches `first_patch`…`last_patch` of `system` (`:arteries` or `:veins`, a key of the notebook's `SYSTEMS`)."
function grow_all(system, first_patch, last_patch)
    @info "Starting a worker and loading the notebook in it: cached stages load in seconds" NOTEBOOK system
    w, pid = start_worker()
    voxels, nterm, cache, radius_um, prefix = remotecall_fetch(w, system) do system
        s = Main.SYSTEMS[system]
        (s.patches.voxels, s.counts, Main.CACHE_DIR, 1000 * s.terminal_mm, s.prefix)
    end
    range = first_patch:min(last_patch, length(voxels))
    log = open(joinpath(cache, "stage2_progress.log"), "a")
    function say(msg)
        line = Dates.format(now(), "yyyy-mm-dd HH:MM:SS") * "  " * msg
        println(line)
        println(log, line)
        flush(log)
        flush(stdout)
    end
    todo = [k for k in range if voxels[k] > 0 && !isfile(joinpath(cache, "$(prefix)_$k.jls"))]
    total = sum(nterm[todo]; init = 0)
    say(@sprintf("stage 2 (%s), patches %d-%d: %d to grow (%d already cached), %.1f million terminals of %.1f µm",
        system, first(range), last(range), length(todo), length(range) - length(todo), total / 1.0e6, radius_um))
    grown = 0
    t0 = time()
    failed = Int[]
    timed_out = Int[]
    for k in todo
        path = joinpath(cache, "$(prefix)_$k.jls")
        limit = time_limit(nterm[k])
        done = RemoteChannel(() -> Channel{Any}(1))   # here, so waiting never needs the worker
        remote_do((k, system, path, done) -> Main.stage2_patch(k, system, path, done), w, k, system, path, done)
        if timedwait(() -> isready(done), limit; pollint = 1.0) === :timed_out
            push!(timed_out, k)
            say(@sprintf("%s %4d: still growing after %.0f s, so it is skipped; restarting the worker", prefix, k, limit))
            stop_worker(w, pid)
            w, pid = start_worker()
            continue
        end
        r = take!(done)
        if !r.ok && r.stage === :save
            say("cannot write $path, so the run stops (free space or change LIMB_CACHE_DIR, then run again): $(r.error)")
            stop_worker(w, pid)
            close(log)
            error("cannot write $path")
        elseif !r.ok
            push!(failed, k)
            say("$prefix $k FAILED: $(r.error)")
            continue
        end
        grown += r.nterm
        left = (time() - t0) / grown * (total - grown) / 3600
        say(@sprintf("%s %4d: %7d terminals in %5.1f s | %5.1f%% done, about %.1f h left",
            prefix, k, r.grown, r.seconds, 100grown / total, left))
    end
    unfinished = [label * join(ks, ", ") for (label, ks) in (("failed: ", failed), ("timed out: ", timed_out)) if !isempty(ks)]
    say(isempty(unfinished) ? "stage 2 ($system) done for patches $(first(range))-$(last(range))" :
        "stage 2 ($system) finished without these patches, $(join(unfinished, "; ")) (run again to retry them)")
    close(log)
    rmprocs(w)
    WORKER_PID[] = 0
    return (; failed, timed_out)
end

const SYSTEM = !isempty(ARGS) && ARGS[1] in ("arteries", "veins") ? Symbol(popfirst!(ARGS)) : :arteries
grow_all(SYSTEM, length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 1, length(ARGS) >= 2 ? parse(Int, ARGS[2]) : typemax(Int))
