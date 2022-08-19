using NewRivaraProductions
using BenchmarkTools
using PkgBenchmark
using Plots

max_threads = Sys.CPU_THREADS

threads = [2^i for i in 0:floor(Int, log2(max_threads))]
benchmarks = Vector{BenchmarkGroup}(undef, length(threads))

bench_results = map(x -> benchmarkpkg(NewRivaraProductions, BenchmarkConfig(id = nothing, juliacmd = `julia -O3`, env = Dict("JULIA_NUM_THREADS" => x))), threads)

# Generating Plots
# Tetrahedra
mean_times_all = map(results -> mean(results.benchmarkgroup["TetrahedralMesh"]["all_refinements"].times), bench_results)
mean_speedup_all = mean_times_all[1] ./ mean_times_all
median_times_all = map(results -> median(results.benchmarkgroup["TetrahedralMesh"]["all_refinements"].times), bench_results)
median_speedup_all = median_times_all[1] ./ median_times_all

p_tetra_all = plot(threads, hcat(mean_speedup_all, median_speedup_all), label=["Mean" "Median"], marker = 3, xlabel="Number of Threads", ylabel="Speedup", title="Refine all tetrahedra")

display(p_tetra_all)
savefig(p_tetra_all, "speedup_tetrahedra_all.svg")

mean_times_last = map(results -> mean(results.benchmarkgroup["TetrahedralMesh"]["last_refinement"].times), bench_results)
mean_speedup_last = mean_times_last[1] ./ mean_times_last
median_times_last = map(results -> median(results.benchmarkgroup["TetrahedralMesh"]["last_refinement"].times), bench_results)
median_speedup_last = median_times_last[1] ./ median_times_last

p_tetra_last = plot(threads, hcat(mean_speedup_last, median_speedup_last), label=["Mean" "Median"], marker = 3, xlabel="Number of Threads", ylabel="Speedup", title="Refine last tetrahedra")

display(p_tetra_last)
savefig(p_tetra_last, "speedup_tetrahedra_last.svg")

# triangles
mean_times_all = map(results -> mean(results.benchmarkgroup["TriangularMesh"]["all_refinements"].times), bench_results)
mean_speedup_all = mean_times_all[1] ./ mean_times_all
median_times_all = map(results -> median(results.benchmarkgroup["TriangularMesh"]["all_refinements"].times), bench_results)
median_speedup_all = median_times_all[1] ./ median_times_all

p_tri_all = plot(threads, hcat(mean_speedup_all, median_speedup_all), label=["Mean" "Median"], marker = 3, xlabel="Number of Threads", ylabel="Speedup", title="Refine all triangles")

display(p_tri_all)
savefig(p_tri_all, "speedup_triangles_all.svg")

mean_times_last = map(results -> mean(results.benchmarkgroup["TriangularMesh"]["last_refinement"].times), bench_results)
mean_speedup_last = mean_times_last[1] ./ mean_times_last
median_times_last = map(results -> median(results.benchmarkgroup["TriangularMesh"]["last_refinement"].times), bench_results)
median_speedup_last = median_times_last[1] ./ median_times_last

p_tri_last = plot(threads, hcat(mean_speedup_last, median_speedup_last), label=["Mean" "Median"], marker = 3, xlabel="Number of Threads", ylabel="Speedup", title="Refine last triangles")

display(p_tri_last)
savefig(p_tri_last, "speedup_triangles_last.svg")


printstyled("\nResults of all the steps for refining a tetrahedral mesh"; bold=true, color=:blue)
for (nthreads, results) in zip(threads, bench_results)
    printstyled("\n\nPara $nthreads threads\n"; color=:cyan)
    display(results.benchmarkgroup["TetrahedralMesh"]["all_refinements"])
end

printstyled("\nResults of the last step of refinenement a tetrahedral mesh"; bold=true, color=:blue)
for (nthreads, results) in zip(threads, bench_results)
    printstyled("\n\nPara $nthreads threads\n"; color=:cyan)
    display(results.benchmarkgroup["TetrahedralMesh"]["last_refinement"])
end

printstyled("\nResults of all the steps for refining a triangular mesh"; bold=true, color=:blue)
for (nthreads, results) in zip(threads, bench_results)
    printstyled("\n\nPara $nthreads threads\n"; color=:cyan)
    display(results.benchmarkgroup["TriangularMesh"]["all_refinements"])
end

printstyled("\nResults of the last step of refinenement a triangular mesh"; bold=true, color=:blue)
for (nthreads, results) in zip(threads, bench_results)
    printstyled("\n\nPara $nthreads threads\n"; color=:cyan)
    display(results.benchmarkgroup["TriangularMesh"]["last_refinement"])
end

printstyled("\nEnd of run_benchmarks\n"; bold=true, color=:green)
