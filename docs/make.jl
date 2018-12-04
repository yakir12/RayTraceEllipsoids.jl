using Documenter, RayTraceEllipsoids

makedocs(
    modules = [RayTraceEllipsoids],
    format = :html,
    sitename = "RayTraceEllipsoids.jl",
    pages = Any["index.md"]
)

deploydocs(
    repo = "github.com/yakir12/RayTraceEllipsoids.jl.git",
    target = "build",
    julia = "1.0",
    deps = nothing,
    make = nothing,
)
