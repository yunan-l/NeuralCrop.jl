function _resolve_catalog_path(base::AbstractString, path::AbstractString)
    expanded = expanduser(path)
    return isabspath(expanded) ? normpath(expanded) : normpath(joinpath(base, expanded))
end

"""Load dataset paths and CFT metadata from a TOML catalog."""
function load_catalog(path::AbstractString)
    raw = TOML.parsefile(path)
    haskey(raw, "datasets") || throw(ArgumentError("catalog must contain a [datasets] section"))
    haskey(raw, "cfts") || throw(ArgumentError("catalog must contain a [cfts] section"))

    base = dirname(abspath(path))
    specs = Dict{Symbol, DatasetSpec}()
    for (name, entry) in raw["datasets"]
        entry isa AbstractDict || throw(ArgumentError("dataset '$name' must be a TOML table"))
        haskey(entry, "path") || throw(ArgumentError("dataset '$name' is missing path"))
        haskey(entry, "variable") || throw(ArgumentError("dataset '$name' is missing variable"))
        specs[Symbol(name)] = DatasetSpec(
            _resolve_catalog_path(base, entry["path"]),
            entry["variable"];
            units = get(entry, "units", ""),
            cft_ids = get(entry, "cft_ids", Int[]),
            rainfed_bands = get(entry, "rainfed_bands", Int[]),
            irrigated_bands = get(entry, "irrigated_bands", Int[]),
        )
    end

    cfts = raw["cfts"]
    haskey(cfts, "ids") || throw(ArgumentError("[cfts] is missing ids"))
    haskey(cfts, "names") || throw(ArgumentError("[cfts] is missing names"))
    registry = CFTRegistry(cfts["ids"], cfts["names"])
    registry_ids = Set(registry.ids)
    for (name, spec) in specs
        all(id -> id in registry_ids, spec.cft_ids) ||
            throw(ArgumentError("dataset '$name' contains CFT ids absent from the global registry"))
        if spec.management_bands !== nothing
            length(spec.management_bands.rainfed) == length(registry.ids) ||
                throw(ArgumentError("dataset '$name' must map every crop CFT"))
        end
    end
    return DatasetCatalog(specs, registry)
end

function dataset(catalog::DatasetCatalog, name::Symbol)
    haskey(catalog.datasets, name) || throw(KeyError(name))
    return catalog.datasets[name]
end

function cft_index(registry::CFTRegistry, id::Integer)
    index = findfirst(==(Int32(id)), registry.ids)
    isnothing(index) && throw(ArgumentError("unknown CFT id $id"))
    return index
end

function cft_index(registry::CFTRegistry, name::AbstractString)
    index = findfirst(==(String(name)), registry.names)
    isnothing(index) && throw(ArgumentError("unknown CFT name '$name'"))
    return index
end

cft_name(registry::CFTRegistry, id::Integer) = registry.names[cft_index(registry, id)]
