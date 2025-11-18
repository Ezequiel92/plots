####################################################################################################
# Data acquisition functions
####################################################################################################

"""
    readGroupCatHeader(path::Union{String,Missing})::GroupCatHeader

Read the header of a group catalog in the HDF5 format.

!!! note

    If each group catalog is made of multiple files, I'll read the header on the first one.

# Arguments

  - `path::Union{String,Missing}`: Path to the group catalog file or folder.

# Returns

  - A [`GroupCatHeader`](@ref).
"""
function readGroupCatHeader(path::Union{String,Missing})::GroupCatHeader

    if ismissing(path)

        (
            !logging[] ||
            @warn("readGroupCatHeader: The group catalog file or folder is missing")
        )

        return GroupCatHeader()

    elseif isfile(path)

        (
            HDF5.ishdf5(path) ||
            throw(ArgumentError("readGroupCatHeader: The file $(path) is not in the \
            HDF5 format, I don't know how to read it"))
        )

        file_path = path

    elseif isdir(path)

        sub_files = glob("$(GC_BASENAME)_*.*.hdf5", path)

        (
            !isempty(sub_files) && all(HDF5.ishdf5, sub_files) ||
            throw(ArgumentError("readGroupCatHeader: The directory $(path) does not contain \
            group catalog sub-files in the HDF5 format"))
        )

        file_path = minimum(sub_files)

    else

        throw(ArgumentError("readGroupCatHeader: $(path) does not exist as a file or folder"))

    end

    header = h5open(file_path, "r") do gc_file

        h = gc_file["Header"]

        # Read which attributes are present
        attrs_present = keys(HDF5.attrs(h))

        missing_attrs = setdiff(
            [
                "BoxSize",
                "HubbleParam",
                "Ngroups_ThisFile",
                "Ngroups_Total",
                "Nsubgroups_ThisFile",
                "Nsubgroups_Total",
                "NumFiles",
                "Omega0",
                "OmegaLambda",
                "Redshift",
                "Time",
            ],
            attrs_present,
        )

        (
            isempty(missing_attrs) ||
            throw(ArgumentError("readGroupCatHeader: The attributes $(missing_attrs) are missing \
            from the header"))
        )

        GroupCatHeader(
            box_size          = read_attribute(h, "BoxSize"),
            h0                = read_attribute(h, "HubbleParam"),
            n_groups_part     = read_attribute(h, "Ngroups_ThisFile"),
            n_groups_total    = read_attribute(h, "Ngroups_Total"),
            n_subgroups_part  = read_attribute(h, "Nsubgroups_ThisFile"),
            n_subgroups_total = read_attribute(h, "Nsubgroups_Total"),
            num_files         = read_attribute(h, "NumFiles"),
            omega_0           = read_attribute(h, "Omega0"),
            omega_l           = read_attribute(h, "OmegaLambda"),
            redshift          = read_attribute(h, "Redshift"),
            time              = read_attribute(h, "Time"),
        )

    end

    return header

end

"""
    readSnapHeader(path::String)::SnapshotHeader

Read the header of a snapshot in the HDF5 format.

!!! note

    If each snapshot is made of multiple files, I'll read the header on the first chunk.

# Arguments

  - `path::String`: Path to the snapshot file or folder.

# Returns

  - A [`SnapshotHeader`](@ref) structure.
"""
function readSnapHeader(path::String)::SnapshotHeader

    if isfile(path)

        (
            HDF5.ishdf5(path) ||
            throw(ArgumentError("readSnapHeader: The file $(path) is not in the HDF5 format, \
            I don't know how to read it"))
        )

        file_path = path

        # Count the number of stellar particles, excluding wind particles
        num_part_stars = countStars(file_path)
        num_total_stars = num_part_stars

    elseif isdir(path)

        sub_files = glob("$(SNAP_BASENAME)_*.*.hdf5", path)

        (
            !isempty(sub_files) && all(HDF5.ishdf5, sub_files) ||
            throw(ArgumentError("readSnapHeader: The directory $(path) does not contain \
            snapshot sub-files in the HDF5 format"))
        )

        file_path = minimum(sub_files)

        # Count the number of stellar particles, excluding wind particles
        num_part_stars = countStars(file_path)
        num_total_stars = sum(countStars, sub_files)

    else

        throw(ArgumentError("readSnapHeader: $(path) does not exist as a file or folder"))

    end

    header = h5open(file_path, "r") do snap_file

        h = snap_file["Header"]

        # Read which attributes are present
        attrs_present = keys(HDF5.attrs(h))

        missing_attrs = setdiff(
            [
                "BoxSize",
                "HubbleParam",
                "MassTable",
                "NumFilesPerSnapshot",
                "NumPart_ThisFile",
                "NumPart_Total",
                "Omega0",
                "OmegaLambda",
                "Redshift",
                "Time",
            ],
            attrs_present,
        )

        (
            isempty(missing_attrs) ||
            throw(ArgumentError("readSnapHeader: The attributes $(missing_attrs) are missing from \
            the header"))
        )

        # Only for the stars edit the number in the header, to exclude wind particles
        num_part = read_attribute(h, "NumPart_ThisFile")
        num_part[PARTICLE_INDEX[:stars] + 1] = num_part_stars
        num_total = read_attribute(h, "NumPart_Total")
        num_total[PARTICLE_INDEX[:stars] + 1] = num_total_stars

        # Check if the length units are in the header, otherwise use the IllustrisTNG values
        if "UnitLength_in_cm" ∈ attrs_present
            l_unit = read_attribute(h, "UnitLength_in_cm") * u"cm"
        else
            l_unit = ILLUSTRIS_L_UNIT
        end

        # Check if the mass units are in the header, otherwise use the IllustrisTNG values
        if "UnitMass_in_g" ∈ attrs_present
            m_unit = read_attribute(h, "UnitMass_in_g") * u"g"
        else
            m_unit = ILLUSTRIS_M_UNIT
        end

        # Check if the velocity units are in the header, otherwise use the IllustrisTNG values
        if "UnitVelocity_in_cm_per_s" ∈ attrs_present
            v_unit = read_attribute(h, "UnitVelocity_in_cm_per_s") * u"cm*s^-1"
        else
            v_unit = ILLUSTRIS_V_UNIT
        end

        SnapshotHeader(
            box_size   = read_attribute(h, "BoxSize"),
            h0         = read_attribute(h, "HubbleParam"),
            mass_table = read_attribute(h, "MassTable"),
            num_files  = read_attribute(h, "NumFilesPerSnapshot"),
            num_part   = num_part,
            num_total  = num_total,
            omega_0    = read_attribute(h, "Omega0"),
            omega_l    = read_attribute(h, "OmegaLambda"),
            redshift   = read_attribute(h, "Redshift"),
            time       = read_attribute(h, "Time"),
            l_unit     = l_unit,
            m_unit     = m_unit,
            v_unit     = v_unit,
        )

    end

    return header

end

"""
    isBlockPresent(block::String, group::HDF5.Group)::Bool

Checks if a given data block exist in a HDF5 group.

# Arguments

  - `block::String`: Target block. The possibilities are the keys of [`QUANTITIES`](@ref).
  - `group::HDF5.Group`: HDF5 group.

# Returns

  - If `block` exist in `group`.
"""
function isBlockPresent(block::String, group::HDF5.Group)::Bool

    (
        block ∈ keys(QUANTITIES) ||
        throw(ArgumentError("isBlockPresent: `block` should be a key of `QUANTITIES`, \
        but I got $(block), see the options in `./src/constants/globals.jl`"))
    )

    return QUANTITIES[block].hdf5_name ∈ keys(group)

end

"""
    isBlockPresent(component::Symbol, block::String, path::String)::Bool

Checks if a given block exist in a snapshot.

!!! note

    If each snapshot is made of multiple files, I'll check only in the first chunk.

# Arguments

  - `component::Symbol`: The cell/particle type of the target block. The possibilities are the keys of [`PARTICLE_INDEX`](@ref).
  - `block::String`: Target block. The possibilities are the keys of [`QUANTITIES`](@ref).
  - `path::String`: Path to the snapshot file or folder.

# Returns

  - If `block` exist in the snapshot.
"""
function isBlockPresent(component::Symbol, block::String, path::String)::Bool

    if isfile(path)

        (
            HDF5.ishdf5(path) ||
            throw(ArgumentError("isBlockPresent: The file $(path) is not in the HDF5 format, \
            I don't know how to read it"))
        )

        file_path = path

    elseif isdir(path)

        sub_files = glob("$(SNAP_BASENAME)_*.*.hdf5", path)

        (
            !isempty(sub_files) && all(HDF5.ishdf5, sub_files) ||
            throw(ArgumentError("isBlockPresent: The directory $(path) does not contain \
            snapshot sub-files in the HDF5 format"))
        )

        file_path = minimum(sub_files)

    else

        throw(ArgumentError("isBlockPresent: $(path) does not exist as a file or folder"))

    end

    response = h5open(file_path, "r") do snapshot

        type_str = PARTICLE_CODE_NAME[component]

        if type_str ∈ keys(snapshot)
            isBlockPresent(block, snapshot[type_str])
        else
            false
        end

    end

    return response

end

"""
    readTime(path::String)::Float64

Read the "Time" field in the header of a snapshot file.

!!! note

    If each snapshot is made of multiple files, I'll read the header on the first chunk.

# Arguments

  - `path::String`: Path to the snapshot file or folder.

# Returns

  - The "Time" field in the header (for cosmological simulations is the scale factor).
"""
function readTime(path::String)::Float64

    if isfile(path)

        (
            HDF5.ishdf5(path) ||
            throw(ArgumentError("readTime: The file $(path) is not in the HDF5 format, \
            I don't know how to read it"))
        )

        file_path = path

    elseif isdir(path)

        sub_files = glob("$(SNAP_BASENAME)_*.*.hdf5", path)

        (
            !isempty(sub_files) && all(HDF5.ishdf5, sub_files) ||
            throw(ArgumentError("readTime: The directory $(path) does not contain \
            snapshot sub-files in the HDF5 format"))
        )

        file_path = minimum(sub_files)

    else

        throw(ArgumentError("readTime: $(path) does not exist as a file or folder"))

    end

    time = h5open(file_path, "r") do file
        read_attribute(file["Header"], "Time")
    end

    return time

end

"""
    readTemperature(file_path::String)::Vector{<:Unitful.Temperature}

Compute the temperature of the gas cells in a snapshot.

# Arguments

  - `file_path::String`: Path to the snapshot file.

# Returns

  - The temperature of the gas cells.
"""
function readTemperature(file_path::String)::Vector{<:Unitful.Temperature}

    if isfile(file_path)

        (
            HDF5.ishdf5(file_path) ||
            throw(ArgumentError("readTemperature: The file $(file_path) is not in the \
            HDF5 format, I don't know how to read it"))
        )

    else

        throw(ArgumentError("readTemperature: $(file_path) does not exist as a file"))

    end

    # List of blocks needed to compute the temperature
    blocks = ["U   ", "NE  "]

    temp_data = h5open(file_path, "r") do snapshot

        group = snapshot[PARTICLE_CODE_NAME[:gas]]

        # Get the indices of the missing blocks
        idx_missing = map(x -> !isBlockPresent(x, group), blocks)

        (
            !any(idx_missing) ||
            throw(ArgumentError("readTemperature: The blocks $(blocks[idx_missing]) \
            are missing, and I need them to compute the temperature"))
        )

        # Compute the unit factor for each block
        units = internalUnits.(blocks, file_path)

        [read(group, QUANTITIES[block].hdf5_name) .* unit for (unit, block) in zip(units, blocks)]

    end

    return computeTemperature(temp_data...)

end

"""
    readGoupCatBlocks(
        file_path::String,
        snapshot_path::String,
        request::Dict{Symbol,Vector{String}},
    )::Dict{Symbol,Dict{String,VecOrMat{<:Number}}}

Read the specified blocks from a group catalog file.

# Arguments

  - `file_path::String`: Path to the group catalog file.
  - `snapshot_path::String`: Path to the corresponding snapshot file or folder. This is needed for unit conversion.
  - `request::Dict{Symbol,Vector{String}}`: The blocks to be read. It must have the shape `group type` -> [`block`, `block`, `block`].

# Returns

  - A dictionary with the following shape: `group type` -> (`block` -> data of `block`).
"""
function readGoupCatBlocks(
    file_path::String,
    snapshot_path::String,
    request::Dict{Symbol,Vector{String}},
)::Dict{Symbol,Dict{String,VecOrMat{<:Number}}}

    if isfile(file_path)

        (
            HDF5.ishdf5(file_path) ||
            throw(ArgumentError("readGoupCatBlocks: The file $(file_path) is not in the \
            HDF5 format, I don't know how to read it"))
        )

    else

        throw(ArgumentError("readGoupCatBlocks: $(file_path) does not exist as a file"))

    end

    # Allocate memory
    output = Dict{Symbol,Dict{String,VecOrMat{<:Number}}}()

    h5open(file_path, "r") do gc_file

        # Read from the request only the group catalog types
        for component in groupcatTypes(request)

            blocks = copy(request[component])

            type_str = titlecase(string(component))

            # Allocate memory
            qty_data = Dict{String,VecOrMat{<:Number}}()

            if type_str ∈ keys(gc_file)

                # Read the HDF5 group
                hdf5_group = gc_file[type_str]

                if isempty(hdf5_group)

                    (
                        !logging[] ||
                        @warn("readGoupCatBlocks: The group catalog type :$(component) \
                        in $(file_path) is empty")
                    )

                    # Return an empty array for every missing block
                    for block in blocks
                        unit = internalUnits(block, snapshot_path)
                        qty_data[block] = typeof(1.0 * unit)[]
                    end

                else

                    for block in blocks

                        unit = internalUnits(block, snapshot_path)

                        if isBlockPresent(block, hdf5_group)

                            qty_data[block] = read(hdf5_group, QUANTITIES[block].hdf5_name) .* unit

                        else

                            (
                                !logging[] ||
                                @warn("readGoupCatBlocks: The block $(block) for the group \
                                catalog type :$(component) in $(file_path) is missing")
                            )

                            # Return an empty array for every missing block
                            qty_data[block] = typeof(1.0 * unit)[]

                        end

                    end

                end

            else

                (
                    !logging[] ||
                    @warn("readGoupCatBlocks: The group catalog type \
                    :$(component) in $(file_path) is missing")
                )

                # Return an empty array for every missing block
                for block in blocks
                    unit = internalUnits(block, snapshot_path)
                    qty_data[block] = typeof(1.0 * unit)[]
                end

            end

            output[component] = qty_data

        end

    end

    return output

end

"""
    readSnapBlocks(
        file_path::String,
        request::Dict{Symbol,Vector{String}},
    )::Dict{Symbol,Dict{String,VecOrMat{<:Number}}}

Read the specified blocks from a snapshot file.

# Arguments

  - `file_path::String`: Path to the snapshot file.
  - `request::Dict{Symbol,Vector{String}}`: The blocks to be read. It must have the shape `cell/particle type` -> [`block`, `block`, `block`].

# Returns

  - A dictionary with the following shape: `cell/particle type` -> (`block` -> data of `block`).
"""
function readSnapBlocks(
    file_path::String,
    request::Dict{Symbol,Vector{String}},
)::Dict{Symbol,Dict{String,VecOrMat{<:Number}}}

    if isfile(file_path)

        (
            HDF5.ishdf5(file_path) ||
            throw(ArgumentError("readSnapBlocks: The file $(file_path) is not in the \
            HDF5 format, I don't know how to read it"))
        )

    else

        throw(ArgumentError("readSnapBlocks: $(file_path) does not exist as a file"))

    end

    # Allocate memory
    output = Dict{Symbol,Dict{String,VecOrMat{<:Number}}}()

    # Read the header
    header = readSnapHeader(file_path)

    h5open(file_path, "r") do snapshot

        # Read from the request only the cell/particle types
        for component in snapshotTypes(request)

            blocks = copy(request[component])

            type_str = PARTICLE_CODE_NAME[component]

            # Allocate memory
            qty_data = Dict{String,VecOrMat{<:Number}}()

            if type_str ∈ keys(snapshot)

                # Read the HDF5 group
                hdf5_group = snapshot[type_str]

                if isempty(hdf5_group)

                    (
                        !logging[] ||
                        @warn("readSnapBlocks: The cell/particle type \
                        :$(component) in $(file_path) is empty")
                    )

                    # Return an empty array for every missing block
                    for block in blocks
                        unit = internalUnits(block, snapshot_path)
                        qty_data[block] = typeof(1.0 * unit)[]
                    end

                else

                    # For the stellar particles, exclude wind particles
                    if component == :stars
                        idxs = findRealStars(file_path)
                    else
                        idxs = (:)
                    end

                    for block in blocks

                        unit = internalUnits(block, file_path)

                        if block == "TEMP"

                            (
                                component == :gas ||
                                throw(ArgumentError("readSnapBlocks: I can't compute the \
                                temperature for cells/particles of type :$(component), \
                                only for :gas"))
                            )

                            qty_data["TEMP"] = readTemperature(file_path)

                        elseif block == "MASS"

                            # Read the mass table from the header
                            mass_table = header.mass_table[PARTICLE_INDEX[component] + 1]

                            if iszero(mass_table)

                                raw = read(hdf5_group, QUANTITIES["MASS"].hdf5_name)
                                qty_data["MASS"] = selectdim(raw, ndims(raw), idxs) .* unit

                            else

                                # All cell/particles have the same mass
                                cp_number = header.num_part[PARTICLE_INDEX[component] + 1]
                                qty_data["MASS"] = fill(mass_table, cp_number) .* unit

                            end

                        elseif isBlockPresent(block, hdf5_group)

                            raw = read(hdf5_group, QUANTITIES[block].hdf5_name)
                            qty_data[block] = selectdim(raw, ndims(raw), idxs) .* unit

                        else

                            (
                                !logging[] ||
                                @warn("readSnapBlocks: The block $(block) for the \
                                cell/particle type :$(component) in $(file_path) is missing")
                            )

                            # Return an empty array for every missing block
                            qty_data[block] = typeof(1.0 * unit)[]

                        end

                    end

                end

            else

                (
                    !logging[] ||
                    @warn("readSnapBlocks: The cell/particle type \
                    :$(component) in $(file_path) is missing")
                )

                for block in blocks
                    unit = internalUnits(block, file_path)
                    qty_data[block] = typeof(1.0 * unit)[]
                end

            end

            output[component] = qty_data

        end

    end

    return output

end

"""
    readGroupCatalog(
        path::Union{String,Missing},
        snapshot_path::String,
        request::Dict{Symbol,Vector{String}},
    )::Dict{Symbol,Dict{String,VecOrMat{<:Number}}}

Read the specified blocks from a group catalog file or folder.

# Arguments

  - `path::Union{String,Missing}`: Path to the group catalog file or folder.
  - `snapshot_path::String`: Path to the corresponding snapshot file or folder. This is needed for unit conversion.
  - `request::Dict{Symbol,Vector{String}}`: Which blocks will be read. It must have the shape `group type` -> [`block`, `block`, `block`].

# Returns

  - A dictionary with the following shape: `group type` -> (`block` -> data of `block`).
"""
function readGroupCatalog(
    path::Union{String,Missing},
    snapshot_path::String,
    request::Dict{Symbol,Vector{String}},
)::Dict{Symbol,Dict{String,VecOrMat{<:Number}}}

    if ismissing(path)

        (
            !logging[] ||
            @warn("readGroupCatalog: The group catalog file or folder is missing")
        )

        return Dict{Symbol,Dict{String,VecOrMat{<:Number}}}()

    elseif isfile(path)

        (
            HDF5.ishdf5(path) ||
            throw(ArgumentError("readGroupCatalog: The file $(path) is not in the HDF5 format, \
            I don't know how to read it"))
        )

        return readGoupCatBlocks(path, snapshot_path, request)

    elseif isdir(path)

        sub_files = glob("$(GC_BASENAME)_*.*.hdf5", path)

        (
            !isempty(sub_files) && all(HDF5.ishdf5, sub_files) ||
            throw(ArgumentError("readGroupCatalog: The directory $(path) does not contain \
            group catalog sub-files in the HDF5 format"))
        )

        # Sort the sub files to concatenate the data in them correctly
        sort!(sub_files)

        # Read the data in each sub file
        data_in_files = [
            readGoupCatBlocks(file, snapshot_path, request) for file in sub_files
        ]

        # Allocate memory
        output = Dict{Symbol,Dict{String,VecOrMat{<:Number}}}()

        for (component, data_blocks) in first(data_in_files)

            qty_data = Dict{String,VecOrMat{<:Number}}()

            for block in keys(data_blocks)

                raw = [
                    data_in_file[component][block] for
                    data_in_file in data_in_files if !isempty(data_in_file[component][block])
                ]

                qty_data[block] = cat(raw...; dims=ndims(first(raw)))

            end

            output[component] = qty_data

        end

        return output

    else

        throw(ArgumentError("readGroupCatalog: $(path) does not exist as a file or folder"))

    end

end

"""
    readSnapshot(
        path::Union{String,Missing},
        request::Dict{Symbol,Vector{String}},
    )::Dict{Symbol,Dict{String,VecOrMat{<:Number}}}

Read the specified blocks from a snapshot file or folder.

# Arguments

  - `path::Union{String,Missing}`: Path to the snapshot file or folder.
  - `request::Dict{Symbol,Vector{String}}`: Which blocks will be read. It must have the shape `cell/particle type` -> [`block`, `block`, `block`].

# Returns

  - A dictionary with the following shape: `cell/particle type` -> (`block` -> data of `block`).
"""
function readSnapshot(
    path::Union{String,Missing},
    request::Dict{Symbol,Vector{String}},
)::Dict{Symbol,Dict{String,VecOrMat{<:Number}}}

    if ismissing(path)

        throw(ArgumentError("readSnapshot: The snapshot file or folder is missing"))

    elseif isfile(path)

        (
            HDF5.ishdf5(path) ||
            throw(ArgumentError("readSnapshot: The file $(path) is not in the HDF5 format, \
            I don't know how to read it"))
        )

        return readSnapBlocks(path, request)

    elseif isdir(path)

        sub_files = glob("$(SNAP_BASENAME)_*.*.hdf5", path)

        (
            !isempty(sub_files) && all(HDF5.ishdf5, sub_files) ||
            throw(ArgumentError("readSnapshot: The directory $(path) does not contain \
            snapshot sub-files in the HDF5 format"))
        )

        # Sort the sub files to concatenate the data in them correctly
        sort!(sub_files)

        # Read the data in each sub file
        data_in_files = [readSnapBlocks(file, request) for file in sub_files]

        # Allocate memory
        output = Dict{Symbol,Dict{String,VecOrMat{<:Number}}}()

        for (component, data_blocks) in first(data_in_files)

            qty_data = Dict{String,VecOrMat{<:Number}}()

            for block in keys(data_blocks)

                raw = [
                    data_in_file[component][block] for
                    data_in_file in data_in_files if !isempty(data_in_file[component][block])
                ]

                if isempty(raw)
                    qty_data[block] = Number[]
                else
                    qty_data[block] = cat(raw...; dims=ndims(first(raw)))
                end

            end

            output[component] = qty_data

        end

        return output

    else

        throw(ArgumentError("readSnapshot: $(path) does not exist as a file or folder"))

    end

end

"""
    getBlock(path::String, component::Symbol, block::String)::VecOrMat{<:Number}

Convenience function to directly get the data associated with one block.

# Arguments

  - `path::String`: Path to the snapshot file or folder.
  - `component::Symbol`: Type of cell/particle. The possibilities are the keys of [`PARTICLE_INDEX`](@ref).
  - `block::String`: Target block. The possibilities are the keys of [`QUANTITIES`](@ref).

# Returns

  - The data for `block`.
"""
function getBlock(path::String, component::Symbol, block::String)::VecOrMat{<:Number}

    return readSnapshot(path, Dict(component => [block]))[component][block]

end

"""
    readSfrFile(
        file_path::String,
        snap_path::String,
    )::Dict{Int32,VecOrMat{<:Number}}

Read the `sfr.txt` file.

# Arguments

  - `file_path::String`: Path to the `sfr.txt` file.
  - `snapshot_path::String`: Path to one snapshot file or folder of the simulation. This is needed for unit conversion.

# Returns

  - A dictionary with the following shape:

      + `1` -> Time or scale factor (internal units).
      + `2` -> Total stellar mass to be formed prior to stochastic sampling (internal units).
      + `3` -> Instantaneous star formation rate of all cells (``\\mathrm{M_\\odot \\, yr^{-1}}``).
      + `4` -> Instantaneous star formation rate of active cells (``\\mathrm{M_\\odot \\, yr^{-1}}``).
      + `5` -> Total mass in stars formed after stochastic sampling (internal units).
      + `6` -> Cumulative stellar mass formed (internal units).
"""
function readSfrFile(
    file_path::String,
    snap_path::String,
)::Dict{Int32,VecOrMat{<:Number}}

    (
        isfile(file_path) ||
        throw(ArgumentError("readSfrFile: $(file_path) does not exist as a file"))
    )

    # Load the data from the `sfr.txt` file
    file_data = readdlm(file_path, Float64)
    n_cols = size(file_data, 2)

    # Check that the data in the file has the correct size
    (
        n_cols <= 6 ||
        throw(ArgumentError("readSfrFile: I don't know how to handle more \
        than 6 columns in `sfr.txt`"))
    )
    (
        !(logging[] && n_cols < 6) ||
        @warn("readSfrFile: I could only find $(n_cols) columns \
        in $(flie_path). I was expecting 6")
    )

    # Load the units for each column
    units = [internalUnits("SFC$(i)", snap_path) for i in 1:n_cols]

    return Dict(i => column .* units[i] for (i, column) in pairs(eachcol(file_data)))

end

"""
    readCpuFile(
        file_path::String,
        targets::Vector{String};
        <keyword arguments>
    )::Dict{String,Matrix{Float64}}

Read the `cpu.txt` file.

For each process in `targets` a matrix with all the CPU usage data is returned.

# Arguments

  - `file_path::String`: Path to the `cpu.txt` file.
  - `targets::Vector{String}`: Target processes.
  - `step::Int=1`: Step used to traverse the rows.

# Returns

  - A dictionary with the following shape:

    `target process` -> matrix with columns:

      + Time step.
      + Simulation time (scale factor for cosmological simulations and physical time for non-cosmological simulations).
      + Clock time in seconds.
      + Clock time as a percentage.
      + Total clock time in seconds.
      + Total clock time as a percentage.
"""
function readCpuFile(
    file_path::String,
    targets::Vector{String};
    step::Int=1,
)::Dict{String,Matrix{Float64}}

    (
        isfile(file_path) ||
        throw(ArgumentError("readCpuFile: $(file_path) does not exist as a file"))
    )

    # Load the data from the `cpu.txt` file
    file_data = eachline(file_path)

    # Set up an auxiliary dictionary
    data_aux = Dict(target => Matrix{Float64}[] for target in targets)

    # Clock time for each sync-point
    time = 0.0
    # Time step
    time_step = 0

    for line in file_data

        # Ignore empty lines
        !isempty(line) || continue

        columns = split(line)
        title = columns[1]

        # Ignore header lines
        !(title == "diff") || continue

        # Use "Step" lines to capture the clock time
        if title == "Step"
            time = parse(Float64, rstrip(columns[4], ','))
            time_step = parse(Int, rstrip(columns[2], ','))
            continue
        end

        if title ∈ targets
            push!(
                data_aux[title],
                [
                    time_step;;                               # Time step
                    time;;                                    # Simulation time for each sync-point
                    parse(Float64, columns[2]);;              # Clock time in seconds
                    parse(Float64, rstrip(columns[3], '%'));; # Clock time as a percentage
                    parse(Float64, columns[4]);;              # Total clock time in seconds
                    parse(Float64, rstrip(columns[5], '%'))   # Total clock time as a percentage
                ],
            )
        end

    end

    (
        !(logging[] && any(isempty, values(data_aux))) ||
        @warn("readCpuFile: I could not find some of the target rows in $(file_path)")
    )

    # Allocate memory
    data_out = Dict{String,Matrix{Float64}}()

    # Try reducing the data size
    for (target, values) in data_aux

        # Ignore empty targets
        !isempty(values) || continue

        l_e = length(values)

        if 1 < step < l_e
            data_out[target] = vcat(values[1:step:end]...)
        else
            data_out[target] = vcat(values...)
            (
                !(logging[] && step > l_e) ||
                @warn("readCpuFile: `step` = $(step) is bigger than the number \
                of time steps in $(file_path)")
            )
        end

    end

    return data_out

end

"""
    getSnapshotPaths(simulation_path::String)::Dict{Symbol,Vector{String}}

Find the path and number of every snapshot in `simulation_path`.

!!! note

    If each snapshot is made of multiple files, the `:paths` field will have paths to folders, each one containing the sub-files of the corresponding snapshot.

# Arguments

  - `simulation_path::String`: Path to the simulation directory, set in the code variable `OutputDir`.

# Returns

  - A dictionary with the following shape:

      + `:numbers` -> The number that characterize each snapshot.
      + `:paths`   -> The full path to each snapshot.
"""
function getSnapshotPaths(simulation_path::String)::Dict{Symbol,Vector{String}}

    (
        isdir(simulation_path) ||
        throw(ArgumentError("getSnapshotPaths: $(simulation_path) does not exist as a directory"))
    )

    # Get the full list of paths to every snapshot in `simulation_path`
    path_list = [
        glob("*/*/$(SNAP_BASENAME)_*", simulation_path)
        glob("*/$(SNAP_BASENAME)_*", simulation_path)
        glob("$(SNAP_BASENAME)_*", simulation_path)
    ]

    # Check for an empty folder
    if isempty(path_list)

        (
            !logging[] ||
            @warn("getSnapshotPaths: I could not find any file named $(SNAP_BASENAME)_*.hdf5 \
            within $(simulation_path), or any of its subfolders")
        )

        return Dict(:numbers => String[], :paths => String[])

    end

    # Get the numbers that characterize each snapshot
    reg = Regex("(?<=$(SNAP_BASENAME)_).*?(?=(?:\\.)|\$)")
    number_list = map(x -> match(reg, x).match, path_list)

    if readSnapHeader(first(path_list)).num_files > 1
        # If there are multiple files per snapshot, get the path to the snapshot directory
        map!(dirname, path_list, path_list)
        # Delete duplicates
        unique!(path_list)
        unique!(number_list)
    end

    return Dict(:numbers => number_list, :paths => sort(normpath.(path_list)))

end

"""
    getGroupCatPaths(simulation_path::String)::Dict{Symbol,Vector{String}}

Find the path and number of every group catalog in `simulation_path`.

!!! note

    If each group catalog is made of multiple files, the `:paths` field will have paths to folders, each one containing the sub-files of the corresponding group catalog.

# Arguments

  - `simulation_path::String`: Path to the simulation directory, set in the code variable `OutputDir`.

# Returns

  - A dictionary with the following shape:

      + `:numbers` -> The number that characterize each group catalog.
      + `:paths`   -> The full path to each group catalog.
"""
function getGroupCatPaths(simulation_path::String)::Dict{Symbol,Vector{String}}

    (
        isdir(simulation_path) ||
        throw(ArgumentError("getGroupCatPaths: $(simulation_path) does not exist as a directory"))
    )

    # Get the full list of paths to every group catalog in `simulation_path`
    path_list = [
        glob("*/*/$(GC_BASENAME)_*", simulation_path)
        glob("*/$(GC_BASENAME)_*", simulation_path)
        glob("$(GC_BASENAME)_*", simulation_path)
    ]

    # Check for an empty folder
    if isempty(path_list)

        (
            !logging[] ||
            @warn("getGroupCatPaths: I could not find any file named $(GC_BASENAME)_*.hdf5 \
            within $(simulation_path), or any of its subfolders")
        )

        return Dict(:numbers => String[], :paths => String[])

    end

    # Get the numbers that characterize each group catalog
    reg = Regex("(?<=$(GC_BASENAME)_).*?(?=(?:\\.)|\$)")
    number_list = map(x -> match(reg, x).match, path_list)

    if readGroupCatHeader(first(path_list)).num_files > 1
        # If there are multiple files per group catalog, get the path to the group catalog directory
        map!(dirname, path_list, path_list)
        # Delete duplicates
        unique!(path_list)
        unique!(number_list)
    end

    return Dict(:numbers => number_list, :paths => sort(normpath.(path_list)))

end

"""
    makeSimulationTable(simulation_path::String)::DataFrame

Construct a dataframe with the path, time stamps and number of each snapshot and group catalog file in `simulation_path`.

# Arguments

  - `simulation_path::String`: Path to the simulation directory, set in the code variable `OutputDir`.

# Returns

  - A dataframe with 8 columns:

      + `:ids`            -> Dataframe index of each snapshot, i.e. if there are 10 snapshots in total it runs from 1 to 10.
      + `:numbers`        -> Number in the file name of each snapshot.
      + `:scale_factors`  -> Scale factor of each snapshot.
      + `:redshifts`      -> Redshift of each snapshot.
      + `:physical_times` -> Physical time since the Big Bang of each snapshot.
      + `:lookback_times` -> Physical time left to reach the last snapshot.
      + `:snapshot_paths` -> Full path to the snapshots.
      + `:groupcat_paths` -> Full path to the group catalog files.
"""
function makeSimulationTable(simulation_path::String)::DataFrame

    # Get the path and number of each snapshot
    snap_source = getSnapshotPaths(simulation_path)
    snapshot_paths = isempty(snap_source[:paths]) ? [missing] : snap_source[:paths]

    # Get the path and number of each group catalog file
    groupcat_source = getGroupCatPaths(simulation_path)
    groupcat_paths  = isempty(groupcat_source[:paths]) ? [missing] : groupcat_source[:paths]

    (
        length(snapshot_paths) >= length(groupcat_paths) ||
        throw(ArgumentError("makeSimulationTable: I found less snapshots \
        ($(length(snapshot_paths))) than group catalogs ($(length(groupcat_paths))) in \
        $(simulation_path), I cannot make the table when not every group catalog has a \
        corresponding snapshot"))
    )

    paths  = [snapshot_paths, groupcat_paths]
    labels = [:snapshot_paths, :groupcat_paths]
    rows   = [[1:length(snapshot_paths);], [1:length(groupcat_paths);]]

    source_table = unstack(flatten(DataFrame(l=labels, p=paths, ids=rows), [:p, :ids]), :l, :p)

    # Add the file name number column
    numbers = snap_source[:numbers]
    insertcols!(source_table, 2, :numbers => isempty(numbers) ? ["000"] : numbers; copycols=false)

    # Get the time stamps of every snapshot
    scale_factors, redshifts, physical_times, lookback_times = computeTimeTicks(snapshot_paths)

    # Add the scale factor column
    insertcols!(source_table, 3, :scale_factors => scale_factors; copycols=false)

    # Add the redshift column
    insertcols!(source_table, 4, :redshifts => redshifts; copycols=false)

    # Add the physical time column
    insertcols!(source_table, 5, :physical_times => physical_times; copycols=false)

    # Add the lookback time column
    insertcols!(source_table, 6, :lookback_times => lookback_times; copycols=false)

    return identity.(DataFrame(source_table))

end

"""
    makeDataDict(
        simulation_path::String,
        slice_n::Int,
        request::Dict{Symbol,Vector{String}},
    )::Dict

Construct a data dictionary for a single snapshot.

# Arguments

  - `simulation_path::String`: Path to the simulation directory, set in the code variable `OutputDir`.
  - `slice_n::Int`: Selects the target snapshot. Starts at 1 and is independent of the number in the file name. If every snapshot is present, the relation is `slice_n` = filename_number + 1.
  - `request::Dict{Symbol,Vector{String}}`: Dictionary with the shape `cell/particle type` -> [`block`, `block`, ...], where the possible types are the keys of [`PARTICLE_INDEX`](@ref), and the possible quantities are the keys of [`QUANTITIES`](@ref).

# Returns

  - A dictionary with the following shape:

      + `:sim_data`          -> ::Simulation (see [`Simulation`](@ref)).
      + `:snap_data`         -> ::Snapshot (see [`Snapshot`](@ref)).
      + `:gc_data`           -> ::GroupCatalog (see [`GroupCatalog`](@ref)).
      + `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `cell/particle type` -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + ...
      + `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + `groupcat type`      -> (`block` -> data of `block`, `block` -> data of `block`, ...).
      + ...
"""
function makeDataDict(
    simulation_path::String,
    slice_n::Int,
    request::Dict{Symbol,Vector{String}},
)::Dict

    # Make a dataframe for every simulation, with the following columns:
    #   - 1. DataFrame index
    #   - 2. Number in the file name
    #   - 3. Scale factor
    #   - 4. Redshift
    #   - 5. Physical time
    #   - 6. Lookback time
    #   - 7. Snapshot path
    #   - 8. Group catalog path
    simulation_table = makeSimulationTable(simulation_path)

    snapshot_numbers = simulation_table[!, :numbers]

    (
        length(snapshot_numbers) >= slice_n ||
        throw(ArgumentError("makeDataDict: The snapshot number $(slice_n) does not exist in  \
        $(simulation_path). There are only $(length(snapshot_numbers)) snapshots. \
        The full simulation table is:\n\n$(simulation_table)"))
    )

    # Select the target snapshot
    snapshot_row = simulation_table[slice_n, :]

    ################################################################################################
    # Compute the metadata for the current snapshot and simulation
    ################################################################################################

    # Get the snapshot file path
    snapshot_path = snapshot_row[:snapshot_paths]
    # Get the group catalog file path
    groupcat_path = snapshot_row[:groupcat_paths]

    # Store the metadata of the current snapshot and simulation
    metadata = Dict(
        :sim_data => Simulation(
            simulation_path,
            1,
            slice_n,
            isCosmological(snapshot_path),
            simulation_table,
        ),
        :snap_data => Snapshot(
            snapshot_path,
            slice_n,
            slice_n,
            snapshot_row[:physical_times],
            snapshot_row[:lookback_times],
            snapshot_row[:scale_factors],
            snapshot_row[:redshifts],
            readSnapHeader(snapshot_path),
        ),

        :gc_data => GroupCatalog(
            groupcat_path,
            readGroupCatHeader(groupcat_path),
        ),
    )

    return merge(
        metadata,
        readSnapshot(snapshot_path, request),
        readGroupCatalog(groupcat_path, snapshot_path, request),
    )

end

"""
    countSnapshot(simulation_path::String)::Int

Count the number of snapshots in `simulation_path`.

!!! note

    This function count the number of snapshots, no the number of snapshot files. So if each snapshot is made of more than one files, the count will not change.

# Arguments

  - `simulation_path::String`: Path to the simulation directory, set in the code variable `OutputDir`.

# Returns

  - The number of snapshots.
"""
function countSnapshot(simulation_path::String)::Int

    (
        isdir(simulation_path) ||
        throw(ArgumentError("countSnapshot: $(simulation_path) does not exist as a directory"))
    )

    # Get the full list of paths to every snapshot in `simulation_path`
    path_list = [
        glob("*/*/$(SNAP_BASENAME)_*", simulation_path)
        glob("*/$(SNAP_BASENAME)_*", simulation_path)
        glob("$(SNAP_BASENAME)_*", simulation_path)
    ]

    # Check for an empty folder
    if isempty(path_list)

        (
            !logging[] ||
            @warn("countSnapshot: I could not find any file named $(SNAP_BASENAME)_*.hdf5 \
            within $(simulation_path), or any of its subfolders")
        )

        return 0

    end

    if readSnapHeader(first(path_list)).num_files > 1
        # If there are multiple files per snapshot, get the path to the snapshot directory
        map!(dirname, path_list, path_list)
        # Delete duplicates
        unique!(path_list)
    end

    return length(path_list)

end

"""
    mergeRequests(requests::Dict{Symbol,Vector{String}}...)::Dict{Symbol,Vector{String}}

Merge several request dictionaries, ignoring duplicates.

# Arguments

  - `requests`: The request dictionaries for [`readSnapshot`](@ref).

# Returns

  - A new dictionary with all the requests.
"""
function mergeRequests(requests::Dict{Symbol,Vector{String}}...)::Dict{Symbol,Vector{String}}

    return Dict(
        type => union([get(request, type, String[]) for request in requests]...) for
        type in union(keys.(requests)...)
    )

end

"""
    addRequest(
        request::Dict{Symbol,Vector{String}},
        addition::Dict{Symbol,Vector{String}},
    )::Dict{Symbol,Vector{String}}

Add the blocks in `addition` to `request`, only for the types already present in `request`.

# Arguments

  - `request::Dict{Symbol,Vector{String}}`: The request dictionary for [`readSnapshot`](@ref).
  - `addition::Dict{Symbol,Vector{String}}`: Request dictionary with the blocks to be added, only for the types already present in `request`.

# Returns

  - A new dictionary with all the requests.
"""
function addRequest(
    request::Dict{Symbol,Vector{String}},
    addition::Dict{Symbol,Vector{String}},
)::Dict{Symbol,Vector{String}}

    return Dict(type => blocks ∪ get(addition, type, String[]) for (type, blocks) in request)

end

"""
    isSubfindActive(path::String)::Bool

Check if there is information about the halos and subhalos in the group catalog file.

# Arguments

  - `path::String`: Path to the group catalog file or folder.

# Returns

  - If there are halo and subhalo information in the group catalog file.
"""
function isSubfindActive(path::String)::Bool

    if isfile(path)

        (
            HDF5.ishdf5(path) ||
            throw(ArgumentError("isSubfindActive: The file $(path) is not in the HDF5 format, \
            I don't know how to read it"))
        )

        file_path = path

    elseif isdir(path)

        sub_files = glob("$(GC_BASENAME)_*.*.hdf5", path)

        (
            !isempty(sub_files) && all(HDF5.ishdf5, sub_files) ||
            throw(ArgumentError("isSubfindActive: The directory $(path) does not contain \
            group catalog sub-files in the HDF5 format"))
        )

        file_path = minimum(sub_files)

    else

        throw(ArgumentError("isSubfindActive: $(path) does not exist as a file or folder"))

    end

    subfind_active = h5open(file_path, "r") do gc_file

        (
            all(in(keys(gc_file)), ["Group", "Subhalo"]) &&
            all(!isempty, [gc_file["Group"], gc_file["Subhalo"]])
        )

    end

    return subfind_active

end

"""
    findRealStars(path::String)::Vector{Bool}

Find which stellar particles are real stars and not wind particles.

# Arguments

  - `path::String`: Path to the snapshot file or folder.

# Returns

  - A boolean vector with true for stars and false for wind particles.
"""
function findRealStars(path::String)::Vector{Bool}

    if isfile(path)

        (
            HDF5.ishdf5(path) ||
            throw(ArgumentError("findRealStars: The file $(path) is not in the \
            HDF5 format, I don't know how to read it"))
        )

        time_of_birth = h5open(path, "r") do snapshot
            if PARTICLE_CODE_NAME[:stars] ∉ keys(snapshot)
                Float64[]
            else
                read(snapshot[PARTICLE_CODE_NAME[:stars]], QUANTITIES["GAGE"].hdf5_name)
            end
        end

        return isempty(time_of_birth) ? Bool[] : map(isPositive, time_of_birth)

    elseif isdir(path)

        sub_files = glob("$(SNAP_BASENAME)_*.*.hdf5", path)

        (
            !isempty(sub_files) && all(HDF5.ishdf5, sub_files) ||
            throw(ArgumentError("findRealStars: The directory $(path) does not contain \
            snapshot sub-files in the HDF5 format"))
        )

        return vcat([findRealStars(sub_file) for sub_file in sub_files]...)

    else

        throw(ArgumentError("findRealStars: $(path) does not exist as a file or folder"))

    end

end

"""
    countStars(path::String)::Int

Count the number of stars in a snapshot, excluding wind particles.

# Arguments

  - `path::String`: Path to the snapshot file or folder.

# Returns

  - The number of stars.
"""
countStars(path::String)::Int = count(findRealStars(path))

"""
    findQtyExtrema(
        simulation_path::String,
        slice_n::Int,
        component::Symbol,
        block::String;
        <keyword arguments>
    )::NTuple{2,<:Number}

Compute the minimum and maximum values of `block`.

# Arguments

  - `simulation_path::String`: Path to the simulation directory, set in the code variable `OutputDir`.
  - `slice_n::Int`: Selects which snapshot to plot, starts at 1 and is independent of the number in the file name. If every snapshot is present, the relation is `slice_n` = filename_number + 1. If set to a negative number, the values in the whole simulation will be compared.
  - `component::Symbol`: Cell/particle type. The possibilities are the keys of [`PARTICLE_INDEX`](@ref).
  - `block::String`: Target block. The possibilities are the keys of [`QUANTITIES`](@ref).
  - `f::Function=identity`: A function with the signature:

    `f(data) -> values`

    where

      + `data::VecOrMat{<:Number}`: Data returned by [`getBlock`](@ref).
      + `values::Vector{<:Number}`: A vector with the values to be compared.

# Returns

  - Tuple with the minimum and maximum values.
"""
function findQtyExtrema(
    simulation_path::String,
    slice_n::Int,
    component::Symbol,
    block::String;
    f::Function=identity,
)::NTuple{2,<:Number}

    (
        isdir(simulation_path) ||
        throw(ArgumentError("findQtyExtrema: $(simulation_path) does not exist as a directory"))
    )

    simulation_table = makeSimulationTable(simulation_path)

    if slice_n > 0

        # Get the number in the filename
        snap_n = safeSelect(simulation_table[!, :numbers], slice_n)

        # Check that after slicing there is one snapshot left
        (
            !isempty(snap_n) ||
            throw(ArgumentError("findQtyExtrema: There are no snapshots with `slice_n` = \
            $(slice_n), the contents of $(simulation_path) are: \n$(simulation_table)"))
        )

        # Find the target row and snapshot path
        snapshot_row = filter(:numbers => ==(lpad(snap_n, 3, "0")), simulation_table)
        snapshot_path = snapshot_row[1, :snapshot_paths]

        (
            !ismissing(snapshot_path) ||
            throw(ArgumentError("findQtyExtrema: The snapshot number $(slice_n) seems \
            to be missing"))

        )

        values = f(getBlock(snapshot_path, component, block))

        return extrema(values)

    end

    snapshot_paths = filter!(!ismissing, snapshot_row[!, :snapshot_paths])

    (
        !isempty(snapshot_paths) ||
        throw(ArgumentError("findQtyExtrema: I could not find any snapshots in $(simulation_path)"))
    )

    values = [f(getBlock(snapshot_path, component, block)) for snapshot_path in snapshot_paths]

    return extrema(Iterators.flatten(values))

end

"""
    isCosmological(path::String)::Bool

Check if the snapshot in `path` comes from a cosmological simulation.

!!! note

    If each snapshot is made of multiple files, I'll read the first chunk to check if the simulation is cosmological.

# Arguments

  - `path::String`: Path to the snapshot file or folder.

# Returns

  - If the simulation is cosmological

      + `false` -> Newtonian simulation    (`ComovingIntegrationOn` = 0, `Redshift` = 0.0).
      + `true`  -> Cosmological simulation (`ComovingIntegrationOn` = 1, `Redshift` != 0.0).
"""
function isCosmological(path::String)::Bool

    if isfile(path)

        (
            HDF5.ishdf5(path) ||
            throw(ArgumentError("isCosmological: The file $(path) is not in the HDF5 format, \
            I don't know how to read it"))
        )

        file_path = path

    elseif isdir(path)

        sub_files = glob("$(SNAP_BASENAME)_*.*.hdf5", path)

        (
            !isempty(sub_files) && all(HDF5.ishdf5, sub_files) ||
            throw(ArgumentError("isCosmological: The directory $(path) does not contain \
            snapshot sub-files in the HDF5 format"))
        )

        file_path = minimum(sub_files)

    else

        throw(ArgumentError("isCosmological: $(path) does not exist as a file or folder"))

    end

    cosmological = h5open(file_path, "r") do snapshot

        if "Parameters" ∈ keys(snapshot)
            # If the param.txt is saved in the snapshot metadata, read `ComovingIntegrationOn`
            read_attribute(snapshot["Parameters"], "ComovingIntegrationOn")
        else
            # Otherwise, use the readshift in the header
            !iszero(read_attribute(snapshot["Header"], "Redshift"))
        end

    end

    return cosmological

end

"""
    internalUnits(quantity::String, path::String)::Union{Unitful.Quantity,Unitful.Units}

Get the factor to convert a plain number into a [Unitful](https://github.com/PainterQubits/Unitful.jl) quantity, using the correct internal code units.

# Arguments

  - `quantity::String`: Target quantity. The options are the keys of [`QUANTITIES`](@ref).
  - `path::String`: Path to the snapshot file or folder.

# Returns

  - A [Unitful](https://github.com/PainterQubits/Unitful.jl) quantity or unit.
"""
function internalUnits(quantity::String, path::String)::Union{Unitful.Quantity,Unitful.Units}

    (
        quantity ∈ keys(QUANTITIES) ||
        throw(ArgumentError("internalUnits: `quantity` should be one of the keys of \
        `QUANTITIES` but I got $(quantity), see the options in `./src/constants/globals.jl`"))
    )

    header = readSnapHeader(path)
    cosmological = isCosmological(path)

    a0 = cosmological ? header.time : 1.0
    h0 = cosmological ? header.h0 : 1.0

    # Set up the struct for unit conversion
    IU = InternalUnits(; l_unit=header.l_unit, m_unit=header.m_unit, v_unit=header.v_unit, a0, h0)

    dimensions = QUANTITIES[quantity].dimensions
    unit = QUANTITIES[quantity].unit

    if unit == :internal
        if dimensions == Unitful.𝐌

            # From internal units to M⊙
            return IU.m_cosmo

        elseif dimensions == Unitful.𝐋

            if !PHYSICAL_UNITS && !cosmological
                @warn(
                    "internalUnits: You have set the unit system to use comoving lengths \
                    (`PHYSICAL_UNITS` = $(PHYSICAL_UNITS)), but the simulation is not \
                    cosmological. I'll keep the lengths physical. Check `PHYSICAL_UNITS` \
                    in `constants/globals.jl`",
                    maxlog=1,
                )
            end

            # From internal units to kpc
            if !PHYSICAL_UNITS && cosmological
                return IU.x_comoving
            else
                return IU.x_cosmo
            end

        elseif dimensions == Unitful.𝐓

            # From internal units to Myr, for non-cosmological simulations,
            # and to a dimensionless quantity for cosmological simulations
            return cosmological ? Unitful.NoUnits : IU.t_cosmo

        elseif dimensions == Unitful.𝐌 * Unitful.𝐋^-3

            # From internal units to g * cm^-3
            return IU.rho_cgs

        elseif dimensions == Unitful.𝐋^2 * Unitful.𝐓^-2

            # From internal units to erg * g^-1
            return IU.U_cgs

        elseif dimensions == Unitful.𝐋 * Unitful.𝐓^-1

            # From internal units to km * s^-1
            return IU.v_cosmo

        elseif dimensions == Unitful.𝐌 * Unitful.𝐋^-1 * Unitful.𝐓^-2

            # From internal units to Pa
            return IU.P_Pa

        else

            error("internalUnits: I don't know the internal units of a quantity \
            with dimensions $(dimensions)")

        end

    elseif unit == :gvel

        # Special case for "G_Vel" (velocity of the group)
        # See the TNG documentation https://www.tng-project.org/data/docs/specifications/
        return IU.v_cosmo / a0^1.5

    else

        return unit

    end

end

"""
    snapshotTypes(data_dict::Dict)::Vector{Symbol}

Find which cell/particle types are part of the keys of `data_dict`.

# Arguments

  - `data_dict::Dict`: A dictionary.

# Returns

  - A vector with the cell/particle types.
"""
snapshotTypes(data_dict::Dict)::Vector{Symbol} = collect(keys(PARTICLE_INDEX) ∩ keys(data_dict))

"""
    snapshotTypes(path::String)::Vector{Symbol}

Find which cell/particle types are part of the snapshot in `path`.

!!! note

    If each snapshot is made of multiple files, I'll check the first chunk.

# Arguments

  - `path::String`: Path to the snapshot file or folder.

# Returns

  - A vector with the cell/particle types.
"""
function snapshotTypes(path::String)::Vector{Symbol}

    if isfile(path)

        (
            HDF5.ishdf5(path) ||
            throw(ArgumentError("snapshotTypes: The file $(path) is not in the HDF5 format, \
            I don't know how to read it"))
        )

        file_path = path

    elseif isdir(path)

        sub_files = glob("$(SNAP_BASENAME)_*.*.hdf5", path)

        (
            !isempty(sub_files) && all(HDF5.ishdf5, sub_files) ||
            throw(ArgumentError("snapshotTypes: The directory $(path) does not contain \
            snapshot sub-files in the HDF5 format"))
        )

        file_path = minimum(sub_files)

    else

        throw(ArgumentError("snapshotTypes: $(path) does not exist as a file or folder"))

    end

    snapshot_types = h5open(file_path, "r") do snapshot
        collect(keys(PARTICLE_TYPE) ∩ keys(snapshot))
    end

    return snapshot_types

end

"""
    groupcatTypes(data_dict::Dict)::Vector{Symbol}

Find which group catalog data types are part of the keys of `data_dict`.

# Arguments

  - `data_dict::Dict`: A dictionary.

# Returns

  - A vector with the group catalog data types.
"""
groupcatTypes(data_dict::Dict)::Vector{Symbol} = [:group, :subhalo] ∩ keys(data_dict)

"""
    groupcatTypes(path::String)::Vector{Symbol}

Find which group catalog data types are part of the snapshot in `path`.

!!! note

    If each snapshot is made of multiple files, I'll check the first chunk.

# Arguments

  - `path::String`: Path to the snapshot file or folder.

# Returns

  - A vector with the group catalog data types.
"""
function groupcatTypes(path::String)::Vector{Symbol}

    if isfile(path)

        (
            HDF5.ishdf5(path) ||
            throw(ArgumentError("groupcatTypes: The file $(path) is not in the HDF5 format, \
            I don't know how to read it"))
        )

        file_path = path

    elseif isdir(path)

        sub_files = glob("$(SNAP_BASENAME)_*.*.hdf5", path)

        (
            !isempty(sub_files) && all(HDF5.ishdf5, sub_files) ||
            throw(ArgumentError("groupcatTypes: The directory $(path) does not contain \
            snapshot sub-files in the HDF5 format"))
        )

        file_path = minimum(sub_files)

    else

        throw(ArgumentError("groupcatTypes: $(path) does not exist as a file or folder"))

    end

    groupcat_types = h5open(file_path, "r") do snapshot
        [:group, :subhalo] ∩ keys(snapshot)
    end

    return groupcat_types

end
