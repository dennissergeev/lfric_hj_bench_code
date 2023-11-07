#!/usr/bin/env python
"""Process LFRic output by interpolating selected fields to a common grid."""
# Standard library
from functools import partial
import os
from pathlib import Path
from time import time
import warnings

# External libraries
from aeolus.const import add_planet_conf_to_cubes, init_const
from aeolus.coord import get_cube_rel_days
from aeolus.io import create_dummy_cube, save_cubelist
from aeolus.lfric import (
    add_equally_spaced_height_coord,
    add_um_height_coord,
    fix_time_coord,
    load_lfric_raw,
    simple_regrid_lfric,
)
from aeolus.log import create_logger
from aeolus.subset import unique_cubes
import click

# Local modules
import paths

# Global definitions and styles
warnings.filterwarnings("ignore")

os.environ["PROJ_IGNORE_CELESTIAL_BODY"] = "YES"


def chunker(seq, size):
    """Iterate over a sequence in chunks."""
    return (seq[pos : pos + size] for pos in range(0, len(seq), size))


@click.command()
@click.option("-i", "--inpdir", type=str, help="Input directory")
@click.option("-o", "--outdir", type=str, help="Output directory")
@click.option(
    "-l", "--label", type=str, required=True, help="Simulation label"
)
@click.option(
    "-p", "--planet", type=str, required=True, help="Planet configuration"
)
@click.option(
    "-c", "--c_num", type=str, default="C48", help="Cubed Sphere Mesh Number"
)
@click.option(
    "--ref_cube",
    type=str,
    default="air_potential_temperature",
    help="Reference cube, to which coordinates all data will be regridded",
)
@click.option(
    "--level_height",
    type=str,
    default="uniform",
    help="Type of the vertical level height coordinate",
)
@click.option(
    "--model_top",
    type=float,
    default=40e3,
    help="If level_height=uniform, set the model top height.",
)
@click.option(
    "--time_prof",
    default="inst",
    help="Type of the time output",
    show_default=True,
    type=click.Choice(["inst", "mean"]),
)
@click.option(
    "--file_chunk_size",
    type=int,
    default=1000,
    help="Number of files to process in one go.",
)
def main(
    inpdir,
    outdir,
    label,
    planet,
    c_num,
    ref_cube,
    level_height,
    model_top,
    time_prof,
    file_chunk_size,
):
    """Main entry point of the script."""
    t0 = time()
    L = create_logger(Path(__file__))
    L.info(f"{planet=}")
    L.info(f"{label=}")
    # Global attributes of the output file
    gl_attrs = {
        "name": label,
        "planet": planet,
        "processed": "True",
    }
    # Create a dummy cube with a target grid
    tgt_cube = create_dummy_cube(nlat=90, nlon=144, pm180=True)

    # Height coordinate
    if level_height == "uniform":
        add_levs = partial(
            add_equally_spaced_height_coord, model_top_height=model_top
        )
    elif level_height.startswith("um_"):
        add_levs = partial(
            add_um_height_coord,
            path_to_levels_file=paths.vert
            / level_height.replace("um_", "vertlevs_"),
        )
    else:
        raise ValueError(f"{level_height=} is not valid.")

    def combi_callback(cube, field, filename):
        [
            fix_time_coord(cube, field, filename),
            add_levs(cube, field, filename),
        ]

    # Input directory
    if inpdir:
        inpdir = Path(inpdir)
    else:
        inpdir = paths.results_raw_lfric / label / c_num
    L.info(f"{inpdir=}")

    # File names
    if time_prof == "inst":
        file_stem = "lfric*diag"
        drop_coord = ["forecast_reference_time"]
    elif time_prof == "mean":
        file_stem = "lfric_averages"
        drop_coord = []
    else:
        raise ValueError(f"{time_prof=} is not valid.")

    # Create a subdirectory for processed data
    if outdir:
        outdir = Path(outdir)
    else:
        outdir = paths.results_proc_lfric / label / c_num
        L.info(f"{outdir=}")
    outdir.mkdir(parents=True, exist_ok=True)

    # Compile the list of files
    fnames = sorted(
        inpdir.glob(f"*/*{c_num}*/{file_stem}.nc"),
        key=lambda x: int(x.parent.parent.name),
    )
    if len(fnames) == 0:
        L.critical("No files found!")
        return
    elif len(fnames) == 1:
        L.info(f"fnames({len(fnames)}) = {fnames[0]}")
    else:
        L.info(f"fnames({len(fnames)}) = {fnames[0]} ... {fnames[-1]}")

    # Load the data in chunks
    for _ichunk, file_batch in enumerate(chunker(fnames, file_chunk_size)):
        L.info(
            f"Loading {file_chunk_size} files ="
            f" {file_batch[0]} ... {file_batch[-1]}"
        )
        # Make a label depending on the chunk
        chunk_label = f"_{int(file_batch[0].parent.parent.name):03d}"
        chunk_label += f"-{int(file_batch[-1].parent.parent.name):03d}"
        gl_attrs["source_files"] = ", ".join(
            [str(Path(*f.parts[-3:])) for f in file_batch]
        )
        # Load the data
        cl_raw = load_lfric_raw(
            file_batch, callback=combi_callback, drop_coord=drop_coord
        )
        if len(cl_raw) == 0:
            L.critical("Files are empty!")
            continue
        cubes_to_regrid = unique_cubes(cl_raw)
        if len(w_cubes := cubes_to_regrid.extract("w_in_wth")) == 2:
            cubes_to_regrid.remove(w_cubes[-1])
        for cube in cubes_to_regrid:
            if cube.units == "ms-1":
                cube.units = "m s-1"
        # L.info(f"{cubes_to_regrid=}")
        # Regrid all cubes
        cl_proc = simple_regrid_lfric(
            cubes_to_regrid, tgt_cube=tgt_cube, ref_cube_constr=ref_cube
        )
        const = init_const(planet, directory=paths.const)
        add_planet_conf_to_cubes(cl_proc, const=const)

        # Write the data to a netCDF file
        days = 0 + get_cube_rel_days(cl_proc[0]).astype(int)
        day_str = f"days{days[0]}"
        if len(days) > 1:
            day_str += f"_{days[-1]}"
        fname_out = (
            outdir
            / f"{label}_{c_num}_{time_prof}_{chunk_label}_regr.nc".lower()
        )
        save_cubelist(cl_proc, fname_out, **gl_attrs)
        L.success(f"Saved to {fname_out}")
    L.info(f"Execution time: {time() - t0:.1f}s")


if __name__ == "__main__":
    main()
