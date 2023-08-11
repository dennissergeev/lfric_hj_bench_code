#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Standard library
from pathlib import Path
import sys

# Third party
import ants
import numpy as np
import iris
from aeolus.io import create_dummy_cube, load_vert_lev

# The ancil regridding script in a separate repo
sys.path.append(
    str(
        Path.home()
        / "sandbox"
        / "ancil_contrib_trunk"
        / "Models"
        / "LFRic"
        / "ProtoGAL"
        / "bin"
    )
)
from ancil_regrid_um_to_mesh import main as regrid  # noqa

# Local
import paths  # noqa

# Create a dummy cube on a UM grid
dummy_cube = create_dummy_cube(n_res=96, pm180=True)
# Load vertical levels file
thlev = load_vert_lev(paths.vert / "vertlevs_L38_29t_9s_40km")
level_height = iris.coords.AuxCoord(
    points=thlev, units="m", var_name="level_height", attributes={"positive": "up"}
)
model_level_number = iris.coords.DimCoord(
    points=np.arange(len(thlev), dtype=int),
    standard_name="model_level_number",
    var_name="model_level_number",
    attributes={"positive": "up"},
    units="1",
)
time_coord = iris.coords.AuxCoord(
    points=[0],
    var_name="time",
    units="hours since 1970-01-01 00:00:00",
)

# Define gas MMRs
gases = {
    "h2": {
        "data": 0.05,
        "long_name": "mass_mixing_ratio_of_molecular_hydrogen",
        "units": "kg kg-1",
    },
    "nh3": {
        "data": 0.026,
        "long_name": "mass_mixing_ratio_of_ammonia",
        "units": "kg kg-1",
    },
    "co2": {
        "data": 0.95,
        "long_name": "mass_mixing_ratio_of_carbon_dioxide",
        "units": "kg kg-1",
    },
    "n2": {
        "data": 0.026,
        "long_name": "mass_mixing_ratio_of_nitrogen",
        "units": "kg kg-1",
    },
    "so2": {
        "data": 0.024,
        "long_name": "mass_mixing_ratio_of_sulphur_dioxide",
        "units": "kg kg-1",
    },
}


# Assemble the data set
dset_final = iris.cube.CubeList()
for gas_var_name, gas_prop in gases.items():
    cl = iris.cube.CubeList()
    for lh, mln in zip(level_height, model_level_number):
        cube = dummy_cube.copy(data=np.ones_like(dummy_cube.data) * gas_prop["data"])
        cube.add_aux_coord(lh, data_dims=())
        cube.add_aux_coord(mln, data_dims=())
        # cube.add_aux_coord(gc, data_dims=())
        cl.append(cube)
    cube3d = cl.merge_cube()
    cube3d.var_name = gas_var_name
    cube3d.long_name = gas_prop["long_name"]
    cube3d.units = gas_prop["units"]
    dset_final.append(cube3d)
# cube.attributes = {'STASH': iris.fileformats.pp.STASH.from_msi("m01s99i999")}
# cube.add_aux_coord(time_coord, data_dims=())
# lsm.attributes["grid_staggering"] = 6

# Define the output file names
anc_um = paths.ancil / f"gas_mmr_{'_'.join(gases.keys())}_n96e_l38.nc"
anc_lfric = paths.ancil / f"gas_mmr_{'_'.join(gases.keys())}_C12_l38.nc"

# Save the data on UM grid
ants.io.save.netcdf(dset_final, str(anc_um))

# Regrid to LFRic's cubed sphere grid
# tgt_grid = (
#     paths.mo_ancil
#     / "basic-gal/Puku/C12/n96e/orography/gmted_ramp2/qrparm.orog.ugrid.nc"
# )
tgt_grid = paths.ancil / "lfric" / "qrparm.orog.ugrid.nc"
regrid(str(anc_um), str(anc_lfric), target_path=str(tgt_grid))
