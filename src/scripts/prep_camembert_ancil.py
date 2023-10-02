#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Standard library
from pathlib import Path
import sys
from typing import Tuple

# Third party
import ants
import numpy as np
import iris
from iris.coords import AuxCoord
from iris.util import promote_aux_coord_to_dim_coord  # , new_axis
import aeolus
from aeolus.io import create_dummy_cube, load_vert_lev
from aeolus.const import init_const
import pandas as pd

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


CAMEMBERT_GITHUB_RAW = (
    "https://raw.githubusercontent.com/projectcuisines/camembert/main/"
)


def dataframe_to_cubelist(df: pd.DataFrame) -> iris.cube.CubeList():
    """Convert dataframe to iris cubelist."""
    try:
        from iris.pandas import as_cubes

        cl = as_cubes(df, aux_coord_cols=["Height_m"])
    except ImportError:
        cl = iris.cube.CubeList()
        for name in df:
            if name.startswith("MMR_"):
                cube = iris.pandas.as_cube(df[name])
                cube.add_aux_coord(
                    iris.coords.AuxCoord(
                        df["Height_m"].values, long_name="Height_m"
                    ),
                    0,
                )
                cube.coord("index").rename("unknown")
                cube.rename(name)
                cl.append(cube)
    return cl


def load_init_cond(file_path: Path) -> pd.DataFrame:
    """Load the initial profiles into a pandas dataframe."""
    return pd.read_csv(file_path, delimiter=r"\s+", comment="#")


def prep_profiles(
    df: pd.DataFrame, const: aeolus.const.const.ConstContainer
) -> pd.DataFrame:
    """Calculate theta and height and append them to the dataframe."""
    # Discard pressures below the reference pressure
    p_zero = const.reference_surface_pressure.data
    df.loc[-1, "Pressure_Pa"] = p_zero
    df = df.sort_values("Pressure_Pa").reset_index(drop=True).interpolate()
    df = df[df["Pressure_Pa"] <= p_zero]
    df = df.reindex(index=df.index[::-1]).reset_index(drop=True)
    # Calculate additional variables: Exner function, theta, height
    kappa = (
        const.dry_air_gas_constant.data / const.dry_air_spec_heat_press.data
    )
    df["Exner"] = (df["Pressure_Pa"] / p_zero) ** kappa
    df["Theta_K"] = df["Temperature_K"] / df["Exner"]
    pres = df["Pressure_Pa"].values
    temp = df["Temperature_K"].values
    # Calculate height using the hydrostatic equation and variable gravity
    z = np.zeros(df.shape[0])
    for i in range(1, z.shape[0]):
        g_local = (
            const.gravity.data
            * (const.radius.data / (const.radius.data + z[i - 1])) ** 2
        )
        scale_height = const.dry_air_gas_constant.data * temp[i - 1] / g_local
        z[i] = z[i - 1] - scale_height * np.log(pres[i] / pres[i - 1])
    df["Height_m"] = z
    return df


def add_mass_mixing_ratios(
    df: pd.DataFrame, const: aeolus.const.const.ConstContainer
) -> pd.DataFrame:
    """Convert VMR to MMR and add them to the dataframe."""
    # data from Socrates: src/radiance_core/gas_list_pcf.F90
    MOLAR_WEIGHTS = {
        "H2O": 18.0153,
        "CO2": 44.01,
        "CO": 28.0106,
        "CH4": 16.043,
        "NH3": 17.0306,
        "H2": 2.01588,
        "He": 4.002602,
    }
    for gas, gas_mw in MOLAR_WEIGHTS.items():
        df[f"MMR_{gas}"] = (
            df[f"VMR_{gas}"]
            * gas_mw
            * 1e-3
            / const.dry_air_molecular_weight.data
        )
    return df


def make_vert_coord(vert_lev_file_path: Path) -> Tuple[AuxCoord]:
    """Create a pair of level height and level number coordinates."""
    # Load vertical levels file
    thlev = load_vert_lev(vert_lev_file_path)
    # Create iris objects
    level_height = AuxCoord(
        points=thlev,
        units="m",
        var_name="level_height",
        attributes={"positive": "up"},
    )
    model_level_number = iris.coords.DimCoord(
        points=np.arange(len(thlev), dtype=int),
        standard_name="model_level_number",
        var_name="model_level_number",
        attributes={"positive": "up"},
        units="1",
    )
    return level_height, model_level_number


def interpolate_mmr_in_height(
    cube_list: iris.cube.CubeList, target_theta_lev: np.array
) -> iris.cube.CubeList:
    """Interpolate MMR cubes in level height."""
    dset_out = iris.cube.CubeList()
    for cube in cube_list:
        if not cube.name().startswith("MMR"):
            # Skip non-MMR cubes
            continue

        cube.var_name = cube.name().lower().lstrip("mmr_")
        cube.long_name = f"{cube.var_name}_mmr"
        cube.units = "kg kg-1"
        cube.remove_coord("unknown")
        promote_aux_coord_to_dim_coord(cube, "Height_m")
        cube.coord("Height_m").units = "m"
        cube.coord("Height_m").rename("level_height")
        dset_out.append(
            cube.interpolate(
                [("level_height", target_theta_lev)], iris.analysis.Linear()
            )
        )
    return dset_out


def broadcast_to_3d(
    cube_list: iris.cube.CubeList,
    level_height: AuxCoord,
    model_level_number: AuxCoord,
) -> iris.cube.CubeList:
    """Broadcast vertical profiles to a horizontal UM grid."""
    # Create a dummy cube on a UM grid
    dummy_cube = create_dummy_cube(n_res=96, pm180=True)
    # dummy_cube.data[:, :] = 1.0
    # dummy_cube.add_aux_coord(level_height[0], data_dims=())
    # dummy_cube = new_axis(dummy_cube, "level_height")

    # # Multiply the dummy cube (1, NLAT, NLON) by the gas MMR cube (NLEV)
    # Commented out: it doesn't work with an old iris version required by ants
    # dset_out = iris.cube.CubeList()
    # for cube in cube_list:
    #     cube3d = dummy_cube * cube
    #     cube3d.add_aux_coord(model_level_number, data_dims=(0,))
    #     cube3d.var_name = cube.var_name
    #     cube3d.long_name = cube.long_name
    #     cube3d.units = cube.units
    #     dset_out.append(cube3d)
    dset_out = iris.cube.CubeList()
    for gas_mmr_cube in cube_list:
        _cl = iris.cube.CubeList()
        for gas_mmr_val, lh, mln in zip(
            gas_mmr_cube.data, level_height, model_level_number
        ):
            _cube = dummy_cube.copy(
                data=np.ones_like(dummy_cube.data) * gas_mmr_val
            )
            _cube.add_aux_coord(lh, data_dims=())
            _cube.add_aux_coord(mln, data_dims=())
            _cl.append(_cube)
        cube3d = _cl.merge_cube()
        cube3d.var_name = gas_mmr_cube.var_name
        cube3d.long_name = gas_mmr_cube.long_name
        cube3d.units = gas_mmr_cube.units
        dset_out.append(cube3d)
    return dset_out


def main(planet="gj1214b", c_num="C12"):
    """Main entry point."""
    # Load the raw data
    filename = f"CAMEMBERT_{planet[:-1].upper()}{planet[-1]}_IC.dat"
    df = load_init_cond(f"{CAMEMBERT_GITHUB_RAW}/InitialConditions/{filename}")

    # Planetary constants
    const = init_const(f"camembert_{planet}", directory=paths.const)

    # Calculate additional thermodynamic variables
    df = prep_profiles(df, const)
    print(f"theta={np.round(df.Theta_K.values, 2)}")
    print(f"height={np.round(df.Height_m.values, 2)}")

    # Convert VMR to MMR
    df = add_mass_mixing_ratios(df, const)
    max_mmr_sum = (
        df[
            [
                "MMR_H2O",
                "MMR_CO2",
                "MMR_CO",
                "MMR_CH4",
                "MMR_NH3",
                "MMR_H2",
                "MMR_He",
            ]
        ]
        .sum(axis=1)
        .max()
    )
    print(f"max total MMR: {max_mmr_sum}")

    # Create the target vertical grid
    level_height, model_level_number = make_vert_coord(
        paths.vert / "Lev_66_4.0e6_uniform"
    )
    nlev = level_height.shape[0]

    # Interpolate in height
    dset = interpolate_mmr_in_height(
        dataframe_to_cubelist(df),
        target_theta_lev=level_height.points,
    )

    dset = broadcast_to_3d(
        dset, level_height=level_height, model_level_number=model_level_number
    )

    # Add zero values
    all_gases = [
        "cfc11",
        "cfc113",
        "cfc12",
        "ch4",
        "co",
        "co2",
        "dms",
        "dmso",
        "h2",
        "h2o",
        "h2o2",
        "h2o2_limit",
        "h2so4",
        "hcfc22",
        "hcn",
        "he",
        "hfc134a",
        "ho2",
        "monoterpene",
        "n2",
        "n2o",
        "nh3",
        "no3",
        "o2",
        "o3",
        "oh",
        "secondary_organic",
        "so2",
    ]
    dset_final = iris.cube.CubeList()
    for gas_name in all_gases:
        try:
            dset_final.append(dset.extract_strict(f"{gas_name}_mmr"))
        except iris.exceptions.ConstraintMismatchError:
            print(gas_name)
            zero_cube = dset[0].copy() * 0.0
            zero_cube.var_name = gas_name
            zero_cube.long_name = f"{gas_name}_mmr"
            dset_final.append(zero_cube)
    # dset_final = dset

    # Define the output file names
    anc_um = (
        paths.ancil
        / "um"
        / f"camembert_case3_{planet}_gas_mmr_n96e_l{nlev-1}.nc"
    )
    anc_lfric = (
        paths.ancil
        / "lfric"
        / c_num
        / f"camembert_case3_{planet}_gas_mmr_{c_num}_l{nlev-1}.nc"
    )

    # Save the data on UM grid
    ants.io.save.netcdf(dset_final, str(anc_um))

    # Regrid to LFRic's cubed sphere grid
    tgt_grid = paths.ancil / "lfric" / c_num / "qrparm.orog.ugrid.nc"
    regrid(str(anc_um), str(anc_lfric), target_path=str(tgt_grid))
    print(f"Result: {anc_lfric}")


if __name__ == "__main__":
    main(c_num="C48")
