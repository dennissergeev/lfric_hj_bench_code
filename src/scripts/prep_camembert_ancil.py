#!/usr/bin/env python
"""Prepare CAMEMBERT input data for the UM and LFRic."""
# Standard library
from pathlib import Path
import sys
from warnings import filterwarnings

# Third party
import aeolus
from aeolus.const import init_const
from aeolus.io import create_dummy_cube, load_vert_lev
import click
import iris
from iris.coords import AuxCoord
from iris.util import promote_aux_coord_to_dim_coord
import numpy as np
import pandas as pd

filterwarnings("ignore", category=UserWarning)  # noqa

# Third party: ants and the ancil regridding script in a separate repo
import ants  # noqa: E402

sys.path.append(
    str(
        Path.home()
        / "sandbox"
        / "contrib_2.0.0"
        / "Models"
        / "LFRic"
        / "ProtoGAL"
        / "bin"
    )
)
# Local
from ancil_regrid_um_to_mesh import main as regrid  # noqa: E402
import paths  # noqa: E402

CAMEMBERT_GITHUB_RAW = (
    "https://raw.githubusercontent.com/projectcuisines/camembert/main/"
)


def load_socrates_gases(file_path: Path, remove_air=True) -> dict:
    """Load a dictionary of gas absorbers available in SOCRATES."""
    import yaml  # noqa

    with file_path.open("r") as f:
        s_gases = yaml.safe_load(f)
        if remove_air:
            _ = s_gases.pop("air")
    return dict(sorted(s_gases.items(), key=lambda x: x[1]["socrates_ip"]))


def dataframe_to_cubelist(df: pd.DataFrame) -> iris.cube.CubeList():
    """Convert dataframe to iris cubelist."""
    try:
        from iris.pandas import as_cubes

        cl = as_cubes(df, aux_coord_cols=["height_m"])
    except ImportError:
        cl = iris.cube.CubeList()
        for name in df:
            if name.startswith("mmr_"):
                cube = iris.pandas.as_cube(df[name])
                cube.add_aux_coord(
                    iris.coords.AuxCoord(
                        df["height_m"].values, long_name="height_m"
                    ),
                    0,
                )
                cube.coord("index").rename("unknown")
                cube.rename(name)
                cl.append(cube)
    return cl


def load_init_cond(file_path: Path) -> pd.DataFrame:
    """Load the initial profiles into a pandas dataframe."""
    df = pd.read_csv(file_path, delimiter=r"\s+", comment="#")
    df = df.rename(columns={i: i.lower() for i in df.columns})
    return df


def prep_profiles(
    df: pd.DataFrame, const: aeolus.const.const.ConstContainer
) -> pd.DataFrame:
    """Calculate theta and height and append them to the dataframe."""
    # Discard pressures below the reference pressure
    p_zero = const.reference_surface_pressure.data
    df.loc[-1, "pressure_pa"] = p_zero
    df = df.sort_values("pressure_pa").reset_index(drop=True).interpolate()
    df = df[df["pressure_pa"] <= p_zero]
    df = df.reindex(index=df.index[::-1]).reset_index(drop=True)
    # Calculate additional variables: exner function, theta, height
    kappa = (
        const.dry_air_gas_constant.data / const.dry_air_spec_heat_press.data
    )
    df["exner"] = (df["pressure_pa"] / p_zero) ** kappa
    df["theta_k"] = df["temperature_k"] / df["exner"]
    pres = df["pressure_pa"].values
    temp = df["temperature_k"].values
    # Calculate height using the hydrostatic equation and variable gravity
    z = np.zeros(df.shape[0])
    for i in range(1, z.shape[0]):
        g_local = (
            const.gravity.data
            * (const.radius.data / (const.radius.data + z[i - 1])) ** 2
        )
        scale_height = const.dry_air_gas_constant.data * temp[i - 1] / g_local
        z[i] = z[i - 1] - scale_height * np.log(pres[i] / pres[i - 1])
    df["height_m"] = z
    return df


def add_mass_mixing_ratios(
    df: pd.DataFrame, gas_dict: dict, const: aeolus.const.const.ConstContainer
) -> pd.DataFrame:
    """Convert VMR to MMR and add them to the dataframe."""
    for gas, gas_prop in gas_dict.items():
        if f"vmr_{gas}" in df.columns:
            df[f"mmr_{gas}"] = (
                df[f"vmr_{gas}"]
                * gas_prop["molar_weight"]
                * 1e-3
                / const.dry_air_molecular_weight.data
            )
    return df


def make_vert_coord(vert_lev_file_path: Path) -> tuple[AuxCoord]:
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
        if not cube.name().startswith("mmr_"):
            # Skip non-MMR cubes
            continue

        cube.var_name = cube.name().replace("mmr_", "")
        cube.long_name = f"mass_fraction_of_{cube.var_name}_in_air"
        cube.units = "kg kg-1"
        cube.remove_coord("unknown")
        promote_aux_coord_to_dim_coord(cube, "height_m")
        cube.coord("height_m").units = "m"
        cube.coord("height_m").rename("level_height")
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
    n_res: int,
) -> iris.cube.CubeList:
    """Broadcast vertical profiles to a horizontal UM grid."""
    # Create a dummy cube on a UM grid
    dummy_cube = create_dummy_cube(n_res=n_res, pm180=True)
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


@click.command()
@click.argument(
    "planet",
    type=click.Choice(["gj1214b", "k2_18b"]),
)
@click.argument(
    "c_num",
    type=click.Choice([f"C{_c}" for _c in [6, 12, 24, 36, 48, 64]]),
)
@click.option(
    "-l",
    "--levels",
    default="Lev_66_4.0e6_uniform",
    type=click.STRING,
)
@click.option(
    "-z",
    "--add_missing_as_zeros",
    is_flag=True,
    help="Add missing gases with zero values.",
)
@click.option(
    "-n",
    "--n_res",
    default=96,
    type=click.INT,
)
@click.option(
    "-v",
    "--verbose",
    is_flag=True,
    help="Show diagnostic messages.",
)
def main(
    planet: str,
    c_num: str,
    levels: str,
    add_missing_as_zeros: bool,
    n_res: int,
    verbose: bool,
) -> None:
    """
    Prepare CAMEMBERT input data for the UM and LFRic.

    1. Load the data from the GitHub repository

    2. Calculate potential temperature and exner function

    3. Integrate the hydrostatic balance equation to get height

    4. Convert VMR to MMR, opionally set missing gases to 0

    5. Interpolate in the vertical to the specified vertical levels set

    6. Broadcast in the horizontal to the UM grid and save to ancillary file

    7. Regrid from the UM's lat-lon grid to the LFRic's cubed-sphere grid
    """
    # Load the raw data
    filename = (
        f"CAMEMBERT_{planet.replace('_', '-')[:-1].upper()}{planet[-1]}_IC.dat"
    )
    df = load_init_cond(f"{CAMEMBERT_GITHUB_RAW}/InitialConditions/{filename}")

    # Planetary constants
    const = init_const(f"camembert_{planet}", directory=paths.const)

    # Calculate additional thermodynamic variables
    df = prep_profiles(df, const)
    click.secho(f"\nsize={df.shape}", fg="blue")
    click.secho(
        f"\ntheta={','.join(np.round(df.theta_k.values, 2).astype(str))}",
        fg="blue",
    )
    click.secho(
        f"\nheight={','.join(np.round(df.height_m.values, 2).astype(str))}",
        fg="blue",
    )

    # Convert VMR to MMR
    s_gases = load_socrates_gases(paths.scripts / "gases.yaml")
    df = add_mass_mixing_ratios(df, s_gases, const)
    max_mmr_sum = (
        df[[i for i in df.columns if i.startswith("mmr_")]].sum(axis=1).max()
    )
    click.secho(f"\nmax total MMR: {max_mmr_sum}", fg="blue")

    # Create the target vertical grid
    level_height, model_level_number = make_vert_coord(paths.vert / levels)
    nlev = level_height.shape[0]

    # Interpolate in height
    dset = interpolate_mmr_in_height(
        dataframe_to_cubelist(df),
        target_theta_lev=level_height.points,
    )

    # Broadcast horizontally
    dset = broadcast_to_3d(
        dset,
        level_height=level_height,
        model_level_number=model_level_number,
        n_res=n_res,
    )

    # Add zero values
    if add_missing_as_zeros:
        note = "all_gases"
        all_gases = {**s_gases}
        dset_final = iris.cube.CubeList()
        for gas_name in all_gases:
            try:
                dset_final.append(dset.extract_strict(f"{gas_name}_mmr"))
            except iris.exceptions.ConstraintMismatchError:
                if verbose:
                    click.secho(f"{gas_name:<10} MMR = 0.0", fg="yellow")
                zero_cube = dset[0].copy() * 0.0
                zero_cube.var_name = gas_name
                zero_cube.long_name = f"{gas_name}_mmr"
                dset_final.append(zero_cube)
    else:
        note = "nonzero_only"
        dset_final = dset

    # Define the output file names
    anc_um = (
        paths.ancil
        / "um"
        / f"camembert_case3_{planet}_gas_mmr_n{n_res}e_l{nlev-1}_{note}.nc"
    )
    anc_lfric = (
        paths.ancil
        / "lfric"
        / c_num
        / f"camembert_case3_{planet}_gas_mmr_{c_num}_l{nlev-1}_{note}.nc"
    )

    # Save the data on UM grid
    ants.io.save.netcdf(dset_final, str(anc_um))
    click.secho(f"\nUM ancil: {anc_um}", fg="green")

    # Regrid to LFRic's cubed sphere grid
    tgt_grid = paths.ancil / "lfric" / c_num / "qrparm.orog.ugrid.nc"
    regrid(str(anc_um), str(anc_lfric), target_path=str(tgt_grid))
    click.secho(f"\nLFRic ancil: {anc_lfric}", fg="green")


if __name__ == "__main__":
    main()
