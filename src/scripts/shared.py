#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Objects used by many scripts in the project."""
# External modules
from aeolus.model import lfric, um

# Local modules
import paths


MODELS = {
    "um": {
        "model": um,
        "data_proc_path": paths.results_proc_um,
        "data_raw_path": paths.results_raw_um,
        "title": "UM",
        "kw_plt": {
            "linestyle": "--",
            "linewidth": 0.75,
            "dash_capstyle": "round",
        },
    },
    "lfric": {
        "model": lfric,
        "data_proc_path": paths.results_proc_lfric,
        "data_raw_path": paths.results_raw_lfric,
        "title": "LFRic-Atmosphere",
        "kw_plt": {"linestyle": "-", "linewidth": 1.25},
    },
}


TF_CASES = {
    "shj": {
        "title": "Shallow Hot Jupiter",
        "short_title": "SHJ",
        "planet": "shj",
        "kw_plt": {"color": "C0"},
        "timestep": 1200,
        "time_mean_period": 1000,
        "proc_fname_suffix": "sigma_p",
    },
    # "dhj": {
    #     "title": "Deep Hot Jupiter",
    #     "short_title": "DHJ",
    #     "planet": "dhj",
    #     "kw_plt": {"color": "C0"},
    #     "timestep": 1200,
    #     "time_mean_period": 1000,
    #     "proc_fname_suffix": "sigma_p",
    # },
}