# atmosphere mode
Atmosphere_Mode_Args = {

    "number_photon": 300000,
    "atm_type": 1,
    "solar_angle": 45.0,
    "solar_elevation": 0.0,
    "zenith_angle": [100.0, 120.0, 140.0, 150.0, 160.0],
    "azimuth_angle": [0.0, 180.0],
    "integration_mode": 2,
    "atmosphere_type": 2,   # Mid latitude summer
    "aerosol_type": 2,  # Continental average
    "AOT": 0.3,
    "cloud_type": 0,    # cloud free
    "surface_type": 2,  # 3-D canopy
    "calculation_mode": 3,  # BRF only
    "BRF_zenith_angles": [0.0, 10.0, 20.0, 30.0, 40.0, 50.0],
    "BRF_azimuth_angles": [0.0, 180],
    "tree_species": 1,
    "optical_parameters": [0.1, 0.05, 0.15, 0.1, 0.1, 0.05],
    "leaf_area_density": [0.5],
    "forest_floor_LAI": 0.5,
    "branch_area_density": [0.5],
    "sbar": [0.25],
}

