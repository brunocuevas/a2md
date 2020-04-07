SET_500_MP2_ENERGY_DATASET = dict(
    idx_list_file=LIBRARY_PATH + '/params/set500.json',
    mp2_energy=LIBRARY_PATH + '/params/set500_mp2.json',
    files_path='F:/aamdlib/a2mdnet/'
)

COMPLETE_SET_ENERGY_DATASET = dict(
    idx_list_file=LIBRARY_PATH + '/params/complete_set.json',
    mp2_energy=LIBRARY_PATH + '/params/complete_set_mp2.json',
    files_path='F:/aamdlib/optimized_mol2/'
)

QM9_DATASET = dict(
    idx_list_file=LIBRARY_PATH + r'/params/qm_subset_%d.json',
    u0=LIBRARY_PATH + r'/params/qm9_energy.json',
    files_path=r'F:/aamdlib/qm9/'
)

COMPLETE_SET_DENSITY_DATASET = dict(
    geometry_path=r'F:/aamdlib/complete_set/mol2/',
    density_path=r'F:/aamdlib/complete_set/ppp/sym/',
    validation_idx=LIBRARY_PATH + r'/params/complete_set_density_validation.json',
    curated_validation_idx=LIBRARY_PATH + r'/params/complete_set_density_validation_curated.json',
    training_idx=LIBRARY_PATH + r'/params/complete_set_density_training.json',
    curated_training_idx=LIBRARY_PATH + r'/params/complete_set_density_training_curated.json',
    all_idx=LIBRARY_PATH + r'/params/complete_set_density_all.json',
    test_idx=LIBRARY_PATH + r'/params/complete_set_density_test.json',
)

ANISET_DATASET = dict(
    geometry_path=r'F:/aniset/mol2/',
    density_path=r'F:/aniset/ppp/',
    density_values_path=r'F:/aniset/npy/',
    validation_idx=LIBRARY_PATH + r'/params/aniset_validation.json',
    training_idx=LIBRARY_PATH + r'/params/aniset_training.json',
    all_idx=LIBRARY_PATH + r'/params/aniset_all.json',
    test_idx=LIBRARY_PATH + r'/params/aniset_test.json',
    energy=LIBRARY_PATH + r'/params/aniset_all_energies.json'
)

FDASET = dict(
    geometry_path=r'F:/fdaset/mol2/',
    density_path=r'F:/fdaset/ppp/',
    density_values_path=r'F:/fdaset/npy/',
    all_idx=LIBRARY_PATH + r'/params/fdaset.json',
    energy=LIBRARY_PATH + r'/params/fdaset_all_energies.json'
)

VALANISET = dict(
    geometry_path=r'F:/valaniset/mol2/',
    density_path=r'F:/valaniset/ppp/',
    density_values_path=r'F:/valaniset/npy/',
    all_idx=LIBRARY_PATH + r'/params/valaniset_all.json',
    energy=LIBRARY_PATH + r'/params/valaniset_all_energies.json'
)