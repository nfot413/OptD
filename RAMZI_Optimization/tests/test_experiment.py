from dataclasses import replace

from src.config import DEFAULT_CONFIG
from src.experiment import resolve_project_path, run_experiment


def test_resolve_project_path_uses_subproject_root():
    path = resolve_project_path("data/fsr1.xlsx")
    assert path.name == "fsr1.xlsx"
    assert path.parent.name == "data"


def test_run_experiment_with_small_optimizer(tmp_path):
    config = replace(
        DEFAULT_CONFIG,
        frequency=replace(DEFAULT_CONFIG.frequency, num_points=31),
        optimizer=replace(DEFAULT_CONFIG.optimizer, max_iterations=1, population_size=2),
        output=replace(DEFAULT_CONFIG.output, results_dir=str(tmp_path)),
    )

    result = run_experiment(config, log=lambda _: None)

    assert result.spectrum_csv.exists()
    assert result.spectrum_plot.exists()
    assert result.best_loss >= 0
    assert result.best_params.shape == (3,)
