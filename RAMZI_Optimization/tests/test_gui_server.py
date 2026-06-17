from src.gui_server import _optimize_payload, _simulate_payload


def test_gui_simulate_payload_returns_fixed_ramzi_spectrum():
    result = _simulate_payload({"frequency": {"num_points": 51}})

    assert len(result["frequency_hz"]) == 51
    assert len(result["wavelength_nm"]) == 51
    assert len(result["target_db"]) == 51
    assert len(result["simulated_db"]) == 51
    assert result["mse"] >= 0
    assert "fai1" in result["parameters"]


def test_gui_square_target_width_changes_target():
    narrow = _simulate_payload(
        {
            "frequency": {"num_points": 101},
            "target": {"passband_width_hz": 4e9},
        }
    )
    wide = _simulate_payload(
        {
            "frequency": {"num_points": 101},
            "target": {"passband_width_hz": 12e9},
        }
    )

    assert narrow["target_db"] != wide["target_db"]


def test_gui_optimize_payload_returns_seven_parameters():
    result = _optimize_payload(
        {
            "frequency": {"num_points": 41},
            "optimizer": {"max_iterations": 1, "population_size": 2, "polish_result": False},
        }
    )

    assert len(result["variables"]) == 7
    assert len(result["parameters"]) == 7
    assert result["best_loss"] >= 0
