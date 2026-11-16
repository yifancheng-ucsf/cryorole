from cryorole.gui.tasks import LatestWinsRunner


def test_discard_suppresses_the_active_result(qtbot):
    results, errors = [], []
    runner = LatestWinsRunner(results.append, errors.append)
    runner.submit(lambda: "stale")
    runner.discard()
    qtbot.waitUntil(lambda: not runner.busy)
    assert results == errors == []

    runner.submit(lambda: "current")
    qtbot.waitUntil(lambda: not runner.busy)
    assert results == ["current"]
