from __future__ import annotations

import json
import subprocess

import pytest

from cryorole.animation.video import (
    FFmpegExecutionError,
    FFprobeExecutionError,
    encode_mp4,
    validate_mp4,
)


def test_ffmpeg_uses_argument_list_and_required_h264_policy(tmp_path, monkeypatch):
    executable = tmp_path / "ffmpeg"
    executable.write_bytes(b"binary")
    frames = tmp_path / "frames"
    frames.mkdir()
    movie = tmp_path / "animation.mp4"
    seen = {}

    def fake_run(arguments, **kwargs):
        seen["arguments"] = arguments
        seen["kwargs"] = kwargs
        movie.write_bytes(b"mp4")
        return subprocess.CompletedProcess(arguments, 0, "stdout", "stderr")

    monkeypatch.setattr(subprocess, "run", fake_run)
    result = encode_mp4(
        ffmpeg_bin=executable,
        composite_dir=frames,
        output_path=movie,
        fps=10.0,
        crf=18,
        log_path=tmp_path / "ffmpeg.log",
    )
    args = seen["arguments"]
    assert isinstance(args, list)
    assert seen["kwargs"]["shell"] is False
    assert ["-framerate", "10", "-i", str((frames / "frame_%06d.png").resolve())] == args[2:6]
    assert ["-c:v", "libx264"] == args[6:8]
    assert ["-crf", "18"] == args[8:10]
    assert ["-pix_fmt", "yuv420p"] == args[10:12]
    assert ["-r", "10"] == args[12:14]
    assert ["-movflags", "+faststart"] == args[14:16]
    assert result.return_code == 0
    assert result.output_path == movie


def test_ffmpeg_failure_preserves_frames_and_reports_return_code(tmp_path, monkeypatch):
    executable = tmp_path / "ffmpeg"
    executable.write_bytes(b"binary")
    frames = tmp_path / "frames"
    frames.mkdir()
    source = frames / "frame_000000.png"
    source.write_bytes(b"source")

    def fake_run(arguments, **kwargs):
        return subprocess.CompletedProcess(arguments, 7, "", "encode failed")

    monkeypatch.setattr(subprocess, "run", fake_run)
    with pytest.raises(FFmpegExecutionError, match="return code 7"):
        encode_mp4(
            ffmpeg_bin=executable,
            composite_dir=frames,
            output_path=tmp_path / "animation.mp4",
            fps=10,
            crf=18,
            log_path=tmp_path / "ffmpeg.log",
        )
    assert source.read_bytes() == b"source"


def test_ffmpeg_zero_return_without_movie_fails(tmp_path, monkeypatch):
    executable = tmp_path / "ffmpeg"
    executable.write_bytes(b"binary")
    frames = tmp_path / "frames"
    frames.mkdir()
    monkeypatch.setattr(
        subprocess,
        "run",
        lambda arguments, **kwargs: subprocess.CompletedProcess(arguments, 0, "", ""),
    )
    with pytest.raises(RuntimeError, match="missing or empty"):
        encode_mp4(
            ffmpeg_bin=executable,
            composite_dir=frames,
            output_path=tmp_path / "animation.mp4",
            fps=10,
            crf=18,
            log_path=tmp_path / "ffmpeg.log",
        )


def test_ffprobe_validates_codec_format_dimensions_fps_and_frame_count(
    tmp_path,
    monkeypatch,
):
    executable = tmp_path / "ffprobe"
    executable.write_bytes(b"binary")
    movie = tmp_path / "animation.mp4"
    movie.write_bytes(b"mp4")
    payload = {
        "streams": [
            {
                "codec_name": "h264",
                "pix_fmt": "yuv420p",
                "width": 1920,
                "height": 1080,
                "avg_frame_rate": "10/1",
                "r_frame_rate": "10/1",
                "nb_frames": "3",
                "duration": "0.3",
            }
        ],
        "format": {"duration": "0.3"},
    }
    seen = {}

    def fake_run(arguments, **kwargs):
        seen["arguments"] = arguments
        seen["kwargs"] = kwargs
        return subprocess.CompletedProcess(arguments, 0, json.dumps(payload), "")

    monkeypatch.setattr(subprocess, "run", fake_run)
    result = validate_mp4(
        ffprobe_bin=executable,
        movie_path=movie,
        expected_dimensions=(1920, 1080),
        expected_fps=10,
        expected_frame_count=3,
        log_path=tmp_path / "ffprobe.log",
    )
    assert result.validation_method == "reported_frame_count"
    assert result.frame_count == 3
    assert result.codec == "h264"
    assert isinstance(seen["arguments"], list)
    assert seen["kwargs"]["shell"] is False


def test_ffprobe_rejects_missing_movie_and_semantic_mismatch(tmp_path, monkeypatch):
    executable = tmp_path / "ffprobe"
    executable.write_bytes(b"binary")
    with pytest.raises(RuntimeError, match="missing or empty"):
        validate_mp4(
            ffprobe_bin=executable,
            movie_path=tmp_path / "missing.mp4",
            expected_dimensions=(1920, 1080),
            expected_fps=10,
            expected_frame_count=3,
            log_path=tmp_path / "ffprobe.log",
        )


def test_ffprobe_nonzero_return_fails(tmp_path, monkeypatch):
    executable = tmp_path / "ffprobe"
    executable.write_bytes(b"binary")
    movie = tmp_path / "animation.mp4"
    movie.write_bytes(b"mp4")
    monkeypatch.setattr(
        subprocess,
        "run",
        lambda arguments, **kwargs: subprocess.CompletedProcess(
            arguments, 6, "", "probe failed"
        ),
    )
    with pytest.raises(FFprobeExecutionError, match="return code 6"):
        validate_mp4(
            ffprobe_bin=executable,
            movie_path=movie,
            expected_dimensions=(1920, 1080),
            expected_fps=10,
            expected_frame_count=3,
            log_path=tmp_path / "ffprobe.log",
        )

    movie = tmp_path / "animation.mp4"
    movie.write_bytes(b"mp4")
    payload = {
        "streams": [{
            "codec_name": "hevc",
            "pix_fmt": "yuv420p",
            "width": 1920,
            "height": 1080,
            "avg_frame_rate": "10/1",
            "nb_frames": "3",
        }],
        "format": {},
    }
    monkeypatch.setattr(
        subprocess,
        "run",
        lambda arguments, **kwargs: subprocess.CompletedProcess(
            arguments, 0, json.dumps(payload), ""
        ),
    )
    with pytest.raises(RuntimeError, match="codec"):
        validate_mp4(
            ffprobe_bin=executable,
            movie_path=movie,
            expected_dimensions=(1920, 1080),
            expected_fps=10,
            expected_frame_count=3,
            log_path=tmp_path / "ffprobe.log",
        )


@pytest.mark.parametrize(
    ("updates", "message"),
    [
        ({"pix_fmt": "yuv444p"}, "pixel format"),
        ({"width": 1280}, "dimensions"),
        ({"avg_frame_rate": "12/1"}, "frame rate"),
        ({"nb_frames": "4"}, "frame count"),
    ],
)
def test_ffprobe_rejects_each_video_contract_mismatch(
    tmp_path,
    monkeypatch,
    updates,
    message,
):
    executable = tmp_path / "ffprobe"
    executable.write_bytes(b"binary")
    movie = tmp_path / "animation.mp4"
    movie.write_bytes(b"mp4")
    stream = {
        "codec_name": "h264",
        "pix_fmt": "yuv420p",
        "width": 1920,
        "height": 1080,
        "avg_frame_rate": "10/1",
        "nb_frames": "3",
        "duration": "0.3",
    }
    stream.update(updates)
    payload = {"streams": [stream], "format": {"duration": "0.3"}}
    monkeypatch.setattr(
        subprocess,
        "run",
        lambda arguments, **kwargs: subprocess.CompletedProcess(
            arguments, 0, json.dumps(payload), ""
        ),
    )
    with pytest.raises(RuntimeError, match=message):
        validate_mp4(
            ffprobe_bin=executable,
            movie_path=movie,
            expected_dimensions=(1920, 1080),
            expected_fps=10,
            expected_frame_count=3,
            log_path=tmp_path / "ffprobe.log",
        )


def test_ffprobe_duration_fallback_has_explicit_tolerance(tmp_path, monkeypatch):
    executable = tmp_path / "ffprobe"
    executable.write_bytes(b"binary")
    movie = tmp_path / "animation.mp4"
    movie.write_bytes(b"mp4")
    payload = {
        "streams": [{
            "codec_name": "h264",
            "pix_fmt": "yuv420p",
            "width": 1920,
            "height": 1080,
            "avg_frame_rate": "10/1",
            "nb_frames": "N/A",
            "duration": "0.3",
        }],
        "format": {"duration": "0.3"},
    }
    monkeypatch.setattr(
        subprocess,
        "run",
        lambda arguments, **kwargs: subprocess.CompletedProcess(
            arguments, 0, json.dumps(payload), ""
        ),
    )
    result = validate_mp4(
        ffprobe_bin=executable,
        movie_path=movie,
        expected_dimensions=(1920, 1080),
        expected_fps=10,
        expected_frame_count=3,
        log_path=tmp_path / "ffprobe.log",
    )
    assert result.validation_method == "duration_times_fps"
    assert result.frame_count == 3


def test_ffprobe_duration_fallback_rejects_count_mismatch(tmp_path, monkeypatch):
    executable = tmp_path / "ffprobe"
    executable.write_bytes(b"binary")
    movie = tmp_path / "animation.mp4"
    movie.write_bytes(b"mp4")
    payload = {
        "streams": [{
            "codec_name": "h264",
            "pix_fmt": "yuv420p",
            "width": 1920,
            "height": 1080,
            "avg_frame_rate": "10/1",
            "nb_frames": "N/A",
            "duration": "0.5",
        }],
        "format": {"duration": "0.5"},
    }
    monkeypatch.setattr(
        subprocess,
        "run",
        lambda arguments, **kwargs: subprocess.CompletedProcess(
            arguments, 0, json.dumps(payload), ""
        ),
    )
    with pytest.raises(RuntimeError, match="duration-derived"):
        validate_mp4(
            ffprobe_bin=executable,
            movie_path=movie,
            expected_dimensions=(1920, 1080),
            expected_fps=10,
            expected_frame_count=3,
            log_path=tmp_path / "ffprobe.log",
        )
