from __future__ import annotations

import shutil

import pytest
from PIL import Image

from cryorole.animation.video import encode_mp4, validate_mp4


@pytest.mark.ffmpeg
def test_real_ffmpeg_h264_round_trip(tmp_path):
    ffmpeg = shutil.which("ffmpeg")
    ffprobe = shutil.which("ffprobe")
    if not ffmpeg or not ffprobe:
        pytest.skip("FFmpeg and FFprobe are required for real encoding validation")
    frames = tmp_path / "frames"
    frames.mkdir()
    for index, color in enumerate(("red", "green", "blue")):
        Image.new("RGB", (320, 180), color).save(
            frames / f"frame_{index:06d}.png"
        )
    movie = tmp_path / "animation.mp4"
    encode_mp4(
        ffmpeg_bin=ffmpeg,
        composite_dir=frames,
        output_path=movie,
        fps=10,
        crf=18,
        log_path=tmp_path / "ffmpeg.log",
    )
    validation = validate_mp4(
        ffprobe_bin=ffprobe,
        movie_path=movie,
        expected_dimensions=(320, 180),
        expected_fps=10,
        expected_frame_count=3,
        log_path=tmp_path / "ffprobe.log",
    )
    assert validation.codec == "h264"
    assert validation.pixel_format == "yuv420p"
