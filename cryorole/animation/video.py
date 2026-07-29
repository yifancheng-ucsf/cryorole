"""FFmpeg encoding and FFprobe semantic validation for Phase 4 movies."""

from __future__ import annotations

import json
import math
import subprocess
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path


@dataclass(frozen=True)
class FFmpegResult:
    output_path: Path
    arguments: tuple[str, ...]
    return_code: int
    log_path: Path


@dataclass(frozen=True)
class MovieValidation:
    movie_path: Path
    codec: str
    pixel_format: str
    dimensions: tuple[int, int]
    fps: float
    frame_count: int
    duration_sec: float | None
    validation_method: str
    frame_count_tolerance: float
    probe_payload: dict[str, object]
    arguments: tuple[str, ...]
    return_code: int
    log_path: Path


class FFmpegExecutionError(RuntimeError):
    """FFmpeg process failure with stable return-code provenance."""

    def __init__(
        self,
        return_code: int,
        log_path: Path,
        *,
        arguments: list[str] | None = None,
    ) -> None:
        self.return_code = int(return_code)
        self.log_path = log_path
        self.arguments = tuple(arguments or ())
        super().__init__(
            f"FFmpeg encoding failed with return code {return_code}; see {log_path}"
        )


class FFprobeExecutionError(RuntimeError):
    """FFprobe process failure with stable return-code provenance."""

    def __init__(
        self,
        return_code: int,
        log_path: Path,
        *,
        arguments: list[str] | None = None,
    ) -> None:
        self.return_code = int(return_code)
        self.log_path = log_path
        self.arguments = tuple(arguments or ())
        super().__init__(
            f"FFprobe validation failed with return code {return_code}; see {log_path}"
        )


class FFmpegOutputError(RuntimeError):
    """FFmpeg returned zero but did not produce a usable movie artifact."""

    def __init__(self, message: str, *, arguments: list[str], log_path: Path) -> None:
        self.return_code = 0
        self.arguments = tuple(arguments)
        self.log_path = log_path
        super().__init__(message)


class MovieValidationError(RuntimeError):
    """FFprobe returned zero but the movie violates the requested contract."""

    def __init__(self, message: str, *, arguments: list[str], log_path: Path) -> None:
        self.return_code = 0
        self.arguments = tuple(arguments)
        self.log_path = log_path
        super().__init__(message)


def encode_mp4(
    *,
    ffmpeg_bin: str | Path,
    composite_dir: str | Path,
    output_path: str | Path,
    fps: float,
    crf: int,
    log_path: str | Path,
) -> FFmpegResult:
    """Encode contiguous composite PNGs as constant-rate H.264/yuv420p."""

    executable = _validated_executable(ffmpeg_bin, "FFmpeg")
    if not math.isfinite(fps) or fps <= 0:
        raise ValueError("FFmpeg fps must be finite and positive")
    if not isinstance(crf, int) or not 0 <= crf <= 51:
        raise ValueError("FFmpeg CRF must be an integer in [0, 51]")
    frames = Path(composite_dir)
    if not frames.is_dir():
        raise ValueError(f"Composite frame directory does not exist: {frames}")
    movie = Path(output_path)
    movie.parent.mkdir(parents=True, exist_ok=True)
    fps_text = f"{fps:g}"
    arguments = [
        str(executable),
        "-y",
        "-framerate",
        fps_text,
        "-i",
        str((frames / "frame_%06d.png").resolve()),
        "-c:v",
        "libx264",
        "-crf",
        str(crf),
        "-pix_fmt",
        "yuv420p",
        "-r",
        fps_text,
        "-movflags",
        "+faststart",
        str(movie.resolve()),
    ]
    result = subprocess.run(
        arguments,
        shell=False,
        capture_output=True,
        text=True,
        check=False,
    )
    log = _write_process_log(
        log_path,
        arguments=arguments,
        return_code=result.returncode,
        stdout=result.stdout,
        stderr=result.stderr,
    )
    if result.returncode != 0:
        raise FFmpegExecutionError(result.returncode, log, arguments=arguments)
    if not movie.is_file() or movie.stat().st_size <= 0:
        raise FFmpegOutputError(
            f"FFmpeg returned success but MP4 is missing or empty: {movie}",
            arguments=arguments,
            log_path=log,
        )
    return FFmpegResult(
        output_path=movie,
        arguments=tuple(arguments),
        return_code=result.returncode,
        log_path=log,
    )


def validate_mp4(
    *,
    ffprobe_bin: str | Path,
    movie_path: str | Path,
    expected_dimensions: tuple[int, int],
    expected_fps: float,
    expected_frame_count: int,
    log_path: str | Path,
) -> MovieValidation:
    """Validate codec, pixel format, dimensions, rate, and frame count."""

    executable = _validated_executable(ffprobe_bin, "FFprobe")
    movie = Path(movie_path)
    if not movie.is_file() or movie.stat().st_size <= 0:
        raise RuntimeError(f"MP4 is missing or empty: {movie}")
    arguments = [
        str(executable),
        "-v",
        "error",
        "-select_streams",
        "v:0",
        "-show_entries",
        (
            "stream=codec_name,pix_fmt,width,height,avg_frame_rate,r_frame_rate,"
            "nb_frames,duration:format=duration"
        ),
        "-of",
        "json",
        str(movie.resolve()),
    ]
    result = subprocess.run(
        arguments,
        shell=False,
        capture_output=True,
        text=True,
        check=False,
    )
    log = _write_process_log(
        log_path,
        arguments=arguments,
        return_code=result.returncode,
        stdout=result.stdout,
        stderr=result.stderr,
    )
    if result.returncode != 0:
        raise FFprobeExecutionError(result.returncode, log, arguments=arguments)
    try:
        payload = json.loads(result.stdout)
    except json.JSONDecodeError as exc:
        raise MovieValidationError(
            "FFprobe returned invalid JSON",
            arguments=arguments,
            log_path=log,
        ) from exc
    streams = payload.get("streams") if isinstance(payload, dict) else None
    if not isinstance(streams, list) or len(streams) != 1 or not isinstance(streams[0], dict):
        raise MovieValidationError(
            "FFprobe must report exactly one selected video stream",
            arguments=arguments,
            log_path=log,
        )
    stream = streams[0]
    codec = str(stream.get("codec_name", ""))
    pixel_format = str(stream.get("pix_fmt", ""))
    try:
        dimensions = (int(stream.get("width", 0)), int(stream.get("height", 0)))
    except (TypeError, ValueError) as exc:
        raise MovieValidationError(
            "FFprobe reported invalid video dimensions",
            arguments=arguments,
            log_path=log,
        ) from exc
    if codec != "h264":
        raise MovieValidationError(
            f"MP4 codec mismatch: expected h264, found {codec!r}",
            arguments=arguments,
            log_path=log,
        )
    if pixel_format != "yuv420p":
        raise MovieValidationError(
            f"MP4 pixel format mismatch: expected yuv420p, found {pixel_format!r}",
            arguments=arguments,
            log_path=log,
        )
    if dimensions != tuple(expected_dimensions):
        raise MovieValidationError(
            f"MP4 dimensions mismatch: {dimensions} != {tuple(expected_dimensions)}",
            arguments=arguments,
            log_path=log,
        )
    try:
        fps = _resolved_stream_rate(stream)
    except RuntimeError as exc:
        raise MovieValidationError(
            str(exc),
            arguments=arguments,
            log_path=log,
        ) from exc
    if not math.isclose(fps, float(expected_fps), rel_tol=1e-6, abs_tol=1e-6):
        raise MovieValidationError(
            f"MP4 frame rate mismatch: {fps} != {expected_fps}",
            arguments=arguments,
            log_path=log,
        )

    reported_count = _positive_int_or_none(stream.get("nb_frames"))
    duration = _duration_seconds(stream, payload)
    if reported_count is not None:
        frame_count = reported_count
        validation_method = "reported_frame_count"
        frame_count_tolerance = 0.0
        if frame_count != expected_frame_count:
            raise MovieValidationError(
                f"MP4 frame count mismatch: {frame_count} != {expected_frame_count}",
                arguments=arguments,
                log_path=log,
            )
    else:
        if duration is None:
            raise MovieValidationError(
                "FFprobe reported neither frame count nor usable duration",
                arguments=arguments,
                log_path=log,
            )
        estimated = duration * fps
        if abs(estimated - expected_frame_count) > 0.5:
            raise MovieValidationError(
                "MP4 duration-derived frame count mismatch: "
                f"{estimated} != {expected_frame_count}",
                arguments=arguments,
                log_path=log,
            )
        frame_count = int(round(estimated))
        validation_method = "duration_times_fps"
        frame_count_tolerance = 0.5
    return MovieValidation(
        movie_path=movie,
        codec=codec,
        pixel_format=pixel_format,
        dimensions=dimensions,
        fps=fps,
        frame_count=frame_count,
        duration_sec=duration,
        validation_method=validation_method,
        frame_count_tolerance=frame_count_tolerance,
        probe_payload=payload,
        arguments=tuple(arguments),
        return_code=result.returncode,
        log_path=log,
    )


def _validated_executable(value: str | Path, label: str) -> Path:
    path = Path(value)
    if not path.is_file():
        raise ValueError(f"{label} executable does not exist: {path}")
    return path


def _write_process_log(
    path: str | Path,
    *,
    arguments: list[str],
    return_code: int,
    stdout: str,
    stderr: str,
) -> Path:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(
            {
                "arguments": arguments,
                "return_code": return_code,
                "stdout": stdout,
                "stderr": stderr,
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    return output


def _parse_rate(value: object) -> float:
    try:
        rate = float(Fraction(str(value)))
    except (ValueError, ZeroDivisionError) as exc:
        raise RuntimeError(f"FFprobe reported invalid frame rate: {value!r}") from exc
    if not math.isfinite(rate) or rate <= 0:
        raise RuntimeError(f"FFprobe reported invalid frame rate: {value!r}")
    return rate


def _resolved_stream_rate(stream: dict[str, object]) -> float:
    errors = []
    for field in ("avg_frame_rate", "r_frame_rate"):
        try:
            return _parse_rate(stream.get(field))
        except RuntimeError as exc:
            errors.append(str(exc))
    raise RuntimeError("FFprobe reported no valid frame rate: " + "; ".join(errors))


def _positive_int_or_none(value: object) -> int | None:
    try:
        parsed = int(str(value))
    except (TypeError, ValueError):
        return None
    return parsed if parsed > 0 else None


def _duration_seconds(stream: dict[str, object], payload: dict[str, object]) -> float | None:
    candidates = [stream.get("duration")]
    format_payload = payload.get("format")
    if isinstance(format_payload, dict):
        candidates.append(format_payload.get("duration"))
    for value in candidates:
        try:
            duration = float(value)
        except (TypeError, ValueError):
            continue
        if math.isfinite(duration) and duration > 0:
            return duration
    return None
