# Migration from cryoROLE 0.x

This page is for users familiar with the original cryoROLE script workflow.
New users should start with `README.md` or `docs/quick_start.md`.

## Command Mapping

| cryoROLE 0.x command | cryoROLE 2.0 command | Main change |
|---|---|---|
| `orientation_analysis` | `cryorole run` | Creates an audited run bundle instead of only flat output files. |
| `landscape_projection` | `cryorole visualize` | Visualization is display-only and stays inside the run bundle. |
| `point_select` | `cryorole select` | Selection is a recorded scientific decision artifact. |
| `particle_backtrack` | `cryorole export` | Export uses recorded source-row provenance from the run. |

## Relative Orientation Run

Old style:

```bash
orientation_analysis --s1 ref.star --s2 mov.star --o ref_against_mov
```

New style:

```bash
cryorole run --ref ref.star --mov mov.star
```

cryoROLE 2.0 records input normalization, identity matching, relative
orientation, SLD density, reports, and manifests in a run bundle.

## RND and SLD

The old `RND` field and the new `SLD` fields share the idea of local density in
relative-orientation space, but they should not be treated as identical columns.

cryoROLE 2.0 writes explicit SLD fields, including:

```text
sld_raw
sld_display
sld_unfloored
sld_was_floored
sld_local_k_mean
sld_effective_local_k_mean
```

`sld_raw` is the scientific density field used by default. `sld_display` is for
display policy and color scaling.

## Row Order and Matching

The 0.x scripts assumed that two STAR files had the same particle count and the
same particle order.

In cryoROLE 2.0, row-order matching must be explicit:

```bash
cryorole run --ref aligned_ref.star --mov aligned_mov.star --row-aligned
```

Without `--row-aligned`, cryoROLE uses auditable identity matching:

- CryoSPARC `.cs`: `uid`
- RELION STAR: `_rlnTomoParticleName`, then `_rlnImageName` / `rlnImageName`
  when safe

If safe default matching is not possible, use `cryorole align` or manually
prepare aligned metadata before `cryorole run --row-aligned`.

## Selection

Old `point_select` selected points inside a sphere in Euler coordinate space.

New `cryorole select` defaults to SO(3) geodesic radius selection when using
Euler centers:

```bash
cryorole select \
  --run-dir cryorole_outputs \
  --selection-id state_1 \
  --space canonical \
  -c 13 0 14 \
  -r 6
```

This avoids treating Euler coordinates as a simple Euclidean space for the
default scientific radius selection.

## Backtracking and Export

Old `particle_backtrack` used the CSV `ID` column to map selected rows back to a
STAR file.

cryoROLE 2.0 records source-row provenance during `run`:

```text
ref_source_row_id
mov_source_row_id
```

`cryorole export` uses this provenance to subset the original ref and/or mov
metadata. It does not reselect particles, rematch particles, or rewrite source
poses with display/canonical coordinates.

## Autoalign and Canonicalization

The old `--autoalign` / inverse-mean-rotation workflow and 2.0
`cryorole canonicalize` both address coordinate-frame interpretation, but they
are not documented as identical algorithms.

Use:

```bash
cryorole canonicalize --run-dir cryorole_outputs
```

Canonicalization writes a separate canonical landscape and canonical frame
artifacts without overwriting the raw landscape.

## CryoSPARC Inputs

Older workflows often converted CryoSPARC `.cs` files to STAR with pyem before
running cryoROLE.

cryoROLE 2.0 treats `.cs` as a native input format:

```bash
cryorole run --ref ref_domain.cs --mov mov_domain.cs
```

The default `.cs` identity key is `uid`, and the native pose field is
`alignments3D/pose`.
