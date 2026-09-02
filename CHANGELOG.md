# Changelog

Entries begin at 0.2.0. For 0.1.0 see the GitHub release and its Zenodo
record; the concept DOI
[10.5281/zenodo.22135973](https://doi.org/10.5281/zenodo.22135973) always
resolves to the latest version.

## 0.2.0 — 2026-09-02

Ports the admissible reference class and the exact upper envelope in base R.
**No existing output changes** and no new dependency: `pbeta`,
`findInterval` and `seq` are all base R.

- New `R/reference.R`, exported through `NAMESPACE`: `admissibility`,
  `interior_q`, `prereq_index_family`, `pi_envelope`, `uniform_reference`,
  `beta_reference`, `point_mass_reference`, `attaining_reference`.
- Checked against the Python golden vectors at 1e-12 over the thirty new
  quantities: twenty-eight bit-identical, largest absolute difference
  1.110e-16. The two non-identical are `PI_hi` on the fixtures where `A1`
  was already non-identical.
- Test suite 155 passing, against 67 at 0.1.0.
- **The scope remains narrower than the Python and JavaScript packages**:
  the ceiling fit, the plausible-value handling and the study runner are
  still out. The README scope table records what is and is not covered.
