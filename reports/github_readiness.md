# Final local GitHub publication audit

Audit baseline: commit `5b75de8` on `main`. This preparation changes repository metadata and documentation only; scientific source, configurations, locked dependencies, numerical data, and final figure exports retain their existing contents.

## Git and release contents

- The existing eight-commit history was inspected, along with tracked files, ignore rules, branch configuration, and Git object integrity (`git fsck --full`).
- Default/local branch is `main`; no remote is configured. No history rewriting was required or performed.
- README, banner, package, configurations, lockfile, tests, CI, CFF citation, model/reproduction documentation, scientific audit, compact numerical data, and main/alternative figures are tracked.
- Both stopping regimes contain 200 replicate rows and 40 parameter cells with five replicates each. All 200 trajectory NPZs validate. The saved 250-substitution principal metrics still match exactly.
- Publication-preparation changes are committed locally. Push commands appear only as instructions in `docs/publishing.md`; the agent did not execute them.

## Security and privacy review

The pre-existing reachable history contained 335 unique blobs (15,744,547 bytes in total). All were scanned for credential/token patterns, private keys, credential-bearing URLs, personal email addresses, and user-specific home paths. A further 2,571 ZIP members were inspected, including Word XML/properties and uncompressed NPZ members. Commit identities and messages were reviewed separately.

No credential or sensitive-path findings required removal or history rewriting. One email-shaped byte sequence occurred in a numerical floating-point array and was confirmed to be a binary false positive. No actual personal email address was found in file content. Commit identities use the existing synthetic `Codex <codex@localhost>` identity. Author names, institutional affiliations/addresses, and the Word document's co-author metadata are intentional scholarly provenance. No unrelated research data, notebooks, environments, keys, or temporary checkpoints are tracked.

The largest historical blob is 1,114,259 bytes (approximately 1.06 MiB); no blob exceeds 50 MiB. There are no oversized GitHub files or LFS requirements for the included material. The reviewed cleanup diff introduces no credentials. Pattern scanning and metadata review reduce exposure risk but do not constitute a guarantee that every possible secret format is detectable.

New ignore rules cover local audit directories, environment secrets, private-key extensions, editor files, and build products. Text is checked out with LF line endings, preserving scientific-source fingerprints across Git platform defaults. Existing binary attributes for NPZ, PNG, PDF, and DOCX remain intact.

## Verification actually executed

- `make install`, `make test`, `make smoke`, `make validate`, Ruff lint, and Ruff format checks pass in the publication repository. Fourteen tests pass; complete-grid validation reports 200/200 and zero principal endpoint difference.
- A temporary isolated checkout inside the ignored local audit directory was installed from `uv.lock` with Python 3.12.7. It passed all fourteen tests, lint/format checks, and a fresh two-substitution smoke run.
- The isolated checkout's `make reproduce` command completed with **all 200 trajectories cached**. It performed no full-grid simulation. Reproduction outputs were written only in that temporary checkout; the repository's valid scientific results and artwork were not regenerated.
- All 207 checksummed publication input/endpoint/trajectory files match the existing artifact manifest; the scientific-source fingerprint matches. Final figure hashes, 600-dpi metadata, vector exports, and README/document relative links pass the artifact checker.
- The SVG banner was rendered locally with `rsvg-convert` and visually inspected. The displayed figure previews were also inspected. No placeholder DOI, journal acceptance, preprint URL, or CI success badge is present.
- All three badge image URLs returned HTTP 200 SVG responses. The referenced checkout/setup-uv action tags and pinned uv release were verified with unauthenticated read-only HTTP requests.
- `cffconvert` 2.0.0 validates `CITATION.cff` against schema 1.2.0. Names, affiliations, manuscript title/date, and unpublished status are retained without inventing publication identifiers.
- GitHub Actions structure is valid and mirrors the locally exercised commands. Its token permissions are explicitly read-only and checkout credential persistence is disabled. No hosted Actions run is claimed before the user's first push.

## Boundaries and outstanding decisions

No push, GitHub repository creation/modification, authentication, upload, visibility change, remote addition, or original-source modification occurred. No scientific model change or result migration was made in this preparation.

The repository is technically prepared for the user's manual upload. Two substantive author decisions remain: licensing terms, and sign-off on the inherited scientific discrepancies already documented in `scientific_fidelity.md`. The existing rights notice and citation were not altered to imply permissions or publication status. The optional 20-seed expansion remains unexecuted and is not required for the completed five-seed release.

Manual creation/push instructions, a suggested About description, and repository topics are in `docs/publishing.md`.
