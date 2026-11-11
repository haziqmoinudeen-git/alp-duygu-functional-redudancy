# Publishing this repository

The repository is prepared locally on `main`. No remote is configured by this preparation, and nothing is pushed or uploaded. The scientific model and completed simulation outputs are unchanged.

## Decisions for the authors

- Select and approve code/data licensing terms. The existing `LICENSE` grants no additional permissions; no open-source licence has been inferred.
- Review the inherited discrepancies in the [scientific fidelity audit](../reports/scientific_fidelity.md), especially support-expansion hazard caching and the Methods wording. This preparation does not resolve those scientific decisions.
- Add a manuscript or preprint URL/DOI only when one actually exists. The current unpublished citation is valid without one.

## Suggested GitHub About text

**Repository name:** `enzyme-specialisation`

**Description:** Reproducible origin–fixation simulations of enzyme specialisation: functional redundancy, loop-local mutation, and evolutionary exposure.

**Topics:** `enzyme-specialisation`, `evolutionary-biology`, `computational-biology`, `enzyme-promiscuity`, `functional-redundancy`, `origin-fixation`, `stochastic-simulation`, `reproducible-research`, `python`

Leave the website field empty until a real paper, preprint, or project URL is available. Add a repository-specific CI badge only after the workflow has actually run on GitHub.

## Create the empty GitHub repository yourself

1. Open [GitHub's new-repository page](https://github.com/new) in your browser and choose your account or organisation as owner.
2. Set the repository name to `enzyme-specialisation` and paste the description above.
3. Choose visibility after resolving the author decisions above.
4. Leave **README**, **.gitignore**, and **licence** initialisation disabled. These files are already committed locally; GitHub should create an empty repository.
5. Create the repository. Then run the following yourself in Bash. The owner prompt accepts either your username or organisation name.

```bash
cd ~/phlora/enzyme-specialisation
read -r -p "GitHub username or organisation: " github_owner
git remote add origin "https://github.com/${github_owner}/enzyme-specialisation.git"
git push --set-upstream origin main
```

The preparation agent has not executed these commands. Authentication, if required when you push, is handled by you through your usual Git setup. No credentials belong in the remote URL or repository.

## After your first push

- Confirm `main` is the default branch and inspect the README on GitHub.
- Paste the topics above into the repository's About settings.
- Check the first `scientific-tests` Actions run. Local equivalent checks pass; a remote CI result cannot be claimed before that run.
- Citation metadata, figures, compact data, and source provenance are already included. No Git LFS setup is required for the current file sizes.

## Local checks

```bash
make install
make test
make smoke
make validate
uv run ruff check src tests scripts
uv run ruff format --check src tests scripts
uv run python scripts/check_release_artifacts.py
```

`make reproduce` is the full documented reproduction command. With the supplied, validated trajectories it skips completed simulations. It was also verified in an isolated local checkout during publication preparation. No full-grid simulation was restarted for that check.
