install:
	uv sync --locked --extra dev

test:
	uv run pytest -q

smoke:
	uv run python -m enzyme_specialisation simulate --config configs/smoke_test.yaml

reproduce:
	uv run python scripts/reproduce_paper.py

figures:
	uv run python -m enzyme_specialisation plot

validate:
	uv run python scripts/validate_outputs.py

extended:
	uv run python -m enzyme_specialisation simulate --config configs/extended.yaml
