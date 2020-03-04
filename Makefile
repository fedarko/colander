# NOTE: based on https://github.com/fedarko/qeeseburger/blob/master/Makefile
.PHONY: test stylecheck style

test:
	python3 -B -m pytest colander/tests --cov colander

stylecheck:
	@# both of these ignored things are things that black and flake8 disagree on
	flake8 --ignore=E203,W503 colander
	black --check -l 79 colander

style:
	black -l 79 colander
