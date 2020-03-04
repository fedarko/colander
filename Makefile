# NOTE: based on https://github.com/fedarko/qeeseburger/blob/master/Makefile
.PHONY: test stylecheck style

test:
	python3 -B -m pytest colander/tests --cov colander

stylecheck:
	flake8 *
	black --check -l 79 *

style:
	black -l 79 *
