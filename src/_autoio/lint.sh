#!/usr/bin/env bash
pylint --rcfile=.pylintrc $(git ls-files '*.py')
