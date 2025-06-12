#!/usr/bin/env bash
pylint --rcfile=.pylintrc $(git ls-files 'ratefit/*.py')
pylint --rcfile=.pylintrc $(git ls-files 'thermfit/*.py')
pylint --rcfile=.pylintrc $(git ls-files 'mechanalyzer/*.py')
