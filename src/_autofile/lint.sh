#!/usr/bin/env bash
pylint --rcfile=.pylintrc $(git ls-files 'autofile/*.py')
