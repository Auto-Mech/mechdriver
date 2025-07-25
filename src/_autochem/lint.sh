#!/usr/bin/env bash

uv run lefthook run pre-commit --all-files
uv run lefthook run pre-push --all-files
