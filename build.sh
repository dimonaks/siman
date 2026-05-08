#!/usr/bin/env bash
set -e

rm -rf build dist siman.egg-info

python -m pip install --upgrade pip setuptools wheel build twine

python -m build

python -m twine check dist/*

python -m twine upload --repository siman dist/*
