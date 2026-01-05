#!/bin/bash --login
set -e

mamba activate $HOME/app/env
exec "$@"