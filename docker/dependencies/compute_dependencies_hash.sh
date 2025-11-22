#!/bin/bash

# Exit the script at the first failure
set -e

find dependencies/kokkos dependencies/kokkos-kernels docker/dependencies -type f -exec sha1sum {} \; | sort | sha1sum | cut -c -8
