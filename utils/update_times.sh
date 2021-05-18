#!/bin/sh

# This is a primitive tool that walks a directory tree recursively and sets the last modification to the current time.
# It can be useful to avoid some overflow deletion on temporary storages.
# Input argument: directory name

for i in "$@"; do
    find "$i" -type f -exec touch {} \;
done
