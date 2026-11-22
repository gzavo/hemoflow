#!/usr/bin/env bash
# Sets up palabos/: upstream Palabos at the pinned commit, trimmed to what HemoFlow builds,
# with the HemoFlow patches of patches/palabos/ applied in order (0001, 0002, ...).
#
#   scripts/setup_palabos.sh           create palabos/ (stops if it exists and differs)
#   scripts/setup_palabos.sh --force   replace an existing palabos/ (local changes there are lost)
#   scripts/setup_palabos.sh --check   compare palabos/ with what the script would create
#
# palabos/HEMOFLOW_PATCHES records the commit and the patches; CMake stops when they no longer
# match this script and patches/palabos/. Needs git and access to gitlab.com.
#
# Another Palabos version: change PALABOS_COMMIT, update or drop the patches that no longer
# apply, run with --force and repeat the regression tests (doc/dev/TODO.md, section 0).
# A new patch: patches/palabos/NNNN-name.patch, a diff relative to the Palabos root
# (paths a/src/...), with a description of the change and why it is needed above the diff.
set -euo pipefail

PALABOS_URL=https://gitlab.com/unigespc/palabos.git
PALABOS_COMMIT=4127697e90169bbef982295f1d1c933cf6e90caa   # upstream master, 3 April 2026 (after 2.3.0)
KEEP=(src externalLibraries dco COPYING README.md CHANGELOG.md)

root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
dest=$root/palabos
shopt -s nullglob
patches=("$root"/patches/palabos/*.patch)   # In name order; none once upstream has them all
shopt -u nullglob

die() { echo "setup_palabos.sh: $*" >&2; exit 1; }

mode=create
case "${1:-}" in
    "") ;;
    --force) mode=force ;;
    --check) mode=check ;;
    -h|--help) sed -n '2,15p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) die "unknown option '$1' (see --help)" ;;
esac
command -v git >/dev/null || die "git is needed"

sha256() {
    if command -v sha256sum >/dev/null; then sha256sum "$1"; else shasum -a 256 "$1"; fi | cut -d' ' -f1
}

# The record that CMakeLists.txt compares (same format there)
stamp() {
    echo "commit $PALABOS_COMMIT"
    for p in ${patches[@]+"${patches[@]}"}; do   # (empty-array form for bash 3)
        echo "patch $(sha256 "$p") $(basename "$p")"
    done
}

work=$(mktemp -d "${TMPDIR:-/tmp}/hemoflow-palabos.XXXXXX")
trap 'rm -rf "$work"' EXIT

echo "Fetching Palabos $PALABOS_COMMIT from $PALABOS_URL"
git init -q "$work/upstream"
git -C "$work/upstream" fetch -q --depth 1 "$PALABOS_URL" "$PALABOS_COMMIT" ||
    die "could not fetch commit $PALABOS_COMMIT"
[ "$(git -C "$work/upstream" rev-parse FETCH_HEAD)" = "$PALABOS_COMMIT" ] || die "fetched the wrong commit"
git -C "$work/upstream" -c advice.detachedHead=false checkout -q FETCH_HEAD

for p in ${patches[@]+"${patches[@]}"}; do   # (empty-array form for bash 3)
    echo "Applying $(basename "$p")"
    git -C "$work/upstream" apply "$p" || die "$(basename "$p") does not apply to $PALABOS_COMMIT"
done

mkdir "$work/palabos"
for k in "${KEEP[@]}"; do
    cp -R "$work/upstream/$k" "$work/palabos/"
done
stamp > "$work/palabos/HEMOFLOW_PATCHES"

if [ "$mode" = check ]; then
    [ -d "$dest" ] || die "palabos/ does not exist; run scripts/setup_palabos.sh"
    if diff -r -q "$work/palabos" "$dest" >"$work/differences"; then
        echo "palabos/ is Palabos $PALABOS_COMMIT with the ${#patches[@]} patches of patches/palabos/."
        exit 0
    fi
    sed -e "s|$work/palabos|(pinned + patches)|g" -e "s|$dest|palabos|g" "$work/differences"
    echo "palabos/ differs from Palabos $PALABOS_COMMIT with the patches of patches/palabos/." >&2
    exit 1
fi

if [ -e "$dest" ]; then
    if diff -r -q "$work/palabos" "$dest" >/dev/null; then
        echo "palabos/ is up to date."
        exit 0
    fi
    [ "$mode" = force ] ||
        die "palabos/ exists and differs from what this script creates (--check lists the differences). --force replaces it; local changes in palabos/ are then lost."
    rm -rf "$dest"
fi
mv "$work/palabos" "$dest"
echo "palabos/ is now Palabos $PALABOS_COMMIT with the ${#patches[@]} patches of patches/palabos/."
