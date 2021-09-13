#!/bin/bash -e

die () {
  echo -e "$*"
  exit 1
}

[[ $1 ]] || die "Usage: ${0} <version>"
[[ $1 = 1.[0-9].[0-9] ]] || die "Version has the wrong format"

cd "$(git rev-parse --show-toplevel)"
sed -i "s/VERSION 1.[0-9].[0-9]/VERSION $1/" CMakeLists.txt
sed -i "/PROJECT_NUMBER         =/s/=.*/= $1/" doc/doxygen/Doxyfile
sed -i "1s/Libecpint.*/Libecpint $1/" README.md

git add -u
git commit -m "Bump version to $1"
git tag v"$1"
