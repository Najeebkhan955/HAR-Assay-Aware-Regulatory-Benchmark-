#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
url='https://github.com/Najeebkhan955/HAR-Research-.git'
command -v git >/dev/null || { echo 'Install Git first'; exit 1; }
[ -d .git ] || git init -b main
git add -A
[ -z "$(git status --porcelain)" ] || git commit -m 'Archive complete available HAR research workflow, data and results'
if git remote get-url origin >/dev/null 2>&1; then git remote set-url origin "$url"; else git remote add origin "$url"; fi
git push -u origin HEAD:main
