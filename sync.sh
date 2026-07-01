#!/bin/bash
set -euo pipefail

BRANCH=$(git rev-parse --abbrev-ref HEAD)

# 1. Commit local changes if any
if ! git diff --quiet && git diff --cached --quiet; then
  echo "=== Committing local changes ==="
  git add -A && git commit -m "sync $(date +%Y%m%d-%H%M)"
elif ! git diff --cached --quiet; then
  echo "=== Staged changes already present ==="
else
  echo "=== No local changes to commit ==="
fi

# 2. Push to GitHub
echo ""
echo "=== Push to GitHub ($BRANCH) ==="
git push origin "$BRANCH"

# 3. Pull on mgt (stash local changes if needed)
echo ""
echo "=== Pull on mgt ==="
ssh mgt "cd ~/github/extrempy && git stash && git pull && git stash pop"

# 4. Pull on cmt.calc4
echo ""
echo "=== Pull on cmt.calc4 ==="
ssh cmt.calc4 "cd /home/deeph/work/zengqy/github/extrempy && git stash && git pull && git stash pop"

echo ""
echo "=== Done ==="
