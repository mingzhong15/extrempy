#!/bin/bash
set -euo pipefail

BRANCH=$(git rev-parse --abbrev-ref HEAD)

echo "=== Push to GitHub ($BRANCH) ==="
git push origin "$BRANCH"

echo ""
echo "=== Pull on mgt ==="
ssh mgt "cd ~/github/extrempy && git pull"

echo ""
echo "=== Done ==="
