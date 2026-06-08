#!/bin/bash
# =============================================================================
# Install THIS package's declared dependencies from its DESCRIPTION.
# Invoked by the SessionStart hook in .claude/settings.json.
# =============================================================================
set -euo pipefail

# CLAUDE_CODE_REMOTE is "true" only inside a cloud session.
if [ "${CLAUDE_CODE_REMOTE:-}" != "true" ]; then
  exit 0
fi

cd "$CLAUDE_PROJECT_DIR"

# Reads DESCRIPTION (Imports/Suggests/Remotes) and installs everything
Rscript -e 'pak::local_install_deps(dependencies = TRUE)'

exit 0