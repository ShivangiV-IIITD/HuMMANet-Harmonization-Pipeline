#!/usr/bin/env bash
set -euo pipefail

PKG_ROOT="$SRC_DIR"
INSTALL_ROOT="$PREFIX/share/humanet-final-pipeline"

mkdir -p "$INSTALL_ROOT"
mkdir -p "$PREFIX/bin"

cp "$PKG_ROOT"/*.py "$INSTALL_ROOT"/
cp "$PKG_ROOT"/*.R "$INSTALL_ROOT"/
cp "$PKG_ROOT"/*.Rmd "$INSTALL_ROOT"/
cp "$PKG_ROOT"/.gitignore "$INSTALL_ROOT"/

mkdir -p "$INSTALL_ROOT/Databases"
cp "$PKG_ROOT"/Databases/Humannet_Library_V1_ungrouped.csv "$INSTALL_ROOT/Databases/"

cat > "$PREFIX/bin/humanet-pipeline" <<EOF
#!/usr/bin/env bash
exec "\$CONDA_PREFIX/bin/python" "$INSTALL_ROOT/pipeline.py" "\$@"
EOF

cat > "$PREFIX/bin/humanet-bootstrap-pubchem-sqlite" <<EOF
#!/usr/bin/env bash
exec "\$CONDA_PREFIX/bin/python" "$INSTALL_ROOT/bootstrap_pubchem_sqlite.py" "\$@"
EOF

cat > "$PREFIX/bin/humanet-download-pubchem-source-files" <<EOF
#!/usr/bin/env bash
exec "\$CONDA_PREFIX/bin/python" "$INSTALL_ROOT/download_pubchem_source_files.py" "\$@"
EOF

cat > "$PREFIX/bin/humanet-build-pubchem-sqlite" <<EOF
#!/usr/bin/env bash
exec "\$CONDA_PREFIX/bin/python" "$INSTALL_ROOT/build_pubchem_sqlite.py" "\$@"
EOF

cat > "$PREFIX/bin/humanet-semicolon-stage3-mapper" <<EOF
#!/usr/bin/env bash
exec "\$CONDA_PREFIX/bin/python" "$INSTALL_ROOT/semicolon_fuzzy_mapper.py" "\$@"
EOF

chmod +x "$PREFIX/bin/humanet-pipeline"
chmod +x "$PREFIX/bin/humanet-bootstrap-pubchem-sqlite"
chmod +x "$PREFIX/bin/humanet-download-pubchem-source-files"
chmod +x "$PREFIX/bin/humanet-build-pubchem-sqlite"
chmod +x "$PREFIX/bin/humanet-semicolon-stage3-mapper"
