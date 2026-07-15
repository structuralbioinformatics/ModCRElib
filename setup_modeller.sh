#!/bin/bash

# 1. Ensure a Conda environment is actively running
if [ -z "$CONDA_PREFIX" ]; then
    echo "❌ Error: No active Conda environment detected!"
    echo "Please activate your environment first (e.g., conda activate modcrelib_env)"
    exit 1
fi

# 2. Dynamically extract the environment root and the Python version
TOPDIR="$CONDA_PREFIX"
PYTHON_VERSION=$(python -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')")

# 3. Establish the standardized Modeller directory structure
TARGET_DIR="${TOPDIR}/lib/modeller-10.8/bin"
mkdir -p "$TARGET_DIR"

echo "⚙️  Configuring environment pathways..."
echo "   - Conda Root: $TOPDIR"
echo "   - Python Target: $PYTHON_VERSION"
echo "   - Output Path: $TARGET_DIR/modpy.sh"

# 4. Generate the modpy.sh script injection block
cat << EOF > "$TARGET_DIR/modpy.sh"
#!/bin/sh

# 1. Force the dynamic script to utilize your Conda environment root
TOPDIR="${TOPDIR}"

# 2. Tell the script where the compiled C-libraries (.dylib files) live in Conda
LLP="\${TOPDIR}/lib"

if test -z "\${LD_LIBRARY_PATH}"; then
  LD_LIBRARY_PATH=\${LLP}
else
  LD_LIBRARY_PATH=\${LD_LIBRARY_PATH}:\${LLP}
fi

if test -z "\${DYLD_LIBRARY_PATH}"; then
  DYLD_LIBRARY_PATH=\${LLP}
else
  DYLD_LIBRARY_PATH=\${DYLD_LIBRARY_PATH}:\${LLP}
fi

if test -z "\${LIBPATH}"; then
  LIBPATH=\${LLP}
else
  LIBPATH=\${LIBPATH}:\${LLP}
fi

# 3. Direct Python to where Conda unpacked MODELLER's Python modules.
PP="\${TOPDIR}/lib/python${PYTHON_VERSION}/site-packages"

if test -z "\${PYTHONPATH}"; then
  PYTHONPATH=\${PP}
else
  ORIGPYPATH="\${PYTHONPATH}"
  PYTHONPATH=\${PYTHONPATH}:\${PP}
fi

export LD_LIBRARY_PATH DYLD_LIBRARY_PATH PYTHONPATH LIBPATH ORIGPYPATH

# 4. Drop the redundant 'python' argument passed by the master program
shift

# 5. Force execution through Rosetta 2 using the Conda Intel Python binary
exec arch -x86_64 "\${TOPDIR}/bin/python" "\$@"
EOF

# 5. Unlock executable permissions flags on the newly created script
chmod +x "$TARGET_DIR/modpy.sh"

echo "✅ Generation Complete! modpy.sh is successfully deployed and ready for use."

