#!/bin/bash
set -e

# Build CenSim.app for macOS
# Usage: ./packaging/build_app.sh

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
BUILD_DIR="$PROJECT_DIR/packaging/build"
APP_NAME="CenSim"
APP_BUNDLE="$BUILD_DIR/$APP_NAME.app"

echo "=== Building $APP_NAME.app ==="
echo "Project dir: $PROJECT_DIR"

# Clean previous build
rm -rf "$BUILD_DIR"
mkdir -p "$BUILD_DIR"

# Step 1: Build C binary (Release)
echo ""
echo "=== Step 1: Building C binary ==="
cd "$PROJECT_DIR"
premake5 gmake2
make -C build config=release clean
make -C build config=release

# Step 2: Create app bundle structure
echo ""
echo "=== Step 2: Creating app bundle structure ==="
mkdir -p "$APP_BUNDLE/Contents/MacOS"
mkdir -p "$APP_BUNDLE/Contents/Resources"
mkdir -p "$APP_BUNDLE/Contents/Frameworks"

# Step 3: Copy main binary
echo "Copying main binary..."
cp "$PROJECT_DIR/bin/Release/censim" "$APP_BUNDLE/Contents/MacOS/"

# Step 4: Bundle raylib dylib
echo "Bundling raylib..."
RAYLIB_PATH=$(otool -L "$APP_BUNDLE/Contents/MacOS/censim" | grep raylib | awk '{print $1}')
if [ -n "$RAYLIB_PATH" ] && [ -f "$RAYLIB_PATH" ]; then
    cp "$RAYLIB_PATH" "$APP_BUNDLE/Contents/Frameworks/"
    RAYLIB_NAME=$(basename "$RAYLIB_PATH")

    # Fix the binary to use bundled raylib
    install_name_tool -change "$RAYLIB_PATH" "@executable_path/../Frameworks/$RAYLIB_NAME" "$APP_BUNDLE/Contents/MacOS/censim"
    echo "  Bundled: $RAYLIB_NAME"
else
    echo "  Warning: Could not find raylib dylib"
fi

# Step 5: Build Python bundle with PyInstaller
echo ""
echo "=== Step 3: Building Python bundle ==="
VENV_DIR="$BUILD_DIR/venv"
python3 -m venv "$VENV_DIR"
source "$VENV_DIR/bin/activate"

pip install --upgrade pip
pip install pyinstaller
pip install -r "$PROJECT_DIR/scripts/requirements.txt"

# Create PyInstaller bundle
cd "$BUILD_DIR"
pyinstaller --onedir --name visualize_umap \
    --distpath "$APP_BUNDLE/Contents/Resources" \
    --workpath "$BUILD_DIR/pyinstaller_work" \
    --specpath "$BUILD_DIR" \
    --noconfirm \
    "$PROJECT_DIR/scripts/visualize_umap.py"

deactivate

# Step 6: Create app icon
echo ""
echo "=== Step 4: Creating app icon ==="
ICON_SOURCE="$SCRIPT_DIR/icon.png"
if [ -f "$ICON_SOURCE" ]; then
    ICONSET_DIR="$BUILD_DIR/CenSim.iconset"
    mkdir -p "$ICONSET_DIR"

    sips -z 16 16 "$ICON_SOURCE" --out "$ICONSET_DIR/icon_16x16.png" >/dev/null
    sips -z 32 32 "$ICON_SOURCE" --out "$ICONSET_DIR/icon_16x16@2x.png" >/dev/null
    sips -z 32 32 "$ICON_SOURCE" --out "$ICONSET_DIR/icon_32x32.png" >/dev/null
    sips -z 64 64 "$ICON_SOURCE" --out "$ICONSET_DIR/icon_32x32@2x.png" >/dev/null
    sips -z 128 128 "$ICON_SOURCE" --out "$ICONSET_DIR/icon_128x128.png" >/dev/null
    sips -z 256 256 "$ICON_SOURCE" --out "$ICONSET_DIR/icon_128x128@2x.png" >/dev/null
    sips -z 256 256 "$ICON_SOURCE" --out "$ICONSET_DIR/icon_256x256.png" >/dev/null
    sips -z 512 512 "$ICON_SOURCE" --out "$ICONSET_DIR/icon_256x256@2x.png" >/dev/null
    sips -z 512 512 "$ICON_SOURCE" --out "$ICONSET_DIR/icon_512x512.png" >/dev/null
    sips -z 1024 1024 "$ICON_SOURCE" --out "$ICONSET_DIR/icon_512x512@2x.png" >/dev/null

    iconutil -c icns "$ICONSET_DIR" -o "$APP_BUNDLE/Contents/Resources/AppIcon.icns"
    echo "  Created AppIcon.icns"
else
    echo "  Warning: No icon.png found in packaging/"
fi

# Step 7: Create Info.plist
echo ""
echo "=== Step 5: Creating Info.plist ==="
cat > "$APP_BUNDLE/Contents/Info.plist" << 'EOF'
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0">
<dict>
    <key>CFBundleName</key>
    <string>CenSim</string>
    <key>CFBundleDisplayName</key>
    <string>CenSim</string>
    <key>CFBundleIdentifier</key>
    <string>com.censim.app</string>
    <key>CFBundleVersion</key>
    <string>1.0.0</string>
    <key>CFBundleShortVersionString</key>
    <string>1.0</string>
    <key>CFBundleExecutable</key>
    <string>censim</string>
    <key>CFBundlePackageType</key>
    <string>APPL</string>
    <key>CFBundleIconFile</key>
    <string>AppIcon</string>
    <key>LSMinimumSystemVersion</key>
    <string>11.0</string>
    <key>NSHighResolutionCapable</key>
    <true/>
    <key>NSSupportsAutomaticGraphicsSwitching</key>
    <true/>
</dict>
</plist>
EOF

# Step 8: Create launcher script that sets up paths
echo ""
echo "=== Step 6: Creating launcher ==="
mv "$APP_BUNDLE/Contents/MacOS/censim" "$APP_BUNDLE/Contents/MacOS/censim-bin"
cat > "$APP_BUNDLE/Contents/MacOS/censim" << 'EOF'
#!/bin/bash
# Launcher script for CenSim
DIR="$(cd "$(dirname "$0")" && pwd)"
export CENSIM_RESOURCES="$DIR/../Resources"
exec "$DIR/censim-bin" "$@"
EOF
chmod +x "$APP_BUNDLE/Contents/MacOS/censim"

echo ""
echo "=== Build complete ==="
echo "App bundle: $APP_BUNDLE"
du -sh "$APP_BUNDLE"
echo ""
echo "To test: open \"$APP_BUNDLE\""
echo "To create DMG: ./packaging/create_dmg.sh"
