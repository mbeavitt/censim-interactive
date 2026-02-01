#!/bin/bash
set -e

# Create DMG for CenSim.app
# Usage: ./packaging/create_dmg.sh

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BUILD_DIR="$SCRIPT_DIR/build"
APP_BUNDLE="$BUILD_DIR/CenSim.app"
DMG_NAME="CenSim-1.0.0"
DMG_PATH="$BUILD_DIR/$DMG_NAME.dmg"

if [ ! -d "$APP_BUNDLE" ]; then
    echo "Error: $APP_BUNDLE not found. Run build_app.sh first."
    exit 1
fi

echo "=== Creating DMG ==="

# Remove old DMG if exists
rm -f "$DMG_PATH"

# Check if create-dmg is available (brew install create-dmg)
if command -v create-dmg &> /dev/null; then
    echo "Using create-dmg for prettier DMG..."
    create-dmg \
        --volname "CenSim" \
        --volicon "$APP_BUNDLE/Contents/Resources/AppIcon.icns" 2>/dev/null || true \
        --window-pos 200 120 \
        --window-size 600 400 \
        --icon-size 100 \
        --icon "CenSim.app" 150 190 \
        --hide-extension "CenSim.app" \
        --app-drop-link 450 185 \
        "$DMG_PATH" \
        "$APP_BUNDLE"
else
    echo "Using hdiutil (install create-dmg for prettier result)..."

    # Create temporary directory for DMG contents
    DMG_TEMP="$BUILD_DIR/dmg_temp"
    rm -rf "$DMG_TEMP"
    mkdir -p "$DMG_TEMP"

    # Copy app bundle
    cp -R "$APP_BUNDLE" "$DMG_TEMP/"

    # Create Applications symlink
    ln -s /Applications "$DMG_TEMP/Applications"

    # Create DMG
    hdiutil create -volname "CenSim" -srcfolder "$DMG_TEMP" -ov -format UDZO "$DMG_PATH"

    # Cleanup
    rm -rf "$DMG_TEMP"
fi

echo ""
echo "=== DMG created ==="
echo "DMG: $DMG_PATH"
echo ""
echo "To test: open \"$DMG_PATH\""
