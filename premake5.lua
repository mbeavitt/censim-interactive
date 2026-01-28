-- Centromere Evolution Simulator
-- Build with: premake5 gmake2 && make config=release

workspace "CenSim"
    configurations { "Debug", "Release" }
    architecture "arm64"
    location "build"

project "censim"
    kind "ConsoleApp"
    language "C"
    cdialect "C99"
    targetdir "bin/%{cfg.buildcfg}"
    objdir "build/obj/%{cfg.buildcfg}"

    files {
        "src/**.h",
        "src/**.c"
    }

    -- raylib via Homebrew
    includedirs {
        "/opt/homebrew/include"
    }
    libdirs {
        "/opt/homebrew/lib"
    }

    -- macOS frameworks required by raylib
    links {
        "raylib",
        "OpenGL.framework",
        "Cocoa.framework",
        "IOKit.framework",
        "CoreVideo.framework"
    }

    filter "configurations:Debug"
        defines { "DEBUG" }
        symbols "On"
        optimize "Debug"

    filter "configurations:Release"
        defines { "NDEBUG" }
        optimize "Speed"
        symbols "Off"

    filter "system:macosx"
        defines { "_MACOSX" }
