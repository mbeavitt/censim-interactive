-- Centromere Evolution Simulator
-- Build with: premake5 gmake2 && make config=release

workspace "CenSim"
    configurations { "Debug", "Release" }
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

    filter "configurations:Debug"
        defines { "DEBUG" }
        symbols "On"
        optimize "Debug"

    filter "configurations:Release"
        defines { "NDEBUG" }
        optimize "Speed"
        symbols "Off"

    -- macOS specific
    filter "system:macosx"
        architecture "arm64"
        defines { "_MACOSX" }
        includedirs { "/opt/homebrew/include" }
        libdirs { "/opt/homebrew/lib" }
        links {
            "raylib",
            "OpenGL.framework",
            "Cocoa.framework",
            "IOKit.framework",
            "CoreVideo.framework"
        }

    -- Linux specific
    filter "system:linux"
        architecture "x86_64"
        defines { "_LINUX" }
        includedirs { "/usr/local/include" }
        libdirs { "/usr/local/lib" }
        links {
            "raylib",
            "GL",
            "m",
            "pthread",
            "dl",
            "rt",
            "X11"
        }
