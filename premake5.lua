-- Centromere Evolution Simulator
-- Build with: premake5 gmake && make config=release

workspace "CenSim"
    configurations { "Debug", "Release" }
    location "build"

project "censim"
    kind "ConsoleApp"
    language "C"
    cdialect "gnu99"
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

-- Single-threaded benchmark/profiling harness (no raylib, no pthreads). Always
-- built with symbols + frame pointers so it's ready for sampling profilers
-- (perf, samply, gperftools, valgrind). Build: premake5 gmake2 && make config=release
-- Profile e.g.: CPUPROFILE=cpu.prof LD_PRELOAD=libprofiler.so ./bin/Release/bench --scan-only 20
-- Tip: clang gives notably better AVX2 codegen here -- build with: make config=release CC=clang
-- (the whole app can be built that way too; CI/release stays on gcc for a portable AppImage)
project "bench"
    kind "ConsoleApp"
    language "C"
    cdialect "gnu99"
    targetdir "bin/%{cfg.buildcfg}"
    objdir "build/obj/%{cfg.buildcfg}"
    includedirs { "src" }
    files { "bench/bench.c", "src/simulation.c", "src/hor.c", "src/hist.c" }
    symbols "On"
    buildoptions { "-fno-omit-frame-pointer" }

    filter "configurations:Debug"
        defines { "DEBUG" }
        optimize "Debug"

    filter "configurations:Release"
        defines { "NDEBUG" }
        optimize "Speed"

    filter "system:macosx"
        architecture "arm64"

    filter "system:linux"
        architecture "x86_64"
        links { "m" }
