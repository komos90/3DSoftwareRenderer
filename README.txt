Dependencies:
    SDL2(2.0.7)

To Build on Windows:
    Download the above dependencies into the libs folder.
    run winbuild.bat
    Copy the required DLLs into the bin folder:
        SDL2.dll

    Copy (or make a symlink to) the res folder into the bin folder.

To Build on Linux:
	Install libsdl2-dev and clang using your distro's package manager
	Give build.sh execution permissions
	Run build.sh
	Copy (or make soft link to) the res folder into the bin folder.

Current key mappings:

	WASD:         Move forwards, left, back, and right.
	Left Key:     Rotate camera left
	Right Key:    Rotate camera right
    1:            Enable/Disable wirefrmes
    2:            Enable/Disable surface drawing
    Mouse:        rotate camera
    Alt-F4 or Esc:Closes the game

