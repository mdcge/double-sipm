# Getting started

1. [Install nix](https://nixos.org/download.html):
   ```
   sh <(curl -L https://nixos.org/nix/install) --daemon
   ```
2. Clone this repository: `git clone https://github.com/mdcge/double-sipm.git`.
3. `cd` into it.
4. Type `nix develop --extra-experimental-features 'nix-command flakes'` (this may take some time as Geant4 will be downloaded and compiled).
5. Type `just double-sipm/` (the trailing `/` is important).

   This will compile and run the example and should take a few seconds. At this point, an image of a small simulated detector should appear. 
   
If all of this works, you now have an environment in which you can execute the Geant4 examples, and develop and run your own Geant4 code.

# Caveats

1. **Linux**: this is developed and tested on Linux. If the above instructions fail to work for you on Linux, this is considered a bug: please report it.

2. **Windows**: use WSL2, where it should work just as well as on Linux.

3. **MacOS**: the support for Geant4 in Nix (on MacOS) was marked as broken in the past. This appears to longer be the case, but we have not verified to what extent it works.

# Ergonomics

## Automatic environment switching with `direnv`

As is stands you have to write `nix develop` in order to activate the environment necessary to run and develop this code. What is more, `nix develop` places you in a minimal bash shell in which any personal configurations you may be used to will be missing.

Both of these problems can be fixed with [direnv](https://direnv.net/) which:
  * automatically enables the environment when you enter the directory
  * places you in your usual shell with all your preset configurations

To use `direnv`:

1. Make sure that it is [installed](https://direnv.net/docs/installation.html) on your system.

2. Don't forget to [hook](https://direnv.net/docs/hook.html) it into your shell.

Depending on which shell you are using, this will involve adding one of the following lines to the end of your shell configuration file:
``` shell
eval "$(direnv hook bash)"  # in ~/.bashrc
eval "$(direnv hook zsh)"   # in ~/.zshrc
eval `direnv hook tcsh`     # in ~/.cshrc
```

# Geant 4 configuration

Various configuration options of Geant4 itself can be changed by editing `flake.nix` here: 

``` nix
(geant4.override {
  enableMultiThreading = false;
  enableInventor       = false;
  enableQt             = true;
  enableXM             = false;
  enableOpenGLX11      = true;
  enablePython         = false;
  enableRaytracerX11   = false;
})
```

If you change the Geant4 configuration (if you are using `direnv`, it will notice the change, and automatically switch to the new configuration (recompiling Geant4, if this is a configuration not seen before) at your next shell prompt).

Be sure to expunge any examples you had compiled with a differently-configured Geant4, otherwise you may get mysterious problems.

# Activating the environment

If you have [direnv](https://direnv.net/) installed and configured, it will automatically activate the development environment when you enter this directory (though you will have to approve it the first time with `direnv allow`) and disable it when you leave the directory.

Without `direnv` you can manually enable the environment with `nix develop`, and disable it by quitting the shell started by `nix develop`.


# Something to try, that should work

Once the environment has been activated, try typing

``` shell
just run B1
```

This should copy the sources of the most basic example that is distributed with Geant4, into the local directory, configure it, compile it and execute it.

If all goes well, an image of a detector should appear. Try typing `/run/beamOn 10` in the `Session` box, and some simulated events should appear on the image.

You should be able to modify the source code (for example increase the value of `env_sizeZ` in `B1/src/DetectorConstruction.cc` (change it from 30 to 130, to make the change obvious)) and run your modified version by repeating the earlier command `just run B1`.

# Running other examples

Many of the other examples can be run in the same way: `just run <example-name>`. 

+ Some of them will fail because necessary libraries have not (yet) been provided in this flake. 

+ Others will fail because the internal organization of the example differs from that of the simple ones, which is assumed by `just run`.

  In many of these cases it should be fairly easy to figure out how to compile and execute the example by hand. The procedure tends to me something like
  
  1. Create a `build` subirectory the example's top-level directory
  2. `cd` into the newly-created `build` directory
  3. `cmake ..`
  4. `make -j N` (where `N` specifies the number of processors you want to use during compilation)
  5. Find the executable which was produced by the previous step, and execute it by preceding its name with `./` In the case of the B1 example, this would be `./exampleB1`. 

Fixes to these problems may be provided eventually. Don't hold your breath.
