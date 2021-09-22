## JF Nested Dielectric

Arnold Shader
Copyright (c) 2014, Psyop Media Company, LLC and Jonah Friedman

By Jonah Friedman
jfriedman@psyop.tv

Binaries builds by Jonah Friedman provided for Windows 64, for 4.0.x and 4.1.x.
4.1 binaries will work in Arnold 4.2 as well. Linux builds by Vladimir Jankijevic and
Benoit Leveau. OS X build by Jan Walter.

Contributions, questions, and suggestions are also welcomed.

Based on the 2002 paper by Charles M. Schmidt and Brian Budge, Simple Nested Dielectrics
in Ray Traced Images available here:
	http://graphics.idav.ucdavis.edu/~bcbudge/deep/research/nested_dielectrics.pdf

See license.txt for license.

### Use

Full documentation: https://github.com/Psyop/jf-nested-dielectric/raw/master/documentation/jf_nested_dielectric_manual.pdf

### Installation

Download the entire repo as a zip, using the "download as
zip" button on the right sidebar of the github page. Then move the appropriate binaries
and helper files (found under source) into the right directories of Softimage or Maya.
Note: If you instead download individual files from github, you may accidentally
download HTML files instead of the actual raw files, which of course won't work.

### Building from Source

To build:

```bash
git clone https://github.com/Psyop/jf-nested-dielectric
cd jf-nested-dielectric
mkdir build
cd build
cmake .. -DCMAKE_PREFIX_PATH=/path/to/install_root
make install
```

### Contributors

 - Benoit Leveau
 - Vladimir Jankijevic
 - Jan Walter
 - Chad Dombrova 
 - Daniel Hennies
 - Vahan Sosoyan  

### Special Thanks To

 - Andy Jones
 - Tony Barbieri
 - Todd Akita
 - Nisa Foster


### Release Notes

#### 1.0.1

- Fixes a bug causing a crash when there are two perfectly overlapping pieces
  of geometry, such as duplicate geometry.

#### 1.0.2

- Fixed a bug where shadow rays were evaluating wrong if they originated from
  inside nested dielectric, for example a standard shader embedded in some nested
  dielectric glass.
- Fixed a crash that occurred with perfectly overlapping geometry inside of a nested
  surface, and fixed a crash with huge numbers of overlapping geometries.

#### 1.0.3

- Fixed a new bug in Arnold 4.2 where Cook-Torrance speculars would render
  incorrectly	(due to the inclusion of Anisotropic Cook-Torrance.)

- Note: Arnold 4.1 is binary compatible with Arnold 4.2, so there are no seperate
  binaries for 4.2.

#### 1.0.4

- Added specular light contributions to direct speculars. Direct refraction is not affected. 

#### 1.0.5

- Fixed issue reported by Vahan Sosoyan where a glass of liquid could look wrong in reflection. 
- Added Cinema 4D helper files by Daniel Hennies

To build on Linux:
Make sure you've added $ARNOLD_PATH to your environment variables
Make sure $ARNOLD_PATH/bin is in your PATH:
