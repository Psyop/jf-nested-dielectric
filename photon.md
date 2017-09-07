
### JF Photon

Disclaimer: These are the instructions I gave to very brave alpha testers. Use with care, this will happily fill your hard drive or network storage with photons. Contributions are welcome, especially ones that improve workflow. 

### Overview:

This does caustics using photon mapping. To use it you need a two step render, using perhaps two render layers. Put JF Nested Dielectric on whatever you want to be made of glass/water/whatever, and then put JF_Photon on whatever you want to catch photons. 

Make a layer that makes photons (Photon generation pass), and one that reads them after they're made. Render the pass that makes photons using the camera as a "light", then you can switch to different passes that read those photons and render them into your scene. 

## General:

* Put a [Frame] token in your file name and set the frame to -1 to do auto frame. 

## Photon generation pass:
* Make a layer where the camera is where you want the light to be (think of the camera field of view as a spot light). 
* Set up shaders:
  * Put the JF photon shader on everything that you want photons to fall on. 
  * Leave JFND on anything you want to cast caustics. 
  * Put a black hole on everything else
* Set up render options
  * Turn off filtering (use box filter). Otherwise, the outline of the buckets will appear in your photon map
  * Hide all lights
* Leave merging on. 10x reduction in size on disk, and load time is usual with reduction, sometimes higher (prism render is 20-30x!). Write radius determines how small a detail gets preserved.
* Use AA samples to effect overall photons, glossy samples for reflected photons, and refraction samples for refracted. These values are normalized so more photons != brighter.
* Ray depth in this pass controls ray depth of photons
* Caustics settings on JFND do not work yet, it just uses whatever the settings of the shader are
* The alternate lighting tools has some features to make kick render more suitable to generating photons

## Photon reader pass:
  * Plug it into emission.
* Sometimes it seems some photons end up in indirect diffuse, use a front back switch to fix this.
* Has a diffuse color, so if you have a diffuse map you want to use on the photons, plug it in here.
* Octree visualization is cool but not useful for normal photon effects.
* Exposure controls overall brightness. 
* Good idea to set JFND to cast black shadows and no caustics in this layer. 