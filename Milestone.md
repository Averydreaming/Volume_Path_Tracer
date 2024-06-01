# Volume Path Tracer

#### Group Member: Ruan Hang, Yijin Guo, Yuhao Wang

## Current Progress

### Environment Setup

- We have tried for several structures and finally choose the \href{https://github.com/BachiLi/lajolla_public}{lajolla} as the basic renderer.

- Since lajolla finally output the rendered images as .exr file format, we visualize the images by \href{https://github.com/wkjarosz/hdrview}{hdrview}. We plan to change the output format to .png.

### Model for paricipanting media

The radiative transfer equation is made of four components: absorption, emission, in-scattering, and out-scattering. Given a ray inside the volume parametrized by distance $p(t)$, the radiance along the ray is modeled by the **radiative transfer equation**:

$$
\frac{d}{dt}L(p(t),\omega)=-(\sigma_a(p(t))+\sigma_s(p(t)))L(p(t),\omega)+L_e(p(t),\omega)+\sigma_s(p(t))\int_{S^2}\rho(p(t),\omega,\omega')L(p(t),\omega')d\omega'
$$

where $L$ is the radiance, $\sigma_a$ is the absorption coefficient, $\sigma_s$ is the scattering coefficient, $L_e$ is the (volumetric) emission, $\rho$ is the *phase function* that is akin to BSDFs in surface rendering, and $S^2$ is the spherical domain.


### Codes

We split the development into 6 steps and build 6 volumetric path tracers. Each has more features than the previous ones. Up to now, we have finished the first two steps. This two renderer is based on the assumptions as follows:

#### Single monochromatic absorption-only homogeneous volume

- There is only a single, homogeneous volume:  are constant.
- The volume does not scatter light: .
- The surfaces in the scene only emit lights and do not reflect/transmit lights.
- The volume is monochromatic: three color channels of  have the same values.

#### Single monochromatic homogeneous volume with absorption and single-scattering, no surface lighting

- There is only a single, homogeneous volume:  and are constant.
- The surfaces in the scene only emit lights and do not reflect/transmit lights.
- The volume is monochromatic: three color channels of  and  have the same values.
- Light only scatters (changes direction) once in the volume.

The codes especially for the two renderers can be found in `src/vol_path_tracing.h`. 

## Preliminary Reuslts 

After completing the first two render, we can obtain the corresponding rendering results. They are as follows:

![Example Image1](docs/test1.exr)

![Example Image1](docs/test2.exr)

The rendering scene file can be found in `/scenes/volpath_test`.


## Work Plan

### Introduction
 
There are still 4 steps remaining for our volume path tracer.

- Multiple monochromatic homogeneous volumes with absorption and multiple-scattering using only phase function sampling, no surface lighting

- Multiple monochromatic homogeneous volumes with absorption and multiple-scattering with both phase function sampling and next event estimation, no surface lighting

- Multiple monochromatic homogeneous volumes with absorption and multiple-scattering with both phase function sampling and next event estimation, with surface lighting

- Multiple chromatic heterogeneous volumes with absorption and multiple-scattering with both phase function sampling and next event estimation, with surface lighting

Furthermore, we plan to generate the rendering scene file of smoke and render beautiful images. We consider to generate animations if there's enough time remaining. 

### Specific Plan

**June 1 - June 9**

- Do something for the output format. For example, use the .png instead of .exr.

- Understand all the remaining steps and complete the codes for the following 3 steps.

**June 10 - June 16**

- Complete the codes for the last step.

- Generate the scene file of smoke and get the rendered results.

- Try to generate the videos if there's enough time.

- Finish the final report and prepare for the final presentation.