# orthonormal_to_perspective
This code converts a depth map under orthonormal projection into perspective one.
It is based on [normal_integration](https://github.com/yqueau/normal_integration)

## 1. Requirements

This code has following dependency:

0) MATLAB (code was tested on R2019a)

## 2. Input

- A depth map `z_orth` under orthonormal projection
- The 3x3 intrinsic paramter matrix `K` of the image size.
- A binary `mask` describing the object of interest in the image.

## 3. Usage
Download the demo dataset by `download.sh` in `data` folder. Then run the `demo_ortho2Persp.m` to get the depth map `z_persp` under perspective projection.

## 4. Reference
[1] "Variational Methods for Normal Integration", Quéau et al., Journal of Mathematical Imaging and Vision 60(4), pp 609--632, 2018. (Arxiv preprint: https://arxiv.org/abs/1709.05965)
