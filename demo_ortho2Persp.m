%%
clear
clc

addpath('utils');
%% load dataset
load('data/o2p_demo.mat');
% z_orth: the depth map under orthonormal projection.
% mask: the mask of the object
% K: the intrinsic matrix

N = Depth2Normals_o2p(z_orth, size(mask), mask, 'CNC');
z_persp = normals2DepthPersp(N, mask, K);