
function keyDescriptor_m = meisam_mirror_SIFT_descriptor(keyDescriptor)
% Use lookup table to generate the SIFT descriptor for the mirrored
% version of the image patch.  This will only work for SIFT, and requires
% the precise ordering returned by Lowes SIFT code (version 4).
%
% INPUT:  keyDescriptor    - 128 element SIFT descriptor vector of some
%                            image patch P.
% OUTPUT: keyDescriptor_m  - 128 element SIFT descriptor vector of some
%                            image patch P reflected about the y-axis.

ind(1:8:25) = [97:8:121];
ind(33:8:33+24) = [65:8:65+24];
ind(65:8:65+24) = [33:8:33+24];
ind(97:8:121) = [1:8:25];

ind(2:2+6) = [98+6:-1:98];
ind(10:10+6) = [106+6:-1:106];
ind(18:18+6) = [114+6:-1:114];
ind(26:26+6) = [122+6:-1:122];

ind(34:34+6) = [66+6:-1:66];
ind(42:42+6) = [74+6:-1:74];
ind(50:50+6) = [82+6:-1:82];
ind(58:58+6) = [90+6:-1:90];

ind(66:66+6) = [34+6:-1:34];
ind(74:74+6) = [42+6:-1:42];
ind(82:82+6) = [50+6:-1:50];
ind(90:90+6) = [58+6:-1:58];

ind(98:98+6) = [2+6:-1:2];
ind(106:106+6) = [10+6:-1:10];
ind(114:114+6) = [18+6:-1:18];
ind(122:122+6) = [26+6:-1:26];

keyDescriptor_m = keyDescriptor(:,ind);

end