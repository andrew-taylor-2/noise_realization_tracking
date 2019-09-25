## Purpose
This function has been created to offer an alternative method of "probabilistic" tractography based on creating and tracking from multiple realizations of the noise distribution.
The pipeline as it is will create the "simulated" datasets, denoise again, run MRtrix's 5ttgen, get diffusion tensors and FA, get kurtosis orientation distribution functions, track from these, combine the tract files, and generate a connectome based off of an MNI-space atlas provided.

## How to run 
Arguments:
	`dwi_folder`: This folder should contain a single DWI in the form of a .nii file and valid bvectors and bvalues files with a matching number of elements. The extensions for the bvecs and bvals files can be .bvec/.bvecs and .bval/.bvals.
	`t1_folder`: A folder containing a T1 .nii file of the chosen subject.
	`out_dir`: A directory in which to place the outputs
	`nsim`: The number of noise realizations to write. The time for this function to run is linearly related to `nsim`
	`atlas4connectome`: An atlas *in MNI space* that will be used to generate a connectome. This function will use SPM's "Old Normalise" to warp the atlas into subject space based on "MNI152_T1_1mm.nii.gz" included with FSL.

## Dependencies
- SPM12 must be on your Matlab path
- d2n2s, which can be cloned with `git clone https://github.com/andrew-taylor-2/d2n2s.git`, must be on your Matlab path
- FSL must be installed on your system and MRtrix must be able to access it
- MRtrix must be installed on your system and callable from the command line (this may necessitate starting Matlab from a shell
- kODF_nii_preprocess and its dependencies must be in your Matlab path (this will be added to the repository soon, but can also be found at github.com/neurolabusc/nii_preprocess on the DKI_pipeline branch)

## Notes
- Currently, this function only works on Mac and Unix, but it can be made to work on Windows very easily
- Using an 88x88x50 DWI with 10 b0s, 64 b1000s, and 64 b2000s; a 192x256x256 T1; and an `nsim` of 5 on my iMac Pro with 8 logical cores, this took 45 minutes to complete

